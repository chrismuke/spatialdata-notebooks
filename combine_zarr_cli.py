#!/usr/bin/env python3
"""
Combine Zarr Files CLI Tool

This tool combines multiple SpatialData zarr files into a single zarr file,
with spatial translation to prevent overlapping. By default, zarr files are
arranged vertically (stacked down), with optional border spacing.

Based on SpatialData concatenate and transformation functionality.
"""

import logging
import sys
from pathlib import Path
from datetime import datetime
from typing import List, Tuple
import json

import click
import numpy as np
import spatialdata as sd
from spatialdata.transformations import Translation, get_transformation, set_transformation


def setup_logging(level: str = "INFO") -> None:
    """Set up logging configuration."""
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )


def get_spatial_bounds(sdata: sd.SpatialData, coordinate_system: str = "global") -> Tuple[float, float, float, float]:
    """
    Get spatial bounds (min_x, min_y, max_x, max_y) of a SpatialData object.
    
    Args:
        sdata: SpatialData object
        coordinate_system: Coordinate system to use for bounds calculation
        
    Returns:
        Tuple of (min_x, min_y, max_x, max_y)
    """
    min_x, min_y = float('inf'), float('inf')
    max_x, max_y = float('-inf'), float('-inf')
    
    # Check shapes for bounds
    for shape_name, shape_data in sdata.shapes.items():
        if coordinate_system in sdata.coordinate_systems:
            bounds = shape_data.total_bounds  # [minx, miny, maxx, maxy]
            min_x = min(min_x, bounds[0])
            min_y = min(min_y, bounds[1])
            max_x = max(max_x, bounds[2])
            max_y = max(max_y, bounds[3])
    
    # Check images for bounds
    for image_name, image_data in sdata.images.items():
        # Get image dimensions - assuming last two dimensions are spatial (y, x)
        if hasattr(image_data, 'sizes'):
            if 'y' in image_data.sizes and 'x' in image_data.sizes:
                height, width = image_data.sizes['y'], image_data.sizes['x']
                min_x = min(min_x, 0)
                min_y = min(min_y, 0)
                max_x = max(max_x, width)
                max_y = max(max_y, height)
    
    # Check labels for bounds
    for label_name, label_data in sdata.labels.items():
        if hasattr(label_data, 'sizes'):
            if 'y' in label_data.sizes and 'x' in label_data.sizes:
                height, width = label_data.sizes['y'], label_data.sizes['x']
                min_x = min(min_x, 0)
                min_y = min(min_y, 0)
                max_x = max(max_x, width)
                max_y = max(max_y, height)
    
    # Check points for bounds  
    for point_name, point_data in sdata.points.items():
        if 'x' in point_data.columns and 'y' in point_data.columns:
            try:
                # Force computation for dask arrays
                points_min_x = point_data['x'].min().compute() if hasattr(point_data['x'].min(), 'compute') else point_data['x'].min()
                points_min_y = point_data['y'].min().compute() if hasattr(point_data['y'].min(), 'compute') else point_data['y'].min()
                points_max_x = point_data['x'].max().compute() if hasattr(point_data['x'].max(), 'compute') else point_data['x'].max()
                points_max_y = point_data['y'].max().compute() if hasattr(point_data['y'].max(), 'compute') else point_data['y'].max()
                
                min_x = min(min_x, points_min_x)
                min_y = min(min_y, points_min_y)
                max_x = max(max_x, points_max_x)
                max_y = max(max_y, points_max_y)
                
                logging.debug(f"Points {point_name} bounds: x=({points_min_x:.1f}, {points_max_x:.1f}), y=({points_min_y:.1f}, {points_max_y:.1f})")
            except Exception as e:
                logging.warning(f"Could not compute bounds for points {point_name}: {e}")
                # Skip this points dataset for bounds calculation
    
    # Handle case where no spatial elements were found
    if min_x == float('inf'):
        min_x, min_y, max_x, max_y = 0, 0, 1000, 1000
        logging.warning("No spatial bounds found, using default bounds (0, 0, 1000, 1000)")
    
    return min_x, min_y, max_x, max_y


def apply_translation_to_sdata(sdata: sd.SpatialData, translation_x: float, translation_y: float, 
                              coordinate_system: str = "global") -> sd.SpatialData:
    """
    Apply translation by DIRECTLY modifying coordinates for all spatial elements.
    This ensures actual coordinate separation that works reliably in napari.
    
    Args:
        sdata: SpatialData object to translate
        translation_x: Translation in x direction
        translation_y: Translation in y direction
        coordinate_system: Target coordinate system
        
    Returns:
        SpatialData object with translation applied
    """
    if translation_x == 0 and translation_y == 0:
        return sdata
        
    logging.info(f"Applying direct coordinate translation ({translation_x:.1f}, {translation_y:.1f}) to all elements")
    
    # Handle points by directly modifying coordinates
    for point_name, point_data in sdata.points.items():
        logging.info(f"Directly translating coordinates for points: {point_name}")
        
        if 'x' in point_data.columns and 'y' in point_data.columns:
            try:
                # Modify coordinates directly to ensure actual position changes
                point_data['x'] = point_data['x'] + translation_x
                point_data['y'] = point_data['y'] + translation_y
                
                logging.info(f"Successfully translated coordinates for {point_name}")
            except Exception as e:
                logging.warning(f"Failed to translate coordinates for {point_name}: {e}")
    
    # Handle shapes by directly modifying geometry coordinates
    for shape_name, shape_data in sdata.shapes.items():
        logging.info(f"Directly translating geometry for shapes: {shape_name}")
        try:
            # For GeoPandas GeoDataFrame, translate the geometry
            if hasattr(shape_data, 'geometry'):
                # Apply translation to all geometries
                from shapely.affinity import translate
                shape_data['geometry'] = shape_data['geometry'].apply(
                    lambda geom: translate(geom, xoff=translation_x, yoff=translation_y)
                )
                logging.info(f"Successfully translated geometry for {shape_name}")
        except Exception as e:
            logging.warning(f"Failed to translate geometry for {shape_name}: {e}")
    
    # For images and labels, apply transformation metadata (can't modify pixel grids directly)
    # But first, let's try to understand their coordinate relationship
    translation = Translation([translation_x, translation_y], axes=("x", "y"))
    
    for element_type in ['images', 'labels']:
        element_dict = getattr(sdata, element_type)
        for element_name, element_data in element_dict.items():
            logging.info(f"Applying transformation to {element_type}: {element_name}")
            try:
                # Get existing transformation and compose with translation
                existing_transform = get_transformation(element_data, to_coordinate_system=coordinate_system)
                
                # Create a sequence of existing transform + translation
                from spatialdata.transformations import Sequence
                new_transform = Sequence([existing_transform, translation])
                
                # Set the new transformation
                set_transformation(element_data, new_transform, to_coordinate_system=coordinate_system)
                logging.info(f"Successfully applied transformation to {element_name}")
            except Exception as e:
                logging.warning(f"Failed to apply transformation to {element_name}: {e}")
    
    return sdata


def calculate_layout_positions(zarr_files: List[Path], border: float = 0.0, 
                              layout: str = "vertical") -> List[Tuple[float, float]]:
    """
    Calculate translation positions for each zarr file to avoid overlap.
    
    Args:
        zarr_files: List of zarr file paths
        border: Border spacing between files
        layout: Layout arrangement ("vertical", "horizontal", "grid")
        
    Returns:
        List of (translation_x, translation_y) tuples for each file
    """
    positions = []
    current_x, current_y = 0.0, 0.0
    
    # First pass: calculate bounds for each file
    bounds_list = []
    for zarr_file in zarr_files:
        try:
            sdata = sd.read_zarr(zarr_file)
            bounds = get_spatial_bounds(sdata)
            bounds_list.append(bounds)
            logging.info(f"Bounds for {zarr_file.name}: {bounds}")
        except Exception as e:
            logging.error(f"Failed to read {zarr_file}: {e}")
            # Use default bounds
            bounds_list.append((0, 0, 1000, 1000))
    
    # Calculate positions based on layout
    if layout == "vertical":
        # For centering, calculate the overall X range and find the widest sample
        all_min_x = min(bounds[0] for bounds in bounds_list)
        all_max_x = max(bounds[2] for bounds in bounds_list)
        overall_width = all_max_x - all_min_x
        
        # Find the center point of the overall bounding box
        overall_center_x = (all_min_x + all_max_x) / 2
        
        # For proper vertical stacking with no overlap
        placement_y = 0.0  # Start placing files at Y=0
        
        for i, (min_x, min_y, max_x, max_y) in enumerate(bounds_list):
            # For vertical stacking with centering and no overlap
            file_width = max_x - min_x
            file_height = max_y - min_y
            file_center_x = (min_x + max_x) / 2
            
            # Center each file relative to the overall center
            translation_x = overall_center_x - file_center_x
            
            # Place file at current placement position, ensuring no overlap
            translation_y = placement_y - min_y  # Move to placement position
            
            positions.append((translation_x, translation_y))
            
            # Calculate where the next file should be placed (bottom of current file + border)
            file_bottom_after_translation = placement_y + file_height
            placement_y = file_bottom_after_translation + border
            
            logging.info(f"File {i+1} ({zarr_files[i].name}): translate by ({translation_x:.1f}, {translation_y:.1f})")
            logging.info(f"  Original bounds: x=({min_x:.1f}, {max_x:.1f}), y=({min_y:.1f}, {max_y:.1f})")
            logging.info(f"  File center: {file_center_x:.1f}, Overall center: {overall_center_x:.1f}")
            logging.info(f"  After translation: placed at Y={placement_y - border - file_height:.1f} to {placement_y - border:.1f}")
            logging.info(f"  Next file will start at Y={placement_y:.1f} (border={border:.1f})")
    
    elif layout == "horizontal":
        for i, (min_x, min_y, max_x, max_y) in enumerate(bounds_list):
            # For horizontal stacking, translate so files are side by side
            translation_x = current_x - min_x  # Stack horizontally
            translation_y = -min_y  # Align top edges
            
            positions.append((translation_x, translation_y))
            
            # Update current_x for next file (right edge of current file + border)
            width = max_x - min_x
            current_x += width + border
            
            logging.info(f"File {i+1} ({zarr_files[i].name}): translate by ({translation_x:.1f}, {translation_y:.1f})")
    
    elif layout == "grid":
        # Simple grid layout - calculate grid size
        n_files = len(zarr_files)
        grid_cols = int(np.ceil(np.sqrt(n_files)))
        grid_rows = int(np.ceil(n_files / grid_cols))
        
        # Calculate max dimensions for grid spacing
        max_width = max(bounds[2] - bounds[0] for bounds in bounds_list)
        max_height = max(bounds[3] - bounds[1] for bounds in bounds_list)
        
        for i, (min_x, min_y, max_x, max_y) in enumerate(bounds_list):
            row = i // grid_cols
            col = i % grid_cols
            
            translation_x = col * (max_width + border) - min_x
            translation_y = row * (max_height + border) - min_y
            
            positions.append((translation_x, translation_y))
            
            logging.info(f"File {i+1} ({zarr_files[i].name}): grid position ({row}, {col}), translate by ({translation_x:.1f}, {translation_y:.1f})")
    
    return positions


def combine_zarr_files(
    zarr_files: List[Path],
    output_path: Path,
    border: float = 0.0,
    layout: str = "vertical",
    coordinate_system: str = "global",
    overwrite: bool = False
) -> None:
    """
    Combine multiple zarr files into a single zarr file with spatial translation.
    
    Args:
        zarr_files: List of input zarr file paths
        output_path: Output zarr file path
        border: Border spacing between files
        layout: Layout arrangement ("vertical", "horizontal", "grid")
        coordinate_system: Coordinate system for alignment
        overwrite: Whether to overwrite existing output
    """
    if output_path.exists() and not overwrite:
        raise ValueError(f"Output file exists: {output_path}. Use --overwrite to overwrite.")
    
    logging.info(f"Combining {len(zarr_files)} zarr files with {layout} layout and {border} border")
    
    # Calculate translation positions
    positions = calculate_layout_positions(zarr_files, border, layout)
    
    # Load and translate each SpatialData object
    translated_sdatas = []
    suffixes = {}
    
    for i, (zarr_file, (tx, ty)) in enumerate(zip(zarr_files, positions)):
        logging.info(f"Processing file {i+1}/{len(zarr_files)}: {zarr_file.name}")
        
        try:
            # Load SpatialData object
            sdata = sd.read_zarr(zarr_file)
            
            # Apply translation if needed
            if tx != 0 or ty != 0:
                logging.info(f"Applying translation ({tx:.1f}, {ty:.1f}) to {zarr_file.name}")
                sdata = apply_translation_to_sdata(sdata, tx, ty, coordinate_system)
            
            translated_sdatas.append(sdata)
            
            # Create suffix for unique naming
            file_stem = zarr_file.stem
            if file_stem.endswith('.zarr'):
                file_stem = file_stem[:-5]  # Remove .zarr suffix
            suffixes[f"sample_{i+1}_{file_stem}"] = sdata
            
        except Exception as e:
            logging.error(f"Failed to process {zarr_file}: {e}")
            raise
    
    # Concatenate all SpatialData objects
    logging.info("Concatenating SpatialData objects...")
    try:
        # Try concatenation with tables first
        try:
            combined_sdata = sd.concatenate(
                suffixes,
                concatenate_tables=True,
                obs_names_make_unique=True
            )
        except ValueError as e:
            if "region_key" in str(e):
                logging.warning(f"Table concatenation failed due to region key mismatch: {e}")
                logging.warning("Concatenating without tables to preserve spatial data including transcripts")
                # Concatenate without tables to preserve spatial data (points, images, etc.)
                combined_sdata = sd.concatenate(
                    suffixes,
                    concatenate_tables=False,
                    obs_names_make_unique=True
                )
            else:
                raise
        
        logging.info(f"Combined SpatialData object created with {len(combined_sdata.coordinate_systems)} coordinate systems")
        
        # Save combined object
        logging.info(f"Saving combined zarr file to: {output_path}")
        combined_sdata.write(output_path)
        
        # Ensure the file is self-contained for napari compatibility
        logging.info("Verifying self-contained status for napari compatibility...")
        try:
            test_reload = sd.read_zarr(output_path)
            if test_reload.is_self_contained():
                logging.info("Combined zarr file is self-contained and ready for napari")
            else:
                logging.warning("Combined zarr file is not self-contained - may have issues with napari")
        except Exception as e:
            logging.warning(f"Could not verify self-contained status: {e}")
        
        # Print summary
        logging.info("Combination complete!")
        logging.info(f"Combined object contains:")
        logging.info(f"  - Images: {len(combined_sdata.images)}")
        logging.info(f"  - Labels: {len(combined_sdata.labels)}")
        logging.info(f"  - Shapes: {len(combined_sdata.shapes)}")
        logging.info(f"  - Points: {len(combined_sdata.points)}")
        logging.info(f"  - Tables: {len(combined_sdata.tables)}")
        
    except Exception as e:
        logging.error(f"Failed to concatenate SpatialData objects: {e}")
        raise


def save_combination_metadata(output_path: Path, zarr_files: List[Path], 
                            border: float, layout: str, positions: List[Tuple[float, float]]) -> None:
    """Save metadata about the combination process."""
    metadata_path = output_path.parent / f"{output_path.stem}_combination_metadata.json"
    
    metadata = {
        "timestamp": datetime.now().isoformat(),
        "input_files": [str(f) for f in zarr_files],
        "output_file": str(output_path),
        "layout": layout,
        "border": border,
        "translations": [{"file": str(f), "translation_x": tx, "translation_y": ty} 
                        for f, (tx, ty) in zip(zarr_files, positions)],
        "tool_version": "1.0.0"
    }
    
    with open(metadata_path, 'w') as f:
        json.dump(metadata, f, indent=2)
    
    logging.info(f"Combination metadata saved: {metadata_path}")


@click.command()
@click.argument("zarr_files", nargs=-1, required=True, type=click.Path(exists=True, path_type=Path))
@click.argument("output_path", type=click.Path(path_type=Path))
@click.option(
    "--border",
    type=float,
    default=0.0,
    help="Border spacing between zarr files in coordinate units.",
    show_default=True
)
@click.option(
    "--layout", 
    type=click.Choice(["vertical", "horizontal", "grid"]),
    default="vertical",
    help="Layout arrangement for combining files.",
    show_default=True
)
@click.option(
    "--coordinate-system",
    default="global",
    help="Coordinate system to use for alignment.",
    show_default=True
)
@click.option(
    "--overwrite",
    is_flag=True,
    help="Overwrite output file if it exists."
)
@click.option(
    "--log-level",
    type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR"], case_sensitive=False),
    default="INFO",
    help="Logging level.",
    show_default=True
)
@click.option(
    "--save-metadata",
    is_flag=True,
    help="Save combination metadata to JSON file."
)
def main(
    zarr_files: List[Path],
    output_path: Path,
    border: float,
    layout: str,
    coordinate_system: str,
    overwrite: bool,
    log_level: str,
    save_metadata: bool
):
    """
    Combine multiple SpatialData zarr files into a single zarr file.
    
    Applies spatial translations to prevent overlap. By default, files are
    arranged vertically (stacked down) with optional border spacing.
    
    ZARR_FILES: Input zarr files to combine (minimum 2 files required)
    
    OUTPUT_PATH: Path for the combined output zarr file
    
    Examples:
    
        # Combine 3 files vertically with 100 unit border
        combine_zarr_cli.py file1.zarr file2.zarr file3.zarr combined.zarr --border 100
        
        # Arrange files horizontally
        combine_zarr_cli.py *.zarr combined.zarr --layout horizontal
        
        # Grid layout with metadata saving
        combine_zarr_cli.py *.zarr combined.zarr --layout grid --save-metadata
    """
    # Setup logging
    setup_logging(log_level.upper())
    
    # Validate inputs
    zarr_files = list(zarr_files)
    if len(zarr_files) < 2:
        click.echo("Error: At least 2 zarr files are required for combination.", err=True)
        sys.exit(1)
    
    if output_path.exists() and not overwrite:
        click.echo(f"Error: Output file already exists: {output_path}. Use --overwrite to overwrite.", err=True)
        sys.exit(1)
    
    if border < 0:
        click.echo("Error: Border spacing must be non-negative.", err=True)
        sys.exit(1)
    
    # Log input parameters
    logging.info(f"Combining {len(zarr_files)} zarr files:")
    for i, zarr_file in enumerate(zarr_files, 1):
        logging.info(f"  {i}. {zarr_file}")
    logging.info(f"Output: {output_path}")
    logging.info(f"Layout: {layout}, Border: {border}, Coordinate system: {coordinate_system}")
    
    try:
        # Calculate positions for metadata if requested
        positions = []
        if save_metadata:
            positions = calculate_layout_positions(zarr_files, border, layout)
        
        # Combine zarr files
        combine_zarr_files(
            zarr_files=zarr_files,
            output_path=output_path,
            border=border,
            layout=layout,
            coordinate_system=coordinate_system,
            overwrite=overwrite
        )
        
        # Save metadata if requested
        if save_metadata:
            if not positions:  # Calculate positions if not done already
                positions = calculate_layout_positions(zarr_files, border, layout)
            save_combination_metadata(output_path, zarr_files, border, layout, positions)
        
        # Standard output
        print(f"Successfully combined {len(zarr_files)} zarr files into: {output_path}")
        
    except Exception as e:
        click.echo(f"❌ Error during combination: {e}", err=True)
        sys.exit(1)


if __name__ == "__main__":
    main()