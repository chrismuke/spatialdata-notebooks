#!/usr/bin/env python3
"""
Convert napari-roi-manager ROIs to roi_umap_analysis format

This script converts ROI files created with napari-roi-manager into the format
expected by roi_umap_analysis.py and integrates them with SpatialData.
"""

import json
import numpy as np
from pathlib import Path
from shapely.geometry import Polygon
from typing import Dict, List, Tuple
import spatialdata as sd
import click

def load_napari_roi_manager_json(json_path: str) -> Dict:
    """Load ROIs from napari-roi-manager JSON format"""
    with open(json_path, 'r') as f:
        roi_data = json.load(f)
    
    print(f"Loaded napari-roi-manager file with {len(roi_data['data'])} ROIs")
    return roi_data

def convert_to_standard_format(napari_roi_data: Dict) -> Dict:
    """Convert napari-roi-manager format to standard ROI format"""
    
    standard_rois = {}
    
    # Extract data from napari-roi-manager format
    roi_coordinates = napari_roi_data['data']
    shape_types = napari_roi_data['shape_type']
    roi_names = napari_roi_data.get('names', [f"ROI_{i+1}" for i in range(len(roi_coordinates))])
    
    print("Converting ROIs to standard format...")
    
    for i, (coords, shape_type, name) in enumerate(zip(roi_coordinates, shape_types, roi_names)):
        if shape_type == 'polygon':
            # Convert coordinates to numpy array
            coords_array = np.array(coords)
            
            # napari-roi-manager typically uses x,y format already
            # But let's verify and ensure proper format
            if coords_array.shape[1] != 2:
                print(f"Warning: ROI {name} has unexpected coordinate format")
                continue
            
            # Create polygon to calculate area and bounds
            try:
                polygon = Polygon(coords_array)
                
                # Ensure polygon is valid
                if not polygon.is_valid:
                    print(f"Warning: ROI {name} is not a valid polygon, attempting to fix...")
                    polygon = polygon.buffer(0)  # Try to fix invalid polygon
                
                if polygon.is_valid:
                    standard_rois[name] = {
                        'coordinates': coords_array.tolist(),
                        'area': float(polygon.area),
                        'bounds': list(polygon.bounds),  # [minx, miny, maxx, maxy]
                        'perimeter': float(polygon.length),
                        'shape_type': shape_type,
                        'created_with': 'napari-roi-manager',
                        'original_index': i
                    }
                    print(f"  ✓ Converted {name}: {len(coords)} points, area = {polygon.area:,.0f}")
                else:
                    print(f"  ✗ Skipped {name}: invalid polygon after repair attempt")
                    
            except Exception as e:
                print(f"  ✗ Error processing {name}: {e}")
        else:
            print(f"  ⚠ Skipped {name}: unsupported shape type '{shape_type}'")
    
    print(f"Successfully converted {len(standard_rois)} out of {len(roi_coordinates)} ROIs")
    return standard_rois

def save_standard_roi_json(roi_data: Dict, output_path: str):
    """Save ROIs in standard format"""
    # Create clean version for saving (remove extra metadata)
    clean_roi_data = {}
    for name, roi_info in roi_data.items():
        clean_roi_data[name] = {
            'coordinates': roi_info['coordinates'],
            'area': roi_info['area'],
            'bounds': roi_info['bounds']
        }
    
    with open(output_path, 'w') as f:
        json.dump(clean_roi_data, f, indent=2)
    
    print(f"Saved standard ROI format to: {output_path}")

def validate_rois_with_spatial_data(roi_data: Dict, zarr_path: str) -> Dict:
    """Validate ROIs against spatial data bounds"""
    print(f"Validating ROIs against spatial data: {zarr_path}")
    
    # Load spatial data to get bounds
    sdata = sd.read_zarr(zarr_path)
    
    # Get spatial data bounds
    data_min_x = data_min_y = float('inf')
    data_max_x = data_max_y = float('-inf')
    
    # Check points bounds
    for point_name, point_data in sdata.points.items():
        if 'x' in point_data.columns and 'y' in point_data.columns:
            # Handle potential dask arrays
            x_vals = point_data['x']
            y_vals = point_data['y']
            
            # Compute if dask
            if hasattr(x_vals.min(), 'compute'):
                x_min, x_max = x_vals.min().compute(), x_vals.max().compute()
                y_min, y_max = y_vals.min().compute(), y_vals.max().compute()
            else:
                x_min, x_max = x_vals.min(), x_vals.max()
                y_min, y_max = y_vals.min(), y_vals.max()
            
            data_min_x = min(data_min_x, x_min)
            data_max_x = max(data_max_x, x_max)
            data_min_y = min(data_min_y, y_min)
            data_max_y = max(data_max_y, y_max)
    
    # Check shapes bounds
    for shape_name, shape_data in sdata.shapes.items():
        bounds = shape_data.total_bounds  # [minx, miny, maxx, maxy]
        data_min_x = min(data_min_x, bounds[0])
        data_min_y = min(data_min_y, bounds[1])
        data_max_x = max(data_max_x, bounds[2])
        data_max_y = max(data_max_y, bounds[3])
    
    print(f"Spatial data bounds: x=({data_min_x:.1f}, {data_max_x:.1f}), y=({data_min_y:.1f}, {data_max_y:.1f})")
    
    # Validate each ROI
    valid_rois = {}
    for name, roi_info in roi_data.items():
        roi_bounds = roi_info['bounds']  # [minx, miny, maxx, maxy]
        
        # Check if ROI overlaps with data bounds
        roi_in_bounds = (
            roi_bounds[0] < data_max_x and roi_bounds[2] > data_min_x and  # x overlap
            roi_bounds[1] < data_max_y and roi_bounds[3] > data_min_y      # y overlap
        )
        
        if roi_in_bounds:
            valid_rois[name] = roi_info
            print(f"  ✓ {name}: within data bounds")
        else:
            print(f"  ⚠ {name}: outside data bounds - ROI bounds {roi_bounds}")
            print(f"    Consider checking coordinate systems or ROI placement")
    
    if len(valid_rois) != len(roi_data):
        print(f"Warning: {len(roi_data) - len(valid_rois)} ROIs are outside spatial data bounds")
    
    return valid_rois

def preview_rois(roi_data: Dict):
    """Print a preview of the ROI data"""
    print("\n=== ROI Preview ===")
    for name, roi_info in roi_data.items():
        coords = roi_info['coordinates']
        area = roi_info['area']
        bounds = roi_info['bounds']
        
        print(f"{name}:")
        print(f"  Points: {len(coords)}")
        print(f"  Area: {area:,.0f}")
        print(f"  Bounds: x=({bounds[0]:.1f}, {bounds[2]:.1f}), y=({bounds[1]:.1f}, {bounds[3]:.1f})")
        
        # Show first few coordinates
        if len(coords) <= 4:
            print(f"  Coordinates: {coords}")
        else:
            print(f"  Coordinates: {coords[:2]} ... {coords[-2:]} ({len(coords)} total)")
        print()

@click.command()
@click.argument('roi_json_path', type=click.Path(exists=True))
@click.argument('zarr_path', type=click.Path(exists=True))
@click.option('--output', '-o', default='converted_rois.json', 
              help='Output path for converted ROI JSON file')
@click.option('--validate/--no-validate', default=True,
              help='Validate ROIs against spatial data bounds')
@click.option('--preview/--no-preview', default=True,
              help='Show preview of converted ROIs')
@click.option('--run-analysis', is_flag=True,
              help='Automatically run ROI analysis after conversion')
def main(roi_json_path: str, zarr_path: str, output: str, validate: bool, 
         preview: bool, run_analysis: bool):
    """
    Convert napari-roi-manager ROIs to roi_umap_analysis format
    
    ROI_JSON_PATH: Path to the napari-roi-manager ROI JSON file
    ZARR_PATH: Path to the SpatialData zarr file
    
    Examples:
        # Basic conversion
        python convert_napari_roi_manager.py rois.json combined_direct_coords_annotated.zarr
        
        # Convert and run analysis
        python convert_napari_roi_manager.py rois.json combined_direct_coords_annotated.zarr --run-analysis
        
        # Custom output path
        python convert_napari_roi_manager.py rois.json combined_direct_coords_annotated.zarr -o my_converted_rois.json
    """
    
    print("=== napari-roi-manager to roi_umap_analysis Converter ===")
    print(f"Input ROI file: {roi_json_path}")
    print(f"Input zarr file: {zarr_path}")
    print(f"Output file: {output}")
    print()
    
    # Step 1: Load napari-roi-manager ROIs
    try:
        napari_roi_data = load_napari_roi_manager_json(roi_json_path)
    except Exception as e:
        print(f"Error loading ROI file: {e}")
        return
    
    # Step 2: Convert to standard format
    try:
        standard_rois = convert_to_standard_format(napari_roi_data)
        if not standard_rois:
            print("No ROIs were successfully converted")
            return
    except Exception as e:
        print(f"Error converting ROIs: {e}")
        return
    
    # Step 3: Validate against spatial data (optional)
    if validate:
        try:
            standard_rois = validate_rois_with_spatial_data(standard_rois, zarr_path)
            if not standard_rois:
                print("No ROIs passed validation")
                return
        except Exception as e:
            print(f"Warning: Could not validate ROIs against spatial data: {e}")
    
    # Step 4: Preview ROIs (optional)
    if preview:
        preview_rois(standard_rois)
    
    # Step 5: Save converted ROIs
    try:
        save_standard_roi_json(standard_rois, output)
    except Exception as e:
        print(f"Error saving converted ROIs: {e}")
        return
    
    # Step 6: Run analysis (optional)
    if run_analysis:
        print("\n=== Running ROI Analysis ===")
        import subprocess
        try:
            cmd = [
                'uv', 'run', 'python', 'roi_umap_analysis.py',
                zarr_path,
                '--roi-file', output,
                '--output-dir', 'napari_roi_manager_results'
            ]
            
            print(f"Running: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode == 0:
                print("✅ ROI analysis completed successfully!")
                print("Check the 'napari_roi_manager_results' directory for results")
            else:
                print("❌ ROI analysis failed:")
                print(result.stderr)
        except Exception as e:
            print(f"Error running analysis: {e}")
    
    print(f"\n✅ Conversion complete!")
    print(f"Your ROIs are now ready for use with roi_umap_analysis.py:")
    print(f"python roi_umap_analysis.py {zarr_path} --roi-file {output}")

if __name__ == "__main__":
    main()