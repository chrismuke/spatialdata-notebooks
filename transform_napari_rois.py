#!/usr/bin/env python3
"""
Transform napari-roi-manager ROIs to spatial data coordinate system

The large coordinates in your napari-roi-manager ROIs indicate they were drawn
at full resolution (scale0), while the spatial analysis uses the global coordinate
system. This script applies the appropriate coordinate transformation.
"""

import json
import numpy as np
from shapely.geometry import Polygon
import click

def calculate_scale_factor(roi_data, spatial_bounds):
    """Calculate the scale factor between ROI and spatial coordinates"""
    
    # Get ROI coordinate ranges
    all_x_coords = []
    all_y_coords = []
    
    for coords in roi_data['data']:
        for x, y in coords:
            all_x_coords.append(x)
            all_y_coords.append(y)
    
    roi_min_x, roi_max_x = min(all_x_coords), max(all_x_coords)
    roi_min_y, roi_max_y = min(all_y_coords), max(all_y_coords)
    
    # Spatial data bounds
    spatial_min_x, spatial_min_y, spatial_max_x, spatial_max_y = spatial_bounds
    
    # Calculate scale factors
    roi_width = roi_max_x - roi_min_x
    roi_height = roi_max_y - roi_min_y
    spatial_width = spatial_max_x - spatial_min_x
    spatial_height = spatial_max_y - spatial_min_y
    
    scale_x = roi_width / spatial_width if spatial_width > 0 else 1
    scale_y = roi_height / spatial_height if spatial_height > 0 else 1
    
    # Use average scale factor
    scale_factor = (scale_x + scale_y) / 2
    
    print(f"ROI bounds: x=({roi_min_x:.1f}, {roi_max_x:.1f}), y=({roi_min_y:.1f}, {roi_max_y:.1f})")
    print(f"Spatial bounds: x=({spatial_min_x:.1f}, {spatial_max_x:.1f}), y=({spatial_min_y:.1f}, {spatial_max_y:.1f})")
    print(f"Scale factors: x={scale_x:.1f}, y={scale_y:.1f}")
    print(f"Average scale factor: {scale_factor:.1f}")
    
    return scale_factor

def transform_coordinates(roi_data, scale_factor, offset_x=0, offset_y=0):
    """Transform ROI coordinates by scaling and offsetting"""
    
    transformed_rois = {}
    
    for i, (coords, name) in enumerate(zip(roi_data['data'], roi_data['names'])):
        # Scale down coordinates
        scaled_coords = [[x / scale_factor + offset_x, y / scale_factor + offset_y] 
                        for x, y in coords]
        
        # Create polygon and calculate properties
        polygon = Polygon(scaled_coords)
        
        if polygon.is_valid:
            transformed_rois[name] = {
                'coordinates': scaled_coords,
                'area': float(polygon.area),
                'bounds': list(polygon.bounds),
                'original_scale_factor': scale_factor,
                'created_with': 'napari-roi-manager',
                'transformed': True
            }
            print(f"✓ Transformed {name}: {len(scaled_coords)} points")
        else:
            print(f"✗ Invalid polygon after transformation: {name}")
    
    return transformed_rois

@click.command()
@click.argument('roi_json_path', type=click.Path(exists=True))
@click.option('--scale-factor', type=float, help='Manual scale factor (auto-detect if not provided)')
@click.option('--spatial-bounds', type=str, default='1738.5,0.0,6413.1,5026.7',
              help='Spatial data bounds as "min_x,min_y,max_x,max_y"')
@click.option('--output', '-o', default='transformed_rois.json',
              help='Output path for transformed ROI JSON file')
@click.option('--offset-x', type=float, default=0.0,
              help='X offset to apply after scaling')
@click.option('--offset-y', type=float, default=0.0,
              help='Y offset to apply after scaling')
def main(roi_json_path, scale_factor, spatial_bounds, output, offset_x, offset_y):
    """
    Transform napari-roi-manager ROI coordinates to spatial data coordinate system
    
    This fixes the coordinate system mismatch where ROIs were drawn at full resolution
    but analysis needs them in the global coordinate system.
    
    ROI_JSON_PATH: Path to the napari-roi-manager ROI JSON file
    
    Examples:
        # Auto-detect scale factor
        python transform_napari_rois.py rois.json
        
        # Manual scale factor
        python transform_napari_rois.py rois.json --scale-factor 18.0
        
        # Custom spatial bounds
        python transform_napari_rois.py rois.json --spatial-bounds "0,0,10000,10000"
    """
    
    print("=== napari-roi-manager Coordinate Transformation ===")
    
    # Load napari-roi-manager ROIs
    with open(roi_json_path, 'r') as f:
        roi_data = json.load(f)
    
    print(f"Loaded {len(roi_data['data'])} ROIs from {roi_json_path}")
    
    # Parse spatial bounds
    bounds_parts = spatial_bounds.split(',')
    if len(bounds_parts) != 4:
        print("Error: spatial-bounds must be in format 'min_x,min_y,max_x,max_y'")
        return
    
    spatial_bounds_vals = [float(x) for x in bounds_parts]
    
    # Calculate or use provided scale factor
    if scale_factor is None:
        scale_factor = calculate_scale_factor(roi_data, spatial_bounds_vals)
        print(f"Auto-detected scale factor: {scale_factor:.1f}")
    else:
        print(f"Using manual scale factor: {scale_factor:.1f}")
    
    # Transform coordinates
    print(f"Transforming coordinates with scale={scale_factor:.1f}, offset=({offset_x:.1f}, {offset_y:.1f})")
    transformed_rois = transform_coordinates(roi_data, scale_factor, offset_x, offset_y)
    
    if not transformed_rois:
        print("No ROIs were successfully transformed")
        return
    
    # Save transformed ROIs in standard format
    clean_rois = {}
    for name, roi_info in transformed_rois.items():
        clean_rois[name] = {
            'coordinates': roi_info['coordinates'],
            'area': roi_info['area'],
            'bounds': roi_info['bounds']
        }
    
    with open(output, 'w') as f:
        json.dump(clean_rois, f, indent=2)
    
    print(f"✅ Saved {len(clean_rois)} transformed ROIs to {output}")
    
    # Show results preview
    print("\n=== Transformed ROI Preview ===")
    for name, roi_info in clean_rois.items():
        bounds = roi_info['bounds']
        print(f"{name}: bounds = ({bounds[0]:.1f}, {bounds[1]:.1f}, {bounds[2]:.1f}, {bounds[3]:.1f})")
    
    print(f"\n✅ ROIs are now ready for analysis:")
    print(f"python roi_umap_analysis.py combined_direct_coords_annotated.zarr --roi-file {output}")

if __name__ == "__main__":
    main()