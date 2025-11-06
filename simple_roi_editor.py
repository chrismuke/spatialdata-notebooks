#!/usr/bin/env python3
"""
Simple ROI Editor for napari

A simplified script to create and save ROIs from napari without the complex analysis.
This focuses just on the napari interaction and JSON saving.
"""

import napari
import spatialdata as sd
import numpy as np
import json
from pathlib import Path

def simple_roi_editor(zarr_path: str, output_json: str = "my_rois.json"):
    """
    Simple ROI editor using napari
    
    Args:
        zarr_path: Path to spatial data zarr file
        output_json: Path to save ROI JSON file
    """
    
    print("=== Simple ROI Editor ===")
    print(f"Loading spatial data from: {zarr_path}")
    
    # Load spatial data
    sdata = sd.read_zarr(zarr_path)
    
    # Create napari viewer
    viewer = napari.Viewer(title="ROI Editor")
    
    # Add image data (sample one image to avoid clutter)
    if sdata.images:
        img_name, img_data = next(iter(sdata.images.items()))
        viewer.add_image(img_data, name=img_name, colormap='gray', contrast_limits=[0, 1000])
        print(f"Added image: {img_name}")
    
    # Add a sample of transcripts for context (reduce number for performance)
    if sdata.points:
        point_name, point_data = next(iter(sdata.points.items()))
        if 'x' in point_data.columns and 'y' in point_data.columns:
            # Sample 5000 random points for visualization
            n_sample = min(5000, len(point_data))
            indices = np.random.choice(len(point_data), n_sample, replace=False)
            coords = point_data.iloc[indices][['y', 'x']].values  # napari uses y,x order
            viewer.add_points(coords, name=f"{point_name}_sample", 
                            size=0.5, face_color='red', opacity=0.3)
            print(f"Added {n_sample} sample points from: {point_name}")
    
    # Add empty shapes layer for ROI drawing
    roi_layer = viewer.add_shapes(
        name='My_ROIs',
        shape_type='polygon',
        edge_width=3,
        edge_color='cyan',
        face_color=[0, 0, 1, 0.1],  # Semi-transparent blue
        text={
            'string': 'ROI_{index}',
            'size': 12,
            'color': 'cyan'
        }
    )
    
    print("\n=== Instructions ===")
    print("1. Select the 'My_ROIs' layer in napari")
    print("2. Make sure 'polygon' tool is selected")
    print("3. Click to draw polygon vertices around regions of interest")
    print("4. Double-click to close each polygon")
    print("5. Draw as many ROIs as needed")
    print("6. Close napari when finished")
    print("\nTips:")
    print("- Hold Shift while drawing for more precise control")
    print("- Use mouse wheel to zoom in/out")
    print("- Use middle mouse button to pan")
    
    # Run napari (this will block until the viewer is closed)
    napari.run()
    
    # After napari is closed, save the ROIs
    if len(roi_layer.data) > 0:
        print(f"\nSaving {len(roi_layer.data)} ROIs...")
        roi_data = {}
        
        for i, shape_coords in enumerate(roi_layer.data):
            # Convert from napari y,x format to x,y format
            coords_xy = shape_coords[:, [1, 0]]
            
            # Calculate area and bounds
            from shapely.geometry import Polygon
            polygon = Polygon(coords_xy)
            
            roi_name = f"ROI_{i+1}"
            roi_data[roi_name] = {
                'coordinates': coords_xy.tolist(),
                'area': float(polygon.area),
                'bounds': list(polygon.bounds),
                'created_with': 'simple_roi_editor'
            }
        
        # Save to JSON
        with open(output_json, 'w') as f:
            json.dump(roi_data, f, indent=2)
        
        print(f"✅ Successfully saved {len(roi_data)} ROIs to: {output_json}")
        
        # Print summary
        print("\n=== ROI Summary ===")
        for name, info in roi_data.items():
            print(f"{name}: Area = {info['area']:,.0f}, Bounds = {info['bounds']}")
        
        return roi_data
    else:
        print("❌ No ROIs were drawn. JSON file not created.")
        return {}

def load_and_view_rois(json_path: str):
    """Load and display existing ROIs"""
    print(f"Loading ROIs from: {json_path}")
    
    with open(json_path, 'r') as f:
        roi_data = json.load(f)
    
    print(f"Found {len(roi_data)} ROIs:")
    for name, info in roi_data.items():
        print(f"  {name}: {len(info['coordinates'])} points, Area = {info['area']:,.0f}")
    
    return roi_data

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 2:
        print("Usage:")
        print("  Create ROIs: python simple_roi_editor.py <zarr_path> [output.json]")
        print("  View ROIs:   python simple_roi_editor.py --view <roi.json>")
        print("")
        print("Examples:")
        print("  python simple_roi_editor.py combined_direct_coords_annotated.zarr my_rois.json")
        print("  python simple_roi_editor.py --view my_rois.json")
        sys.exit(1)
    
    if sys.argv[1] == "--view":
        if len(sys.argv) < 3:
            print("Error: Please provide JSON file path")
            sys.exit(1)
        load_and_view_rois(sys.argv[2])
    else:
        zarr_path = sys.argv[1]
        output_json = sys.argv[2] if len(sys.argv) > 2 else "my_rois.json"
        
        if not Path(zarr_path).exists():
            print(f"Error: Zarr file not found: {zarr_path}")
            sys.exit(1)
        
        simple_roi_editor(zarr_path, output_json)