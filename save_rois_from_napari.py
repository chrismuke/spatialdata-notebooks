#!/usr/bin/env python3
"""
Save ROI JSON from napari

This script shows different methods to save ROI data from napari:
1. Interactive napari session with automatic JSON saving
2. Manual extraction from napari shapes layer
3. Programmatic ROI creation and saving
"""

import napari
import spatialdata as sd
import numpy as np
import json
from pathlib import Path
from shapely.geometry import Polygon
from typing import Dict, List

def method1_interactive_napari_with_json_save(zarr_path: str, output_json: str):
    """
    Method 1: Interactive napari with automatic JSON saving
    """
    print("=== Method 1: Interactive napari with JSON saving ===")
    print("Instructions:")
    print("1. napari will open with your spatial data")
    print("2. Use the 'ROIs' shapes layer to draw polygons")
    print("3. Draw around areas of interest")
    print("4. Close napari when done")
    print("5. ROIs will be automatically saved to JSON")
    
    # Load spatial data
    sdata = sd.read_zarr(zarr_path)
    
    # Create napari viewer
    viewer = napari.Viewer()
    
    try:
        # Add spatial data layers
        for img_name, img_data in sdata.images.items():
            viewer.add_image(img_data, name=img_name)
        
        # Add points (transcripts) - sample subset for performance
        for point_name, point_data in sdata.points.items():
            if 'x' in point_data.columns and 'y' in point_data.columns:
                # Sample points for visualization (napari can be slow with too many points)
                n_points = min(10000, len(point_data))
                indices = np.random.choice(len(point_data), n_points, replace=False)
                coords = point_data.iloc[indices][['y', 'x']].values  # napari expects y,x
                viewer.add_points(coords, name=f"{point_name}_sample", size=0.5, opacity=0.3)
        
        # Add shapes (cell boundaries) - sample subset
        for shape_name, shape_data in sdata.shapes.items():
            if 'cell_boundaries' in shape_name:
                shapes_list = []
                # Sample shapes for performance
                n_shapes = min(1000, len(shape_data))
                for i, geom in enumerate(shape_data.geometry.iloc[:n_shapes]):
                    if geom.geom_type == 'Polygon':
                        coords = np.array(geom.exterior.coords)
                        coords = coords[:, [1, 0]]  # swap x,y to y,x for napari
                        shapes_list.append(coords)
                
                if shapes_list:
                    viewer.add_shapes(shapes_list, name=f"{shape_name}_sample", 
                                    shape_type='polygon', edge_width=0.3, 
                                    edge_color='white', face_color='transparent')
        
        # Add empty ROI layer for drawing
        roi_layer = viewer.add_shapes(name='ROIs', shape_type='polygon', 
                                    edge_width=3, edge_color='red', face_color='transparent')
        
        print(f"napari opened. Draw your ROIs and close when done.")
        
        # Show viewer (blocks until closed)
        napari.run()
        
        # Extract and save ROIs after napari is closed
        if len(roi_layer.data) > 0:
            roi_data = extract_roi_data_from_layer(roi_layer)
            save_roi_json(roi_data, output_json)
            print(f"✅ Saved {len(roi_data)} ROIs to {output_json}")
            return roi_data
        else:
            print("❌ No ROIs were drawn")
            return {}
            
    except Exception as e:
        print(f"Error: {e}")
        return {}

def method2_manual_extraction(shapes_layer, output_json: str):
    """
    Method 2: Manual extraction from existing napari shapes layer
    """
    print("=== Method 2: Manual extraction from shapes layer ===")
    
    roi_data = extract_roi_data_from_layer(shapes_layer)
    save_roi_json(roi_data, output_json)
    print(f"✅ Extracted and saved {len(roi_data)} ROIs to {output_json}")
    return roi_data

def method3_programmatic_creation(output_json: str):
    """
    Method 3: Programmatically create ROIs and save to JSON
    """
    print("=== Method 3: Programmatic ROI creation ===")
    
    # Create example ROIs programmatically
    roi_data = {
        "Upper_Left_Region": {
            "coordinates": [
                [1000, 0], [3000, 0], [3000, 1200], [1000, 1200]
            ],
            "area": 2400000,
            "bounds": [1000, 0, 3000, 1200],
            "description": "Upper left tissue region"
        },
        "Central_Region": {
            "coordinates": [
                [2500, 1000], [5500, 1000], [5500, 3000], [2500, 3000]
            ],
            "area": 6000000,
            "bounds": [2500, 1000, 5500, 3000],
            "description": "Central tissue region"
        },
        "Lower_Right_Region": {
            "coordinates": [
                [4000, 3500], [6500, 3500], [6500, 5000], [4000, 5000]
            ],
            "area": 3750000,
            "bounds": [4000, 3500, 6500, 5000],
            "description": "Lower right tissue region"
        }
    }
    
    save_roi_json(roi_data, output_json)
    print(f"✅ Created and saved {len(roi_data)} programmatic ROIs to {output_json}")
    return roi_data

def extract_roi_data_from_layer(roi_layer) -> Dict:
    """Extract ROI data from napari shapes layer"""
    roi_data = {}
    
    for i, shape_coords in enumerate(roi_layer.data):
        # Convert from napari y,x format to x,y format
        coords_xy = shape_coords[:, [1, 0]]
        
        # Create polygon
        polygon = Polygon(coords_xy)
        
        # Generate ROI name
        roi_name = f"ROI_{i+1}"
        if hasattr(roi_layer, 'text') and len(roi_layer.text) > i and roi_layer.text[i]:
            roi_name = roi_layer.text[i]
        
        roi_data[roi_name] = {
            'coordinates': coords_xy.tolist(),
            'area': polygon.area,
            'bounds': list(polygon.bounds),
            'perimeter': polygon.length,
            'created_method': 'napari_interactive'
        }
    
    return roi_data

def save_roi_json(roi_data: Dict, output_path: str):
    """Save ROI data to JSON file with pretty formatting"""
    roi_json = {}
    
    for name, roi_info in roi_data.items():
        # Ensure all data is JSON serializable
        roi_json[name] = {
            'coordinates': roi_info['coordinates'],
            'area': float(roi_info['area']),
            'bounds': [float(x) for x in roi_info['bounds']],
        }
        
        # Add optional fields if present
        if 'perimeter' in roi_info:
            roi_json[name]['perimeter'] = float(roi_info['perimeter'])
        if 'description' in roi_info:
            roi_json[name]['description'] = roi_info['description']
        if 'created_method' in roi_info:
            roi_json[name]['created_method'] = roi_info['created_method']
    
    # Save with pretty formatting
    with open(output_path, 'w') as f:
        json.dump(roi_json, f, indent=2, sort_keys=True)

def load_roi_json(json_path: str) -> Dict:
    """Load ROI data from JSON file"""
    with open(json_path, 'r') as f:
        return json.load(f)

def validate_roi_json(json_path: str) -> bool:
    """Validate ROI JSON file format"""
    try:
        roi_data = load_roi_json(json_path)
        
        for roi_name, roi_info in roi_data.items():
            # Check required fields
            required_fields = ['coordinates', 'area', 'bounds']
            for field in required_fields:
                if field not in roi_info:
                    print(f"❌ Missing required field '{field}' in ROI '{roi_name}'")
                    return False
            
            # Check coordinate format
            coords = roi_info['coordinates']
            if not isinstance(coords, list) or len(coords) < 3:
                print(f"❌ Invalid coordinates format in ROI '{roi_name}'")
                return False
            
            # Check bounds format
            bounds = roi_info['bounds']
            if not isinstance(bounds, list) or len(bounds) != 4:
                print(f"❌ Invalid bounds format in ROI '{roi_name}'")
                return False
        
        print(f"✅ ROI JSON file is valid with {len(roi_data)} ROIs")
        return True
        
    except Exception as e:
        print(f"❌ Error validating ROI JSON: {e}")
        return False

def demo_all_methods():
    """Demonstrate all methods for saving ROI JSON"""
    print("ROI JSON Saving Methods Demo")
    print("=" * 40)
    
    # Method 3: Programmatic (always works)
    method3_programmatic_creation("programmatic_rois.json")
    
    # Validate the created file
    validate_roi_json("programmatic_rois.json")
    
    # Show the content
    print("\n=== Example ROI JSON Content ===")
    with open("programmatic_rois.json", 'r') as f:
        content = f.read()
        print(content[:500] + "..." if len(content) > 500 else content)
    
    print("\n=== Usage Instructions ===")
    print("1. For interactive ROI definition:")
    print("   python save_rois_from_napari.py --interactive combined_direct_coords_annotated.zarr my_rois.json")
    print()
    print("2. For programmatic ROI creation:")
    print("   python save_rois_from_napari.py --programmatic my_rois.json")
    print()
    print("3. To use with ROI analysis:")
    print("   python roi_umap_analysis.py combined_direct_coords_annotated.zarr --roi-file my_rois.json")

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Save ROI JSON from napari")
    parser.add_argument("--interactive", nargs=2, metavar=("ZARR_PATH", "OUTPUT_JSON"),
                       help="Interactive napari ROI definition")
    parser.add_argument("--programmatic", metavar="OUTPUT_JSON", 
                       help="Create programmatic ROIs")
    parser.add_argument("--validate", metavar="JSON_PATH",
                       help="Validate existing ROI JSON file")
    parser.add_argument("--demo", action="store_true",
                       help="Run demonstration of all methods")
    
    args = parser.parse_args()
    
    if args.interactive:
        zarr_path, output_json = args.interactive
        method1_interactive_napari_with_json_save(zarr_path, output_json)
    elif args.programmatic:
        method3_programmatic_creation(args.programmatic)
    elif args.validate:
        validate_roi_json(args.validate)
    elif args.demo:
        demo_all_methods()
    else:
        demo_all_methods()