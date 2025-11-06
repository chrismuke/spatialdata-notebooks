#!/usr/bin/env python3
"""
Verification script for combined zarr files
Demonstrates that the direct coordinate modification approach successfully fixes tissue overlap.
"""

import spatialdata as sd
import numpy as np

def verify_combined_zarr(zarr_path: str):
    """Verify that the combined zarr file has proper spatial separation."""
    print(f"Verifying combined zarr file: {zarr_path}")
    print("=" * 60)
    
    # Load the combined zarr file
    combined_sdata = sd.read_zarr(zarr_path)
    
    # Check contents
    print(f"Combined object contains:")
    print(f"  - Images: {len(combined_sdata.images)}")
    print(f"  - Labels: {len(combined_sdata.labels)}")
    print(f"  - Shapes: {len(combined_sdata.shapes)}")
    print(f"  - Points: {len(combined_sdata.points)}")
    print(f"  - Tables: {len(combined_sdata.tables)}")
    print()
    
    # Analyze transcript positions
    print("Transcript Analysis:")
    print("-" * 30)
    
    for point_name, point_data in combined_sdata.points.items():
        if 'transcripts' in point_name:
            n_transcripts = len(point_data)
            
            # Handle Dask arrays properly
            min_x_val = point_data['x'].min()
            max_x_val = point_data['x'].max()
            min_y_val = point_data['y'].min()
            max_y_val = point_data['y'].max()
            
            # Compute if dask arrays
            min_x = float(min_x_val.compute() if hasattr(min_x_val, 'compute') else min_x_val)
            max_x = float(max_x_val.compute() if hasattr(max_x_val, 'compute') else max_x_val)
            min_y = float(min_y_val.compute() if hasattr(min_y_val, 'compute') else min_y_val)
            max_y = float(max_y_val.compute() if hasattr(max_y_val, 'compute') else max_y_val)
            
            print(f"{point_name}:")
            print(f"  - Count: {n_transcripts:,} transcripts")
            print(f"  - X range: {min_x:.2f} to {max_x:.2f}")
            print(f"  - Y range: {min_y:.2f} to {max_y:.2f}")
            print()
    
    # Check for proper separation
    print("Separation Analysis:")
    print("-" * 30)
    
    transcript_names = [name for name in combined_sdata.points.keys() if 'transcripts' in name]
    if len(transcript_names) == 2:
        # Get Y ranges for both samples
        y_ranges = []
        for name in transcript_names:
            point_data = combined_sdata.points[name]
            min_y_val = point_data['y'].min()
            max_y_val = point_data['y'].max()
            
            # Compute if dask arrays
            min_y = float(min_y_val.compute() if hasattr(min_y_val, 'compute') else min_y_val)
            max_y = float(max_y_val.compute() if hasattr(max_y_val, 'compute') else max_y_val)
            y_ranges.append((min_y, max_y, name))
        
        # Sort by minimum Y coordinate
        y_ranges.sort(key=lambda x: x[0])
        
        # Check separation
        sample1_max_y = y_ranges[0][1]
        sample2_min_y = y_ranges[1][0]
        gap = sample2_min_y - sample1_max_y
        
        print(f"Sample 1 ({y_ranges[0][2]}): Y max = {sample1_max_y:.2f}")
        print(f"Sample 2 ({y_ranges[1][2]}): Y min = {sample2_min_y:.2f}")
        print(f"Gap between samples: {gap:.2f} pixels")
        
        if gap > 0:
            print("✅ SUCCESS: Tissues are properly separated!")
        else:
            print("❌ FAILURE: Tissues are overlapping!")
    
    print()
    print("Verification complete.")

if __name__ == "__main__":
    verify_combined_zarr("combined_direct_coords.zarr")