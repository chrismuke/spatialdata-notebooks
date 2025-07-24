#!/usr/bin/env python3
"""Minimal script to read Xenium zarr file with spatialdata."""

import spatialdata as sd

def main():
    zarr_path = "/Users/chrism/datasets/10x_xenium_mouse_datasets/zarr/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP.zarr"
    
    # Read the SpatialData object
    sdata = sd.read_zarr(zarr_path)
    
    # Display basic information
    print(f"Loaded SpatialData object from: {zarr_path}")
    print(f"Elements: {list(sdata.elements.keys())}")
    print(f"Coordinate systems: {list(sdata.coordinate_systems.keys())}")
    
    return sdata

if __name__ == "__main__":
    sdata = main()