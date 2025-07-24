#!/usr/bin/env python3
"""Inspect Xenium zarr file structure."""

import zarr
import os

def inspect_zarr_structure(zarr_path):
    """Inspect the structure of a zarr store."""
    if not os.path.exists(zarr_path):
        print(f"Error: Path {zarr_path} does not exist")
        return
    
    root = zarr.open(zarr_path, mode='r')
    
    print(f"Zarr store: {zarr_path}")
    print(f"Root keys: {list(root.keys())}")
    
    def print_group_structure(group, prefix=""):
        for key in group.keys():
            item = group[key]
            if isinstance(item, zarr.Group):
                print(f"{prefix}{key}/ (Group)")
                print_group_structure(item, prefix + "  ")
            else:
                print(f"{prefix}{key} (Array) - shape: {item.shape}, dtype: {item.dtype}")
    
    print_group_structure(root)

def main():
    zarr_path = "/Users/chrism/datasets/10x_xenium_mouse_datasets/zarr/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP.zarr"
    inspect_zarr_structure(zarr_path)

if __name__ == "__main__":
    main()