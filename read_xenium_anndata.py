#!/usr/bin/env python3
"""Read SpatialData object from Xenium zarr file with manual reconstruction."""

import zarr
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
import spatialdata as sd
import xarray as xr
from xarray.core.datatree import DataTree
from spatialdata.models import Image2DModel, Labels2DModel, PointsModel, TableModel

def read_anndata_from_zarr_group(table_group):
    """Read AnnData from a zarr table group."""
    # Read sparse matrix data
    X_group = table_group['X']
    data = X_group['data'][:]
    indices = X_group['indices'][:]
    indptr = X_group['indptr'][:]
    
    # Create sparse matrix
    X = sparse.csr_matrix((data, indices, indptr))
    
    # Read observations (cells)
    obs_group = table_group['obs']
    obs_data = {}
    for key in obs_group.keys():
        if key == '_index':
            obs_index = obs_group[key][:].astype(str)
        elif key == 'region':
            # Handle categorical data
            categories = obs_group[key]['categories'][:].astype(str)
            codes = obs_group[key]['codes'][:]
            obs_data[key] = pd.Categorical.from_codes(codes, categories=categories)
        else:
            obs_data[key] = obs_group[key][:]
    
    obs_df = pd.DataFrame(obs_data, index=obs_index)
    
    # Read variables (genes)
    var_group = table_group['var']
    var_data = {}
    for key in var_group.keys():
        if key == '_index':
            var_index = var_group[key][:].astype(str)
        else:
            var_data[key] = var_group[key][:].astype(str)
    
    # Create var DataFrame with gene symbols as index and metadata as columns
    var_df = pd.DataFrame(var_data, index=var_index)
    
    # Ensure the columns are in the expected order
    expected_columns = ['gene_ids', 'feature_types', 'genome']
    var_df = var_df.reindex(columns=[col for col in expected_columns if col in var_df.columns])
    
    # Read spatial coordinates
    spatial_coords = table_group['obsm']['spatial'][:]
    obsm = {'spatial': spatial_coords}
    
    # Create AnnData object
    adata = ad.AnnData(X=X, obs=obs_df, var=var_df, obsm=obsm)
    
    return adata

def read_spatialdata_from_zarr(zarr_path):
    """Manually reconstruct SpatialData object from zarr store."""
    root = zarr.open(zarr_path, mode='r')
    
    # Initialize containers for SpatialData elements
    images = {}
    labels = {}
    points = {}
    tables = {}
    
    # Read images
    if 'images' in root:
        for img_name in root['images'].keys():
            img_group = root['images'][img_name]
            
            # Get the highest resolution level (0)
            level_data = img_group['0'][:]
            dims = ['c', 'y', 'x'] if len(level_data.shape) == 3 else ['y', 'x']
            
            # Create proper coordinates for the level
            coords = {}
            for i, dim in enumerate(dims):
                coords[dim] = np.arange(level_data.shape[i])
            
            img_xr = xr.DataArray(
                level_data,
                dims=dims,
                coords=coords
            )
            
            # Calculate scale factors from the existing levels
            scale_factors = []
            prev_shape = level_data.shape
            for level in sorted(img_group.keys(), key=int)[1:]:  # Skip level 0
                curr_data = img_group[level]
                curr_shape = curr_data.shape
                # Calculate the scale factor for each dimension (use y dimension)
                scale_factor = prev_shape[-2] / curr_shape[-2]  # y dimension
                scale_factors.append(int(round(scale_factor)))
                prev_shape = curr_shape
            
            # Apply SpatialData model for multiscale images
            images[img_name] = Image2DModel.parse(
                img_xr, 
                dims=dims, 
                scale_factors=scale_factors if scale_factors else None
            )
    
    # Read labels
    if 'labels' in root:
        for label_name in root['labels'].keys():
            label_group = root['labels'][label_name]
            
            # Get the highest resolution level (0)
            level_data = label_group['0'][:]
            dims = ['y', 'x']
            
            # Create proper coordinates for the level
            coords = {}
            for i, dim in enumerate(dims):
                coords[dim] = np.arange(level_data.shape[i])
            
            label_xr = xr.DataArray(
                level_data,
                dims=dims,
                coords=coords
            )
            
            # Calculate scale factors from the existing levels
            scale_factors = []
            prev_shape = level_data.shape
            for level in sorted(label_group.keys(), key=int)[1:]:  # Skip level 0
                curr_data = label_group[level]
                curr_shape = curr_data.shape
                # Calculate the scale factor for each dimension (use y dimension)
                scale_factor = prev_shape[0] / curr_shape[0]  # y dimension
                scale_factors.append(int(round(scale_factor)))
                prev_shape = curr_shape
            
            # Apply SpatialData model for multiscale labels
            labels[label_name] = Labels2DModel.parse(
                label_xr, 
                dims=dims, 
                scale_factors=scale_factors if scale_factors else None
            )
    
    # Read points (transcripts)
    if 'points' in root:
        for point_name in root['points'].keys():
            point_group = root['points'][point_name]
            
            # Check if the group has any arrays
            if len(list(point_group.keys())) > 0:
                # Read point coordinates and attributes
                point_data = {}
                for attr in point_group.keys():
                    if isinstance(point_group[attr], zarr.Array):
                        point_data[attr] = point_group[attr][:]
                
                if point_data:  # Only create if we have data
                    # Create DataFrame
                    points_df = pd.DataFrame(point_data)
                    
                    # Apply SpatialData model - need to specify coordinate columns
                    # Assume x, y coordinates exist or find them
                    coord_cols = {}
                    if 'x' in points_df.columns and 'y' in points_df.columns:
                        coord_cols = {'x': 'x', 'y': 'y'}
                    elif 'X' in points_df.columns and 'Y' in points_df.columns:
                        coord_cols = {'x': 'X', 'y': 'Y'}
                    
                    if coord_cols:
                        points[point_name] = PointsModel.parse(points_df, coordinates=coord_cols)
                    else:
                        print(f"Warning: Could not find x,y coordinates for {point_name}, skipping")
            else:
                print(f"Warning: Points group '{point_name}' is empty, skipping")
    
    # Read tables
    if 'tables' in root:
        for table_name in root['tables'].keys():
            table_group = root['tables'][table_name]
            adata = read_anndata_from_zarr_group(table_group)
            
            # Set up proper spatial annotations for the table
            # Check if the table has spatial information in uns
            if 'spatialdata_attrs' in table_group.get('uns', {}):
                attrs_group = table_group['uns']['spatialdata_attrs']
                
                # Read spatial annotations - these are scalar values
                instance_key = str(attrs_group['instance_key'][()]) if 'instance_key' in attrs_group else 'cell_id'
                region = str(attrs_group['region'][()]) if 'region' in attrs_group else 'cell_labels'
                region_key = str(attrs_group['region_key'][()]) if 'region_key' in attrs_group else 'region'
                
                # Set the spatial annotations
                adata.uns['spatialdata_attrs'] = {
                    'instance_key': instance_key,
                    'region': region,
                    'region_key': region_key
                }
                
                # Ensure the region column exists in obs and points to the correct labels
                if region_key in adata.obs.columns:
                    # Update region values to match available labels
                    if 'cell_labels' in labels:
                        adata.obs[region_key] = pd.Categorical([region] * len(adata.obs), categories=[region])
            else:
                # Set default spatial annotations if not found
                instance_key = 'cell_id' if 'cell_id' in adata.obs.columns else None
                region = 'cell_labels' if 'cell_labels' in labels else None
                region_key = 'region'
                
                if region:
                    adata.uns['spatialdata_attrs'] = {
                        'instance_key': instance_key,
                        'region': region,
                        'region_key': region_key
                    }
                    
                    # Add region column if it doesn't exist
                    if region_key not in adata.obs.columns:
                        adata.obs[region_key] = pd.Categorical([region] * len(adata.obs), categories=[region])
            
            # Apply SpatialData model - it will use the spatialdata_attrs from uns
            tables[table_name] = TableModel.parse(adata)
    
    # Create SpatialData object
    sdata = sd.SpatialData(
        images=images,
        labels=labels,
        points=points,
        tables=tables
    )
    
    return sdata

def main():
    zarr_path = "/Users/chrism/datasets/10x_xenium_mouse_datasets/zarr/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP.zarr"
    
    print(f"Reading SpatialData from: {zarr_path}")
    sdata = read_spatialdata_from_zarr(zarr_path)
    
    print(f"SpatialData object, with associated Zarr store: {zarr_path}")
    print(sdata)
    
    return sdata

if __name__ == "__main__":
    sdata = main()