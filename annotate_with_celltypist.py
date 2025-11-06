#!/usr/bin/env python3
"""
Simple CellTypist annotation for spatial data
"""

import spatialdata as sd
import celltypist
from celltypist import models
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad

def annotate_spatial_data(zarr_path, output_path, model_name='Adult_Human_Skin.pkl'):
    """Annotate spatial data using CellTypist"""
    
    print(f"Loading spatial data from: {zarr_path}")
    sdata = sd.read_zarr(zarr_path)
    
    # Get the table
    table = sdata.tables['table']
    print(f"Data shape: {table.shape}")
    
    # Convert to AnnData for CellTypist
    adata = ad.AnnData(X=table.X, obs=table.obs, var=table.var)
    
    # Normalize the data
    print("Normalizing data...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # Load the model
    print(f"Loading model: {model_name}")
    model = models.Model.load(model=model_name)
    
    # Predict cell types
    print("Predicting cell types...")
    predictions = celltypist.annotate(adata, model=model, majority_voting=True)
    
    # Add predictions back to the table
    predicted_labels = predictions.predicted_labels
    
    # Save predictions
    print(f"Saving predictions to: {output_path}")
    predicted_labels.to_csv(output_path)
    
    # Check available columns
    print("Available columns in predictions:")
    print(predicted_labels.columns.tolist())
    
    # Also save the predicted labels back to the spatial data
    spatial_output = zarr_path.replace('.zarr', '_annotated.zarr')
    table.obs['cell_type_predicted'] = predicted_labels['predicted_labels'].values
    
    # Add confidence scores if available
    if 'conf_score' in predicted_labels.columns:
        table.obs['cell_type_conf_score'] = predicted_labels['conf_score'].values
    elif 'majority_voting' in predicted_labels.columns:
        table.obs['cell_type_conf_score'] = predicted_labels['majority_voting'].values
    
    # Save the annotated spatial data
    print(f"Saving annotated spatial data to: {spatial_output}")
    sdata.write(spatial_output)
    
    print("Annotation complete!")
    return predicted_labels

if __name__ == "__main__":
    # Annotate the combined data
    predictions = annotate_spatial_data(
        'combined_direct_coords.zarr',
        'combined_celltypist_predictions.csv',
        'Adult_Human_Skin.pkl'
    )
    
    print("Cell type distribution:")
    print(predictions['predicted_labels'].value_counts())