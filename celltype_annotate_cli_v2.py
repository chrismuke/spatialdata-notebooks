#!/usr/bin/env python3
"""
Cell Type Annotation CLI Tool - Version 2

This tool performs cell type annotation on Xenium spatial transcriptomics data
using a single-cell RNA-seq reference dataset and CellTypist.

Creates organized output directories with comprehensive reporting.

Based on the celltype_annotation_celltypist_mouse.ipynb notebook.
"""

import logging
import sys
from pathlib import Path
from datetime import datetime
import json
import base64

import anndata
import celltypist
import click
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import spatialdata as sd
import spatialdata_plot  # noqa
import squidpy as sq
from copy import deepcopy
import openpyxl  # For Excel export


def create_output_directory(xenium_path: Path, reference_path: Path, min_clusters: int = None, max_clusters: int = None, base_results_dir: str = "results") -> Path:
    """
    Create a structured output directory based on input files and parameters.
    
    Creates directories with format: xenium_filename___reference_filename___[clusters]___timestamp
    
    Args:
        xenium_path: Path to Xenium data
        reference_path: Path to reference data
        min_clusters: Minimum clusters parameter
        max_clusters: Maximum clusters parameter
        base_results_dir: Base directory for results (default: "results")
        
    Returns:
        Path to created output directory
        
    Examples:
        xenium.zarr + ref.h5ad -> results/xenium___ref___20250724_120110/
        With clusters -> results/xenium___ref___clusters_min5_max10___20250724_120110/
    """
    # Use exact filenames without extensions
    xenium_name = xenium_path.stem  # e.g. "combined_direct_coords_annotated" from "combined_direct_coords_annotated.zarr"
    ref_name = reference_path.stem  # e.g. "reference_data" from "reference_data.h5ad"
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Create folder name: xenium___reference___[clusters]___timestamp
    folder_parts = [xenium_name, ref_name]
    
    if min_clusters is not None or max_clusters is not None:
        cluster_str = "clusters"
        if min_clusters is not None:
            cluster_str += f"_min{min_clusters}"
        if max_clusters is not None:
            cluster_str += f"_max{max_clusters}"
        folder_parts.append(cluster_str)
    
    folder_parts.append(timestamp)
    
    output_dir = Path(base_results_dir) / "___".join(folder_parts)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create subdirectories
    (output_dir / "plots").mkdir(exist_ok=True)
    (output_dir / "data").mkdir(exist_ok=True)
    (output_dir / "logs").mkdir(exist_ok=True)
    
    return output_dir


def setup_logging(output_dir: Path, level: str = "INFO") -> None:
    """
    Set up logging configuration with file output.
    
    Args:
        output_dir: Output directory for log files
        level: Logging level
    """
    log_file = output_dir / "logs" / "celltype_annotation.log"
    
    # Create formatters
    file_formatter = logging.Formatter(
        "%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )
    console_formatter = logging.Formatter("%(message)s")
    
    # Setup file handler
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(getattr(logging, level.upper()))
    file_handler.setFormatter(file_formatter)
    
    # Setup console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(getattr(logging, level.upper()))
    console_handler.setFormatter(console_formatter)
    
    # Configure root logger
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        handlers=[file_handler, console_handler],
        force=True
    )


def harmonize_gene_names(adata_xenium: anndata.AnnData, adata_ref: anndata.AnnData) -> tuple[anndata.AnnData, anndata.AnnData]:
    """
    Harmonize gene names between Xenium and reference datasets.
    
    Args:
        adata_xenium: Xenium spatial data
        adata_ref: Single-cell reference data
        
    Returns:
        Tuple of harmonized datasets with matching gene names
    """
    logging.info("Harmonizing gene names between datasets...")
    
    # Create standardized gene IDs for Xenium data
    if 'gene_ids' in adata_xenium.var.columns:
        adata_xenium.var['gene_ids_stripped'] = adata_xenium.var['gene_ids'].str.split('.').str[0]
    else:
        logging.warning("No 'gene_ids' column found in Xenium data, using var_names")
        adata_xenium.var['gene_ids_stripped'] = adata_xenium.var_names.str.split('.').str[0]
    
    # Determine reference gene IDs
    if 'feature_id' in adata_ref.var.columns:
        ref_ids_series = adata_ref.var['feature_id']
    else:
        ref_ids_series = pd.Series(adata_ref.var.index)
    
    ref_ids_stripped = ref_ids_series.str.split('.').str[0]
    
    # Find common genes
    intersecting_genes = set(adata_xenium.var['gene_ids_stripped']).intersection(set(ref_ids_stripped))
    logging.info(f"Found {len(intersecting_genes)} common genes after standardizing IDs")
    
    if len(intersecting_genes) < 50:
        logging.warning("⚠️ Very few common genes found. Annotation quality may be poor.")
    
    # Filter reference data
    adata_ref_mask = ref_ids_stripped.isin(intersecting_genes)
    adata_ref = adata_ref[:, adata_ref_mask].copy()
    
    # Create mapping from gene ID to symbol
    gene_id_to_name_map = pd.Series(
        adata_xenium.var.index.values, 
        index=adata_xenium.var['gene_ids_stripped']
    ).drop_duplicates()
    
    # Update reference var index
    if 'feature_id' in adata_ref.var.columns:
        filtered_ref_ids_stripped = adata_ref.var['feature_id'].str.split('.').str[0]
    else:
        filtered_ref_ids_stripped = pd.Series(adata_ref.var.index).str.split('.').str[0]
    
    adata_ref.var.index = filtered_ref_ids_stripped.map(gene_id_to_name_map)
    adata_ref.var.index.name = 'gene_symbol'
    
    # Clean up temporary column
    del adata_xenium.var['gene_ids_stripped']
    
    # Final harmonization - keep only common genes
    common_genes = adata_ref.var_names.intersection(adata_xenium.var_names)
    logging.info(f"Final common genes: {len(common_genes)}")
    
    adata_ref = adata_ref[:, common_genes].copy()
    adata_xenium = adata_xenium[:, common_genes].copy()
    
    return adata_xenium, adata_ref


def normalize_data(adata: anndata.AnnData, target_sum: float = 1e4) -> None:
    """Normalize and log-transform expression data."""
    logging.info(f"Normalizing data to {target_sum} counts per cell")
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)


def create_model_save_path(reference_path: Path, label_column: str = "cell_type") -> Path:
    """
    Create a unique model save path in the same directory as the reference data.

    Args:
        reference_path: Path to reference h5ad file
        label_column: Column name used for training labels

    Returns:
        Path for saving the model file with unique name

    Examples:
        reference_data.h5ad -> reference_data_cell_type_model_20250724_120110.pkl
    """
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    ref_name = reference_path.stem  # e.g. "reference_data" from "reference_data.h5ad"
    model_filename = f"{ref_name}_{label_column}_model_{timestamp}.pkl"
    return reference_path.parent / model_filename


def save_celltypist_model(model: celltypist.models.Model, save_path: Path) -> None:
    """
    Save a trained CellTypist model to disk.

    Args:
        model: Trained CellTypist model
        save_path: Path where to save the model
    """
    logging.info(f"Saving CellTypist model to: {save_path}")
    model.write(str(save_path))
    logging.info("Model saved successfully")


def load_celltypist_model(model_path: Path) -> celltypist.models.Model:
    """
    Load a CellTypist model from disk.

    Args:
        model_path: Path to saved model file

    Returns:
        Loaded CellTypist model
    """
    logging.info(f"Loading CellTypist model from: {model_path}")
    if not model_path.exists():
        raise FileNotFoundError(f"Model file not found: {model_path}")

    model = celltypist.models.Model.load(str(model_path))
    logging.info("Model loaded successfully")
    return model


def train_celltypist_model(adata_ref: anndata.AnnData, label_column: str, feature_selection: bool = False,
                          save_model_path: Path = None) -> celltypist.models.Model:
    """
    Train CellTypist model on reference data.

    Args:
        adata_ref: Reference dataset for training
        label_column: Column name containing cell type labels
        feature_selection: Whether to use feature selection
        save_model_path: Path to save the trained model (optional)

    Returns:
        Trained CellTypist model
    """
    logging.info(f"Training CellTypist model using label column: {label_column}")

    if label_column not in adata_ref.obs.columns:
        raise ValueError(f"Label column '{label_column}' not found in reference data")

    model = celltypist.train(adata_ref, labels=label_column, feature_selection=feature_selection)
    logging.info("Model training completed")

    # Save model if path provided
    if save_model_path is not None:
        save_celltypist_model(model, save_model_path)

    return model


def predict_cell_types(adata_query: anndata.AnnData, model: celltypist.models.Model) -> pd.DataFrame:
    """Predict cell types using trained model."""
    logging.info("Predicting cell types on spatial data")
    predictions = celltypist.annotate(adata_query, model=model)
    return predictions.predicted_labels


def optimize_leiden_resolution(adata: anndata.AnnData, min_clusters: int = None, max_clusters: int = None) -> float:
    """
    Optimize Leiden clustering resolution to achieve desired cluster count.
    
    Args:
        adata: AnnData object with computed neighbors
        min_clusters: Minimum number of clusters desired
        max_clusters: Maximum number of clusters desired
        
    Returns:
        Optimized resolution value
    """
    if min_clusters is None and max_clusters is None:
        return 0.5
    
    resolutions = [0.2, 0.4, 0.6, 0.8, 1.0, 1.5]
    best_resolution = 0.5
    best_score = float('inf')
    
    logging.info(f"Optimizing clustering resolution for {min_clusters}-{max_clusters} clusters...")
    
    for resolution in resolutions:
        sc.tl.leiden(adata, resolution=resolution, key_added='leiden_temp')
        n_clusters = len(adata.obs['leiden_temp'].unique())
        
        # Calculate score based on distance from target range
        if min_clusters is not None and max_clusters is not None:
            if min_clusters <= n_clusters <= max_clusters:
                score = 0  # Perfect fit
                best_resolution = resolution
                best_score = score
                del adata.obs['leiden_temp']
                break
            elif n_clusters < min_clusters:
                score = min_clusters - n_clusters
            else:  # n_clusters > max_clusters
                score = n_clusters - max_clusters
        elif min_clusters is not None:
            if n_clusters >= min_clusters:
                score = 0
                best_resolution = resolution
                best_score = score
                del adata.obs['leiden_temp']
                break
            else:
                score = min_clusters - n_clusters
        elif max_clusters is not None:
            if n_clusters <= max_clusters:
                score = 0
                best_resolution = resolution
                best_score = score
                del adata.obs['leiden_temp']
                break
            else:
                score = n_clusters - max_clusters
        
        if score < best_score:
            best_score = score
            best_resolution = resolution
        
        del adata.obs['leiden_temp']
    
    logging.info(f"Optimized resolution: {best_resolution} (targeting {min_clusters}-{max_clusters} clusters)")
    return best_resolution


def perform_clustering_analysis(adata: anndata.AnnData, min_clusters: int = None, max_clusters: int = None) -> anndata.AnnData:
    """
    Perform clustering analysis on the data.
    
    Args:
        adata: AnnData object to analyze
        min_clusters: Minimum number of clusters desired
        max_clusters: Maximum number of clusters desired
        
    Returns:
        Processed AnnData object with clustering results
    """
    # Make a copy to avoid modifying the original
    adata_analysis = adata.copy()
    
    # Filter cells and genes
    sc.pp.filter_cells(adata_analysis, min_counts=10)
    sc.pp.filter_genes(adata_analysis, min_cells=5)
    
    # Calculate spatial neighbors
    sq.gr.spatial_neighbors(adata_analysis)
    
    # Store raw counts
    adata_analysis.layers["counts"] = adata_analysis.X.copy()
    
    # Normalize and log transform
    sc.pp.normalize_total(adata_analysis, inplace=True)
    sc.pp.log1p(adata_analysis)
    
    # PCA and neighbors
    sc.pp.pca(adata_analysis)
    sc.pp.neighbors(adata_analysis)
    
    # UMAP
    sc.tl.umap(adata_analysis)
    
    # Optimize Leiden clustering resolution
    optimal_resolution = optimize_leiden_resolution(adata_analysis, min_clusters, max_clusters)
    sc.tl.leiden(adata_analysis, resolution=optimal_resolution)
    
    n_clusters = len(adata_analysis.obs['leiden'].unique())
    logging.info(f"Final clustering: {n_clusters} clusters with resolution {optimal_resolution}")
    
    return adata_analysis


def generate_visualizations(
    sdata: sd.SpatialData,
    xenium_path: Path,
    reference_path: Path,
    output_dir: Path,
    table_name: str = "table",
    prediction_column: str = "cell_type_predicted",
    min_clusters: int = None,
    max_clusters: int = None,
    show_unknown_cells: bool = False
) -> None:
    """
    Generate and save visualization plots.
    
    Args:
        sdata: Annotated SpatialData object
        xenium_path: Path to original Xenium data (for filename generation)
        reference_path: Path to reference data (for filename generation)
        output_dir: Base output directory
        table_name: Name of table in SpatialData object
        prediction_column: Name of prediction column
        min_clusters: Minimum clusters parameter
        max_clusters: Maximum clusters parameter
    """
    plots_dir = output_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate descriptive filename prefix
    xenium_name = xenium_path.stem
    ref_name = reference_path.parent.name
    prefix = f"{xenium_name}_annotated_with_{ref_name}"
    
    logging.info(f"Generating visualizations with prefix: {prefix}")
    
    # Get the main adata object for analysis
    adata_main = sdata.tables[table_name].copy()
    
    # Perform clustering analysis to get properly prepared adata for spatial plotting
    adata_clustered = None
    try:
        adata_clustered = perform_clustering_analysis(adata_main, min_clusters, max_clusters)
    except Exception as e:
        logging.warning(f"Failed to perform clustering analysis: {e}")
    
    # 1. Spatial cell type map (create this using adata_clustered which has proper spatial setup)
    if adata_clustered is not None:
        try:
            logging.info("Creating spatial cell type map...")
            
            # Prepare data for plotting - filter unknown cells if requested
            adata_for_plot = adata_clustered.copy()
            
            if not show_unknown_cells:
                # Define common patterns for unknown/unassigned cell types
                unknown_patterns = ['unknown', 'unassigned', 'unlabeled', 'ambiguous', 'unclear', 'na', 'none', 'nan']
                
                # Filter out cells with unknown cell types
                cell_types = adata_for_plot.obs[prediction_column].astype(str).str.lower()
                is_known = ~cell_types.isin(unknown_patterns)
                
                # Also filter out actual NaN values and empty strings
                is_known = is_known & (adata_for_plot.obs[prediction_column].notna()) & (adata_for_plot.obs[prediction_column] != '')
                
                original_count = len(adata_for_plot)
                adata_for_plot = adata_for_plot[is_known]
                filtered_count = len(adata_for_plot)
                
                logging.info(f"Filtered spatial plot: showing {filtered_count:,} cells with known types (hiding {original_count - filtered_count:,} unknown cells)")
                title_suffix = " (Known Cell Types Only)"
            else:
                logging.info(f"Showing all {len(adata_for_plot):,} cells including unknown types")
                title_suffix = " (All Cells)"
            
            if len(adata_for_plot) > 0:
                # Use the clustered adata which has proper spatial neighbors calculated
                # This ensures the spatial plot works correctly (same approach as leiden clusters)
                fig, ax = plt.subplots(figsize=(12, 12))
                sq.pl.spatial_scatter(
                    adata_for_plot,
                    library_id="spatial", 
                    shape=None,
                    color=[prediction_column],
                    ax=ax,
                    frameon=False
                )
                ax.set_title(f'{xenium_name}: Predicted Cell Types{title_suffix}')
                
                plt.tight_layout()
                spatial_map_path = plots_dir / f"{prefix}_spatial_celltype_map.png"
                plt.savefig(spatial_map_path, dpi=300, bbox_inches='tight')
                plt.close()
                logging.info(f"Saved spatial cell type map: {spatial_map_path}")
            else:
                logging.warning("No cells to plot after filtering - all cells have unknown types")
            
        except Exception as e:
            logging.warning(f"Failed to create spatial cell type map: {e}")
            logging.warning(f"Error details: {str(e)}")
    else:
        logging.warning("Cannot create spatial cell type map: clustering analysis failed")
    
    # 2. Cell type distribution bar plot
    try:
        logging.info("Creating cell type distribution plot...")
        cell_type_counts = sdata.tables[table_name].obs[prediction_column].value_counts()
        
        fig, ax = plt.subplots(figsize=(12, 8))
        bars = ax.bar(range(len(cell_type_counts)), cell_type_counts.values)
        ax.set_xlabel('Cell Type')
        ax.set_ylabel('Number of Cells')
        ax.set_title(f'{xenium_name}: Cell Type Distribution')
        ax.set_xticks(range(len(cell_type_counts)))
        ax.set_xticklabels(cell_type_counts.index, rotation=45, ha='right')
        
        # Add count labels on bars
        for bar, count in zip(bars, cell_type_counts.values):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01*max(cell_type_counts.values),
                   str(count), ha='center', va='bottom', fontsize=8)
        
        plt.tight_layout()
        distribution_path = plots_dir / f"{prefix}_celltype_distribution.png"
        plt.savefig(distribution_path, dpi=300, bbox_inches='tight')
        plt.close()
        logging.info(f"Saved cell type distribution: {distribution_path}")
        
    except Exception as e:
        logging.warning(f"Failed to create cell type distribution plot: {e}")
    
    # 3. Create UMAP plots using existing clustered data
    if adata_clustered is not None:
        try:
            # Create UMAP plots
            logging.info("Creating UMAP plots...")
            fig, axes = plt.subplots(1, 2, figsize=(16, 6))
            
            # UMAP colored by total counts
            sc.pl.umap(adata_clustered, color="total_counts", ax=axes[0], show=False, frameon=False)
            axes[0].set_title('UMAP: Total Counts')
            
            # UMAP colored by leiden clusters
            sc.pl.umap(adata_clustered, color="leiden", ax=axes[1], show=False, frameon=False)
            axes[1].set_title('UMAP: Leiden Clusters')
            
            plt.tight_layout()
            umap_path = plots_dir / f"{prefix}_umap_analysis.png"
            plt.savefig(umap_path, dpi=300, bbox_inches='tight')
            plt.close()
            logging.info(f"Saved UMAP analysis: {umap_path}")
            
        except Exception as e:
            logging.warning(f"Failed to create UMAP plots: {e}")
    
    # 4. Spatial scatter plot with leiden clusters
    if adata_clustered is not None:
        try:
            logging.info("Creating spatial scatter plot with clusters...")
            
            # Remove leiden_colors if it exists to avoid conflicts
            if 'leiden_colors' in adata_clustered.uns:
                adata_clustered.uns.pop('leiden_colors')
            
            fig, ax = plt.subplots(figsize=(12, 12))
            sq.pl.spatial_scatter(
                adata_clustered,
                library_id="spatial",
                shape=None,
                color=["leiden"],
                ax=ax,
                frameon=False
            )
            ax.set_title('Spatial Distribution: Leiden Clusters')
            
            plt.tight_layout()
            spatial_clusters_path = plots_dir / f"{prefix}_spatial_leiden_clusters.png"
            plt.savefig(spatial_clusters_path, dpi=300, bbox_inches='tight')
            plt.close()
            logging.info(f"Saved spatial clusters plot: {spatial_clusters_path}")
            
        except Exception as e:
            logging.warning(f"Failed to create spatial scatter plot: {e}")
    
    # 5. Marker genes analysis
    if adata_clustered is not None:
        try:
            logging.info("Performing marker gene analysis...")
            
            # Calculate marker genes
            sc.tl.rank_genes_groups(adata_clustered, 'leiden', method='wilcoxon')
            
            # Create marker genes plot
            logging.info("Creating marker genes plot...")
            fig, ax = plt.subplots(figsize=(15, 10))
            sc.pl.rank_genes_groups(adata_clustered, n_genes=10, sharey=False, ax=ax, show=False)
            ax.set_title('Top Marker Genes per Leiden Cluster')
            
            plt.tight_layout()
            marker_genes_path = plots_dir / f"{prefix}_marker_genes.png"
            plt.savefig(marker_genes_path, dpi=300, bbox_inches='tight')
            plt.close()
            logging.info(f"Saved marker genes plot: {marker_genes_path}")
            
        except Exception as e:
            logging.warning(f"Failed to create marker genes plot: {e}")


def generate_html_report(
    output_dir: Path,
    xenium_path: Path,
    reference_path: Path,
    sdata: sd.SpatialData,
    table_name: str = "table",
    prediction_column: str = "cell_type_predicted",
    run_parameters: dict = None,
    zarr_saved: bool = True
) -> None:
    """
    Generate a comprehensive HTML report with all results.
    
    Args:
        output_dir: Output directory
        xenium_path: Path to Xenium data
        reference_path: Path to reference data
        sdata: Annotated SpatialData object
        table_name: Name of table in SpatialData object
        prediction_column: Name of prediction column
        run_parameters: Dictionary of run parameters
        zarr_saved: Whether the zarr file was saved
    """
    plots_dir = output_dir / "plots"
    report_path = output_dir / "annotation_report.html"
    
    # Get cell type statistics
    adata = sdata.tables[table_name]
    cell_type_counts = adata.obs[prediction_column].value_counts()
    
    # Helper function to encode images as base64
    def image_to_base64(image_path: Path) -> str:
        if image_path.exists():
            with open(image_path, "rb") as img_file:
                return base64.b64encode(img_file.read()).decode('utf-8')
        return ""
    
    # Get image data
    prefix = f"{xenium_path.stem}_annotated_with_{reference_path.parent.name}"
    images = {
        "spatial_celltype_map": plots_dir / f"{prefix}_spatial_celltype_map.png",
        "celltype_distribution": plots_dir / f"{prefix}_celltype_distribution.png",
        "umap_analysis": plots_dir / f"{prefix}_umap_analysis.png",
        "spatial_leiden_clusters": plots_dir / f"{prefix}_spatial_leiden_clusters.png",
        "marker_genes": plots_dir / f"{prefix}_marker_genes.png"
    }
    
    # Create HTML content
    html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Cell Type Annotation Report</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            margin: 20px;
            line-height: 1.6;
        }}
        .header {{
            background-color: #f4f4f4;
            padding: 20px;
            border-radius: 5px;
            margin-bottom: 20px;
        }}
        .section {{
            margin-bottom: 30px;
        }}
        .image-container {{
            text-align: center;
            margin: 20px 0;
        }}
        .image-container img {{
            max-width: 100%;
            height: auto;
            border: 1px solid #ddd;
            border-radius: 5px;
        }}
        table {{
            border-collapse: collapse;
            width: 100%;
            margin: 10px 0;
        }}
        th, td {{
            border: 1px solid #ddd;
            padding: 8px;
            text-align: left;
        }}
        th {{
            background-color: #f2f2f2;
        }}
        .parameter-table {{
            max-width: 600px;
        }}
        .stats-table {{
            max-width: 500px;
        }}
        h1 {{
            color: #333;
        }}
        h2 {{
            color: #666;
            border-bottom: 2px solid #eee;
            padding-bottom: 5px;
        }}
        .timestamp {{
            color: #888;
            font-size: 0.9em;
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>Cell Type Annotation Report</h1>
        <p class="timestamp">Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        <p><strong>Xenium Data:</strong> {xenium_path.name}</p>
        <p><strong>Reference Data:</strong> {reference_path.name}</p>
    </div>

    <div class="section">
        <h2>Run Parameters</h2>
        <table class="parameter-table">"""
    
    # Add parameters table
    if run_parameters:
        for key, value in run_parameters.items():
            html_content += f"<tr><td>{key}</td><td>{value}</td></tr>"
    
    html_content += f"""
        </table>
    </div>

    <div class="section">
        <h2>Dataset Overview</h2>
        <table class="stats-table">
            <tr><td>Total Cells</td><td>{adata.n_obs:,}</td></tr>
            <tr><td>Total Genes</td><td>{adata.n_vars:,}</td></tr>
            <tr><td>Cell Types Predicted</td><td>{len(cell_type_counts)}</td></tr>
        </table>
    </div>

    <div class="section">
        <h2>Cell Type Distribution</h2>
        <table>
            <tr><th>Cell Type</th><th>Count</th><th>Percentage</th></tr>"""
    
    # Add cell type counts
    total_cells = cell_type_counts.sum()
    for cell_type, count in cell_type_counts.items():
        percentage = (count / total_cells) * 100
        html_content += f"<tr><td>{cell_type}</td><td>{count:,}</td><td>{percentage:.1f}%</td></tr>"
    
    html_content += """
        </table>
    </div>

    <div class="section">
        <h2>Visualizations</h2>"""
    
    # Add images
    image_titles = {
        "spatial_celltype_map": "Spatial Cell Type Map",
        "celltype_distribution": "Cell Type Distribution",
        "umap_analysis": "UMAP Analysis",
        "spatial_leiden_clusters": "Spatial Leiden Clusters",
        "marker_genes": "Marker Genes Analysis"
    }
    
    for key, title in image_titles.items():
        if key in images and images[key].exists():
            img_base64 = image_to_base64(images[key])
            html_content += f"""
        <h3>{title}</h3>
        <div class="image-container">
            <img src="data:image/png;base64,{img_base64}" alt="{title}">
        </div>"""
    
    html_content += f"""
    </div>

    <div class="section">
        <h2>Files Generated</h2>
        <ul>"""
    
    # Add zarr file info conditionally
    if zarr_saved:
        html_content += f"<li><strong>Annotated Data:</strong> data/{xenium_path.stem}_annotated.zarr</li>"
    else:
        html_content += "<li><strong>Annotated Data:</strong> Not generated (use --save-annotated-zarr to enable)</li>"
    
    html_content += f"""
            <li><strong>Cell Type Summary (Excel):</strong> data/{xenium_path.stem}_celltype_summary.xlsx</li>
            <li><strong>Log File:</strong> logs/celltype_annotation.log</li>
            <li><strong>Plots Directory:</strong> plots/</li>
            <li><strong>Run Metadata:</strong> run_metadata.json</li>
            <li><strong>This Report:</strong> annotation_report.html</li>
        </ul>
    </div>

</body>
</html>"""
    
    # Write HTML file
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    logging.info(f"HTML report generated: {report_path}")


def create_label_lookup_table(sdata: sd.SpatialData, table_name: str = "table", label_name: str = "cell_labels") -> None:
    """
    Create a lookup mechanism for napari to display cell type colors on labels.

    Since labels have arbitrary integer pixel values and table has string cell_ids,
    this function ensures napari-spatialdata can correctly map between them.

    The key insight: napari-spatialdata uses the shapes layer cell_ids as the source
    of truth, so we need to ensure our table's cell_ids match those in the shapes.

    Args:
        sdata: SpatialData object
        table_name: Name of table in SpatialData object
        label_name: Name of label layer to check

    Side Effects:
        Logs information about ID compatibility
    """
    table = sdata.tables[table_name]

    # Check if cell_ids are already integers
    if table.obs['cell_id'].dtype in [np.int32, np.int64]:
        logging.info("✅ Cell IDs are integers - labels should work directly in napari")
        return

    # For string cell_ids, check if there's a shapes layer we can reference
    has_shapes = False
    shape_id_source = None

    for shape_name in ['cell_boundaries', 'cell_circles', 'nucleus_boundaries']:
        if shape_name in sdata.shapes:
            has_shapes = True
            shape_id_source = shape_name
            break

    if has_shapes:
        # Get the shapes dataframe to check IDs
        shapes_df = sdata.shapes[shape_id_source]
        shape_ids = set(shapes_df.index)
        table_ids = set(table.obs['cell_id'])

        if shape_ids == table_ids:
            logging.info(f"✅ Cell IDs match between table and {shape_id_source}")
            logging.info(f"   Labels visualization in napari: Will work via shape reference")
            logging.info(f"   Recommendation: Use {shape_id_source} layer in napari for colored cell types")
        else:
            logging.warning(f"⚠️  Cell ID mismatch between table and {shape_id_source}")
            logging.warning(f"   Table has {len(table_ids)} IDs, shapes have {len(shape_ids)} IDs")
    else:
        logging.warning("⚠️  String cell_ids with no shapes layer - label visualization may not work in napari")
        logging.warning("   Recommendation: Use --attach-to shapes for napari visualization")


def attach_table_to_regions(
    sdata: sd.SpatialData,
    table_name: str = "table",
    attach_to: list[str] = None,
    primary_region: str = None,
    create_int_mapping: bool = True
) -> list[str]:
    """
    Attach table annotations to multiple spatial elements (shapes and/or labels).

    The annotations are LINKED (not copied) - there's only one table that multiple
    spatial elements can reference. This means no storage overhead for attaching
    to multiple regions.

    Args:
        sdata: SpatialData object
        table_name: Name of table in SpatialData object
        attach_to: List of region names to attach to. Can include:
            - 'cell_labels' (raster - fast in napari)
            - 'cell_boundaries' (shapes - slower in napari)
            - 'nucleus_labels' (raster - fast in napari)
            - 'nucleus_boundaries' (shapes - slower in napari)
            - 'cell_circles' (shapes)
            If None, auto-detects and prefers labels over shapes for performance
        primary_region: Primary region for the obs['region'] column. If None,
            uses the first available region.
        create_int_mapping: If True and attaching to labels, creates integer cell_id
            mapping for napari compatibility (default: True)

    Returns:
        List of regions the table was successfully attached to
    """
    # Auto-detect available regions if not specified
    if attach_to is None:
        attach_to = []
        # Prefer shapes first (they work for visualization), then labels
        # Note: The PRIMARY region (first in list) is what napari will use for display
        for shape_name in ['cell_boundaries', 'cell_circles', 'nucleus_boundaries']:
            if shape_name in sdata.shapes:
                attach_to.append(shape_name)
        # Also include labels (faster loading, but won't show annotations in napari)
        for label_name in ['cell_labels', 'nucleus_labels']:
            if label_name in sdata.labels:
                attach_to.append(label_name)

    # Filter to only existing regions
    available_regions = []
    for region in attach_to:
        if region in sdata.shapes or region in sdata.labels:
            available_regions.append(region)
        else:
            logging.warning(f"Region '{region}' not found in spatial data, skipping")

    if not available_regions:
        logging.warning("No valid regions found to attach table to")
        return []

    # Check if we're attaching to any labels
    has_labels = any(r in sdata.labels for r in available_regions)

    # Check label compatibility for napari visualization
    if has_labels and create_int_mapping:
        logging.info("Detected label regions - checking napari visualization compatibility")
        create_label_lookup_table(sdata, table_name)

    # Always use 'cell_id' as instance_key - napari-spatialdata handles the mapping
    instance_key = 'cell_id'

    # Determine primary region for obs column
    if primary_region is None:
        primary_region = available_regions[0]
    elif primary_region not in available_regions:
        logging.warning(f"Primary region '{primary_region}' not in available regions, using '{available_regions[0]}'")
        primary_region = available_regions[0]

    # Set the obs region column to primary region
    sdata.tables[table_name].obs['region'] = primary_region
    sdata.tables[table_name].obs['region'] = sdata.tables[table_name].obs['region'].astype('category')

    # Set spatialdata_attrs to link to all regions
    sdata.tables[table_name].uns['spatialdata_attrs'] = {
        'region': available_regions,  # List of all regions (annotations are linked, not copied)
        'region_key': 'region',
        'instance_key': instance_key
    }

    logging.info(f"✅ Table '{table_name}' linked to {len(available_regions)} region(s): {', '.join(available_regions)}")
    logging.info(f"   Primary region (obs['region']): {primary_region}")
    logging.info(f"   Instance key: {instance_key}")
    logging.info(f"   📊 Annotations are LINKED (not copied) - no storage overhead")

    # Log performance info
    labels_count = sum(1 for r in available_regions if r in sdata.labels)
    shapes_count = sum(1 for r in available_regions if r in sdata.shapes)
    if labels_count > 0:
        logging.info(f"   ⚡ {labels_count} label region(s) - FAST loading in napari")
    if shapes_count > 0:
        logging.info(f"   🐌 {shapes_count} shape region(s) - slower loading in napari")

    # Important napari visualization note
    if len(available_regions) > 1:
        logging.info(f"\n   ⚠️  IMPORTANT for napari visualization:")
        logging.info(f"   Table annotations will only appear in napari for the PRIMARY region: {primary_region}")
        logging.info(f"   This is because obs['region'] can only point to one spatial element at a time.")
        logging.info(f"   To see annotations on a different layer, you'll need to:")
        logging.info(f"   1. Re-run with that layer as primary (e.g., --attach-to {available_regions[0]})")
        logging.info(f"   2. Or manually change obs['region'] values in the zarr file")

    return available_regions


def export_celltype_summary_to_excel(
    sdata: sd.SpatialData,
    output_dir: Path,
    xenium_path: Path,
    reference_path: Path,
    table_name: str = "table",
    prediction_column: str = "cell_type_predicted"
) -> None:
    """
    Export cell type summary to Excel file.

    Creates an Excel file with columns: ID, celltype, cellcount
    Each row represents one cell type with its count.

    Args:
        sdata: SpatialData object with annotations
        output_dir: Output directory
        xenium_path: Path to Xenium data (for filename generation)
        reference_path: Path to reference data (for filename generation)
        table_name: Name of table in SpatialData object
        prediction_column: Name of prediction column
    """
    data_dir = output_dir / "data"
    data_dir.mkdir(parents=True, exist_ok=True)

    # Generate descriptive filename
    xenium_name = xenium_path.stem
    ref_name = reference_path.stem
    excel_filename = f"{xenium_name}_celltype_summary.xlsx"
    excel_path = data_dir / excel_filename

    # Get cell type counts
    cell_type_counts = sdata.tables[table_name].obs[prediction_column].value_counts()

    # Create DataFrame with required columns: ID, celltype, cellcount
    summary_df = pd.DataFrame({
        'ID': range(1, len(cell_type_counts) + 1),
        'celltype': cell_type_counts.index,
        'cellcount': cell_type_counts.values
    })

    # Export to Excel
    logging.info(f"Exporting cell type summary to Excel: {excel_path}")
    summary_df.to_excel(excel_path, index=False, engine='openpyxl')

    # Log summary statistics
    logging.info(f"Excel export complete:")
    logging.info(f"  - Total cell types: {len(summary_df)}")
    logging.info(f"  - Total cells: {summary_df['cellcount'].sum():,}")
    logging.info(f"  - File: {excel_path}")


def save_run_metadata(output_dir: Path, run_parameters: dict) -> None:
    """
    Save run metadata and parameters to JSON file.

    Args:
        output_dir: Output directory
        run_parameters: Dictionary of run parameters
    """
    metadata_path = output_dir / "run_metadata.json"

    metadata = {
        "timestamp": datetime.now().isoformat(),
        "parameters": run_parameters,
        "tool_version": "2.0.0",
        "command": " ".join(sys.argv)
    }

    with open(metadata_path, 'w') as f:
        json.dump(metadata, f, indent=2)

    logging.info(f"Run metadata saved: {metadata_path}")


def annotate_spatial_data(
    xenium_path: Path,
    reference_path: Path,
    output_dir: Path,
    label_column: str = "cell_type",
    target_sum: float = 1e4,
    table_name: str = "table",
    prediction_column: str = "cell_type_predicted",
    feature_selection: bool = False,
    min_common_genes: int = 50,
    generate_plots: bool = True,
    min_clusters: int = None,
    max_clusters: int = None,
    consolidate_data: bool = False,
    overwrite: bool = False,
    show_unknown_cells: bool = False,
    save_annotated_zarr: bool = False,
    load_model: Path = None,
    attach_to: str = "auto"
) -> str:
    """
    Main annotation pipeline.

    Args:
        xenium_path: Path to Xenium .zarr file
        reference_path: Path to reference .h5ad file
        output_dir: Output directory for all results
        label_column: Column name containing cell type labels in reference
        target_sum: Target sum for normalization
        table_name: Name of table in SpatialData object
        prediction_column: Name for prediction column
        feature_selection: Whether to use feature selection in model training
        min_common_genes: Minimum number of common genes required
        generate_plots: Whether to generate visualizations
        min_clusters: Minimum number of clusters for Leiden clustering
        max_clusters: Maximum number of clusters for Leiden clustering
        consolidate_data: Whether to create self-contained data copy
        overwrite: Whether to overwrite existing output
        show_unknown_cells: Whether to include unknown cells in spatial plots
        save_annotated_zarr: Whether to save the annotated zarr file (default: False)
        load_model: Path to existing model file to load instead of training new one
        attach_to: Comma-separated list of regions to attach annotations to, or special values:
            - "auto" (default): Auto-detect and attach to all available regions (labels + shapes)
            - "labels": Only attach to label regions (cell_labels, nucleus_labels) - fastest in napari
            - "shapes": Only attach to shape regions (cell_boundaries, nucleus_boundaries, cell_circles)
            - Custom: e.g., "cell_labels,cell_boundaries" for specific regions

    Returns:
        Path to saved zarr file or empty string if not saved
    """
    # Load datasets
    logging.info(f"Loading Xenium data from: {xenium_path}")
    sdata = sd.read_zarr(xenium_path)
    
    logging.info(f"Loading reference data from: {reference_path}")
    adata_ref = anndata.read_h5ad(reference_path)
    
    # Extract Xenium table
    if table_name not in sdata.tables:
        raise ValueError(f"Table '{table_name}' not found in spatial data. Available tables: {list(sdata.tables.keys())}")
    
    adata_xenium = sdata.tables[table_name].copy()
    
    logging.info(f"Xenium data: {adata_xenium.n_obs} cells, {adata_xenium.n_vars} genes")
    logging.info(f"Reference data: {adata_ref.n_obs} cells, {adata_ref.n_vars} genes")
    
    # Harmonize gene names
    adata_xenium, adata_ref = harmonize_gene_names(adata_xenium, adata_ref)
    
    # Check minimum gene requirement
    if adata_xenium.n_vars < min_common_genes:
        raise ValueError(f"Only {adata_xenium.n_vars} common genes found, minimum {min_common_genes} required")
    
    # Normalize data
    normalize_data(adata_ref, target_sum)
    normalize_data(adata_xenium, target_sum)

    # Model handling: load existing model or train new one
    if load_model is not None:
        # Load existing model
        logging.info(f"Loading existing CellTypist model from: {load_model}")
        model = load_celltypist_model(load_model)
    else:
        # Train new model and save it by default
        model_save_path = create_model_save_path(reference_path, label_column)
        logging.info(f"Training new CellTypist model and saving to: {model_save_path}")
        model = train_celltypist_model(adata_ref, label_column, feature_selection, model_save_path)
    
    # Predict cell types
    predictions = predict_cell_types(adata_xenium, model)
    
    # Add predictions to spatial data
    logging.info(f"Adding predictions as column: {prediction_column}")
    sdata.tables[table_name].obs[prediction_column] = predictions['predicted_labels']
    sdata.tables[table_name].obs[prediction_column] = sdata.tables[table_name].obs[prediction_column].astype('category')

    # Export cell type summary to Excel
    export_celltype_summary_to_excel(
        sdata=sdata,
        output_dir=output_dir,
        xenium_path=xenium_path,
        reference_path=reference_path,
        table_name=table_name,
        prediction_column=prediction_column
    )

    # Parse attach_to parameter and attach table to regions
    logging.info("Attaching table annotations to spatial regions")
    attach_regions = None

    if attach_to == "auto":
        # Auto-detect all available regions
        attach_regions = None  # Let attach_table_to_regions auto-detect
    elif attach_to == "labels":
        # Only labels (fast in napari)
        attach_regions = ['cell_labels', 'nucleus_labels']
    elif attach_to == "shapes":
        # Only shapes (slower in napari)
        attach_regions = ['cell_boundaries', 'nucleus_boundaries', 'cell_circles']
    else:
        # Custom comma-separated list
        attach_regions = [r.strip() for r in attach_to.split(',')]

    attached_regions = attach_table_to_regions(
        sdata,
        table_name=table_name,
        attach_to=attach_regions
    )
    
    # Save annotated data to data subdirectory if requested
    zarr_output_path = ""
    if save_annotated_zarr:
        zarr_output_path = output_dir / "data" / f"{xenium_path.stem}_annotated.zarr"
        logging.info(f"Saving annotated spatial data to: {zarr_output_path}")
        
        if consolidate_data:
            logging.info("Creating self-contained data copy (this may take some time for large datasets)...")
            # Note: SpatialData write() doesn't have direct consolidate parameter
            # For now, we use the standard write and document the limitation
            sdata.write(zarr_output_path, overwrite=overwrite)
            logging.info("✅ Data saved to results directory")
            logging.info("ℹ️  Note: Large data elements (images, shapes) may still reference original locations")
            logging.info("   This is a current limitation of SpatialData - full consolidation not yet supported")
        else:
            # Standard save - keeps references to original data locations
            sdata.write(zarr_output_path, overwrite=overwrite)
            logging.info("ℹ️  Data saved with references to original locations (not self-contained)")
            logging.info("   Use --consolidate-data flag for more detailed consolidation information")
        zarr_output_path = str(zarr_output_path)
    else:
        logging.info("Skipping zarr file generation (use --save-annotated-zarr to enable)")
    
    # Generate visualizations if requested
    if generate_plots:
        # Suppress scanpy and squidpy plotting settings to avoid conflicts
        sc.settings.verbosity = 1  # Reduce scanpy verbosity
        sc.settings.set_figure_params(dpi=80, facecolor='white')  # Set figure parameters
        
        generate_visualizations(sdata, xenium_path, reference_path, output_dir, table_name, prediction_column, min_clusters, max_clusters, show_unknown_cells)
    
    # Print summary
    cell_type_counts = sdata.tables[table_name].obs[prediction_column].value_counts()
    print(f"Annotation summary:")
    for cell_type, count in cell_type_counts.items():
        print(f"  {cell_type}: {count} cells")
    
    # Save run parameters for HTML report
    run_parameters = {
        "Xenium Data": str(xenium_path),
        "Reference Data": str(reference_path),
        "Label Column": label_column,
        "Target Sum": target_sum,
        "Table Name": table_name,
        "Prediction Column": prediction_column,
        "Feature Selection": feature_selection,
        "Min Common Genes": min_common_genes,
        "Min Clusters": min_clusters if min_clusters is not None else "Not specified",
        "Max Clusters": max_clusters if max_clusters is not None else "Not specified",
        "Generate Plots": generate_plots,
        "Save Annotated Zarr": save_annotated_zarr,
        "Load Model": str(load_model) if load_model is not None else "None (trained new model)",
        "Attach To Regions": f"{attach_to} → {', '.join(attached_regions)}" if attached_regions else attach_to
    }
    
    # Save metadata
    save_run_metadata(output_dir, run_parameters)

    # Generate HTML report (only if plots were generated)
    if generate_plots:
        generate_html_report(
            output_dir, xenium_path, reference_path, sdata,
            table_name, prediction_column, run_parameters, save_annotated_zarr
        )
    else:
        logging.info("Skipping HTML report generation (no plots generated)")

    return zarr_output_path


@click.command()
@click.argument("xenium_data", type=click.Path(exists=True, path_type=Path))
@click.argument("reference_data", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--label-column",
    default="cell_type",
    help="Column name containing cell type labels in reference data.",
    show_default=True
)
@click.option(
    "--target-sum",
    type=float,
    default=1e4,
    help="Target sum for normalization.",
    show_default=True
)
@click.option(
    "--table-name",
    default="table",
    help="Name of table in SpatialData object containing cell data.",
    show_default=True
)
@click.option(
    "--prediction-column",
    default="cell_type_predicted",
    help="Name for the prediction column to add.",
    show_default=True
)
@click.option(
    "--feature-selection",
    is_flag=True,
    help="Use feature selection in model training."
)
@click.option(
    "--min-common-genes",
    type=int,
    default=50,
    help="Minimum number of common genes required.",
    show_default=True
)
@click.option(
    "--log-level",
    type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR"], case_sensitive=False),
    default="INFO",
    help="Logging level.",
    show_default=True
)
@click.option(
    "--overwrite",
    is_flag=True,
    help="Overwrite output directory if it exists."
)
@click.option(
    "--generate-plots",
    is_flag=True,
    default=True,
    help="Generate visualization plots.",
    show_default=True
)
@click.option(
    "--min-clusters",
    type=int,
    help="Minimum number of clusters for Leiden clustering. If not set, clustering is not constrained."
)
@click.option(
    "--max-clusters",
    type=int,
    help="Maximum number of clusters for Leiden clustering. If not set, clustering is not constrained."
)
@click.option(
    "--results-dir",
    default="results",
    help="Base directory for storing results (default: 'results').",
    show_default=True
)
@click.option(
    "--consolidate-data",
    is_flag=True,
    default=False,
    help="Create a self-contained copy of all data in the results directory (increases storage but improves portability).",
    show_default=True
)
@click.option(
    "--show-unknown-cells",
    is_flag=True,
    default=False,
    help="Include cells with unknown/unassigned cell types in spatial plots (default: hide them for cleaner visualization).",
    show_default=True
)
@click.option(
    "--save-annotated-zarr",
    is_flag=True,
    default=False,
    help="Save the annotated zarr file to the results directory (default: False, only generates plots and reports).",
    show_default=True
)
@click.option(
    "--load-model",
    type=click.Path(exists=True, path_type=Path),
    help="Path to existing CellTypist model file (.pkl). If provided, skips model training and uses this model instead."
)
@click.option(
    "--attach-to",
    default="auto",
    help=(
        "Regions to attach cell type annotations to. Options: "
        "'auto' (all available - default), "
        "'labels' (cell_labels, nucleus_labels - fastest in napari), "
        "'shapes' (cell_boundaries, nucleus_boundaries, cell_circles), "
        "or custom comma-separated list (e.g., 'cell_labels,cell_boundaries'). "
        "Annotations are LINKED (not copied), so no storage overhead."
    ),
    show_default=True
)
@click.option(
    "--zarr-only",
    is_flag=True,
    default=False,
    help=(
        "Only generate the annotated zarr file, skip all analysis and visualizations. "
        "Automatically enables --save-annotated-zarr and disables --generate-plots. "
        "Much faster when you only need the annotations."
    ),
    show_default=True
)
def main(
    xenium_data: Path,
    reference_data: Path,
    label_column: str,
    target_sum: float,
    table_name: str,
    prediction_column: str,
    feature_selection: bool,
    min_common_genes: int,
    log_level: str,
    overwrite: bool,
    generate_plots: bool,
    min_clusters: int,
    max_clusters: int,
    results_dir: str,
    consolidate_data: bool,
    show_unknown_cells: bool,
    save_annotated_zarr: bool,
    load_model: Path,
    attach_to: str,
    zarr_only: bool
):
    """
    Annotate cell types in Xenium spatial transcriptomics data using single-cell reference.

    Creates an organized output directory with visualizations, logs, and HTML report.
    Use --save-annotated-zarr to also save the annotated zarr file.

    XENIUM_DATA: Path to Xenium spatial data (.zarr file)

    REFERENCE_DATA: Path to single-cell reference data (.h5ad file)
    """
    # If zarr-only mode, override related flags
    if zarr_only:
        save_annotated_zarr = True
        generate_plots = False

    # Create output directory
    output_dir = create_output_directory(xenium_data, reference_data, min_clusters, max_clusters, results_dir)

    # Setup logging with output directory
    setup_logging(output_dir, log_level.upper())

    # Log zarr-only mode after logging is configured
    if zarr_only:
        logging.info("=" * 80)
        logging.info("ZARR-ONLY MODE ENABLED")
        logging.info("  - Annotated zarr file will be saved")
        logging.info("  - All visualizations and plots will be SKIPPED")
        logging.info("  - HTML report generation will be SKIPPED")
        logging.info("  - This mode is much faster when you only need the annotations")
        logging.info("=" * 80)
    
    # Validate inputs (only check zarr output if we're going to save it)
    if save_annotated_zarr:
        zarr_output_path = output_dir / "data" / f"{xenium_data.stem}_annotated.zarr"
        if zarr_output_path.exists() and not overwrite:
            click.echo(f"Error: Output file already exists: {zarr_output_path}. Use --overwrite to overwrite.", err=True)
            sys.exit(1)
    
    if min_clusters is not None and max_clusters is not None and min_clusters > max_clusters:
        click.echo(f"Error: min-clusters ({min_clusters}) cannot be greater than max-clusters ({max_clusters}).", err=True)
        sys.exit(1)
    
    if min_clusters is not None and min_clusters < 1:
        click.echo(f"Error: min-clusters must be at least 1, got {min_clusters}.", err=True)
        sys.exit(1)
    
    if max_clusters is not None and max_clusters < 1:
        click.echo(f"Error: max-clusters must be at least 1, got {max_clusters}.", err=True)
        sys.exit(1)
    
    try:
        result_path = annotate_spatial_data(
            xenium_path=xenium_data,
            reference_path=reference_data,
            output_dir=output_dir,
            label_column=label_column,
            target_sum=target_sum,
            table_name=table_name,
            prediction_column=prediction_column,
            feature_selection=feature_selection,
            min_common_genes=min_common_genes,
            generate_plots=generate_plots,
            min_clusters=min_clusters,
            max_clusters=max_clusters,
            consolidate_data=consolidate_data,
            overwrite=overwrite,
            show_unknown_cells=show_unknown_cells,
            save_annotated_zarr=save_annotated_zarr,
            load_model=load_model,
            attach_to=attach_to
        )
        
        # Standard output: paths to key outputs
        if result_path:
            print(f"Zarr file: {result_path}")
        else:
            print("Zarr file: Not generated (use --save-annotated-zarr to enable)")
        print(f"Output directory: {output_dir}")
        print(f"HTML report: {output_dir / 'annotation_report.html'}")
        
    except Exception as e:
        click.echo(f"❌ Error during annotation: {e}", err=True)
        sys.exit(1)


if __name__ == "__main__":
    main()