#!/usr/bin/env python3
"""
Cell Type Annotation CLI Tool

This tool performs cell type annotation on Xenium spatial transcriptomics data
using a single-cell RNA-seq reference dataset and CellTypist.

Based on the celltype_annotation_celltypist_mouse.ipynb notebook.
"""

import logging
import sys
from pathlib import Path

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


def setup_logging(level: str = "INFO") -> None:
    """Set up logging configuration."""
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
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


def train_celltypist_model(adata_ref: anndata.AnnData, label_column: str, feature_selection: bool = False) -> celltypist.models.Model:
    """Train CellTypist model on reference data."""
    logging.info(f"Training CellTypist model using label column: {label_column}")
    
    if label_column not in adata_ref.obs.columns:
        raise ValueError(f"Label column '{label_column}' not found in reference data")
    
    model = celltypist.train(adata_ref, labels=label_column, feature_selection=feature_selection)
    logging.info("Model training completed")
    
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
        # Use default resolution
        return 0.5
    
    # Test fewer resolutions for speed
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
                # Clean up and break early if we found a perfect fit
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
        
        # Clean up temporary clustering
        del adata.obs['leiden_temp']
    
    logging.info(f"Optimized resolution: {best_resolution} (targeting {min_clusters}-{max_clusters} clusters)")
    return best_resolution


def perform_clustering_analysis(adata: anndata.AnnData, min_clusters: int = None, max_clusters: int = None) -> anndata.AnnData:
    """
    Perform clustering analysis on the data.
    
    Args:
        adata: AnnData object to analyze
        
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
    max_clusters: int = None
) -> None:
    """
    Generate and save visualization plots.
    
    Args:
        sdata: Annotated SpatialData object
        xenium_path: Path to original Xenium data (for filename generation)
        reference_path: Path to reference data (for filename generation)
        output_dir: Directory to save images
        table_name: Name of table in SpatialData object
        prediction_column: Name of prediction column
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate descriptive filename prefix
    xenium_name = xenium_path.stem
    ref_name = reference_path.parent.name  # Use parent directory name for reference
    prefix = f"{xenium_name}_annotated_with_{ref_name}"
    
    logging.info(f"Generating visualizations with prefix: {prefix}")
    
    # Get the main adata object for analysis
    adata_main = sdata.tables[table_name].copy()
    
    # 1. Spatial cell type map
    if 'cell_boundaries' in sdata.shapes:
        try:
            logging.info("Creating spatial cell type map...")
            fig, ax = plt.subplots(figsize=(20, 20))
            
            sdata.pl.render_shapes(
                element="cell_boundaries",
                color=prediction_column,
                legend_kwargs={"loc": "upper left", "bbox_to_anchor": (1, 1)},
                legend_fontsize=8,
                title=f"{xenium_name}: Predicted Cell Types",
                ax=ax
            )
            
            spatial_map_path = output_dir / f"{prefix}_spatial_celltype_map.png"
            plt.savefig(spatial_map_path, dpi=300, bbox_inches='tight')
            plt.close()
            logging.info(f"Saved spatial cell type map: {spatial_map_path}")
            
        except Exception as e:
            logging.warning(f"Failed to create spatial cell type map: {e}")
    
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
        distribution_path = output_dir / f"{prefix}_celltype_distribution.png"
        plt.savefig(distribution_path, dpi=300, bbox_inches='tight')
        plt.close()
        logging.info(f"Saved cell type distribution: {distribution_path}")
        
    except Exception as e:
        logging.warning(f"Failed to create cell type distribution plot: {e}")
    
    # 3. Summary statistics plot
    try:
        logging.info("Creating annotation summary plot...")
        adata = sdata.tables[table_name]
        
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle(f'{xenium_name}: Annotation Summary', fontsize=16)
        
        # Total cells and genes
        axes[0, 0].text(0.5, 0.7, f'Total Cells: {adata.n_obs:,}', 
                        ha='center', va='center', fontsize=14, transform=axes[0, 0].transAxes)
        axes[0, 0].text(0.5, 0.5, f'Total Genes: {adata.n_vars:,}', 
                        ha='center', va='center', fontsize=14, transform=axes[0, 0].transAxes)
        axes[0, 0].text(0.5, 0.3, f'Cell Types: {len(cell_type_counts)}', 
                        ha='center', va='center', fontsize=14, transform=axes[0, 0].transAxes)
        axes[0, 0].set_title('Dataset Overview')
        axes[0, 0].axis('off')
        
        # Top 5 cell types pie chart
        top5_counts = cell_type_counts.head(5)
        other_count = cell_type_counts.iloc[5:].sum() if len(cell_type_counts) > 5 else 0
        
        if other_count > 0:
            pie_data = list(top5_counts.values) + [other_count]
            pie_labels = list(top5_counts.index) + ['Others']
        else:
            pie_data = list(top5_counts.values)
            pie_labels = list(top5_counts.index)
        
        axes[0, 1].pie(pie_data, labels=pie_labels, autopct='%1.1f%%', startangle=90)
        axes[0, 1].set_title('Top Cell Types Distribution')
        
        # Cell type counts (top 10)
        top10_counts = cell_type_counts.head(10)
        axes[1, 0].barh(range(len(top10_counts)), top10_counts.values)
        axes[1, 0].set_yticks(range(len(top10_counts)))
        axes[1, 0].set_yticklabels(top10_counts.index)
        axes[1, 0].set_xlabel('Number of Cells')
        axes[1, 0].set_title('Top 10 Cell Types')
        axes[1, 0].invert_yaxis()
        
        # Reference dataset info
        ref_info_text = f"Reference: {ref_name}\nXenium Data: {xenium_name}"
        axes[1, 1].text(0.1, 0.5, ref_info_text, ha='left', va='center', 
                        fontsize=12, transform=axes[1, 1].transAxes)
        axes[1, 1].set_title('Dataset Information')
        axes[1, 1].axis('off')
        
        plt.tight_layout()
        summary_path = output_dir / f"{prefix}_annotation_summary.png"
        plt.savefig(summary_path, dpi=300, bbox_inches='tight')
        plt.close()
        logging.info(f"Saved annotation summary: {summary_path}")
        
    except Exception as e:
        logging.warning(f"Failed to create annotation summary plot: {e}")
    
    # 4. UMAP plots with clustering analysis
    try:
        logging.info("Performing clustering analysis for UMAP plots...")
        adata_clustered = perform_clustering_analysis(adata_main, min_clusters, max_clusters)
        
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
        umap_path = output_dir / f"{prefix}_umap_analysis.png"
        plt.savefig(umap_path, dpi=300, bbox_inches='tight')
        plt.close()
        logging.info(f"Saved UMAP analysis: {umap_path}")
        
    except Exception as e:
        logging.warning(f"Failed to create UMAP plots: {e}")
        adata_clustered = None
    
    # 5. Spatial scatter plot with leiden clusters
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
            spatial_clusters_path = output_dir / f"{prefix}_spatial_leiden_clusters.png"
            plt.savefig(spatial_clusters_path, dpi=300, bbox_inches='tight')
            plt.close()
            logging.info(f"Saved spatial clusters plot: {spatial_clusters_path}")
            
        except Exception as e:
            logging.warning(f"Failed to create spatial scatter plot: {e}")
    
    # 6. Marker genes analysis
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
            marker_genes_path = output_dir / f"{prefix}_marker_genes.png"
            plt.savefig(marker_genes_path, dpi=300, bbox_inches='tight')
            plt.close()
            logging.info(f"Saved marker genes plot: {marker_genes_path}")
            
        except Exception as e:
            logging.warning(f"Failed to create marker genes plot: {e}")
    
    # 7. Detailed analysis for endothelial cells (if present)
    try:
        cell_types = adata_main.obs[prediction_column].unique()
        if 'endothelial cell' in cell_types:
            logging.info("Creating detailed endothelial cell analysis...")
            
            # Filter for endothelial cells
            sdata_endothelial = deepcopy(sdata)
            endothelial_mask = sdata_endothelial.tables[table_name].obs[prediction_column] == 'endothelial cell'
            sdata_endothelial.tables[table_name] = sdata_endothelial.tables[table_name][endothelial_mask].copy()
            
            # Perform clustering analysis on endothelial cells
            adata_endo = sdata_endothelial.tables[table_name]
            adata_endo_clustered = perform_clustering_analysis(adata_endo, min_clusters, max_clusters)
            
            # Create endothelial cell UMAP
            fig, axes = plt.subplots(1, 2, figsize=(16, 6))
            
            sc.pl.umap(adata_endo_clustered, color="total_counts", ax=axes[0], show=False, frameon=False)
            axes[0].set_title('Endothelial Cells UMAP: Total Counts')
            
            sc.pl.umap(adata_endo_clustered, color="leiden", ax=axes[1], show=False, frameon=False)
            axes[1].set_title('Endothelial Cells UMAP: Leiden Clusters')
            
            plt.tight_layout()
            endo_umap_path = output_dir / f"{prefix}_endothelial_umap.png"
            plt.savefig(endo_umap_path, dpi=300, bbox_inches='tight')
            plt.close()
            logging.info(f"Saved endothelial UMAP: {endo_umap_path}")
            
            # Create endothelial spatial scatter
            if 'leiden_colors' in adata_endo_clustered.uns:
                adata_endo_clustered.uns.pop('leiden_colors')
            
            fig, ax = plt.subplots(figsize=(12, 12))
            sq.pl.spatial_scatter(
                adata_endo_clustered,
                library_id="spatial",
                shape=None,
                color=["leiden"],
                ax=ax,
                frameon=False
            )
            ax.set_title('Endothelial Cells: Spatial Leiden Clusters')
            
            plt.tight_layout()
            endo_spatial_path = output_dir / f"{prefix}_endothelial_spatial_clusters.png"
            plt.savefig(endo_spatial_path, dpi=300, bbox_inches='tight')
            plt.close()
            logging.info(f"Saved endothelial spatial clusters: {endo_spatial_path}")
            
        else:
            logging.info("No endothelial cells found for detailed analysis")
            
    except Exception as e:
        logging.warning(f"Failed to create endothelial cell analysis: {e}")


def annotate_spatial_data(
    xenium_path: Path,
    reference_path: Path,
    output_path: Path,
    label_column: str = "cell_type",
    target_sum: float = 1e4,
    table_name: str = "table",
    prediction_column: str = "cell_type_predicted",
    feature_selection: bool = False,
    min_common_genes: int = 50,
    generate_plots: bool = True,
    plot_output_dir: Path = None,
    min_clusters: int = None,
    max_clusters: int = None
) -> str:
    """
    Main annotation pipeline.
    
    Args:
        xenium_path: Path to Xenium .zarr file
        reference_path: Path to reference .h5ad file
        output_path: Path to save annotated spatial data
        label_column: Column name containing cell type labels in reference
        target_sum: Target sum for normalization
        table_name: Name of table in SpatialData object
        prediction_column: Name for prediction column
        feature_selection: Whether to use feature selection in model training
        min_common_genes: Minimum number of common genes required
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
    
    # Train model
    model = train_celltypist_model(adata_ref, label_column, feature_selection)
    
    # Predict cell types
    predictions = predict_cell_types(adata_xenium, model)
    
    # Add predictions to spatial data
    logging.info(f"Adding predictions as column: {prediction_column}")
    sdata.tables[table_name].obs[prediction_column] = predictions['predicted_labels']
    sdata.tables[table_name].obs[prediction_column] = sdata.tables[table_name].obs[prediction_column].astype('category')
    
    # Fix spatial data metadata for consistency
    logging.info("Fixing spatial data metadata")
    if 'cell_boundaries' in sdata.shapes:
        sdata.tables[table_name].obs['region'] = 'cell_boundaries'
        sdata.tables[table_name].obs['region'] = sdata.tables[table_name].obs['region'].astype('category')
        
        sdata.tables[table_name].uns['spatialdata_attrs'] = {
            'region': 'cell_boundaries',
            'region_key': 'region',
            'instance_key': 'cell_id'
        }
    
    # Save annotated data
    logging.info(f"Saving annotated spatial data to: {output_path}")
    sdata.write(output_path)
    
    # Generate visualizations if requested
    if generate_plots:
        if plot_output_dir is None:
            plot_output_dir = output_path.parent / "plots"
        
        # Suppress scanpy and squidpy plotting settings to avoid conflicts
        sc.settings.verbosity = 1  # Reduce scanpy verbosity
        sc.settings.set_figure_params(dpi=80, facecolor='white')  # Set figure parameters
        
        generate_visualizations(sdata, xenium_path, reference_path, plot_output_dir, table_name, prediction_column, min_clusters, max_clusters)
    
    # Print summary (now quieter for standard output)
    cell_type_counts = sdata.tables[table_name].obs[prediction_column].value_counts()
    print(f"Annotation summary:")
    for cell_type, count in cell_type_counts.items():
        print(f"  {cell_type}: {count} cells")
    
    return str(output_path)


@click.command()
@click.argument("xenium_data", type=click.Path(exists=True, path_type=Path))
@click.argument("reference_data", type=click.Path(exists=True, path_type=Path))
@click.argument("output_path", type=click.Path(path_type=Path))
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
    help="Overwrite output file if it exists."
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
def main(
    xenium_data: Path,
    reference_data: Path,
    output_path: Path,
    label_column: str,
    target_sum: float,
    table_name: str,
    prediction_column: str,
    feature_selection: bool,
    min_common_genes: int,
    log_level: str,
    overwrite: bool,
    generate_plots: bool,
    plot_output_dir: Path,
    min_clusters: int,
    max_clusters: int
):
    """
    Annotate cell types in Xenium spatial transcriptomics data using single-cell reference.
    
    XENIUM_DATA: Path to Xenium spatial data (.zarr file)
    
    REFERENCE_DATA: Path to single-cell reference data (.h5ad file)
    
    OUTPUT_PATH: Path to save annotated spatial data (.zarr file)
    """
    # Setup logging
    setup_logging(log_level.upper())
    
    # Validate inputs
    if output_path.exists() and not overwrite:
        click.echo(f"Error: Output file already exists: {output_path}. Use --overwrite to overwrite.", err=True)
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
            output_path=output_path,
            label_column=label_column,
            target_sum=target_sum,
            table_name=table_name,
            prediction_column=prediction_column,
            feature_selection=feature_selection,
            min_common_genes=min_common_genes,
            generate_plots=generate_plots,
            plot_output_dir=plot_output_dir,
            min_clusters=min_clusters,
            max_clusters=max_clusters
        )
        
        # Standard output: just the path to the zarr file
        print(result_path)
        
    except Exception as e:
        click.echo(f"❌ Error during annotation: {e}", err=True)
        sys.exit(1)


if __name__ == "__main__":
    main()