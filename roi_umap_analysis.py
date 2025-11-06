#!/usr/bin/env python3
"""
ROI Analysis and UMAP Visualization

This script provides tools to:
1. Define ROIs in napari and save them
2. Extract cells within ROIs from spatial data
3. Perform UMAP clustering on ROI-selected cells
4. Visualize ROI-based clustering results

Usage:
    python roi_umap_analysis.py --help
"""

import napari
import spatialdata as sd
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import json
import anndata as ad
from shapely.geometry import Point, Polygon
import geopandas as gpd
from typing import List, Dict, Optional, Tuple
import click
import warnings
warnings.filterwarnings('ignore')

# Set up scanpy
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white')

class ROIAnalyzer:
    """Class to handle ROI definition and UMAP analysis"""
    
    def __init__(self, zarr_path: str):
        """Initialize with spatial data"""
        self.zarr_path = zarr_path
        self.sdata = sd.read_zarr(zarr_path)
        self.table = self.sdata.tables['table']
        self.roi_data = {}
        
    def launch_napari_for_roi_definition(self, save_path: Optional[str] = None):
        """Launch napari for ROI definition"""
        print("Launching napari for ROI definition...")
        print("Instructions:")
        print("1. Use the 'Shapes' layer to draw ROIs")
        print("2. Draw polygons around areas of interest")
        print("3. Name your ROIs in the shapes layer")
        print("4. Save the ROIs when finished")
        
        # Create napari viewer
        viewer = napari.Viewer()
        
        # Add the spatial data to napari
        try:
            # Add images
            for img_name, img_data in self.sdata.images.items():
                viewer.add_image(img_data, name=img_name)
            
            # Add labels
            for label_name, label_data in self.sdata.labels.items():
                viewer.add_labels(label_data, name=label_name)
                
            # Add points (transcripts)
            for point_name, point_data in self.sdata.points.items():
                if 'x' in point_data.columns and 'y' in point_data.columns:
                    coords = point_data[['y', 'x']].values  # napari expects y,x order
                    viewer.add_points(coords, name=point_name, size=1, opacity=0.5)
            
            # Add shapes (cell boundaries)
            for shape_name, shape_data in self.sdata.shapes.items():
                # Convert geometries to napari format
                shapes_list = []
                for geom in shape_data.geometry:
                    if geom.geom_type == 'Polygon':
                        # Extract exterior coordinates and convert to y,x order
                        coords = np.array(geom.exterior.coords)
                        coords = coords[:, [1, 0]]  # swap x,y to y,x
                        shapes_list.append(coords)
                
                if shapes_list:
                    viewer.add_shapes(shapes_list, name=shape_name, shape_type='polygon', 
                                    edge_width=0.5, edge_color='white', face_color='transparent')
            
            # Add an empty shapes layer for ROI drawing
            roi_layer = viewer.add_shapes(name='ROIs', shape_type='polygon', 
                                        edge_width=2, edge_color='red', face_color='transparent')
            
            print("Napari viewer launched. Draw your ROIs and close the viewer when done.")
            
            # Show the viewer
            napari.run()
            
            # Extract ROI data after viewer is closed
            if len(roi_layer.data) > 0:
                self.roi_data = self._extract_roi_data(roi_layer)
                if save_path:
                    self._save_roi_data(save_path)
                    print(f"ROI data saved to: {save_path}")
                return self.roi_data
            else:
                print("No ROIs were defined.")
                return {}
                
        except Exception as e:
            print(f"Error launching napari: {e}")
            return {}
    
    def _extract_roi_data(self, roi_layer) -> Dict:
        """Extract ROI data from napari shapes layer"""
        roi_data = {}
        
        for i, shape_coords in enumerate(roi_layer.data):
            # Convert from napari y,x format to x,y format
            coords_xy = shape_coords[:, [1, 0]]
            
            # Create polygon
            polygon = Polygon(coords_xy)
            
            # Use layer properties if available, otherwise use index
            roi_name = f"ROI_{i+1}"
            if hasattr(roi_layer, 'text') and len(roi_layer.text) > i:
                roi_name = roi_layer.text[i] or roi_name
            
            roi_data[roi_name] = {
                'polygon': polygon,
                'coordinates': coords_xy.tolist(),
                'area': polygon.area,
                'bounds': polygon.bounds
            }
        
        return roi_data
    
    def load_roi_data(self, roi_path: str) -> Dict:
        """Load ROI data from JSON file"""
        with open(roi_path, 'r') as f:
            data = json.load(f)
        
        # Reconstruct polygons
        roi_data = {}
        for name, roi_info in data.items():
            polygon = Polygon(roi_info['coordinates'])
            roi_data[name] = {
                'polygon': polygon,
                'coordinates': roi_info['coordinates'],
                'area': roi_info['area'],
                'bounds': roi_info['bounds']
            }
        
        self.roi_data = roi_data
        return roi_data
    
    def _save_roi_data(self, save_path: str):
        """Save ROI data to JSON file"""
        # Convert polygons to serializable format
        save_data = {}
        for name, roi_info in self.roi_data.items():
            save_data[name] = {
                'coordinates': roi_info['coordinates'],
                'area': roi_info['area'],
                'bounds': roi_info['bounds']
            }
        
        with open(save_path, 'w') as f:
            json.dump(save_data, f, indent=2)
    
    def extract_cells_in_rois(self) -> Dict[str, List[int]]:
        """Extract cell indices that fall within each ROI"""
        if not self.roi_data:
            print("No ROI data available. Please define ROIs first.")
            return {}
        
        print("Extracting cells within ROIs...")
        
        # Get cell positions from cell_circles or centroids
        cell_positions = None
        cell_indices = None
        
        # Try to get cell positions from shapes (cell_circles)
        for shape_name, shape_data in self.sdata.shapes.items():
            if 'cell_circles' in shape_name or 'circles' in shape_name:
                # Extract centroids of circles
                centroids = []
                indices = []
                for idx, geom in enumerate(shape_data.geometry):
                    if geom.geom_type == 'Point':
                        centroids.append([geom.x, geom.y])
                        indices.append(idx)
                    elif hasattr(geom, 'centroid'):
                        centroids.append([geom.centroid.x, geom.centroid.y])
                        indices.append(idx)
                
                if centroids:
                    cell_positions = np.array(centroids)
                    cell_indices = np.array(indices)
                    break
        
        # If no cell positions found, use table obs coordinates if available
        if cell_positions is None:
            obs_data = self.table.obs
            if 'x' in obs_data.columns and 'y' in obs_data.columns:
                cell_positions = obs_data[['x', 'y']].values
                cell_indices = np.arange(len(obs_data))
            else:
                print("No cell position information found.")
                return {}
        
        # Find cells in each ROI
        roi_cells = {}
        for roi_name, roi_info in self.roi_data.items():
            polygon = roi_info['polygon']
            cells_in_roi = []
            
            for i, (x, y) in enumerate(cell_positions):
                point = Point(x, y)
                if polygon.contains(point):
                    cells_in_roi.append(cell_indices[i])
            
            roi_cells[roi_name] = cells_in_roi
            print(f"{roi_name}: {len(cells_in_roi)} cells")
        
        return roi_cells
    
    def perform_umap_analysis(self, roi_cells: Dict[str, List[int]], 
                            output_dir: str = "roi_analysis_results") -> Dict:
        """Perform UMAP analysis on ROI-selected cells"""
        
        # Create output directory
        Path(output_dir).mkdir(exist_ok=True)
        
        results = {}
        
        for roi_name, cell_indices in roi_cells.items():
            if len(cell_indices) < 10:  # Skip ROIs with too few cells
                print(f"Skipping {roi_name}: too few cells ({len(cell_indices)})")
                continue
            
            print(f"Analyzing {roi_name} with {len(cell_indices)} cells...")
            
            # Extract cells for this ROI
            roi_adata = self.table[cell_indices].copy()
            
            # Basic preprocessing
            sc.pp.normalize_total(roi_adata, target_sum=1e4)
            sc.pp.log1p(roi_adata)
            
            # Feature selection if enough cells
            if len(cell_indices) > 50:
                sc.pp.highly_variable_genes(roi_adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
                roi_adata.raw = roi_adata
                sc.pp.scale(roi_adata, max_value=10)
            
            # PCA
            sc.tl.pca(roi_adata, svd_solver='arpack')
            
            # Compute neighborhood graph
            sc.pp.neighbors(roi_adata, n_neighbors=min(10, len(cell_indices)//2), n_pcs=40)
            
            # UMAP
            sc.tl.umap(roi_adata)
            
            # Leiden clustering
            sc.tl.leiden(roi_adata, resolution=0.5)
            
            # Store results
            results[roi_name] = {
                'adata': roi_adata,
                'n_cells': len(cell_indices),
                'cell_indices': cell_indices
            }
            
            # Save UMAP plot
            self._plot_roi_umap(roi_adata, roi_name, output_dir)
        
        return results
    
    def _plot_roi_umap(self, adata: ad.AnnData, roi_name: str, output_dir: str):
        """Create UMAP plots for ROI"""
        
        # Create figure with subplots
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        
        # UMAP colored by clusters
        sc.pl.umap(adata, color='leiden', ax=axes[0], show=False, frameon=False)
        axes[0].set_title(f'{roi_name} - Leiden Clusters')
        
        # UMAP colored by cell type if available
        if 'cell_type_predicted' in adata.obs.columns:
            sc.pl.umap(adata, color='cell_type_predicted', ax=axes[1], show=False, frameon=False)
            axes[1].set_title(f'{roi_name} - Cell Types')
        else:
            axes[1].text(0.5, 0.5, 'No cell type\nannotations', ha='center', va='center', 
                        transform=axes[1].transAxes)
            axes[1].set_title(f'{roi_name} - Cell Types (N/A)')
        
        # UMAP colored by total UMI count
        sc.pl.umap(adata, color='total_counts', ax=axes[2], show=False, frameon=False)
        axes[2].set_title(f'{roi_name} - Total UMI Count')
        
        plt.tight_layout()
        plt.savefig(f'{output_dir}/{roi_name}_umap_analysis.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # Create cluster composition plot if cell types available
        if 'cell_type_predicted' in adata.obs.columns:
            self._plot_cluster_composition(adata, roi_name, output_dir)
    
    def _plot_cluster_composition(self, adata: ad.AnnData, roi_name: str, output_dir: str):
        """Plot cluster composition by cell type"""
        
        # Create composition DataFrame
        comp_df = pd.crosstab(adata.obs['leiden'], adata.obs['cell_type_predicted'])
        comp_df_pct = comp_df.div(comp_df.sum(axis=1), axis=0) * 100
        
        # Plot
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        # Absolute counts
        sns.heatmap(comp_df.T, annot=True, fmt='d', cmap='Blues', ax=ax1)
        ax1.set_title(f'{roi_name} - Cell Type Counts by Cluster')
        ax1.set_ylabel('Cell Type')
        ax1.set_xlabel('Leiden Cluster')
        
        # Percentages
        sns.heatmap(comp_df_pct.T, annot=True, fmt='.1f', cmap='Reds', ax=ax2)
        ax2.set_title(f'{roi_name} - Cell Type Percentages by Cluster')
        ax2.set_ylabel('Cell Type')
        ax2.set_xlabel('Leiden Cluster')
        
        plt.tight_layout()
        plt.savefig(f'{output_dir}/{roi_name}_cluster_composition.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    def create_combined_analysis(self, results: Dict, output_dir: str = "roi_analysis_results"):
        """Create combined analysis across all ROIs"""
        
        print("Creating combined ROI analysis...")
        
        # Combine data from all ROIs
        combined_data = []
        for roi_name, roi_result in results.items():
            adata = roi_result['adata']
            roi_df = pd.DataFrame({
                'roi': roi_name,
                'cell_id': adata.obs.index,
                'leiden_cluster': adata.obs['leiden'],
                'umap_1': adata.obsm['X_umap'][:, 0],
                'umap_2': adata.obsm['X_umap'][:, 1],
                'n_genes': adata.obs.get('n_genes_by_counts', adata.obs.get('n_genes', 0)),
                'total_counts': adata.obs.get('total_counts', 0)
            })
            
            # Add cell type if available
            if 'cell_type_predicted' in adata.obs.columns:
                roi_df['cell_type'] = adata.obs['cell_type_predicted']
            
            combined_data.append(roi_df)
        
        # Create combined DataFrame
        combined_df = pd.concat(combined_data, ignore_index=True)
        
        # Save combined results
        combined_df.to_csv(f'{output_dir}/combined_roi_analysis.csv', index=False)
        
        # Create comparative plots
        self._plot_roi_comparison(combined_df, output_dir)
        
        return combined_df
    
    def _plot_roi_comparison(self, combined_df: pd.DataFrame, output_dir: str):
        """Create comparative plots across ROIs"""
        
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # ROI cell counts
        roi_counts = combined_df['roi'].value_counts()
        roi_counts.plot(kind='bar', ax=axes[0, 0])
        axes[0, 0].set_title('Number of Cells per ROI')
        axes[0, 0].set_ylabel('Cell Count')
        axes[0, 0].tick_params(axis='x', rotation=45)
        
        # Gene expression distribution by ROI
        sns.boxplot(data=combined_df, x='roi', y='n_genes', ax=axes[0, 1])
        axes[0, 1].set_title('Gene Expression Distribution by ROI')
        axes[0, 1].set_ylabel('Number of Genes')
        axes[0, 1].tick_params(axis='x', rotation=45)
        
        # Total counts distribution by ROI
        sns.boxplot(data=combined_df, x='roi', y='total_counts', ax=axes[1, 0])
        axes[1, 0].set_title('Total UMI Count Distribution by ROI')
        axes[1, 0].set_ylabel('Total UMI Count')
        axes[1, 0].tick_params(axis='x', rotation=45)
        
        # Cell type distribution by ROI (if available)
        if 'cell_type' in combined_df.columns:
            celltype_counts = combined_df.groupby(['roi', 'cell_type']).size().unstack(fill_value=0)
            celltype_counts.plot(kind='bar', stacked=True, ax=axes[1, 1])
            axes[1, 1].set_title('Cell Type Distribution by ROI')
            axes[1, 1].set_ylabel('Cell Count')
            axes[1, 1].tick_params(axis='x', rotation=45)
            axes[1, 1].legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        else:
            axes[1, 1].text(0.5, 0.5, 'No cell type\nannotations available', 
                          ha='center', va='center', transform=axes[1, 1].transAxes)
        
        plt.tight_layout()
        plt.savefig(f'{output_dir}/roi_comparison.png', dpi=300, bbox_inches='tight')
        plt.close()


@click.command()
@click.argument('zarr_path', type=click.Path(exists=True))
@click.option('--roi-file', type=click.Path(), help='Load ROIs from JSON file instead of defining new ones')
@click.option('--output-dir', default='roi_analysis_results', help='Output directory for results')
@click.option('--save-rois', type=click.Path(), help='Save defined ROIs to JSON file')
@click.option('--skip-napari', is_flag=True, help='Skip napari ROI definition (requires --roi-file)')
def main(zarr_path: str, roi_file: Optional[str], output_dir: str, 
         save_rois: Optional[str], skip_napari: bool):
    """
    ROI Analysis and UMAP Visualization
    
    Define ROIs in napari and analyze them with UMAP clustering.
    
    ZARR_PATH: Path to the spatial data zarr file
    
    Examples:
        # Define ROIs in napari and analyze
        python roi_umap_analysis.py combined_direct_coords_annotated.zarr
        
        # Load existing ROIs and analyze
        python roi_umap_analysis.py combined_direct_coords_annotated.zarr --roi-file my_rois.json
        
        # Define ROIs and save them for later use
        python roi_umap_analysis.py combined_direct_coords_annotated.zarr --save-rois my_rois.json
    """
    
    # Initialize analyzer
    analyzer = ROIAnalyzer(zarr_path)
    
    # Load or define ROIs
    if roi_file:
        print(f"Loading ROIs from: {roi_file}")
        analyzer.load_roi_data(roi_file)
    elif not skip_napari:
        print("Launching napari for ROI definition...")
        analyzer.launch_napari_for_roi_definition(save_rois)
    else:
        print("Error: Must provide --roi-file when using --skip-napari")
        return
    
    if not analyzer.roi_data:
        print("No ROI data available. Exiting.")
        return
    
    # Extract cells in ROIs
    roi_cells = analyzer.extract_cells_in_rois()
    
    if not roi_cells:
        print("No cells found in ROIs. Exiting.")
        return
    
    # Perform UMAP analysis
    results = analyzer.perform_umap_analysis(roi_cells, output_dir)
    
    # Create combined analysis
    if results:
        combined_df = analyzer.create_combined_analysis(results, output_dir)
        print(f"\nAnalysis complete! Results saved to: {output_dir}/")
        print(f"Combined results: {output_dir}/combined_roi_analysis.csv")
        
        # Print summary
        print("\nSummary:")
        for roi_name, roi_result in results.items():
            print(f"  {roi_name}: {roi_result['n_cells']} cells analyzed")
    else:
        print("No ROI analysis results generated.")


if __name__ == "__main__":
    main()