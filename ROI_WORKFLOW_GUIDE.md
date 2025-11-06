# ROI Analysis and UMAP Visualization Workflow

This guide explains how to define regions of interest (ROIs) in napari and perform UMAP clustering analysis on spatial transcriptomics data.

## Overview

The ROI workflow allows you to:
1. **Define ROIs interactively** in napari by drawing polygons on spatial data
2. **Extract cells within ROIs** to focus analysis on specific tissue regions
3. **Perform UMAP clustering** on ROI-selected cells to identify cell populations
4. **Visualize and compare** results across different ROIs

## Quick Start

### Option 1: Interactive ROI Definition
```bash
# Define ROIs in napari and analyze
python roi_umap_analysis.py combined_direct_coords_annotated.zarr --save-rois my_rois.json

# Or load existing ROIs
python roi_umap_analysis.py combined_direct_coords_annotated.zarr --roi-file my_rois.json
```

### Option 2: Use Example ROIs
```bash
# Create example ROIs for testing
python -c "
import json
example_rois = {
    'Region_1': {
        'coordinates': [[2000, 0], [6000, 0], [6000, 1500], [2000, 1500]],
        'area': 6000000,
        'bounds': [2000, 0, 6000, 1500]
    }
}
with open('my_rois.json', 'w') as f:
    json.dump(example_rois, f, indent=2)
"

# Run analysis
python roi_umap_analysis.py combined_direct_coords_annotated.zarr --roi-file my_rois.json
```

## Detailed Workflow

### Step 1: Define ROIs in napari

When you run the interactive mode, napari will open with:
- **Images**: Morphology and other image layers
- **Labels**: Cell and nucleus segmentation
- **Points**: Transcript locations
- **Shapes**: Cell boundaries
- **ROIs layer**: Empty layer for drawing ROIs

**Instructions for napari:**
1. Select the "ROIs" shapes layer
2. Choose the polygon tool
3. Draw polygons around areas of interest
4. Draw multiple ROIs as needed
5. Close napari when finished

### Step 2: Cell Extraction

The script automatically:
- Finds cells whose positions fall within each ROI polygon
- Uses cell centroids from `cell_circles` shapes or table coordinates
- Reports the number of cells found in each ROI

### Step 3: UMAP Analysis

For each ROI with sufficient cells (>10), the script performs:
- **Normalization**: Scale to 10,000 counts per cell
- **Log transformation**: Log(counts + 1)
- **Feature selection**: Highly variable genes (if >50 cells)
- **PCA**: Principal component analysis
- **Neighborhood graph**: K-nearest neighbors
- **UMAP**: 2D embedding
- **Leiden clustering**: Community detection

### Step 4: Visualization

The script generates:
- **UMAP plots** colored by clusters, cell types, and UMI counts
- **Cluster composition** heatmaps showing cell type distributions
- **Comparative plots** across all ROIs
- **Statistical comparisons** between ROIs

## Output Files

The analysis creates several output files:

```
roi_analysis_results/
├── combined_roi_analysis.csv           # Combined data from all ROIs
├── roi_comparison.png                  # Comparative plots across ROIs
├── [ROI_NAME]_umap_analysis.png       # UMAP plots for each ROI
├── [ROI_NAME]_cluster_composition.png # Cell type composition
└── [ROI_NAME]_analysis.h5ad           # Single-cell data for each ROI
```

### Key Output Files

1. **`combined_roi_analysis.csv`**: Main results file containing:
   - ROI assignments for each cell
   - UMAP coordinates
   - Cluster assignments
   - Cell type annotations
   - Gene expression metrics

2. **UMAP plots**: Three-panel plots showing:
   - Leiden clusters
   - Cell type annotations
   - Total UMI counts

3. **Cluster composition**: Heatmaps showing:
   - Cell type counts per cluster
   - Percentages by cluster

4. **ROI comparison**: Multi-panel plots showing:
   - Cell counts per ROI
   - Gene expression distributions
   - UMI count distributions
   - Cell type proportions

## Data Analysis Examples

### Load and Explore Results

```python
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load combined results
df = pd.read_csv('roi_analysis_results/combined_roi_analysis.csv')

# Basic statistics
print(f"Total cells: {len(df)}")
print(f"ROIs: {df['roi'].nunique()}")
print(f"Clusters: {df['leiden_cluster'].nunique()}")
print(f"Cell types: {df['cell_type'].nunique()}")

# Cell counts per ROI
print("\nCells per ROI:")
print(df['roi'].value_counts())

# Top cell types
print("\nTop cell types:")
print(df['cell_type'].value_counts().head(10))
```

### Visualize UMAP Results

```python
# Plot UMAP colored by ROI
plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)
for roi in df['roi'].unique():
    roi_data = df[df['roi'] == roi]
    plt.scatter(roi_data['umap_1'], roi_data['umap_2'], 
               label=roi, alpha=0.6, s=10)
plt.xlabel('UMAP 1')
plt.ylabel('UMAP 2')
plt.title('UMAP Colored by ROI')
plt.legend()

plt.subplot(1, 2, 2)
for cluster in df['leiden_cluster'].unique():
    cluster_data = df[df['leiden_cluster'] == cluster]
    plt.scatter(cluster_data['umap_1'], cluster_data['umap_2'], 
               label=f'Cluster {cluster}', alpha=0.6, s=10)
plt.xlabel('UMAP 1')
plt.ylabel('UMAP 2')
plt.title('UMAP Colored by Cluster')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout()
plt.show()
```

### Statistical Analysis

```python
from scipy import stats

# Compare gene expression between ROIs
roi1_data = df[df['roi'] == 'ROI_1']
roi2_data = df[df['roi'] == 'ROI_2']

# T-test for gene expression
t_stat, p_val = stats.ttest_ind(roi1_data['n_genes'], roi2_data['n_genes'])
print(f"Gene expression comparison: t={t_stat:.3f}, p={p_val:.6f}")

# Summary statistics by ROI
summary = df.groupby('roi').agg({
    'n_genes': ['mean', 'std'],
    'total_counts': ['mean', 'std'],
    'leiden_cluster': 'nunique'
}).round(2)
print("\nSummary by ROI:")
print(summary)
```

## Advanced Usage

### Custom ROI Shapes

You can create complex ROI shapes programmatically:

```python
from shapely.geometry import Polygon
import json

# Create custom ROI shapes
custom_rois = {
    "Central_Region": {
        "coordinates": [
            [3000, 1000], [5000, 1000], [6000, 2000], 
            [5000, 3000], [3000, 3000], [2000, 2000]
        ],
        "area": 8000000,
        "bounds": [2000, 1000, 6000, 3000]
    }
}

# Save and use
with open('custom_rois.json', 'w') as f:
    json.dump(custom_rois, f, indent=2)
```

### Integration with Jupyter Notebooks

The `notebooks/roi_umap_workflow.ipynb` notebook provides an interactive version of this workflow with:
- Step-by-step explanations
- Interactive visualizations
- Statistical analysis examples
- Export capabilities

### Command Line Options

```bash
# Full command line options
python roi_umap_analysis.py --help

# Key options:
--roi-file ROI_FILE          # Load ROIs from JSON file
--output-dir OUTPUT_DIR      # Output directory for results
--save-rois SAVE_ROIS        # Save defined ROIs to JSON file
--skip-napari               # Skip napari (requires --roi-file)
```

## Troubleshooting

### Common Issues

1. **No cells found in ROIs**: 
   - Check ROI coordinates are within data bounds
   - Verify spatial data has position information
   - Ensure ROIs are large enough

2. **Too few cells for analysis**:
   - Increase ROI size
   - Lower the minimum cell threshold (edit script)
   - Check data quality

3. **napari display issues**:
   - Ensure all dependencies are installed
   - Check coordinate systems match
   - Verify data is properly loaded

### Memory Considerations

For large datasets:
- Use smaller ROIs to reduce memory usage
- Process ROIs individually if needed
- Consider downsampling for visualization

## Integration with Other Tools

### Export to Seurat (R)

```python
import pandas as pd
import numpy as np

# Load ROI results
df = pd.read_csv('roi_analysis_results/combined_roi_analysis.csv')

# Create metadata for Seurat
metadata = df[['cell_id', 'roi', 'leiden_cluster', 'cell_type']].copy()
metadata.columns = ['cell_id', 'ROI', 'Cluster', 'CellType']

# Save for Seurat
metadata.to_csv('roi_metadata_for_seurat.csv', index=False)
```

### Export to CellxGene

```python
import scanpy as sc

# Load individual ROI analysis
adata = sc.read('roi_analysis_results/ROI_1_analysis.h5ad')

# Add metadata
adata.obs['roi'] = 'ROI_1'
adata.obs['analysis_type'] = 'ROI_based'

# Save for CellxGene
adata.write('roi_for_cellxgene.h5ad')
```

## Best Practices

1. **ROI Definition**:
   - Draw ROIs on areas with good data quality
   - Include sufficient cells (>100 recommended)
   - Consider biological relevance of regions

2. **Analysis Parameters**:
   - Adjust clustering resolution based on expected cell types
   - Use appropriate normalization for your data type
   - Consider batch effects between ROIs

3. **Interpretation**:
   - Compare results across ROIs carefully
   - Validate findings with known biology
   - Consider spatial context in interpretation

4. **Documentation**:
   - Save ROI definitions for reproducibility
   - Document analysis parameters
   - Keep track of biological hypotheses

This workflow provides a comprehensive approach to spatial ROI analysis with UMAP clustering, enabling detailed exploration of tissue heterogeneity and spatial organization.