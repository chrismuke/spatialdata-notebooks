#!/usr/bin/env python3
"""
Interactive ROI Workflow Demo

This script demonstrates how to:
1. Define ROIs in napari
2. Extract cells within ROIs
3. Perform UMAP analysis on ROI-selected cells
4. Visualize the results

Run this script to see the workflow in action.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from roi_umap_analysis import ROIAnalyzer
import warnings
warnings.filterwarnings('ignore')

def demo_roi_workflow():
    """Demonstrate the ROI workflow"""
    
    print("=== ROI Analysis Workflow Demo ===\n")
    
    # Step 1: Initialize the analyzer
    print("Step 1: Loading spatial data...")
    zarr_path = "combined_direct_coords_annotated.zarr"
    analyzer = ROIAnalyzer(zarr_path)
    
    print(f"✓ Loaded spatial data with {analyzer.table.n_obs} cells and {analyzer.table.n_vars} genes")
    
    # Step 2: Check if we have existing ROIs or need to define new ones
    roi_file = "demo_rois.json"
    
    try:
        print("\nStep 2: Loading existing ROIs...")
        analyzer.load_roi_data(roi_file)
        print(f"✓ Loaded {len(analyzer.roi_data)} existing ROIs")
    except FileNotFoundError:
        print("\nStep 2: No existing ROIs found. Instructions for defining ROIs:")
        print("""
        To define ROIs:
        1. Run: python roi_umap_analysis.py combined_direct_coords_annotated.zarr --save-rois demo_rois.json
        2. In napari:
           - Use the 'ROIs' shapes layer to draw polygons
           - Draw around areas of interest (different tissue regions)
           - You can draw multiple ROIs
           - Close napari when done
        3. Re-run this demo script
        """)
        return
    
    # Step 3: Extract cells within ROIs
    print("\nStep 3: Extracting cells within ROIs...")
    roi_cells = analyzer.extract_cells_in_rois()
    
    if not roi_cells:
        print("❌ No cells found in ROIs")
        return
    
    # Step 4: Perform UMAP analysis
    print("\nStep 4: Performing UMAP analysis on ROI-selected cells...")
    results = analyzer.perform_umap_analysis(roi_cells, "demo_results")
    
    if not results:
        print("❌ No analysis results generated")
        return
    
    # Step 5: Create combined analysis
    print("\nStep 5: Creating combined analysis...")
    combined_df = analyzer.create_combined_analysis(results, "demo_results")
    
    # Step 6: Display results
    print("\nStep 6: Analysis Results Summary")
    print("=" * 50)
    
    for roi_name, roi_result in results.items():
        adata = roi_result['adata']
        n_cells = roi_result['n_cells']
        n_clusters = len(adata.obs['leiden'].unique())
        
        print(f"\n{roi_name}:")
        print(f"  - Cells: {n_cells}")
        print(f"  - Clusters: {n_clusters}")
        
        if 'cell_type_predicted' in adata.obs.columns:
            top_celltypes = adata.obs['cell_type_predicted'].value_counts().head(3)
            print(f"  - Top cell types:")
            for celltype, count in top_celltypes.items():
                pct = (count / n_cells) * 100
                print(f"    • {celltype}: {count} cells ({pct:.1f}%)")
    
    # Step 7: Show file outputs
    print(f"\nStep 7: Output Files Generated")
    print("=" * 50)
    print("📁 demo_results/")
    print("  ├── combined_roi_analysis.csv       # Combined data from all ROIs")
    print("  ├── roi_comparison.png              # Comparative plots across ROIs")
    
    for roi_name in results.keys():
        print(f"  ├── {roi_name}_umap_analysis.png      # UMAP plots for {roi_name}")
        print(f"  └── {roi_name}_cluster_composition.png # Cell type composition")
    
    print(f"\n✅ ROI analysis complete! Check the demo_results/ directory for visualizations.")
    
    # Step 8: Quick data preview
    print(f"\nStep 8: Quick Data Preview")
    print("=" * 50)
    
    if not combined_df.empty:
        print("Combined ROI Analysis Summary:")
        print(f"Total cells analyzed: {len(combined_df)}")
        print(f"ROIs analyzed: {combined_df['roi'].nunique()}")
        print(f"Total clusters found: {combined_df['leiden_cluster'].nunique()}")
        
        if 'cell_type' in combined_df.columns:
            print(f"Cell types identified: {combined_df['cell_type'].nunique()}")
            print("\nTop cell types across all ROIs:")
            top_types = combined_df['cell_type'].value_counts().head(5)
            for celltype, count in top_types.items():
                print(f"  • {celltype}: {count} cells")


def create_example_roi_file():
    """Create an example ROI file for demonstration"""
    
    print("Creating example ROI file...")
    
    # Define some example ROIs (you would normally get these from napari)
    example_rois = {
        "ROI_1": {
            "coordinates": [
                [2000, 500],
                [4000, 500], 
                [4000, 1500],
                [2000, 1500]
            ],
            "area": 3000000,
            "bounds": [2000, 500, 4000, 1500]
        },
        "ROI_2": {
            "coordinates": [
                [1000, 3000],
                [3000, 3000],
                [3000, 4500], 
                [1000, 4500]
            ],
            "area": 3000000,
            "bounds": [1000, 3000, 3000, 4500]
        }
    }
    
    import json
    with open("demo_rois.json", "w") as f:
        json.dump(example_rois, f, indent=2)
    
    print("✓ Created demo_rois.json with example ROIs")


if __name__ == "__main__":
    print("Choose an option:")
    print("1. Run ROI workflow demo (requires existing ROIs)")
    print("2. Create example ROI file")
    print("3. Instructions for defining ROIs in napari")
    
    choice = input("\nEnter choice (1/2/3): ").strip()
    
    if choice == "1":
        demo_roi_workflow()
    elif choice == "2":
        create_example_roi_file()
        print("\nNow you can run option 1 to see the workflow in action!")
    elif choice == "3":
        print("""
ROI Definition Instructions:
============================

1. Run the following command:
   python roi_umap_analysis.py combined_direct_coords_annotated.zarr --save-rois my_rois.json

2. Napari will open with your spatial data loaded

3. In napari:
   - You'll see an empty 'ROIs' shapes layer
   - Select the shapes layer
   - Choose 'polygon' tool
   - Draw polygons around areas of interest
   - You can draw multiple ROIs
   - Each ROI will be named automatically (ROI_1, ROI_2, etc.)

4. Close napari when finished

5. Your ROIs will be saved to my_rois.json

6. Run the analysis:
   python roi_umap_analysis.py combined_direct_coords_annotated.zarr --roi-file my_rois.json

Alternative: Use the demo workflow
   python roi_workflow_demo.py
        """)
    else:
        print("Invalid choice. Please run the script again.")