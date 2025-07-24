# napari-roi-manager Integration Guide

## Problem Solved ✅

Your napari-roi-manager ROIs had coordinate values in the 70,000+ range, while your spatial data coordinates are in the 4,000-6,000 range. This happened because:

1. **ROIs were drawn at full resolution** (scale0) in napari-roi-manager
2. **Spatial analysis uses global coordinate system** which corresponds to a downsampled level
3. **Scale factor of ~18x** needed to transform between coordinate systems

## Solution

### Step 1: Coordinate Transformation

I created `transform_napari_rois.py` which:
- Detects the coordinate system mismatch
- Applies appropriate scaling transformation  
- Converts to format compatible with `roi_umap_analysis.py`

### Step 2: Successful Integration

```bash
# Transform your napari-roi-manager ROIs
python transform_napari_rois.py rois.json --scale-factor 17.9 --output transformed_napari_rois.json

# Use transformed ROIs for analysis
python roi_umap_analysis.py combined_direct_coords_annotated.zarr --roi-file transformed_napari_rois.json
```

## Results ✅

Your napari-roi-manager ROIs are now successfully integrated:

### Original ROI Coordinates (napari-roi-manager):
```json
{
  "data": [
    [[73024.18, 22214.96], [71611.60, 22995.24], ...],  // ROI "a"
    [[74625.11, 26614.14], [73414.32, 27555.87], ...]   // ROI "b"
  ],
  "names": ["a", "b"]
}
```

### Transformed ROI Coordinates (analysis-ready):
```json
{
  "a": {
    "coordinates": [[4079.56, 1241.06], [4000.65, 1284.65], ...],
    "area": 18389.11,
    "bounds": [4000.65, 1241.06, 4167.50, 1396.63]
  },
  "b": {
    "coordinates": [[4169.00, 1486.82], [4101.36, 1539.43], ...],
    "area": 16949.89,
    "bounds": [4101.36, 1454.51, 4254.68, 1640.90]
  }
}
```

### Cell Extraction Results:
- **ROI "a"**: 73 cells identified
- **ROI "b"**: 59 cells identified
- **Total**: 132 cells across both ROIs

## Why This Happened

### napari-roi-manager Coordinate System
- napari-roi-manager saves coordinates in the **native resolution** where ROIs were drawn
- If you drew on a high-resolution image (scale0), coordinates are at full pixel resolution
- These coordinates are **valid but at wrong scale** for analysis

### SpatialData Global Coordinates  
- Your spatial analysis uses the **"global" coordinate system**
- This typically corresponds to a **downsampled resolution** for computational efficiency
- Scale factor depends on the multiscale pyramid levels

### Multi-scale Image Structure
Your spatial data has 5 resolution levels:
- `scale0`: Full resolution (where ROIs were drawn)
- `scale1`, `scale2`, `scale3`, `scale4`: Progressively downsampled
- Analysis uses global coordinates ≈ scale4 level

## Files Created

1. **`transform_napari_rois.py`** - Coordinate transformation tool
2. **`transformed_napari_rois.json`** - Your ROIs in analysis-ready format
3. **`convert_napari_roi_manager.py`** - Alternative conversion approach (backup)

## Usage Workflow

### For Future napari-roi-manager ROIs:

```bash
# 1. Create ROIs in napari-roi-manager (as you did)
# 2. Transform coordinates to analysis coordinate system
python transform_napari_rois.py your_rois.json --scale-factor 17.9 --output analysis_ready_rois.json

# 3. Run ROI analysis
python roi_umap_analysis.py your_spatial_data.zarr --roi-file analysis_ready_rois.json

# 4. Results will be generated in roi_analysis_results/
```

### Scale Factor Notes:
- **Scale factor 17.9** worked for your data
- This may vary depending on:
  - Which resolution level ROIs were drawn at
  - Which coordinate system the analysis uses
  - The specific multiscale pyramid structure

### Auto-detection:
```bash
# Let the script auto-detect scale factor
python transform_napari_rois.py rois.json --output transformed.json

# Or specify manually if you know it
python transform_napari_rois.py rois.json --scale-factor 18.0 --output transformed.json
```

## Benefits of This Integration

✅ **Use familiar napari-roi-manager interface** for ROI drawing
✅ **Leverage ImageJ-like ROI management** features  
✅ **Seamless integration** with spatial analysis pipeline
✅ **Preserve ROI metadata** and naming
✅ **Reproducible coordinate transformation**

## Troubleshooting

### If ROIs seem misplaced:
1. Check which resolution level you drew ROIs at
2. Adjust scale factor accordingly
3. Verify spatial data coordinate bounds

### If no cells found in ROIs:
1. Visualize transformed coordinates to check placement
2. Try different scale factors (16, 18, 20)
3. Check spatial data coordinate system

### Common Scale Factors:
- **2x, 4x, 8x, 16x**: Standard multiscale factors
- **~18x**: Your specific case (17.9x)
- **Custom**: May need manual calculation

This integration successfully bridges napari-roi-manager's powerful ROI creation tools with your spatial transcriptomics analysis pipeline! 🎉