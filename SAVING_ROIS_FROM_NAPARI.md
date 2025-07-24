# How to Save ROI JSON from napari

This guide shows you multiple ways to create and save Region of Interest (ROI) data from napari in JSON format for use with the ROI analysis pipeline.

## Quick Start

### Method 1: Use the Simple ROI Editor (Recommended)

```bash
# Open napari to draw ROIs and save to JSON
python simple_roi_editor.py combined_direct_coords_annotated.zarr my_custom_rois.json

# View existing ROIs
python simple_roi_editor.py --view my_custom_rois.json

# Use the ROIs for analysis
python roi_umap_analysis.py combined_direct_coords_annotated.zarr --roi-file my_custom_rois.json
```

### Method 2: Use the Full ROI Analysis Tool

```bash
# Define ROIs in napari and automatically run analysis
python roi_umap_analysis.py combined_direct_coords_annotated.zarr --save-rois my_rois.json
```

### Method 3: Create ROIs Programmatically

```bash
# Create example ROIs without napari
python save_rois_from_napari.py --programmatic my_programmatic_rois.json
```

## Detailed Methods

### Method 1: Simple ROI Editor (Easiest)

The `simple_roi_editor.py` script provides a streamlined napari interface:

**Steps:**
1. Run the script with your zarr file
2. napari opens with your spatial data loaded
3. Select the "My_ROIs" layer
4. Draw polygons around regions of interest
5. Close napari to automatically save JSON

**napari Drawing Instructions:**
- **Select tool**: Make sure "polygon" tool is selected in the shapes layer
- **Draw polygons**: Click to place vertices, double-click to close
- **Navigate**: Mouse wheel to zoom, middle button to pan
- **Precision**: Hold Shift while drawing for more control

**Example:**
```bash
python simple_roi_editor.py combined_direct_coords_annotated.zarr my_rois.json
```

### Method 2: Programmatic ROI Creation

Create ROIs directly in code without napari:

```python
import json

# Define ROI coordinates (x, y format)
roi_data = {
    "Tissue_Region_1": {
        "coordinates": [
            [1000, 500],   # x, y coordinates
            [3000, 500],
            [3000, 2000],
            [1000, 2000]
        ],
        "area": 3000000,
        "bounds": [1000, 500, 3000, 2000],
        "description": "Upper tissue region"
    },
    "Tissue_Region_2": {
        "coordinates": [
            [2000, 2500],
            [4000, 2500],
            [4000, 4000],
            [2000, 4000]
        ],
        "area": 3000000,
        "bounds": [2000, 2500, 4000, 4000],
        "description": "Lower tissue region"
    }
}

# Save to JSON
with open('my_custom_rois.json', 'w') as f:
    json.dump(roi_data, f, indent=2)
```

### Method 3: Advanced napari Integration

For more control over the napari session:

```python
import napari
import spatialdata as sd
import json
from shapely.geometry import Polygon

def create_custom_rois(zarr_path, output_json):
    # Load data
    sdata = sd.read_zarr(zarr_path)
    
    # Create viewer with custom settings
    viewer = napari.Viewer()
    
    # Add your data layers
    # ... (add images, points, etc.)
    
    # Add ROI layer
    roi_layer = viewer.add_shapes(name='ROIs', shape_type='polygon')
    
    # Run napari
    napari.run()
    
    # Extract and save ROIs
    roi_data = {}
    for i, coords in enumerate(roi_layer.data):
        coords_xy = coords[:, [1, 0]]  # Convert y,x to x,y
        polygon = Polygon(coords_xy)
        
        roi_data[f"Custom_ROI_{i+1}"] = {
            'coordinates': coords_xy.tolist(),
            'area': polygon.area,
            'bounds': list(polygon.bounds)
        }
    
    # Save JSON
    with open(output_json, 'w') as f:
        json.dump(roi_data, f, indent=2)
```

## ROI JSON Format

The ROI JSON files follow this structure:

```json
{
  "ROI_Name": {
    "coordinates": [
      [x1, y1],
      [x2, y2],
      [x3, y3],
      [x4, y4]
    ],
    "area": 1000000.0,
    "bounds": [min_x, min_y, max_x, max_y],
    "description": "Optional description"
  }
}
```

### Required Fields:
- **`coordinates`**: List of [x, y] coordinate pairs defining the polygon
- **`area`**: Area of the polygon in coordinate units
- **`bounds`**: Bounding box as [min_x, min_y, max_x, max_y]

### Optional Fields:
- **`description`**: Text description of the ROI
- **`perimeter`**: Perimeter length
- **`created_with`**: Tool used to create the ROI
- **`created_method`**: Method used (e.g., "napari_interactive", "programmatic")

## Working with Coordinate Systems

### Understanding Coordinates

napari uses **Y, X** coordinate order while most spatial data uses **X, Y**. The scripts handle this conversion automatically:

```python
# napari format (Y, X)
napari_coords = [[y1, x1], [y2, x2], [y3, x3]]

# Standard format (X, Y) - used in JSON
standard_coords = [[x1, y1], [x2, y2], [x3, y3]]
```

### Coordinate Validation

Your ROI coordinates should be within the bounds of your spatial data:

```python
# Check if coordinates are within data bounds
def validate_coordinates(coords, spatial_data):
    # Get data bounds
    min_x, min_y, max_x, max_y = get_spatial_bounds(spatial_data)
    
    for x, y in coords:
        if not (min_x <= x <= max_x and min_y <= y <= max_y):
            print(f"Warning: Coordinate ({x}, {y}) is outside data bounds")
```

## Tips for Drawing Good ROIs

### 1. Choose Biologically Meaningful Regions
- **Tissue boundaries**: Draw around different tissue types
- **Structural features**: Include/exclude specific anatomical structures
- **Gradients**: Capture regions with different expression gradients

### 2. Size Considerations
- **Minimum cells**: Ensure ROIs contain enough cells for analysis (>50 recommended)
- **Maximum size**: Very large ROIs may be less specific
- **Multiple small ROIs**: Better than one huge ROI for comparative analysis

### 3. Drawing Precision
- **Zoom in**: Use napari's zoom for precise boundary drawing
- **Follow boundaries**: Align with visible tissue structures
- **Avoid artifacts**: Don't include areas with poor data quality

### 4. Consistent Strategy
- **Same criteria**: Use consistent criteria across all ROIs
- **Document decisions**: Keep notes about why you drew ROIs where you did
- **Save intermediate**: Save ROIs frequently during drawing sessions

## Troubleshooting

### Common Issues

1. **napari doesn't open**:
   ```bash
   # Check if napari is installed
   uv run python -c "import napari; print('napari installed')"
   
   # Try with display backend
   export QT_QPA_PLATFORM=offscreen  # For headless systems
   ```

2. **No ROIs saved**:
   - Make sure you drew polygons (not just clicked points)
   - Check that you double-clicked to close polygons
   - Verify the shapes layer was selected when drawing

3. **Coordinates seem wrong**:
   - Check coordinate system (napari uses Y,X, JSON uses X,Y)
   - Verify data bounds match expected tissue locations
   - Use visualization tools to check ROI placement

4. **JSON file invalid**:
   ```bash
   # Validate JSON format
   python save_rois_from_napari.py --validate my_rois.json
   ```

### Performance Tips

1. **Large datasets**: Sample points/shapes for napari visualization
2. **Memory issues**: Close other applications before opening napari
3. **Slow rendering**: Reduce opacity of point layers

## Integration with Analysis Pipeline

### Using ROIs for Analysis

Once you have your ROI JSON file:

```bash
# Run full analysis
python roi_umap_analysis.py your_data.zarr --roi-file your_rois.json --output-dir results/

# Results will include:
# - UMAP plots for each ROI
# - Cell type composition analysis
# - Cross-ROI comparisons
# - Statistical analysis
```

### Modifying Existing ROIs

You can edit the JSON file directly:

```python
import json

# Load existing ROIs
with open('my_rois.json', 'r') as f:
    rois = json.load(f)

# Add a new ROI
rois['New_Region'] = {
    'coordinates': [[1000, 1000], [2000, 1000], [2000, 2000], [1000, 2000]],
    'area': 1000000,
    'bounds': [1000, 1000, 2000, 2000],
    'description': 'Manually added region'
}

# Save updated ROIs
with open('my_rois_updated.json', 'w') as f:
    json.dump(rois, f, indent=2)
```

### Converting Between Formats

If you have ROIs in other formats:

```python
# From ImageJ ROI files
def convert_imagej_roi_to_json(roi_file, output_json):
    # Use roifile library to read ImageJ ROIs
    from roifile import ImagejRoi
    
    roi = ImagejRoi.fromfile(roi_file)
    coords = roi.coordinates()
    
    roi_data = {
        'ImageJ_ROI': {
            'coordinates': coords.tolist(),
            'area': calculate_area(coords),
            'bounds': calculate_bounds(coords)
        }
    }
    
    with open(output_json, 'w') as f:
        json.dump(roi_data, f, indent=2)

# From QuPath annotations
def convert_qupath_to_json(qupath_geojson, output_json):
    # Convert QuPath GeoJSON to ROI JSON format
    # ... implementation depends on QuPath format
    pass
```

## Example Workflows

### Workflow 1: Comparative Analysis
1. Draw ROIs around different tissue regions
2. Run analysis to compare cell populations
3. Identify region-specific cell types or states

### Workflow 2: Time Series Analysis
1. Draw ROIs at consistent locations across time points
2. Analyze changes in cell composition over time
3. Track cellular dynamics within regions

### Workflow 3: Disease vs Control
1. Draw ROIs in comparable regions across samples
2. Compare healthy vs diseased tissue regions
3. Identify disease-associated cellular changes

This comprehensive guide should help you successfully create, save, and use ROI JSON data from napari for your spatial transcriptomics analysis!