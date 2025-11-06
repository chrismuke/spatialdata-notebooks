# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

This is the spatialdata-notebooks repository containing tutorials and examples for SpatialData, a framework for spatial omics data. The repository consists primarily of Jupyter notebooks demonstrating various spatial data analysis workflows and techniques.

## Development Environment

### Package Management
- **Primary tool**: `uv` (ultra-fast Python package manager)
- **Lock file**: `uv.lock` (tracks exact dependency versions)
- **Configuration**: `pyproject.toml` (project metadata and dependencies)

### Common Commands

#### Environment Setup
```bash
# Install dependencies using uv
uv sync

# Install with optional dev dependencies
uv sync --extra dev

# Install with documentation dependencies
uv sync --extra doc
```

#### Development Tasks
```bash
# Run notebooks using uv
uv run jupyter lab

# Run individual notebooks
uv run jupyter notebook notebooks/examples/intro.ipynb

# Format code using ruff (via pre-commit)
uv run ruff format .

# Lint code using ruff
uv run ruff check . --fix

# Run pre-commit hooks manually
uv run pre-commit run --all-files

# Run custom Python scripts and CLI tools
uv run python script_name.py
uv run python celltype_annotate_cli_v2.py --help

# Example: Cell type annotation with organized output
uv run python celltype_annotate_cli_v2.py xenium_data.zarr reference_data.h5ad --min-clusters 5 --max-clusters 10

# Example: Cell type annotation with custom results directory
uv run python celltype_annotate_cli_v2.py xenium_data.zarr reference_data.h5ad --results-dir /path/to/my/results

# Example: Cell type annotation showing all cells (including unknown types)
uv run python celltype_annotate_cli_v2.py xenium_data.zarr reference_data.h5ad --show-unknown-cells

# Example: Cell type annotation using existing trained model
uv run python celltype_annotate_cli_v2.py xenium_data.zarr reference_data.h5ad --load-model reference_data_cell_type_model_20250724_120110.pkl

# Example: ROI analysis with custom output directory
uv run python roi_umap_analysis.py combined_data.zarr --output-dir /path/to/roi_results

# Example: Combine multiple zarr files with spatial translation
uv run python combine_zarr_cli.py file1.zarr file2.zarr file3.zarr combined.zarr --border 100 --layout vertical
```

#### Documentation Building
```bash
# Build documentation using Sphinx
make html

# Clean build artifacts
make clean

# Alternative Sphinx commands
python3 -msphinx -M html . _build
```

## Code Architecture

### Repository Structure
- **`notebooks/`**: Main content directory
  - **`examples/`**: Core tutorial notebooks for different spatial data technologies and analysis workflows
  - **`developers_resources/`**: Technical notebooks for storage format demonstrations
  - **`paper_reproducibility/`**: Research paper reproduction notebooks
- **`datasets/`**: Sample data in Zarr format (mouse liver dataset)
- **`_static/`**: Static assets for documentation (images, CSS)
- **`results/`**: Generated analysis outputs and HTML reports

### Key Notebook Categories
1. **Technology-specific tutorials**: Visium, Xenium, CosMx, MERFISH, etc. (`technology_*.ipynb`)
2. **Analysis workflows**: Aggregation, spatial queries, transformations
3. **Integration examples**: Squidpy integration, model applications
4. **Advanced topics**: Landmark-based alignment, ROI analysis

### Python Utilities
- **`notebooks/examples/generate_toc.py`**: Generates table of contents for notebooks
- **`notebooks/developers_resources/storage_format/io_utils.py`**: Utilities for Zarr data I/O and consistency checking

## Development Workflow

### Code Quality
- **Linting**: Ruff configured in `pyproject.toml` (line length: 120)
- **Formatting**: Ruff formatter with pre-commit hooks
- **Pre-commit**: Configured in `.pre-commit-config.yaml` with ruff, file checks, and security scanning

### Key Dependencies
- **Core**: `spatialdata`, `spatialdata-plot`, `anndata`, `scanpy`
- **Analysis**: `squidpy`, `sopa`, `celltypist`
- **Notebooks**: `jupyterlab`, `jupyterlab-vim`, `notebook`
- **Visualization**: Standard scientific Python stack

### Testing and Validation
- Notebooks should run without errors
- Data consistency checks for Zarr stores using utilities in `io_utils.py`
- Pre-commit hooks validate Python syntax and file formatting

## Important Notes

### Data Handling
- Large datasets stored as Zarr stores in `datasets/` and generated in notebooks
- Use `io_utils.py` functions for consistent Zarr I/O operations
- Sample data available at specific S3 URLs for examples

### Documentation Integration
- Built with Sphinx using `sphinx-book-theme`
- MyST-NB for notebook integration
- Notebooks excluded from certain builds via `conf.py` patterns

### Environment Considerations
- Requires Python >=3.10
- Uses uv for fast dependency resolution
- Optional dependencies for development and documentation building

### Results Directory Management
Analysis tools create large results folders that can quickly consume disk space:

- **Cell type annotation**: `celltype_annotate_cli_v2.py` supports `--results-dir` to specify custom base directory
  - Creates folders with format: `xenium_filename___reference_filename___[clusters]___timestamp`
  - Example: `combined_data___mouse_ref___clusters_min5_max10___20250724_120110/`
  - Use `--consolidate-data` for more information about data portability
- **ROI analysis**: `roi_umap_analysis.py` supports `--output-dir` to specify output location
- **napari ROI conversion**: `convert_napari_roi_manager.py` supports `--results-dir` for analysis results

### Understanding SpatialData Warnings
You may see these INFO messages during analysis:

```
INFO: The SpatialData object is not self-contained...
INFO: The Zarr backing store has been changed from [original] to [new location]...
```

**What this means:**
- **Not self-contained**: Results reference data in multiple locations (efficient but less portable)
- **Backing store changed**: Main metadata moved to results folder, but large data (images) stay at original location
- **This is normal behavior** - SpatialData avoids copying large files unnecessarily

**Implications:**
- ✅ **Analysis works fine** - all visualizations and results are correct
- ✅ **Storage efficient** - avoids duplicating large image files
- ⚠️ **Less portable** - results depend on original data location
- ⚠️ **Broken if original files move** - results become unusable if source data is relocated

### Cell Type Visualization Options
By default, spatial cell type maps hide cells with unknown/unassigned cell types for cleaner visualization:

**Default behavior (hide unknown cells):**
```bash
uv run python celltype_annotate_cli_v2.py data.zarr ref.h5ad
# Creates clean spatial plots showing only confidently identified cell types
```

**Show all cells including unknown types:**
```bash  
uv run python celltype_annotate_cli_v2.py data.zarr ref.h5ad --show-unknown-cells
# Includes cells labeled as "unknown", "unassigned", etc. in spatial plots
```

**Unknown cell patterns filtered by default:**
- "unknown", "unassigned", "unlabeled"
- "ambiguous", "unclear", "na", "none", "nan"
- Empty strings and NaN values

Examples:
```bash
# Store results in external drive
uv run python celltype_annotate_cli_v2.py data.zarr ref.h5ad --results-dir /Volumes/ExternalDrive/analysis_results

# Store ROI analysis in project-specific directory
uv run python roi_umap_analysis.py data.zarr --output-dir ./project_analysis/roi_results
```

## Testing Files

### Standard Test Files for Development
When testing zarr combination and spatial data processing tools, always use these standardized test files:
- **File 1**: `/Users/chrism/datasets/lv_spatialdat_liu/lv_0046706_007.zarr`
- **File 2**: `/Users/chrism/datasets/lv_spatialdat_liu/lv_0046706_117.zarr`

These files are known to have consistent structure and are representative of typical SpatialData zarr files used in this project.