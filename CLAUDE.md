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

## Testing Files

### Standard Test Files for Development
When testing zarr combination and spatial data processing tools, always use these standardized test files:
- **File 1**: `/Users/chrism/datasets/lv_spatialdat_liu/lv_0046706_007.zarr`
- **File 2**: `/Users/chrism/datasets/lv_spatialdat_liu/lv_0046706_117.zarr`

These files are known to have consistent structure and are representative of typical SpatialData zarr files used in this project.