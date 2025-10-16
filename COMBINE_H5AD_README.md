# Combine H5AD Files CLI Tool

Command-line tool for combining multiple h5ad single-cell RNA-seq files with intelligent duplicate handling and conflict resolution.

## Features

- ✅ **Combine 2+ h5ad files** into a single dataset
- ✅ **Automatic gene harmonization** (uses common genes across all files)
- ✅ **Duplicate cell detection** across files
- ✅ **Cell type conflict resolution** with multiple strategies:
  - Majority vote (when 2+ files agree)
  - Random selection (for ties)
  - Exclusion (optional: remove conflicting cells)
- ✅ **Source tracking** (records which file each cell came from)
- ✅ **Detailed logging** with warnings for conflicts

## Installation

The tool is already part of this repository and uses `uv` for dependency management:

```bash
cd /Users/chrism/git/spatialdata-notebooks
uv sync
```

## Usage

### Basic Usage

Combine two files:
```bash
uv run python combine_h5ad_cli.py file1.h5ad file2.h5ad output.h5ad
```

Combine three or more files:
```bash
uv run python combine_h5ad_cli.py file1.h5ad file2.h5ad file3.h5ad output.h5ad
```

### Options

- `--exclude-conflicts`: Remove cells with cell type conflicts instead of resolving them
- `--cell-type-column TEXT`: Specify cell type column name (default: auto-detect)
- `--overwrite`: Overwrite output file if it exists
- `--log-level [DEBUG|INFO|WARNING|ERROR]`: Set logging verbosity

### Examples

**Example 1: Combine colon and skin datasets**
```bash
uv run python combine_h5ad_cli.py \
  /Users/chrism/datasets/scrnaseq_human_colon/d3aaede6-87e6-4d18-9aa1-d757baf5d4d4.h5ad \
  /Users/chrism/datasets/scrnaseq_human_skin/d7a7bb0b-472d-4d95-a300-7271cea56e66.h5ad \
  combined_colon_skin.h5ad
```

**Example 2: Exclude cells with conflicts**
```bash
uv run python combine_h5ad_cli.py \
  file1.h5ad file2.h5ad file3.h5ad \
  combined.h5ad \
  --exclude-conflicts
```

**Example 3: Use specific cell type column**
```bash
uv run python combine_h5ad_cli.py \
  file1.h5ad file2.h5ad \
  combined.h5ad \
  --cell-type-column "Celltype"
```

## How It Works

### 1. Gene Harmonization

The tool identifies **common genes** across all input files:
- Only genes present in ALL files are retained
- Unique genes are excluded (with detailed reporting)
- All files are subset to the same gene set before combining

**Example output:**
```
Harmonizing genes...
  Common genes across all files: 25,384
  file1.h5ad: 1,905 unique genes (will be excluded)
  file2.h5ad: 3,037 unique genes (will be excluded)
```

### 2. Duplicate Detection

The tool checks for **duplicate cell IDs** across files:
- Cell IDs (obs_names) are compared across all input files
- Duplicates are identified and reported
- If no duplicates: files are simply concatenated

**Example output:**
```
Analyzing duplicate cells...
  Found 150 duplicate cell IDs
```

### 3. Conflict Resolution

For cells with the same ID but different cell types:

#### **Majority Vote (Default)**
- If 2+ files agree on a cell type → use that type
- Example: 3 files with types [A, A, B] → resolved to A

#### **Random Selection (Ties)**
- When files equally disagree → pick randomly
- Example: 2 files with types [A, B] → randomly choose A or B
- ⚠️ Warning is logged for random selections

#### **Exclusion (--exclude-conflicts)**
- Removes ALL cells with conflicts
- Ensures only unambiguous cells remain

**Example output:**
```
⚠️  Cell type conflicts detected for 75 cells
   - Resolved by majority vote: 60
   - Resolved randomly (tie): 15

Conflict examples:
  1. cell_12345
     Cell types: {'Type_A': 2, 'Type_B': 1}
     Resolved to: 'Type_A' (majority)
  2. cell_67890
     Cell types: {'Type_X': 1, 'Type_Y': 1}
     Resolved to: 'Type_X' (random)
```

### 4. Duplicate Removal

After conflict resolution:
- **First occurrence** of each duplicate is kept
- Later occurrences are removed
- Cell type is set to resolved type (if conflicts occurred)

### 5. Metadata Tracking

The combined file includes:
- `source_file`: Original filename for each cell
- `source_index`: File index (0, 1, 2, ...)
- `combine_metadata` in `.uns`:
  - List of source files
  - Combine date/time
  - Number of duplicates/conflicts
  - Conflict resolution settings

## Output Structure

### Combined AnnData Object

```python
import anndata as ad

# Load combined file
adata = ad.read_h5ad('combined.h5ad')

# Check dimensions
print(adata.shape)  # (total_cells, common_genes)

# Check source information
print(adata.obs['source_file'].value_counts())

# Check metadata
print(adata.uns['combine_metadata'])
```

### Example Output
```python
# Shape
(293527, 25384)  # 293,527 cells × 25,384 genes

# Source distribution
source_file
d7a7bb0b-472d-4d95-a300-7271cea56e66    195739
d3aaede6-87e6-4d18-9aa1-d757baf5d4d4     97788

# Metadata
{
  'source_files': ['/path/to/file1.h5ad', '/path/to/file2.h5ad'],
  'combine_date': '2025-10-16T14:24:19.342932',
  'num_files': 2,
  'exclude_conflicts': False,
  'num_duplicates': 0,
  'num_conflicts': 0
}
```

## Cell Type Column Auto-Detection

The tool automatically finds the cell type column by searching for:
1. `cell_type` (most common)
2. `celltype`
3. `Celltype`
4. `cell_label`
5. `annotation`
6. `cluster`
7. Any column containing "type" or "label"

If auto-detection fails, use `--cell-type-column` to specify explicitly.

## Use Cases

### 1. Combining Reference Atlases

Merge multiple tissue-specific reference datasets:
```bash
uv run python combine_h5ad_cli.py \
  colon_atlas.h5ad \
  skin_atlas.h5ad \
  lung_atlas.h5ad \
  multi_tissue_reference.h5ad
```

### 2. Batch Integration

Combine data from multiple experimental batches:
```bash
uv run python combine_h5ad_cli.py \
  batch1.h5ad batch2.h5ad batch3.h5ad \
  integrated.h5ad
```

### 3. Quality Control

Identify technical duplicates across datasets:
```bash
# This will warn about all duplicate cell IDs
uv run python combine_h5ad_cli.py \
  sample1.h5ad sample2.h5ad \
  combined.h5ad
```

### 4. Clean Merging

Exclude ambiguous cells for high-confidence analysis:
```bash
uv run python combine_h5ad_cli.py \
  dataset1.h5ad dataset2.h5ad \
  clean_combined.h5ad \
  --exclude-conflicts
```

## Warnings and Considerations

### ⚠️ Gene Loss
- Genes not present in ALL files will be excluded
- Check log output to see how many genes are lost per file
- Consider aligning gene IDs before combining if many genes are lost

### ⚠️ Duplicate Cell IDs
- Duplicate cell IDs across files may indicate:
  - Technical replicates
  - Data processing errors
  - ID collision from different sources
- Review warnings carefully to understand duplicates

### ⚠️ Cell Type Conflicts
- Conflicts may indicate:
  - Different annotation criteria between studies
  - Cell type ambiguity (transitional states)
  - Annotation errors
- `--exclude-conflicts` removes ~5-10% of data in typical merges

### ⚠️ Batch Effects
- This tool does NOT perform batch correction
- Combined data may have batch effects
- Consider using Scanpy's `sc.pp.combat()` or Harmony after combining

## Test Data

Successfully tested with:
- **Colon dataset**: 97,788 cells × 27,289 genes
- **Skin dataset**: 195,739 cells × 28,421 genes
- **Combined**: 293,527 cells × 25,384 genes (no duplicates)

## Troubleshooting

### Error: "No cell type column found"
**Solution**: Use `--cell-type-column` to specify the correct column name

### Error: "Output file exists"
**Solution**: Use `--overwrite` flag or delete the existing file

### Many genes excluded
**Check**: Are gene IDs consistent across files? (e.g., ENSEMBL vs. gene symbols)
**Solution**: Standardize gene IDs before combining

### High conflict rate
**Check**: Are datasets annotated with same cell type ontology?
**Solution**: Consider harmonizing annotations first or use `--exclude-conflicts`

## Advanced Usage

### Programmatic Access

```python
from pathlib import Path
from combine_h5ad_cli import combine_h5ad_files

combine_h5ad_files(
    file_paths=[
        Path('file1.h5ad'),
        Path('file2.h5ad'),
        Path('file3.h5ad')
    ],
    output_path=Path('combined.h5ad'),
    exclude_conflicts=False,
    cell_type_column='cell_type',
    overwrite=True
)
```

### Analyzing Conflicts Programmatically

```python
import anndata as ad

# Load files
adata1 = ad.read_h5ad('file1.h5ad')
adata2 = ad.read_h5ad('file2.h5ad')

# Find duplicates
duplicates = set(adata1.obs_names) & set(adata2.obs_names)

# Compare cell types
for cell_id in duplicates:
    type1 = adata1.obs.loc[cell_id, 'cell_type']
    type2 = adata2.obs.loc[cell_id, 'cell_type']
    if type1 != type2:
        print(f"{cell_id}: {type1} vs {type2}")
```

## Related Tools

This tool complements:
- **scanpy**: For downstream analysis and batch correction
- **anndata**: For data manipulation
- **scvi-tools**: For integration with deep learning models
- **celltype_annotate_cli_v2.py**: For annotating spatial data with combined references

## Citation

If you use this tool, please cite the underlying packages:
- AnnData: https://anndata.readthedocs.io
- Scanpy: https://scanpy.readthedocs.io

## Support

For issues or questions, please check:
1. This README
2. The example output from `--help`
3. Log files with `--log-level DEBUG`
