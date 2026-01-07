#!/usr/bin/env python3
"""
Combine Multiple H5AD Files CLI Tool

This tool combines two or more h5ad single-cell RNA-seq files into a single dataset.
Handles duplicate cells, gene harmonization, and cell type conflict resolution.
"""

import gc
import io
import logging
import shlex
import subprocess
import sys
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple
import random

import anndata as ad
import click
import numpy as np
import psutil
from collections import Counter


def get_memory_usage() -> Dict[str, float]:
    """Get current memory usage in MB."""
    process = psutil.Process()
    mem_info = process.memory_info()
    return {
        'rss_mb': mem_info.rss / 1024 / 1024,  # Resident Set Size
        'vms_mb': mem_info.vms / 1024 / 1024,  # Virtual Memory Size
        'percent': process.memory_percent()
    }


def log_memory(prefix: str = "") -> None:
    """Log current memory usage."""
    mem = get_memory_usage()
    msg = f"Memory: RSS={mem['rss_mb']:.1f}MB, VMS={mem['vms_mb']:.1f}MB, {mem['percent']:.1f}%"
    if prefix:
        msg = f"{prefix} - {msg}"
    logging.info(msg)


def setup_logging(level: str = "INFO", log_file: Path = None) -> None:
    """Setup logging configuration with optional file output."""
    handlers = [logging.StreamHandler()]

    if log_file:
        file_handler = logging.FileHandler(log_file, mode='w', encoding='utf-8')
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        handlers.append(file_handler)

    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format="%(message)s",
        handlers=handlers,
        force=True  # Override any existing config
    )


def find_cell_type_column(adata: ad.AnnData) -> str:
    """
    Find the column containing cell type information.

    Args:
        adata: AnnData object

    Returns:
        Name of cell type column

    Raises:
        ValueError: If no cell type column found
    """
    # Common cell type column names, in order of preference
    common_names = ['cell_type', 'celltype', 'Celltype', 'cell_label', 'annotation', 'cluster']

    for name in common_names:
        if name in adata.obs.columns:
            return name

    # If not found, look for columns containing 'type' or 'label'
    for col in adata.obs.columns:
        if 'type' in col.lower() or 'label' in col.lower():
            logging.warning(f"Using column '{col}' as cell type column")
            return col

    raise ValueError(f"No cell type column found in obs. Available columns: {list(adata.obs.columns)}")


def load_metadata_only(file_path: Path) -> ad.AnnData:
    """
    Load only metadata (obs, var) without the expression matrix to save memory.

    Args:
        file_path: Path to h5ad file

    Returns:
        Metadata-only AnnData object (no X matrix)
    """
    logging.debug(f"  Loading metadata from {file_path.name}")
    # Load with backed mode to avoid loading X into memory
    adata = ad.read_h5ad(file_path, backed='r')

    # Extract only what we need
    obs_df = adata.obs.copy()
    var_names = adata.var_names.copy()

    # Close the backed file
    adata.file.close()
    del adata
    gc.collect()

    # Create a minimal AnnData with no X matrix
    metadata_adata = ad.AnnData(obs=obs_df, var=var_names.to_frame())

    return metadata_adata


def detect_duplicates(adatas: List[ad.AnnData], file_paths: List[Path]) -> Dict[str, List[int]]:
    """
    Detect duplicate cell IDs across multiple datasets.

    Args:
        adatas: List of AnnData objects (can be metadata-only)
        file_paths: List of file paths for reporting

    Returns:
        Dictionary mapping cell IDs to list of file indices where they appear
    """
    logging.debug("Scanning for duplicate cell IDs...")
    cell_to_files = {}

    for i, adata in enumerate(adatas):
        for cell_id in adata.obs_names:
            if cell_id not in cell_to_files:
                cell_to_files[cell_id] = []
            cell_to_files[cell_id].append(i)

    # Keep only duplicates
    duplicates = {cell_id: file_indices for cell_id, file_indices in cell_to_files.items()
                  if len(file_indices) > 1}

    logging.debug(f"Found {len(duplicates)} duplicate cell IDs")
    return duplicates


def analyze_cell_type_conflicts(
    duplicates: Dict[str, List[int]],
    adatas: List[ad.AnnData],
    cell_type_columns: List[str]
) -> Tuple[Dict[str, Dict], int, int]:
    """
    Analyze cell type conflicts for duplicate cells.

    Args:
        duplicates: Dictionary of duplicate cell IDs to file indices
        adatas: List of AnnData objects
        cell_type_columns: List of cell type column names for each dataset

    Returns:
        Tuple of (conflict_details, num_conflicts, num_resolved)
    """
    conflict_details = {}
    num_conflicts = 0
    num_resolved = 0

    for cell_id, file_indices in duplicates.items():
        # Get cell types for this cell from each file
        cell_types = []
        for file_idx in file_indices:
            ct_col = cell_type_columns[file_idx]
            cell_type = adatas[file_idx].obs.loc[cell_id, ct_col]
            cell_types.append(cell_type)

        # Check if there's a conflict
        unique_types = set(cell_types)
        if len(unique_types) > 1:
            num_conflicts += 1
            # Determine resolution by majority vote
            type_counts = Counter(cell_types)
            most_common = type_counts.most_common()

            # Check if there's a clear winner
            if len(most_common) > 1 and most_common[0][1] == most_common[1][1]:
                # Tie - pick randomly
                resolved_type = random.choice([t for t, c in most_common if c == most_common[0][1]])
                resolution = 'random'
            else:
                # Majority vote
                resolved_type = most_common[0][0]
                resolution = 'majority'
                num_resolved += 1

            conflict_details[cell_id] = {
                'cell_types': cell_types,
                'file_indices': file_indices,
                'type_counts': dict(type_counts),
                'resolved_type': resolved_type,
                'resolution': resolution
            }

    return conflict_details, num_conflicts, num_resolved


def run_git_command(args: List[str]) -> str | None:
    """Run a git command and return its stripped stdout, or None if it fails."""
    try:
        result = subprocess.run(
            ["git", *args],
            check=True,
            capture_output=True,
            text=True
        )
        return result.stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None


def get_git_metadata() -> Dict[str, str]:
    """Collect git metadata to describe repository state for reproducibility."""
    commit = run_git_command(["rev-parse", "HEAD"]) or "unknown"
    describe = run_git_command(["describe", "--tags", "--always", "--dirty"]) or "unknown"
    branch = run_git_command(["rev-parse", "--abbrev-ref", "HEAD"]) or "unknown"
    status = run_git_command(["status", "--short", "--branch"]) or ""
    toplevel = run_git_command(["rev-parse", "--show-toplevel"]) or ""
    dirty_flag = "true" if run_git_command(["status", "--porcelain"]) else "false"
    remote = run_git_command(["remote", "get-url", "origin"]) or ""

    return {
        "commit": commit,
        "describe": describe,
        "branch": branch,
        "status": status,
        "toplevel": toplevel,
        "is_dirty": dirty_flag,
        "remote": remote
    }


def write_provenance_file(
    output_path: Path,
    command_line: str,
    git_info: Dict[str, str],
    cli_output: str
) -> Path:
    """Write companion provenance file capturing command invocation, git metadata, and CLI output."""
    provenance_path = output_path.with_suffix(".txt")
    status_lines = git_info.get("status", "").splitlines()
    status_block = "\n".join(f"  {line}" for line in status_lines) if status_lines else "  <none>"
    cli_output = cli_output.rstrip("\n")
    provenance_path.write_text(
        "command:\n"
        f"  {command_line}\n"
        "git_repository:\n"
        f"  path: {git_info.get('toplevel', '') or '<unknown>'}\n"
        f"  branch: {git_info.get('branch', 'unknown')}\n"
        f"  commit: {git_info.get('commit', 'unknown')}\n"
        f"  describe: {git_info.get('describe', 'unknown')}\n"
        f"  remote_origin: {git_info.get('remote', '') or '<none>'}\n"
        f"  is_dirty: {git_info.get('is_dirty', 'unknown')}\n"
        "git_status:\n"
        f"{status_block}\n"
        "cli_output:\n"
        f"{cli_output}\n",
        encoding="utf-8"
    )
    return provenance_path


def combine_h5ad_files(
    file_paths: List[Path],
    output_path: Path,
    exclude_conflicts: bool = False,
    cell_type_column: str = None,
    overwrite: bool = False
) -> None:
    """
    Combine multiple h5ad files into a single dataset with memory optimization.

    Uses a two-pass approach:
    1. First pass: Load metadata only to detect duplicates and find common genes
    2. Second pass: Load and process full data one file at a time

    Args:
        file_paths: List of paths to h5ad files to combine
        output_path: Path for output combined h5ad file
        exclude_conflicts: If True, exclude cells with cell type conflicts
        cell_type_column: Specific cell type column name (if None, auto-detect)
        overwrite: Whether to overwrite existing output file
    """
    if output_path.exists() and not overwrite:
        raise FileExistsError(f"Output file exists: {output_path}. Use --overwrite to replace.")

    log_memory("Starting combination")

    # ====== FIRST PASS: Load metadata only ======
    logging.info(f"Pass 1: Loading metadata from {len(file_paths)} h5ad files...")
    log_memory("Before metadata loading")

    metadata_adatas = []
    cell_type_columns = []

    for i, fpath in enumerate(file_paths, 1):
        logging.info(f"  [{i}/{len(file_paths)}] {fpath.name}")
        metadata_adata = load_metadata_only(fpath)

        # Find or verify cell type column
        if cell_type_column:
            if cell_type_column not in metadata_adata.obs.columns:
                raise ValueError(f"Column '{cell_type_column}' not found in {fpath.name}")
            ct_col = cell_type_column
        else:
            ct_col = find_cell_type_column(metadata_adata)

        cell_type_columns.append(ct_col)
        logging.info(f"      {metadata_adata.n_obs:,} cells x {metadata_adata.n_vars:,} genes (cell_type column: '{ct_col}')")
        metadata_adatas.append(metadata_adata)

    log_memory("After metadata loading")

    # Detect duplicates using metadata only
    logging.info("\nAnalyzing duplicate cells...")
    duplicates = detect_duplicates(metadata_adatas, file_paths)

    if duplicates:
        logging.info(f"  Found {len(duplicates):,} duplicate cell IDs")

        # Analyze conflicts
        conflict_details, num_conflicts, num_resolved = analyze_cell_type_conflicts(
            duplicates, metadata_adatas, cell_type_columns
        )

        if num_conflicts > 0:
            logging.warning(f"\n⚠️  Cell type conflicts detected for {num_conflicts:,} cells")
            logging.warning(f"   - Resolved by majority vote: {num_resolved:,}")
            logging.warning(f"   - Resolved randomly (tie): {num_conflicts - num_resolved:,}")

            # Show examples
            logging.info("\nConflict examples:")
            for i, (cell_id, details) in enumerate(list(conflict_details.items())[:5], 1):
                logging.info(f"  {i}. {cell_id}")
                logging.info(f"     Cell types: {details['type_counts']}")
                logging.info(f"     Resolved to: '{details['resolved_type']}' ({details['resolution']})")

            if len(conflict_details) > 5:
                logging.info(f"     ... and {len(conflict_details) - 5} more conflicts")
    else:
        logging.info("  No duplicate cell IDs found")
        conflict_details = {}

    log_memory("After duplicate analysis")

    # Find common genes using metadata
    logging.info("\nHarmonizing genes...")
    all_genes = [set(adata.var_names) for adata in metadata_adatas]
    common_genes = set.intersection(*all_genes)
    logging.info(f"  Common genes across all files: {len(common_genes):,}")

    # Report unique genes per file
    for i, (adata, fpath) in enumerate(zip(metadata_adatas, file_paths)):
        unique_genes = set(adata.var_names) - common_genes
        if unique_genes:
            logging.info(f"  {fpath.name}: {len(unique_genes):,} unique genes (will be excluded)")

    common_genes_sorted = sorted(common_genes)

    # Prepare conflict resolution data
    cells_to_remove = set()
    cells_to_update = {}

    if duplicates:
        logging.info("\nPreparing duplicate resolution...")
        for cell_id, file_indices in duplicates.items():
            if cell_id in conflict_details:
                if exclude_conflicts:
                    cells_to_remove.add(cell_id)
                else:
                    resolved_type = conflict_details[cell_id]['resolved_type']
                    cells_to_update[cell_id] = resolved_type

        if exclude_conflicts and cells_to_remove:
            logging.info(f"  Will exclude {len(cells_to_remove):,} cells with conflicts")
        if cells_to_update:
            logging.info(f"  Will update cell types for {len(cells_to_update):,} conflict cells")

    # Free metadata memory
    del metadata_adatas
    gc.collect()
    log_memory("After freeing metadata")

    # ====== SECOND PASS: Load and process files incrementally ======
    logging.info(f"\nPass 2: Loading and processing full data...")
    cells_seen = set()
    processed_adatas = []

    for i, fpath in enumerate(file_paths):
        logging.info(f"  [{i+1}/{len(file_paths)}] Processing {fpath.name}")
        log_memory(f"  Before loading file {i+1}")

        # Load full data
        adata = ad.read_h5ad(fpath)
        logging.debug(f"    Loaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

        # Subset to common genes
        adata = adata[:, common_genes_sorted].copy()
        logging.debug(f"    After gene subsetting: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

        # Standardize cell type column name
        if cell_type_columns[i] != 'cell_type':
            adata.obs['cell_type'] = adata.obs[cell_type_columns[i]]

        # Convert categorical cell_type to string to allow updates
        if adata.obs['cell_type'].dtype.name == 'category':
            logging.debug(f"    Converting categorical cell_type to string")
            adata.obs['cell_type'] = adata.obs['cell_type'].astype(str)

        # Add source information
        adata.obs['source_file'] = fpath.stem
        adata.obs['source_index'] = i

        # Apply conflict resolutions
        if cells_to_update:
            for cell_id, resolved_type in cells_to_update.items():
                if cell_id in adata.obs_names:
                    adata.obs.loc[cell_id, 'cell_type'] = resolved_type

        # Remove conflict cells if needed
        if cells_to_remove:
            mask = ~adata.obs_names.isin(cells_to_remove)
            adata = adata[mask].copy()
            logging.debug(f"    After removing conflicts: {adata.n_obs:,} cells")

        # Remove duplicates (keep first occurrence)
        mask = []
        removed = 0
        for cell_id in adata.obs_names:
            if cell_id not in cells_seen:
                cells_seen.add(cell_id)
                mask.append(True)
            else:
                mask.append(False)
                removed += 1

        if removed > 0:
            logging.info(f"    Removed {removed:,} duplicate cells (keeping first occurrence)")
            adata = adata[np.array(mask)].copy()

        logging.debug(f"    Final: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
        processed_adatas.append(adata)
        log_memory(f"  After processing file {i+1}")

    log_memory("Before concatenation")

    # Concatenate processed datasets
    logging.info("\nCombining datasets...")
    combined = ad.concat(processed_adatas, join='inner', index_unique=None)
    logging.info(f"  Concatenation complete: {combined.n_obs:,} cells x {combined.n_vars:,} genes")

    # Free processed datasets memory
    del processed_adatas
    gc.collect()
    log_memory("After concatenation")

    # Add metadata
    combined.uns['combine_metadata'] = {
        'source_files': [str(p) for p in file_paths],
        'combine_date': datetime.now().isoformat(),
        'num_files': len(file_paths),
        'exclude_conflicts': exclude_conflicts,
        'num_duplicates': len(duplicates),
        'num_conflicts': num_conflicts if duplicates else 0
    }

    # Summary statistics
    logging.info("\nCombined dataset summary:")
    logging.info(f"  Total cells: {combined.n_obs:,}")
    logging.info(f"  Total genes: {combined.n_vars:,}")
    logging.info(f"  Unique cell types: {combined.obs['cell_type'].nunique()}")
    logging.info(f"\nCell type distribution:")

    ct_counts = combined.obs['cell_type'].value_counts()
    for ct, count in ct_counts.items():
        logging.info(f"    {ct:40s}: {count:6,} cells")

    log_memory("Before saving")

    # Save
    logging.info(f"\nSaving combined dataset to: {output_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    combined.write_h5ad(output_path)
    logging.info("✅ Done!")

    log_memory("After saving")
    logging.info(f"\nOutput file size: {output_path.stat().st_size / 1024 / 1024:.1f} MB")


@click.command()
@click.argument("input_files", nargs=-1, required=True, type=click.Path(exists=True, path_type=Path))
@click.argument("output_file", type=click.Path(path_type=Path))
@click.option(
    "--exclude-conflicts",
    is_flag=True,
    help="Exclude cells with cell type conflicts instead of resolving them."
)
@click.option(
    "--cell-type-column",
    type=str,
    help="Specific cell type column name (default: auto-detect 'cell_type' or similar)."
)
@click.option(
    "--overwrite",
    is_flag=True,
    help="Overwrite output file if it exists."
)
@click.option(
    "--log-level",
    type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR"], case_sensitive=False),
    default="INFO",
    help="Logging level.",
    show_default=True
)
def main(
    input_files: Tuple[Path, ...],
    output_file: Path,
    exclude_conflicts: bool,
    cell_type_column: str,
    overwrite: bool,
    log_level: str
):
    """
    Combine multiple h5ad single-cell RNA-seq files with memory optimization.

    Handles duplicate cells and cell type conflicts:
    - Detects duplicate cell IDs across files
    - Resolves cell type conflicts by majority vote
    - For ties, selects randomly or excludes based on --exclude-conflicts
    - Uses two-pass approach to minimize memory usage

    \b
    INPUT_FILES: Two or more h5ad files to combine
    OUTPUT_FILE: Path for combined output h5ad file

    \b
    Examples:
        # Combine two files, resolve conflicts
        uv run python combine_h5ad_cli.py file1.h5ad file2.h5ad combined.h5ad

        # Combine three files, exclude conflicts
        uv run python combine_h5ad_cli.py file1.h5ad file2.h5ad file3.h5ad combined.h5ad --exclude-conflicts

        # Use specific cell type column
        uv run python combine_h5ad_cli.py file1.h5ad file2.h5ad combined.h5ad --cell-type-column celltype

    A log file (.log) and provenance text file (.txt) are written alongside the output
    containing memory usage, full command, captured log output, and repository state
    for reproducibility.
    """
    # Create log file with same name as output but .log extension
    log_file = output_file.with_suffix('.log')
    setup_logging(log_level, log_file)

    logging.info("="*80)
    logging.info("Combine H5AD Files - Memory Optimized Version")
    logging.info("="*80)
    logging.info(f"Log file: {log_file}")
    logging.info(f"Log level: {log_level}")
    logging.info(f"Input files: {len(input_files)}")
    for i, f in enumerate(input_files, 1):
        logging.info(f"  [{i}] {f}")
    logging.info(f"Output file: {output_file}")
    logging.info("="*80)

    if len(input_files) < 2:
        click.echo("Error: At least 2 input files required", err=True)
        sys.exit(1)

    command_line = shlex.join(sys.argv)
    git_info = get_git_metadata()
    capture_stream = io.StringIO()
    capture_handler = logging.StreamHandler(capture_stream)
    capture_handler.setFormatter(logging.Formatter("%(message)s"))
    root_logger = logging.getLogger()
    root_logger.addHandler(capture_handler)

    try:
        combine_h5ad_files(
            file_paths=list(input_files),
            output_path=output_file,
            exclude_conflicts=exclude_conflicts,
            cell_type_column=cell_type_column,
            overwrite=overwrite
        )
        log_contents = capture_stream.getvalue()
        provenance_path = write_provenance_file(output_file, command_line, git_info, log_contents)
        logging.info(f"Provenance written to: {provenance_path}")
        updated_log_contents = capture_stream.getvalue()
        if updated_log_contents != log_contents:
            write_provenance_file(output_file, command_line, git_info, updated_log_contents)
    except Exception as e:
        click.echo(f"❌ Error: {e}", err=True)
        if log_level == "DEBUG":
            raise
        sys.exit(1)
    finally:
        root_logger.removeHandler(capture_handler)
        capture_handler.close()


if __name__ == "__main__":
    main()
