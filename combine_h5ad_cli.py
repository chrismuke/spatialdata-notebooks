#!/usr/bin/env python3
"""
Combine Multiple H5AD Files CLI Tool

This tool combines two or more h5ad single-cell RNA-seq files into a single dataset.
Handles duplicate cells, gene harmonization, and cell type conflict resolution.
"""

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
from collections import Counter


def setup_logging(level: str = "INFO") -> None:
    """Setup logging configuration."""
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format="%(message)s",
        handlers=[logging.StreamHandler()]
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


def detect_duplicates(adatas: List[ad.AnnData], file_paths: List[Path]) -> Dict[str, List[int]]:
    """
    Detect duplicate cell IDs across multiple datasets.

    Args:
        adatas: List of AnnData objects
        file_paths: List of file paths for reporting

    Returns:
        Dictionary mapping cell IDs to list of file indices where they appear
    """
    cell_to_files = {}

    for i, adata in enumerate(adatas):
        for cell_id in adata.obs_names:
            if cell_id not in cell_to_files:
                cell_to_files[cell_id] = []
            cell_to_files[cell_id].append(i)

    # Keep only duplicates
    duplicates = {cell_id: file_indices for cell_id, file_indices in cell_to_files.items()
                  if len(file_indices) > 1}

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
    Combine multiple h5ad files into a single dataset.

    Args:
        file_paths: List of paths to h5ad files to combine
        output_path: Path for output combined h5ad file
        exclude_conflicts: If True, exclude cells with cell type conflicts
        cell_type_column: Specific cell type column name (if None, auto-detect)
        overwrite: Whether to overwrite existing output file
    """
    if output_path.exists() and not overwrite:
        raise FileExistsError(f"Output file exists: {output_path}. Use --overwrite to replace.")

    # Load all files
    logging.info(f"Loading {len(file_paths)} h5ad files...")
    adatas = []
    cell_type_columns = []

    for i, fpath in enumerate(file_paths, 1):
        logging.info(f"  [{i}/{len(file_paths)}] {fpath.name}")
        adata = ad.read_h5ad(fpath)

        # Find or verify cell type column
        if cell_type_column:
            if cell_type_column not in adata.obs.columns:
                raise ValueError(f"Column '{cell_type_column}' not found in {fpath.name}")
            ct_col = cell_type_column
        else:
            ct_col = find_cell_type_column(adata)

        cell_type_columns.append(ct_col)
        logging.info(f"      {adata.n_obs:,} cells x {adata.n_vars:,} genes (cell_type column: '{ct_col}')")
        adatas.append(adata)

    # Detect duplicates
    logging.info("\nAnalyzing duplicate cells...")
    duplicates = detect_duplicates(adatas, file_paths)

    if duplicates:
        logging.info(f"  Found {len(duplicates):,} duplicate cell IDs")

        # Analyze conflicts
        conflict_details, num_conflicts, num_resolved = analyze_cell_type_conflicts(
            duplicates, adatas, cell_type_columns
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

    # Find common genes
    logging.info("\nHarmonizing genes...")
    all_genes = [set(adata.var_names) for adata in adatas]
    common_genes = set.intersection(*all_genes)
    logging.info(f"  Common genes across all files: {len(common_genes):,}")

    # Report unique genes per file
    for i, (adata, fpath) in enumerate(zip(adatas, file_paths)):
        unique_genes = set(adata.var_names) - common_genes
        if unique_genes:
            logging.info(f"  {fpath.name}: {len(unique_genes):,} unique genes (will be excluded)")

    # Subset to common genes
    common_genes_sorted = sorted(common_genes)
    adatas_subset = [adata[:, common_genes_sorted].copy() for adata in adatas]

    # Standardize cell type column name
    logging.info("\nStandardizing cell type annotations...")
    for i, adata in enumerate(adatas_subset):
        if cell_type_columns[i] != 'cell_type':
            adata.obs['cell_type'] = adata.obs[cell_type_columns[i]]

    # Add source file information
    for i, (adata, fpath) in enumerate(zip(adatas_subset, file_paths)):
        adata.obs['source_file'] = fpath.stem
        adata.obs['source_index'] = i

    # Handle duplicates and conflicts
    if duplicates:
        logging.info("\nResolving duplicates...")
        cells_to_remove = set()
        cells_to_update = {}

        for cell_id, file_indices in duplicates.items():
            if cell_id in conflict_details:
                # Has conflict
                if exclude_conflicts:
                    # Remove from all datasets
                    cells_to_remove.add(cell_id)
                else:
                    # Update to resolved type
                    resolved_type = conflict_details[cell_id]['resolved_type']
                    cells_to_update[cell_id] = resolved_type

        # Apply removals
        if exclude_conflicts and cells_to_remove:
            logging.info(f"  Excluding {len(cells_to_remove):,} cells with conflicts")
            for adata in adatas_subset:
                mask = ~adata.obs_names.isin(cells_to_remove)
                adata._inplace_subset_obs(mask)

        # Apply updates
        if cells_to_update:
            logging.info(f"  Updating cell types for {len(cells_to_update):,} conflict cells")
            for adata in adatas_subset:
                for cell_id, resolved_type in cells_to_update.items():
                    if cell_id in adata.obs_names:
                        adata.obs.loc[cell_id, 'cell_type'] = resolved_type

        # Remove duplicate cells (keep first occurrence)
        logging.info("  Removing duplicate entries (keeping first occurrence)")
        cells_seen = set()
        for i, adata in enumerate(adatas_subset):
            mask = []
            for cell_id in adata.obs_names:
                if cell_id not in cells_seen:
                    cells_seen.add(cell_id)
                    mask.append(True)
                else:
                    mask.append(False)
            mask = np.array(mask)
            removed = (~mask).sum()
            if removed > 0:
                logging.info(f"    File {i+1}: Removed {removed:,} duplicate cells")
            adatas_subset[i] = adata[mask].copy()

    # Concatenate
    logging.info("\nCombining datasets...")
    combined = ad.concat(adatas_subset, join='inner', index_unique=None)

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

    # Save
    logging.info(f"\nSaving combined dataset to: {output_path}")
    combined.write_h5ad(output_path)
    logging.info("✅ Done!")


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
    Combine multiple h5ad single-cell RNA-seq files.

    Handles duplicate cells and cell type conflicts:
    - Detects duplicate cell IDs across files
    - Resolves cell type conflicts by majority vote
    - For ties, selects randomly or excludes based on --exclude-conflicts

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

    A provenance text file is written alongside the output containing the full command,
    captured log output, and repository state for reproducibility.
    """
    setup_logging(log_level)

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
