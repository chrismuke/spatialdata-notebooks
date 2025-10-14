# Repository Guidelines

## Project Structure & Module Organization
- Core notebooks live in `notebooks/`; use descriptive snake_case names such as `celltype_annotation_celltypist.ipynb`.
- CLI utilities sit at the repo root (for example `combine_zarr_cli.py`, `run_celltype_batch.py`); keep them focused and documented.
- Documentation sources (`conf.py`, `index.md`, `_static/`, `_build/`) reside at the root; avoid editing `_build/` by hand.
- Sample inputs stay in `datasets/`; tutorial outputs and other generated artifacts belong in `results/` or another dedicated output directory.

## Build, Test, and Development Commands
- `uv sync --extra dev --extra doc` sets up the full development environment.
- `uv run jupyter lab` launches the primary authoring environment for notebooks.
- `uv run python combine_zarr_cli.py --help` or similar `uv run` invocations validate CLIs without altering your shell PYTHONPATH.
- `uv run pre-commit run --all-files` runs Ruff formatting/checking and hygiene hooks.
- `make html` builds the Sphinx documentation; follow with `make clean` when regenerating from scratch.

## Coding Style & Naming Conventions
- Target Python 3.10+ with 4-space indentation and Ruff-enforced formatting (line length 120).
- Keep imports ordered and minimize module scope; prefer small, composable helpers.
- Use snake_case for files, functions, and variables; reserve `_cli.py` suffixes for entry points.
- Stick to ASCII unless a file already uses extended characters for scientific notation or labels.

## Testing Guidelines
- Treat `uv run pre-commit run --all-files` as the baseline gate.
- Ensure `make html` completes without warnings before publishing documentation changes.
- Execute representative notebooks end-to-end; clear stale outputs and set deterministic seeds where feasible.
- For CLIs, confirm `--help` succeeds and exercise tiny fixtures from `datasets/`.

## Commit & Pull Request Guidelines
- Write imperative commit subjects with clear scope (e.g., “Add zarr combination CLI”, “Fix spatial colocalization bug”); avoid generic “update” phrasing.
- Pull requests should summarize intent, link tracking issues, note environment deltas, and attach screenshots or key notebook outputs when relevant.
- Verify pre-commit and docs builds locally before opening or updating a PR.

## Security & Data Handling
- Never commit credentials or large datasets; rely on paths under `datasets/` and document expected inputs explicitly.
- INFO logs about SpatialData backing stores are expected—call out only warnings or regressions that affect reproducibility.
