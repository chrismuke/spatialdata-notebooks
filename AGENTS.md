# Repository Guidelines

## Project Structure & Module Organization
- Notebooks: `notebooks/` (tutorials, examples, paper reproducibility). Prefer descriptive, snake_case filenames (e.g., `celltype_annotation_celltypist.ipynb`).
- Scripts & CLIs: Python utilities at repo root (e.g., `combine_zarr_cli.py`, `run_celltype_batch.py`).
- Docs: Sphinx config in repo root (`conf.py`, `index.md`, `Makefile`), static assets in `_static/`, build output in `_build/`.
- Data & results: Example inputs in `datasets/`, demo outputs in `results/` and similar. Avoid committing large binaries.

## Build, Test, and Development Commands
- Environment (primary): `uv sync` (add-ons: `--extra dev`, `--extra doc`).
- Run tools via uv: `uv run jupyter lab`, `uv run python combine_zarr_cli.py --help`.
- Lint/format: `uv run pre-commit run --all-files` (Ruff check/format and hygiene hooks).
- Docs build: `make html` → open `_build/html/index.html`; clean with `make clean`.
- Quick CLI examples:
  - `uv run python celltype_annotate_cli_v2.py data.zarr ref.h5ad --show-unknown-cells`
  - `uv run python combine_zarr_cli.py a.zarr b.zarr out.zarr --border 100 --layout vertical`

## Coding Style & Naming Conventions
- Language: Python ≥ 3.10.
- Formatting/linting: Ruff (line length 120) via pre-commit; `uv run ruff check . --fix`, `uv run ruff format .`.
- Indentation: 4 spaces; keep imports tidy; prefer small, focused modules/CLIs.
- Naming: snake_case for modules, functions, variables; descriptive notebook filenames; `_cli.py` suffix for command-line tools.

## Testing Guidelines
- No formal pytest suite; treat Sphinx build and notebooks as smoke tests.
- Before a PR: (1) `uv run pre-commit run --all-files`, (2) `make html` succeeds, (3) key notebooks run cleanly, (4) CLIs respond to `--help` and run on a tiny sample.
- Keep validation data under `datasets/`; do not commit large artifacts.

## Commit & Pull Request Guidelines
- Commits: imperative mood, clear scope (e.g., "Add zarr combination CLI", "Fix spatial colocalization bug"). Avoid vague "update".
- PRs: include description, linked issues, and screenshots of notebook outputs when relevant; note non-default env details.
- Quality gate: local checks only — run pre-commit and docs build; ensure notebooks are reproducible (minimal outputs, deterministic seeds where applicable).

## Results & Storage Notes
- SpatialData favors storage efficiency; INFO messages about non self-contained objects/backing store changes are expected.
- Results can be large; prefer `--results-dir`/`--output-dir` flags for tools (e.g., `celltype_annotate_cli_v2.py`, `roi_umap_analysis.py`).
- Default plots may hide unknown cell types; use `--show-unknown-cells` to include them.

## Security & Data Handling
- Never commit credentials or private data; pre-commit scans for secrets.
- Do not commit large datasets; keep Zarr stores under `datasets/` locally and reference paths instead of binaries.
