# core-unit-bioinformatics/workflow-nxf-somatic-variantcomp: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
Dates follow ISO 8601 (`YYYY-MM-DD`).

## v0.2.0-dev - [2026-10-07]

### `Added`

- Tab-separated samplesheet (`sample`, `technology`, `library_type`, `tumor_id`, `normal_id`, `caller`, `vcf`, `vcf_index`, `tumor_bam`, `tumor_bai`, `normal_bam`, `normal_bai`, `target_bed`) with `n/a` for empty fields
- Relative samplesheet paths are resolved against the samplesheet's directory
- Cross-row samplesheet checks (duplicates, truth technology per sample, one BED and library type per technology, BAM/index consistency, file existence) reported all at once
- Parameters `--truth_technology` and `--comparison_mode`
- `pipeline_info/samplesheet_validated.tsv` with resolved paths
- Synthetic edge-case test dataset generator `tests/synthetic/make_synthetic_data.py`; `-profile test` uses its output from the test-data repo
- nf-test function tests for the samplesheet checks

## v0.1.0-dev - [2026-10-07]

Initial development version, created with the [nf-core](https://nf-co.re/) template (nf-core/tools 4.1.0).

### `Added`

- nf-core pipeline template (internal name `variantcomp`, repo `workflow-nxf-somatic-variantcomp`), default branch `main`
- Development conda environment `envs/variantcomp-dev.yaml`
- CUBI metadata: `CITATION.cff`, `pyproject.toml` (version + ruff/pytest config)
- Design and roadmap: `docs/pipeline-plan.md`, developer instructions: `CLAUDE.md`

### `Fixed`

### `Dependencies`

### `Deprecated`
