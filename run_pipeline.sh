#!/usr/bin/env bash
# Run the CUT&RUN Snakemake workflow with per-rule conda envs.
#
# The workflow is the standard-layout workflow/Snakefile: a single run builds the
# primary pipeline AND the QC report (unified DAG). Pass extra snakemake args
# (cores, a target, -n, ...) straight through.
#
# Usage:
#   ./run_pipeline.sh                          # everything (primary -> QC), 4 cores
#   ./run_pipeline.sh --cores 16               # everything, 16 cores
#   ./run_pipeline.sh --cores 16 cutandrun_all # primary pipeline only
#   ./run_pipeline.sh --cores 16 qc_all        # QC only (after primary)
#   ./run_pipeline.sh -n                        # dry run: check the DAG first
set -euo pipefail

if [ "$#" -eq 0 ]; then set -- --cores 4; fi   # default snakemake args

echo ">> snakemake --use-conda -s workflow/Snakefile $*"
snakemake --use-conda -s workflow/Snakefile "$@"
