# Dry-run harness

This directory is a **dry-run harness**, not an execution dataset. Its `ref/` and `data/`
files are empty placeholders and its `config/` points at them, so

```bash
snakemake -s workflow/Snakefile -d .test -n <target>
```

validates the DAG (rule wiring, wildcard resolution, per-sample narrow/broad and IgG
routing) without needing real genomes or reads. `config/samples.csv` here exercises a
2-replicate narrow condition (→ IDR), a 1-replicate broad condition (→ single), and a
shared IgG control. Running the workflow for real (without `-n`) against these placeholders
will not produce meaningful results.
