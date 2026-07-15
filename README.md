# CUT&RUN Analysis Pipeline

A Snakemake workflow for paired-end CUT&RUN data — from raw reads to peaks, consensus
count matrices, extensive QC, and differential binding. It mirrors the structure of the
sibling ATAC-seq pipeline but is adapted for CUT&RUN: **no spike-in**, per-sample
**narrow/broad** peak mode, per-sample **IgG/Input controls**, **both MACS2 and SEACR**
peak callers (each with its own consensus count matrix), and **IgG-subtracted** signal
tracks.

## Overview

Two stages live in one standard-layout `workflow/Snakefile` and build in a single unified
DAG:

1. **Primary pipeline** (`cutandrun_all`): FastQC → fastp → human-only Bowtie2 alignment →
   unique/properly-paired filter (+ mitochondrial-% QC) → Picard dedup → ENCODE blacklist
   filter → RPGC + IgG-subtracted bigWigs → **MACS2** (narrow/broad, IgG as `-c` control)
   **and SEACR** (stringent/relaxed, IgG control) peak calls → a fixed-width **MACS2
   consensus** (Corces-2018 SPM, majority/IDR/single reproducibility) and an overlap-based
   variable-width **SEACR consensus**, each with a `featureCounts` matrix.
2. **QC pipeline** (`qc_all`): deepTools fragment-size / fingerprint / correlation / PCA /
   GC-bias / TSS enrichment, FRiP (MACS2 **and** SEACR), IDR on replicate pairs, library
   complexity (NRF/PBC1/PBC2), reads-in-annotation, peak summaries, a FastQC-only MultiQC,
   and a self-contained interactive **`cutandrun_qc_report.html`**.

A separate R/Bioconductor notebook, **`cutandrun_Dx.ipynb`**, runs generic DESeq2
differential binding on the consensus matrices.

```
Raw FASTQ → FastQC → fastp
  → Bowtie2 (human genome) → unique + properly-paired filter → mito-% QC → keep chr1-22/X
     → Picard dedup → blacklist filter
        → RPGC bigWig  +  log2(treatment/IgG) bigWig
        → MACS2 peaks (narrow|broad, IgG -c)  → IDR/consensus → featureCounts matrix
        → SEACR peaks (stringent|relaxed, IgG) → overlap consensus → featureCounts matrix
```

## Control vs treatment

The pipeline decides a sample's role from the sample sheet: a row with an **empty
`peak_mode`** is a **control** (IgG/Input) — it is aligned, filtered, deduplicated,
blacklist-filtered and turned into an RPGC track, and it is used as the MACS2 `-c` control,
the SEACR control, and the `bamCompare -b2` for its matched treatments, **but it is never
peak-called**. Every other row is a **treatment** whose `peak_mode` (`narrow`/`broad`)
picks MACS2 narrow-vs-`--broad` and the SEACR stringency, and whose `input_control` names
the control to pair with it.

## Peak callers

- **MACS2** is the backbone: it feeds IDR reproducibility, the fixed-width consensus, the
  `featureCounts` matrix, and (by default) the differential notebook.
- **SEACR** runs in parallel as a CUT&RUN-native alternative, with its own FRiP/peak-count
  QC **and its own count matrix** (`results/consensus_seacr/`) built by an overlap-based
  reproducible variable-width union. Both matrices are available to the differential
  notebook (`differential_counts: macs2 | seacr | both`).

## Configuration

All parameters live in `config/config.yaml`, validated against the single-source-of-truth
schema [`workflow/schemas/config.schema.yaml`](workflow/schemas/config.schema.yaml). The
sample sheet and full parameter reference are in
[`config/README.md`](config/README.md). At minimum set the reference paths (`genome_fasta`,
`blacklist`, `gtf`, `promoter_bed`, `enhancer_bed`) and edit `config/samples.csv`.

### Sample sheet (`config/samples.csv`)

`sample_id,condition,replicate,input_control,peak_mode,notes` — e.g.:

```csv
sample_id,condition,replicate,input_control,peak_mode,notes
GSF2801-ChIPseq-OVCAR3-3D-IP-cJun_S4,cJUN_3D,1,GSF2801-ChIPseq-OVCAR3-3D-IP-IgG_S5,narrow,3D-cJUN
GSF2801-ChIPseq-OVCAR3-3D-IP-IgG_S5,IgG_3D,1,,,3D-Igg
```

`condition` is the reproducibility group (≥3 reps → majority vote, 2 reps → IDR, 1 rep →
single). Give biologically distinct groups distinct `condition` labels. All replicates of a
condition must share one `peak_mode`. See `config/README.md`.

Reads must be `data/<sample_id>_R1_001.fastq.gz` / `_R2_001.fastq.gz`.

## Reference data (`ref/`)

Provide (not shipped):

- `ref/genome.fa` — chr-prefixed UCSC hg38 FASTA
- `ref/gencode.v36.annotation.gtf` — GENCODE annotation (TSS QC)
- `ref/hg38.2bit` — for `computeGCBias`
- `ref/picard.jar` — Picard MarkDuplicates

Shipped: `ref/hg38_blacklist_regions.bed`, `ref/promoter_chr1-22X.bed`,
`ref/enhancer_chr1-22X.bed`. Built automatically: the Bowtie2 index (`ref/genome/`) and
`ref/genome.chrom.sizes`.

```bash
cd ref
curl -O https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz && gunzip hg38.fa.gz && mv hg38.fa genome.fa
curl -O https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit
curl -O http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.annotation.gtf.gz && gunzip gencode.v36.annotation.gtf.gz
curl -L -o picard.jar https://github.com/broadinstitute/picard/releases/latest/download/picard.jar
cd ..
```

## Running

You need Snakemake + conda/mamba as the driver; per-rule tool environments
(`workflow/envs/*.yaml`) are created on the first `--use-conda` run.

```bash
mamba create -n cutandrun -c conda-forge -c bioconda snakemake-minimal pandas
conda activate cutandrun

# Dry run (check the DAG)
snakemake -s workflow/Snakefile -n

# Everything (primary -> QC) in one dependency-ordered run
snakemake -s workflow/Snakefile --use-conda --cores 20

# Or a single stage
snakemake -s workflow/Snakefile --use-conda --cores 20 cutandrun_all   # primary only
snakemake -s workflow/Snakefile --use-conda --cores 20 qc_all          # QC only
```

`./run_pipeline.sh --cores 16 [target]` is a thin wrapper around the same command.

### Differential binding

Generate the notebook (default MACS2 matrix; `seacr` or `both` also accepted), then run it
in the `cutandrun_Dx` env:

```bash
python workflow/scripts/build_diffbind_notebook.py macs2
jupyter nbconvert --to notebook --execute --inplace \
    --ExecutePreprocessor.kernel_name=ir cutandrun_Dx.ipynb
```

It runs DESeq2 (median-of-ratios) for every pairwise treatment-condition contrast, split
into promoter vs distal peaks, with PCA / MA / volcano / ChIPseeker annotation →
`results/diff_region/<matrix>/`.

## Outputs

```
results/
├── fastqc/ fastp/ aligned/ filtered/ dedup/ blacklist_filtered/   # QC + analysis-ready BAMs
├── bigwig/                 # RPGC bigWigs ({sample}.bw)
├── log2ratio_bigwig/       # IgG-subtracted log2(treatment/IgG) bigWigs
├── peaks/                  # MACS2 peaks (*_peaks.narrowPeak | *_peaks.broadPeak)
├── seacr/                  # SEACR peaks (*.stringent.bed | *.relaxed.bed) + bedgraph/
├── consensus/              # MACS2 fixed-width consensus + consensus_counts.txt
├── consensus_seacr/        # SEACR variable-width consensus + consensus_counts.txt
├── deeptools/ FRiP/ idr/ library_complexity/ peak_annotation/     # QC
├── qc/                     # cutandrun_qc_report.html, multiqc_fastqc.html, summaries
└── diff_region/            # DESeq2 differential binding (from cutandrun_Dx.ipynb)
```

## Citing the tools

Snakemake (Köster & Rahmann 2012); Bowtie2 (Langmead & Salzberg 2012); SAMtools (Li et al.
2009); Picard; fastp (Chen et al. 2018); FastQC (Andrews 2010); MACS2 (Zhang et al. 2008);
SEACR (Meers, Tenenbaum & Henikoff 2019); deepTools (Ramírez et al. 2016); BEDTools (Quinlan
& Hall 2010); IDR (Li et al. 2011); featureCounts / Subread (Liao et al. 2014); consensus
peaks (Corces et al. 2018); DESeq2 (Love, Huber & Anders 2014); ChIPseeker (Yu et al. 2015).

## License

MIT.
