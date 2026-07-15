# CUT&RUN Snakemake Workflow — Design

**Date:** 2026-07-14
**Working dir:** `/easley/scratch/projects/amitra/amitra2016502/snakemake_CutandRun_seq`
**Mirrors:** `../snakemake_ATACseq_spikein` (same standard Snakemake Workflow Catalog layout)

## 1. Goal

Produce a CUT&RUN paired-end analysis pipeline that mirrors the structure and code
style of the sibling ATAC-seq spike-in workflow, but adapted for CUT&RUN:

- **No spike-in.** Module A (concatenated human+spike-in alignment, spike-in extract/count,
  normalization factors, spike-in-scaled bigWigs, spike-in QC) is removed entirely.
  Alignment is to a human-only Bowtie2 index.
- **Per-sample peak mode.** Each treatment sample chooses `narrow` or `broad` peaks via a
  `peak_mode` column in the sample sheet.
- **Per-sample IgG/Input control.** Each treatment sample names its matched IgG/Input in an
  `input_control` column; that control BAM is passed to the peak callers.
- **Two peak callers:** MACS2 (backbone) **and** SEACR (CUT&RUN-native), run in parallel.
- **IgG-subtracted signal tracks:** deepTools `bamCompare` log2(treatment/IgG) bigWigs in
  addition to per-sample RPGC bigWigs.
- **Generic DESeq2 differential-binding notebook** (no spike-in, condition-vs-condition),
  replacing the ATAC repo's bespoke NICD3/Notch notebook.

The layout follows the Snakemake Workflow Catalog conventions: `workflow/Snakefile` includes
`rules/common.smk`, `rules/cutandrun.smk`, `rules/qc.smk`; parameters validated against
`workflow/schemas/config.schema.yaml`; per-rule conda envs under `workflow/envs/`.

## 2. Sample sheet (`config/samples.csv`)

Columns (matching the reference format
`gynecoloji/SnakeMake_CutandRun_seq/ref/samples.csv`):

| column          | required | meaning |
|-----------------|----------|---------|
| `sample_id`     | yes      | Sample name. Reads must be `data/<sample_id>_R1_001.fastq.gz` / `_R2_001.fastq.gz`. |
| `condition`     | yes      | Free-text label. For **treatment** rows this is the replicate group used for reproducibility. |
| `replicate`     | yes      | Replicate index within the condition (integer). |
| `input_control` | no       | `sample_id` of the matched IgG/Input control for this row. Empty for control rows. |
| `peak_mode`     | no       | `narrow` or `broad`. Empty marks the row as a **control** (IgG/Input), not peak-called. |
| `notes`         | no       | Free text; ignored by the pipeline. |

Example (from the reference):

```csv
sample_id,condition,replicate,input_control,peak_mode,notes
GSF2801-ChIPseq-OVCAR3-3D-IP-cJun_S4,cJUN,1,GSF2801-ChIPseq-OVCAR3-3D-IP-IgG_S5,narrow,3D-cJUN
GSF2801-ChIPseq-OVCAR3-3D-IP-IgG_S5,Igg,1,,,3D-Igg
GSF2801-ChIPseq-OVCAR3-Control-IP-cJun_S1,cJUN,1,GSF2801-ChIPseq-OVCAR3-Control-IP-IgG_S2,narrow,Ctrl-cJUN
GSF2801-ChIPseq-OVCAR3-Control-IP-IgG_S2,Igg,1,,,Ctrl-Igg
```

### Derived sample sets (in `common.smk`)

- `ALL_SAMPLES` — every row. All get: FastQC → fastp → align → filter → dedup →
  blacklist-filter → RPGC bigWig, and alignment/complexity/correlation QC.
- `TREATMENT_SAMPLES` — rows with non-empty `peak_mode`. Get peak calling (MACS2 + SEACR),
  FRiP, peak summary, consensus, differential counts. Also get an IgG-subtracted bigWig **iff**
  `input_control` is set.
- `CONTROL_SAMPLES` — rows with empty `peak_mode` (the IgG/Input samples). Never peak-called;
  serve as MACS2 `-c` / SEACR control / bamCompare `-b2`.
- `NARROW_SAMPLES` / `BROAD_SAMPLES` — treatment samples partitioned by `peak_mode` (drives the
  MACS2 output extension and SEACR stringency, via `wildcard_constraints` like the ATAC
  `_alt(...)` helper).
- `GROUPS` — `TREATMENT_SAMPLES` grouped by `condition` → list of `sample_id`. Per-group
  reproducibility method derived from replicate count, exactly as ATAC:
  ≥3 → majority vote, ==2 → IDR, ==1 → single (passthrough).

### Validation (fail fast in `common.smk`)

1. Every `input_control` value refers to an existing `CONTROL_SAMPLES` `sample_id`.
2. Within any `condition`, all treatment members share the same `peak_mode` (consensus/IDR are
   built per condition and cannot mix narrow+broad).
3. `peak_mode` values ∈ {narrow, broad, empty}.

## 3. What changes relative to the ATAC pipeline

**Removed (spike-in / Module A):**
- Rules: `build_combined_genome`, `spikein_extract_count`, `compute_spikein_factors`,
  `create_spikein_bigwig`, and QC `spikein_qc`.
- Scripts: `compute_spikein_factors.py`.
- Config keys: `spikein_fasta`, `spikein_prefix`, `combined_index`, `spikein_pct_min`,
  `spikein_pct_max`.
- Spike-in read splitting inside the filter rule (no prefix logic).

**Changed:**
- `build_combined_genome` → `build_genome_index`: `bowtie2-build` on the human FASTA
  (optionally subset to `align_chroms`); also emit `ref/<genome>.chrom.sizes` (for SEACR).
- `bowtie2_align`: align to the human-only index; CUT&RUN-appropriate params
  (`--end-to-end --very-sensitive --no-mixed --no-discordant --no-unal -I 10 -X {max_fragment_length}`,
  default `max_fragment_length: 700`).
- `samtools_sort_filter_index`: same unique + properly-paired + orphan-removal + mito-% idxstats +
  `keep_chroms` restriction, but **no** spike-in prefix filtering.
- `call_peaks`: now split into `call_peaks_macs2_narrow` / `call_peaks_macs2_broad`, each passing
  the matched IgG as `-c` when `input_control` is set.
- `consensus_peaks` + `relaxed_peaks` + `reproducible_idr`: made narrow/broad aware.
- Naming: `atacseq_all` → `cutandrun_all`; `atacseq_qc_report.html` → `cutandrun_qc_report.html`;
  `ATACseq_Dx.ipynb` → `cutandrun_Dx.ipynb`.

**Added:**
- `call_peaks_seacr` (+ `seacr_bedgraph`, `genome_chrom_sizes`).
- `create_log2ratio_bigwig` (deepTools `bamCompare`, treatment vs its IgG).
- `workflow/envs/seacr.yaml`.
- `remove_duplicates` gated by a config toggle `remove_duplicates: true` (CUT&RUN low-input runs
  sometimes keep duplicates). Default true = mirror ATAC.

## 4. Primary pipeline (`workflow/rules/cutandrun.smk`)

Aggregate target `cutandrun_all`. Rules (I/O sketched; log/threads/conda as in ATAC):

1. `fastqc` — unchanged. `data/{s}_R{1,2}_001.fastq.gz` → FastQC html.
2. `fastp` — unchanged. Adapter auto-detect / override via `adapter_r1/r2`.
3. `genome_chrom_sizes` — `samtools faidx {genome_fasta}` → `ref/genome.chrom.sizes`
   (`cut -f1,2 .fai`). Used by SEACR bedgraph + `bedtools genomecov -g`.
4. `build_genome_index` — subset human FASTA to `align_chroms` (or keep all), `bowtie2-build`
   → `{genome_index}.*`, touch `.build.done`.
5. `bowtie2_align` — align trimmed reads to `{genome_index}`; coordinate-sort + index; keep
   `{s}.bowtie2.log`.
6. `samtools_sort_filter_index` — properly-paired (`-f 2 -F 2316`) + unique (`grep -v XS:i:`) +
   orphan removal (`process_sam.py`); idxstats (mito-% QC) on pre-keep BAM; restrict to
   `keep_chroms`; flagstat.
7. `remove_duplicates` — Picard MarkDuplicates → `{s}.dedup.bam` + metrics. Config toggle
   `remove_duplicates` maps to Picard `REMOVE_DUPLICATES=true|false`: when false, duplicates are
   marked but kept (still emits valid metrics), so the `.dedup.bam` path and QC stay stable for
   low-input CUT&RUN runs. Default true = mirror ATAC.
8. `filter_blacklist` — fragment-level ENCODE blacklist removal (BEDPE) → `{s}.nobl.bam`.
9. `create_bigwig` — deepTools `bamCoverage --normalizeUsing RPGC` per sample (ALL samples).
10. `create_log2ratio_bigwig` — for TREATMENT samples with `input_control`:
    `bamCompare -b1 {s}.nobl.bam -b2 {igg}.nobl.bam --operation log2 --normalizeUsing CPM
    --binSize {bin_size} --blackListFileName {blacklist} -o {s}.log2ratio.bw`.
11. `call_peaks_macs2_narrow` (sample ∈ NARROW_SAMPLES) → `results/peaks/{s}_peaks.narrowPeak`.
    `macs2 callpeak -t {s}.nobl.bam [-c {igg}.nobl.bam] -f BAMPE -g {macs2_genome} --nomodel
    -q {macs2_qvalue}`.
12. `call_peaks_macs2_broad` (sample ∈ BROAD_SAMPLES) → `results/peaks/{s}_peaks.broadPeak`.
    Same but `--broad --broad-cutoff {macs2_broad_cutoff}`.
13. `seacr_bedgraph` — fragment bedGraph per sample needed by SEACR (treatment + referenced
    controls): name-sort BAM → `bedtools bamtobed -bedpe` → fragment BED → `bedtools genomecov
    -bg -g genome.chrom.sizes` → `results/seacr/bedgraph/{s}.fragments.bedgraph`.
14. `call_peaks_seacr` (TREATMENT samples with `input_control`) — `SEACR_1.3.sh
    {s}.bedgraph {igg}.bedgraph norm {stringency} results/seacr/{s}` where
    `stringency = stringent` if `peak_mode==narrow` else `relaxed`. Output
    `results/seacr/{s}.{stringency}.bed`.
15. Reproducibility (MACS2 backbone), condition = replicate group:
    - `relaxed_peaks` — MACS2 relaxed calls for 2-rep conditions, narrow/broad aware, top-N.
    - `reproducible_idr` — IDR on the pair, `--input-file-type {narrowPeak|broadPeak}`.
16. `consensus_peaks` — build fixed-width SPM consensus from MACS2 peaks per condition
    (majority / idr / single), narrow **and** broad aware (broad → summit = peak midpoint) →
    `results/consensus/consensus_peaks.{bed,saf}`.
17. `count_fragments_consensus` — `featureCounts -p --countReadPairs` over the MACS2 consensus SAF
    for all TREATMENT samples → `results/consensus/consensus_counts.txt`.
18. `seacr_consensus_peaks` — build a **variable-width reproducible union** from the per-sample
    SEACR BEDs (new `seacr_consensus.py`). Per condition, keep a SEACR peak if it overlaps peaks
    in ≥K replicates (K = `consensus_min_replicates` for ≥3 reps, 2-of-2 for exactly 2 reps,
    passthrough for 1 rep — **overlap-based, no IDR**, since SEACR has no p-value rank). Union the
    reproducible peaks across conditions, `bedtools merge` to a non-overlapping set, drop
    blacklist/off-target chroms → `results/consensus_seacr/consensus_peaks.{bed,saf}`.
19. `count_fragments_seacr_consensus` — `featureCounts -p --countReadPairs` over the SEACR
    consensus SAF for all TREATMENT samples → `results/consensus_seacr/consensus_counts.txt`.
20. `blacklist_stats` — unchanged (all samples).

### Peak-caller routing / relationship

- **MACS2 backbone**: feeds IDR-based reproducibility, the fixed-width SPM consensus, its
  featureCounts matrix, and (by default) the differential notebook. Its narrowPeak/broadPeak
  columns carry the score/summit that the SPM consensus and IDR rank on.
- **SEACR parallel track**: run per treatment sample as a CUT&RUN-native alternative call. It gets
  its own FRiP/peak-count QC **and its own count matrix** — via an overlap-based reproducible
  variable-width union (`results/consensus_seacr/`), not the IDR/fixed-width machinery (SEACR has
  no p-value rank or summit to feed those). Both matrices are available to the differential
  notebook, so differential binding can be run on MACS2 peaks, SEACR peaks, or both.

## 5. QC pipeline (`workflow/rules/qc.smk`)

Aggregate target `qc_all`. Keep the ATAC QC set, minus spike-in QC, plus SEACR:

- `deeptools_bedgraph` (RPGC), `deeptools_fragmentsize`, `deeptools_plotfingerprint`,
  `deeptools_cor_multibam` + scatter/heatmap/pca, `deeptools_gc_bias` — all samples.
- `deeptools_tss` (heatmap/profile) + `deeptools_tss_heatmap_downsample` + `tss_enrichment_score`
  — kept (informative at promoters for CUT&RUN; uses the RPGC bigWigs + GTF).
- `FRiP` — per TREATMENT sample against its MACS2 peak set (correct narrow/broad extension) **and**
  its SEACR peak set.
- `qc_relaxed_peaks` + `idr` — within-condition replicate pairs (narrow/broad aware).
- `calculate_library_complexity` — all samples (NRF/PBC1/PBC2 on pre-dedup filtered BAM).
- `reads_in_annotations` — reads in promoters vs enhancers, all samples.
- `peak_summary` — per TREATMENT sample: peak count/width + FRiP, for MACS2 and SEACR.
- `multiqc_fastqc` — FastQC-only MultiQC.
- `qc_report` — interactive `cutandrun_qc_report.html`, adapted from `build_qc_report.py`:
  drop spike-in %/norm-factor sections; add MACS2-vs-SEACR peak counts and IgG-pairing info; keep
  alignment/mito/dup/blacklist/FRiP/TSS/complexity/fragment-size/GC/correlation/PCA/fingerprint/
  reads-in-annotation/consensus panels.
- **Removed:** `spikein_qc`.

## 6. Differential binding (`cutandrun_Dx.ipynb`)

Generic R/Bioconductor DESeq2 notebook generated by an adapted
`workflow/scripts/build_diffbind_notebook.py`:

- Reads `config/samples.csv` (condition, replicate) + a consensus count matrix. Parameterized by a
  config `differential_counts` key (`macs2` | `seacr` | `both`, default `macs2`) selecting
  `results/consensus/consensus_counts.txt` and/or `results/consensus_seacr/consensus_counts.txt`;
  when `both`, the notebook runs the full analysis once per matrix.
- DESeq2 with **default median-of-ratios** normalization (no spike-in), `~ condition` (add
  `~ replicate + condition` paired design automatically when replicates are balanced).
- Contrasts: all pairwise treatment-condition comparisons (IgG excluded).
- Split consensus regions into **promoter vs distal** via `promoter_bed`; test separately.
- Per contrast: results table + significant subset (`padj<0.05 & |log2FC|>1`), MA + volcano
  plots, sample PCA (VST), ChIPseeker nearest-gene annotation → `results/diff_region/`.
- Bespoke Notch/Gviz positive-control section is **dropped** (dataset-specific). An optional Gviz
  track section over a config-provided `browser_genes` list may be added later; out of scope for v1.
- Runs in a `cutandrun_Dx` conda env (renamed from `ATACseq_Dx.yaml`, R/Bioconductor unchanged).

## 7. Config schema (`workflow/schemas/config.schema.yaml`)

Keep: `samples_table`, `adapter_r1/r2`, `align_chroms`, `keep_chroms`, `blacklist`,
`effective_genome_size`, `bin_size`, `consensus_window`, `consensus_min_replicates`,
`idr_threshold`, `idr_relaxed_pvalue`, `idr_top_n_peaks`, `keep_chroms_regex`, `macs2_genome`,
`gtf`, `promoter_bed`, `enhancer_bed`.

Rename: `human_fasta` → `genome_fasta`; `combined_index` → `genome_index`.

Remove: `spikein_fasta`, `spikein_prefix`, `spikein_pct_min`, `spikein_pct_max`.

Add:
- `max_fragment_length` (int, default 700) — Bowtie2 `-X`.
- `remove_duplicates` (bool, default true).
- `macs2_qvalue` (number, default 0.05).
- `macs2_broad_cutoff` (number, default 0.1).
- `seacr_norm` (enum `norm`/`non`, default `norm`).
- `seacr_narrow_stringency` / `seacr_broad_stringency` (enum `stringent`/`relaxed`,
  defaults `stringent` / `relaxed`).
- `differential_counts` (enum `macs2`/`seacr`/`both`, default `macs2`) — which consensus count
  matrix the differential notebook analyzes.

`peak_types` is no longer a global config knob — narrow/broad is per-sample in the sheet. QC
targets are built from `NARROW_SAMPLES`/`BROAD_SAMPLES` explicitly.

## 8. Conda envs (`workflow/envs/`)

- `snakemake.yaml` — reuse as-is (bowtie2, samtools, fastp, fastqc, bedtools, featureCounts,
  multiqc, pandas).
- `macs2.yaml`, `deeptools.yaml`, `bedtools.yaml`, `idr.yaml` — reuse as-is (deeptools already
  provides `bamCompare`).
- `seacr.yaml` — **new**: `seacr=1.3`, `bedtools`, `r-base`, `samtools`.
- `cutandrun_Dx.yaml` — renamed from `ATACSeq_Dx.yaml`; R/Bioconductor stack unchanged
  (DESeq2, ChIPseeker, GenomicRanges, ggplot2, subread).

## 9. Reference data (`ref/`)

- `genome.fa` (hg38, chr-prefixed UCSC) — user provides.
- `hg38_blacklist_regions.bed` — ship (copy from ATAC repo).
- `gtf`, `promoter_bed`, `enhancer_bed`, `hg38.2bit`, `picard.jar` — same as ATAC (documented
  download/generate steps in README).
- Built by pipeline: `genome_index` (Bowtie2), `genome.chrom.sizes`.

## 10. Scripts (`workflow/scripts/`)

- Reuse unchanged: `process_sam.py`, `blacklist-stats-script.py`, `tss_score.py`,
  `downsample_tss_matrix.py`, `build_promoter_beds.py`.
- Adapt: `consensus_peaks.py` (narrow+broad peak loading), `build_qc_report.py` (drop spike-in,
  add SEACR), `build_diffbind_notebook.py` (generic DESeq2, MACS2/SEACR matrix selectable),
  `diffbind_helpers.R` (drop spike-in).
- New: `seacr_consensus.py` (overlap-based reproducible variable-width union of SEACR peaks).
- Drop: `compute_spikein_factors.py`.

## 11. Testing / validation

- `snakemake -n` (dry-run) on a small sample sheet mirroring the reference (2 cJUN treatments +
  2 IgG controls, 1 narrow + validate a broad case) to confirm the DAG builds and all target
  paths resolve, with no reference genome required.
- Rulegraph render (`--rulegraph`) sanity check.
- Optional `.test/` mini dataset (as in the ATAC repo) for CI-style execution — out of scope for
  the first cut unless requested.

## 12. Out of scope (v1)

- SEACR feeding the **IDR** reproducibility step (SEACR consensus uses overlap-based
  reproducibility instead; SEACR still gets its own consensus + count matrix, in scope).
- Gviz browser-track section in the notebook.
- Non-human genomes (parameterized but only hg38 documented/defaulted).
- Executable `.test/` dataset + GitHub Actions CI + Docker images (can follow the ATAC repo
  later if wanted).

## 13. Notes

- Git is not initialized in this working dir; the spec is written to `docs/` but not committed.
- `condition` is the reproducibility group. Biologically distinct groups (e.g. cJUN in two
  contexts) must be given distinct `condition` labels to avoid being treated as replicates —
  documented in `config/README.md`.
