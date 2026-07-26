# CUT&RUN: ENCODE QC metrics + differential-binding (diffopen) stage — Design

**Date:** 2026-07-25
**Working dir:** `/easley/scratch/projects/amitra/amitra2016502/snakemake_CutandRun_seq`
**Mirrors:** `../snakemake_ATACseq_spikein` (its `rules/diffopen.smk`, `scripts/diffopen*.R`,
`scripts/build_diffopen_report.py`, `scripts/build_constitutive_ctcf.py`).

Builds on the completed CUT&RUN pipeline. Two independent additions:

- **Part 1 — ENCODE-standard QC metrics:** phantompeakqualtools cross-correlation
  (NSC / RSC / estimated fragment length / QualityTag) and deepTools fingerprint quality
  metrics (JS distance vs IgG, % genome enriched), both wired into the interactive QC report
  with ENCODE thresholds.
- **Part 2 — Differential-binding stage:** mirror the ATAC config-driven "differential
  openness" Snakemake subsystem, adapted for CUT&RUN (no spike-in), run on **both** the MACS2
  and SEACR consensus matrices, with `none` / `anchor` / `rnastable` normalizations, plus the
  downstream gene-annotation → GO-enrichment → per-mode scaled bigWigs → Gviz tracks → HTML
  report chain. Retires the placeholder `cutandrun_Dx.ipynb` notebook.

---

## Part 1 — ENCODE QC metrics (`workflow/rules/qc.smk`)

The existing QC already covers FRiP, NRF/PBC1/PBC2, fragment size, TSS, fingerprint,
correlation/PCA, GC, reads-in-annotation. The two ENCODE ChIP/CUT&RUN signal-quality metrics
it lacks:

### 1a. Cross-correlation (NSC / RSC) — phantompeakqualtools

- New env `workflow/envs/phantompeak.yaml`: `phantompeakqualtools=1.2.2` (bundles `run_spp.R`
  + `r-spp`), `samtools`, `r-base`.
- Rule `cross_correlation` (per sample, ALL samples): run `run_spp.R -c=<{sample}.nobl.bam>
  -out=<{sample}.spp.out> -savp=<{sample}.spp.pdf>` → tab file whose columns 9/10/11 are
  **NSC / RSC / QualityTag** (col 3 = estimated fragment length, comma-separated top-3).
- Rule `cross_correlation_summary`: parse every `{sample}.spp.out` → `xcor_summary.tsv`
  (sample, est_frag_len, NSC, RSC, quality_tag) + a MultiQC-style `_mqc.txt`.
- Directory: `results/xcor/`.

### 1b. Fingerprint quality metrics — deepTools

- Extend the existing `deeptools_plotfingerprint` rule (do not add a new rule) with
  `--outQualityMetrics {DEEPTOOLS_DIR}/fingerprint_quality_metrics.tab` and, when ≥1 control
  (IgG) sample exists, `--JSDsample <first-control>.nobl.bam` (deepTools takes ONE JSD
  reference; using the first control is the documented simplification — noted in the report).
  The quality-metrics table gains AUC, Synthetic AUC, X-intercept, Elbow Point, and — with
  `--JSDsample` — **JS Distance** and **% genome enriched** per sample.

### 1c. QC report integration (`workflow/scripts/build_qc_report.py`)

- Add a **Cross-correlation (NSC/RSC)** section from `results/xcor/xcor_summary.tsv`, flagged
  against ENCODE thresholds: **NSC ≥ 1.05** (warn < 1.05, fail < 1.0) and **RSC ≥ 0.8**
  (warn < 0.8, fail < 0.5).
- Add a **Fingerprint quality** section from `fingerprint_quality_metrics.tab` (JS distance,
  % genome enriched, elbow point); no hard flag (informational), with a note that JS distance
  is relative to the first IgG.
- Both loaders wrapped so a missing file skips the panel (report still renders).

### 1d. `qc_all`

Append `results/xcor/xcor_summary.tsv` and `.../fingerprint_quality_metrics.tab` to `qc_all`
(the report rule already gates on the sections it renders).

---

## Part 2 — Differential-binding (`diffopen`) stage

Mirror the ATAC `rules/diffopen.smk` subsystem. **Opt-in** (not in `rule all`): a differential
test needs ≥2 conditions; request it with `snakemake --use-conda --cores N diffopen_all`. The
ATAC name `diffopen` is kept verbatim (per decision).

### 2a. Two dimensions: caller × mode

- **Callers** (`diffopen_callers`, default `[macs2, seacr]`) select the consensus count matrix:
  `macs2` → `results/consensus/consensus_counts.txt`, `seacr` → `results/consensus_seacr/
  consensus_counts.txt`. Both matrices contain **treatment samples only** (IgG excluded),
  which is what the differential test needs.
- **Modes** (`diffopen_modes`, default `[none, anchor]`; `rnastable` opt-in):
  - `none` — DESeq2 median-of-ratios over all consensus peaks (baseline).
  - `anchor` — median-of-ratios restricted to a fixed **invariant reference BED**
    (`anchor_bed`, default the shipped constitutive-CTCF set), with iterative trimming of
    anchors that move between conditions. Spike-in-free. Generalized from ATAC's `ctcf` mode
    (rename `ctcf`→`anchor` throughout). Meaningful when the target has signal at the anchor
    regions (CTCF/cohesin CUT&RUN, or open-chromatin-correlated marks); documented as such.
  - `rnastable` — median-of-ratios restricted to promoter-class peaks at RNA-seq-stable genes
    (needs `diffopen_rna_table`). Opt-in, spike-in-free, ported unchanged.
  - **Dropped:** `spikein` mode and the `anchor_shape` hybrid (both spike-in-dependent).

Every output lives under `results/diffopen/<caller>/<mode>/…`.

### 2b. Rules (`workflow/rules/diffopen.smk`)

Adapted from ATAC's `diffopen.smk`, adding a `{caller}` wildcard and dropping the two spike-in
rules:

- `diffopen` (wildcards `caller ∈ {macs2,seacr}`, `mode ∈ {none,anchor,rnastable}`): runs
  `diffopen.R`; input counts chosen by caller; mode-specific extra input via
  `_diffopen_extra_input` (anchor→`anchor_bed`, rnastable→rna table + gene-models RDS). Outputs
  the DA table, promoter/enhancer splits, nominal-p05/p01 subsets, `size_factors.tsv`,
  `run_summary.txt`, `MA_plot.png`.
- `diffopen_gene_models` (caller/mode-independent): parse the GTF once → `gene_models.rds`.
- `diffopen_annotate` (per caller×mode): nearest-TSS gene assignment → `genes/`.
- `diffopen_enrich` (per caller×mode): offline GO via clusterProfiler → `enrichment/`.
- `diffopen_bigwig` (per caller×mode×sample): bamCoverage scaled by `1/size_factor` (no RPGC).
- `diffopen_tracks` (per caller×mode): Gviz tracks for top up/down genes.
- `diffopen_report` (per caller): HTML comparing that caller's modes side by side →
  `results/diffopen/<caller>/diffopen_report.html`.
- `diffopen_all`: expand every downstream target over `diffopen_callers × diffopen_modes`.

`build_constitutive_ctcf.py` is copied over (regeneration helper), but the pre-built
`ref/constitutive_ctcf_hg38.bed` is **shipped** (copied from the ATAC `ref/`), so the `anchor`
mode runs without any ENCODE downloads.

### 2c. `diffopen.R` adaptations (`workflow/scripts/diffopen.R`)

Copy from ATAC, then:

1. **`read_design`** — ATAC reads the `type` column and derives pairing from a
   `_<n>_S<lane>` regex. Rewrite for the CUT&RUN sheet: filter to **treatment rows**
   (`peak_mode` non-empty), use the **`condition`** column for the factor (reference =
   `--ref-label`), and use the explicit **`replicate`** column for pairing:
   ```r
   read_design <- function(path, ref_label) {
     s <- read.csv(path, stringsAsFactors = FALSE)
     s$peak_mode <- trimws(ifelse(is.na(s$peak_mode), "", as.character(s$peak_mode)))
     s <- s[s$peak_mode != "", ]                        # treatment samples only
     cond_raw <- s$condition
     lv <- c(ref_label, setdiff(unique(cond_raw), ref_label))
     condition <- factor(cond_raw, levels = lv)
     pair <- factor(as.character(s$replicate))
     data.frame(sample_id = s$sample_id, condition = condition, pair = pair,
                stringsAsFactors = FALSE)
   }
   ```
2. **Remove the `spikein` branch** in `main()` and the `size_factors_spikein` /
   `read_spikein_reads` helpers.
3. **Rename `ctcf` → `anchor`**: the `--mode` token and validation, the `--ctcf` arg →
   `--anchor`, `ctcf_overlap`→`anchor_overlap`, `size_factors_ctcf`→`size_factors_anchor`, and
   the `"CTCF-…"` messages → `"anchor-…"`. The underlying median-of-ratios-over-anchor math is
   unchanged.
4. Fix the `main()` "need >= 2 conditions … 'type' column" message to say `condition`.
5. Keep `none`, `rnastable`, `fit_class`, `median_of_ratios`, `classify_peaks`,
   `size_factors_deseq2`, `write_results`, the MA plot, and all argument plumbing unchanged.

`diffopen_annotate.R`, `diffopen_enrich.R`, `diffopen_tracks.R` are **copied verbatim**
(confirmed spike-in-free). `spikein_anchor_shape.R` is **not** ported.

### 2d. `build_diffopen_report.py`

Copied from ATAC and adapted: it enumerates the modes present under a caller's
`results/diffopen/<caller>/` and renders per-caller. Drop any `anchor_shape`/spike-in-specific
handling; the schema for the ported modes (`none`/`anchor`/`rnastable`) is uniform.

### 2e. Config (`config/config.yaml` + schema)

Remove `differential_counts` (the retired notebook's key). Add:

- `diffopen_callers: [macs2, seacr]` — which consensus matrices to test.
- `diffopen_modes: [none, anchor]` — normalizations for `diffopen_all` (add `rnastable` to
  opt in; requires `diffopen_rna_table`).
- `diffopen_ref_label: "cJUN_Ctrl"` — reference level of the `condition` column (user sets to
  one of their conditions; the shipped value matches the example sheet).
- `anchor_bed: "ref/constitutive_ctcf_hg38.bed"` — invariant reference regions.
- `anchor_trim_k: 2.5`, `anchor_trim_iter: 2`, `anchor_min_anchors: 100`.
- `diffopen_min_genes: 10`, `diffopen_go_ont: "BP"`, `diffopen_track_tier: "p01"`,
  `diffopen_track_top: 5`.
- rnastable (opt-in): `diffopen_rna_table` (unset), `diffopen_rna_tss_window: 2000`,
  `diffopen_rna_gene_col/lfc_col/padj_col/basemean_col`, thresholds, `diffopen_rna_min_anchors`.

Schema validation for `rnastable` mirrors ATAC's fast-fail: if `rnastable` ∈ `diffopen_modes`
and `diffopen_rna_table` is unset, raise at DAG-build time.

### 2f. Retire the notebook

Delete all three notebook-era files — `cutandrun_Dx.ipynb`,
`workflow/scripts/build_diffbind_notebook.py`, and `workflow/scripts/diffbind_helpers.R` (its
`run_deseq2_contrast` is superseded by `diffopen.R`). Remove the notebook from
`test_*`/tests if referenced, and update the README's differential-analysis section to point
at the `diffopen_all` target instead of `build_diffbind_notebook.py`.

---

## Envs (`workflow/envs/`)

- **New** `phantompeak.yaml` — phantompeakqualtools + r-spp + samtools.
- **New** `r-diffopen.yaml` — copied verbatim from ATAC (DESeq2, apeglm, GenomicRanges,
  rtracklayer, Gviz, GenomicFeatures, org.Hs.eg.db, clusterProfiler, ggplot2).
- `deeptools.yaml` reused for `diffopen_bigwig` and the fingerprint metrics.

## Reference data (`ref/`)

- **Ship** `ref/constitutive_ctcf_hg38.bed` (copied from ATAC `ref/`; ~752 KB, tracked with a
  `.gitignore` exception like the other shipped BEDs).

## Scripts summary

- New: (none purely new besides the report/env) — cross-correlation uses `run_spp.R` from the
  env; the summary is a small inline shell/awk in the rule.
- Copied verbatim: `diffopen_annotate.R`, `diffopen_enrich.R`, `diffopen_tracks.R`,
  `build_constitutive_ctcf.py`.
- Copied + adapted: `diffopen.R` (§2c), `build_diffopen_report.py` (§2d), `build_qc_report.py`
  (§1c).
- Removed: `build_diffbind_notebook.py`, `diffbind_helpers.R`, `cutandrun_Dx.ipynb`.

## Testing

- Extend `.test/config/config.yaml` with the diffopen keys and `diffopen_ref_label` set to a
  `.test` condition (`cJUN`), plus placeholder `.test/ref/constitutive_ctcf_hg38.bed`.
- `snakemake -n -d .test` dry-runs: `qc_all` (now incl. xcor + fingerprint metrics) and
  `diffopen_all` (caller×mode expansion resolves; no exceptions).
- Unit tests: the cross-correlation summary parser and the report's new section loaders are
  small; add a pytest smoke test of `build_qc_report.py` with a synthetic `xcor_summary.tsv`
  to confirm the NSC/RSC section renders. `diffopen.R` is validated by dry-run + `Rscript`
  parse (no R execution without real BAMs/DESeq2).
- The `.test` diffopen dry-run needs ≥2 conditions in `.test/config/samples.csv` — the current
  harness has `cJUN` (2 reps) + `H3K27me3` (1 rep) + IgG, i.e. 2 treatment conditions, so
  `diffopen_ref_label: cJUN` gives a valid 2-condition contrast (cJUN vs H3K27me3). (This is a
  DAG-validity fixture, not a biologically meaningful contrast.)

## Out of scope (v1)

- Spike-in `anchor_shape` hybrid (Method 6) and the `spikein` normalization mode.
- Real-data execution of diffopen (needs DESeq2/Gviz envs + real BAMs).
- Per-treatment-vs-own-IgG fingerprint JSD (uses a single representative IgG instead).
- Regenerating the constitutive-CTCF BED (ship the pre-built one; builder included for
  optional refresh).

## Notes

- `condition` remains the reproducibility/contrast group. `diffopen_ref_label` must be one of
  the treatment `condition` labels.
- Running diffopen on both callers × both default modes = 4 DA runs + downstream; it is opt-in,
  so cost is incurred only when explicitly requested.
