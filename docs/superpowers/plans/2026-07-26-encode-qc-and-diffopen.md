# ENCODE QC metrics + diffopen differential-binding — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add ENCODE-standard CUT&RUN QC (phantompeakqualtools NSC/RSC + deepTools fingerprint quality metrics) and a config-driven, opt-in `diffopen` differential-binding stage mirrored from the ATAC workflow, adapted for CUT&RUN (no spike-in), run on both the MACS2 and SEACR consensus matrices.

**Architecture:** Part 1 extends `qc.smk` + `build_qc_report.py`. Part 2 adds `rules/diffopen.smk` (caller×mode wildcards) driven by `diffopen.R` (median-of-ratios / anchor-restricted / RNA-stable normalizations) with downstream annotate→enrich→bigwig→tracks→report rules copied from ATAC. The placeholder notebook is retired.

**Tech Stack:** phantompeakqualtools (run_spp.R / r-spp), deepTools, R/Bioconductor (DESeq2, apeglm, Gviz, clusterProfiler, org.Hs.eg.db), snakemake, pandas, pytest.

## Global Constraints

- **Source of truth for mirrored code:** `/easley/scratch/projects/amitra/amitra2016502/snakemake_ATACseq_spikein` (`$ATAC`). "Copy verbatim" = byte-for-byte unless a change is specified.
- **No spike-in:** never port the `spikein` diffopen mode, `spikein_anchor_shape.R`, `read_spikein_reads`, or `size_factors_spikein`.
- **diffopen is OPT-IN:** it is included in the Snakefile but NOT in `rule all`. Run via `snakemake --use-conda --cores N diffopen_all`.
- **Callers:** `diffopen_callers` (default `[macs2, seacr]`) select the count matrix — `macs2` → `results/consensus/consensus_counts.txt`, `seacr` → `results/consensus_seacr/consensus_counts.txt`. Both contain treatment samples only.
- **Modes:** `diffopen_modes` (default `[none, anchor]`; `rnastable` opt-in, needs `diffopen_rna_table`). Outputs live under `results/diffopen/<caller>/<mode>/`.
- **Sample sheet columns:** `sample_id,condition,replicate,input_control,peak_mode,notes`. Treatment = non-empty `peak_mode`; `condition` is the contrast group; `replicate` is the pairing block.
- **Validation env:** `/users/jiwang1/.conda/envs/crun_smk/bin` on PATH (snakemake 9.23, pandas, jsonschema, pytest, nbformat). `Rscript` is available for R parse checks.
- **Dry-run harness:** `snakemake -s workflow/Snakefile -d .test -n <target>`; expect "Building DAG of jobs…" with no exceptions. `.test/config/samples.csv` has treatment conditions `cJUN` (2 reps) + `H3K27me3` (1 rep) + `IgG` control → 2 conditions, valid for diffopen with `diffopen_ref_label: cJUN`.
- Commit after each task (`git add -A && git commit -m …`); repo is on branch `master`, no remote.

---

## File Structure

```
workflow/
├── envs/
│   ├── phantompeak.yaml          # Task 1 — new (run_spp.R / r-spp)
│   └── r-diffopen.yaml           # Task 4 — copy verbatim from $ATAC
├── rules/
│   ├── qc.smk                    # Task 1,2 — cross_correlation rules + fingerprint metrics
│   ├── common.smk                # Task 8 — diffopen dirs/callers/modes/helpers
│   └── diffopen.smk              # Task 9 — new (caller×mode stage)
├── scripts/
│   ├── build_qc_report.py        # Task 3 — add NSC/RSC + fingerprint-quality sections
│   ├── diffopen.R                # Task 6 — copy + adapt (read_design, drop spikein, ctcf→anchor)
│   ├── diffopen_annotate.R       # Task 4 — copy verbatim
│   ├── diffopen_enrich.R         # Task 4 — copy verbatim
│   ├── diffopen_tracks.R         # Task 4 — copy verbatim
│   ├── build_constitutive_ctcf.py# Task 4 — copy verbatim (optional refresh helper)
│   ├── build_diffopen_report.py  # Task 7 — copy + adapt (MODES, per-caller)
│   ├── build_diffbind_notebook.py# Task 4 — DELETE
│   └── diffbind_helpers.R        # Task 4 — DELETE
├── schemas/config.schema.yaml    # Task 5 — remove differential_counts; add diffopen keys
config/{config.yaml,README.md}    # Task 5,10
ref/constitutive_ctcf_hg38.bed    # Task 4 — ship (copy from $ATAC)
cutandrun_Dx.ipynb                # Task 4 — DELETE
.test/config/{config.yaml,samples.csv}  # Task 5 — diffopen keys; harness unchanged otherwise
.test/ref/constitutive_ctcf_hg38.bed    # Task 4 — placeholder
Snakefile                         # Task 9 — include diffopen.smk (not in rule all)
tests/test_config_schema.py       # Task 5 — still green after key changes
```

---

# PART 1 — ENCODE QC METRICS

## Task 1: Cross-correlation (NSC/RSC) via phantompeakqualtools

**Files:**
- Create: `workflow/envs/phantompeak.yaml`
- Modify: `workflow/rules/qc.smk` (append two rules; extend `qc_all`)

**Interfaces:**
- Produces: `results/xcor/{sample}.spp.out` (per sample), `results/xcor/xcor_summary.tsv`, `results/xcor/xcor_summary_mqc.txt`.

- [ ] **Step 1: Create `workflow/envs/phantompeak.yaml`**

```yaml
name: phantompeak
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - phantompeakqualtools=1.2.2
  - samtools=1.21
  - r-base>=4.0
  - r-spp
```

- [ ] **Step 2: Add `XCOR_DIR` to `common.smk`**

In `workflow/rules/common.smk`, next to the other QC dir constants (after `ANNOT_DIR`), add:
```python
XCOR_DIR       = f"{RESULT_DIR}/xcor"
```

- [ ] **Step 3: Append the cross-correlation rules to `qc.smk`**

Append to `workflow/rules/qc.smk`:
```python

# ── ENCODE cross-correlation (NSC / RSC / est. fragment length) ───────────
rule cross_correlation:
    input:
        bam = os.path.join(RMD_BAM_DIR, "{sample}.nobl.bam")
    output:
        spp = os.path.join(XCOR_DIR, "{sample}.spp.out"),
        pdf = os.path.join(XCOR_DIR, "{sample}.spp.pdf")
    threads: 4
    conda:
        "../envs/phantompeak.yaml"
    log:
        "logs/xcor/{sample}.log"
    shell:
        r"""
        mkdir -p {XCOR_DIR} logs/xcor
        run_spp.R -c={input.bam} -p={threads} \
            -savp={output.pdf} -out={output.spp} > {log} 2>&1
        """


# ── Aggregate NSC/RSC across samples (ENCODE: NSC>=1.05, RSC>=0.8) ─────────
rule cross_correlation_summary:
    input:
        spp = expand(os.path.join(XCOR_DIR, "{s}.spp.out"), s=SAMPLES)
    output:
        tsv = os.path.join(XCOR_DIR, "xcor_summary.tsv"),
        mqc = os.path.join(XCOR_DIR, "xcor_summary_mqc.txt")
    params:
        samples = SAMPLES,
        xcordir = XCOR_DIR
    conda:
        "../envs/snakemake.yaml"
    log:
        "logs/xcor/summary.log"
    shell:
        r"""
        mkdir -p {XCOR_DIR} logs/xcor
        printf "# id: cross_correlation\n# section_name: 'Cross-correlation (NSC/RSC)'\n# description: 'phantompeakqualtools strand cross-correlation; ENCODE targets NSC>=1.05, RSC>=0.8.'\n# plot_type: 'table'\nSample\tEst. frag len\tNSC\tRSC\tQuality tag\n" > {output.mqc}
        echo -e "sample\test_frag_len\tNSC\tRSC\tquality_tag" > {output.tsv}
        for s in {params.samples}; do
            # run_spp.R .out columns: 3=estFragLen(csv) 9=NSC 10=RSC 11=QualityTag
            awk -F'\t' -v s=$s 'BEGIN{{OFS="\t"}}
                {{split($3,a,","); print s, a[1], $9, $10, $11}}' \
                {params.xcordir}/$s.spp.out | tee -a {output.tsv} >> {output.mqc}
        done 2> {log}
        """
```

- [ ] **Step 4: Extend `qc_all` with the xcor summary**

In `qc.smk`'s `rule qc_all` input list, add:
```python
        os.path.join(XCOR_DIR, "xcor_summary.tsv"),
```

- [ ] **Step 5: Dry-run**

Run: `snakemake -s workflow/Snakefile -d .test -n results/xcor/xcor_summary.tsv`
Expected: DAG builds; `cross_correlation` (4 jobs, all samples) + `cross_correlation_summary`; no exceptions.

- [ ] **Step 6: Commit**

```bash
git add -A && git commit -m "feat: ENCODE cross-correlation QC (phantompeakqualtools NSC/RSC)"
```

---

## Task 2: deepTools fingerprint quality metrics

**Files:**
- Modify: `workflow/rules/qc.smk` (extend `deeptools_plotfingerprint`; extend `qc_all`)

**Interfaces:**
- Produces: `results/deeptools/fingerprint_quality_metrics.tab` (added output of the existing rule).

- [ ] **Step 1: Replace the `deeptools_plotfingerprint` rule body**

Replace the existing `rule deeptools_plotfingerprint` in `qc.smk` with (adds `--outQualityMetrics` and a conditional `--JSDsample` on the first control sample):
```python
# 3. Fingerprint (signal-to-noise) + ENCODE quality metrics (JS distance, % enriched)
rule deeptools_plotfingerprint:
    input:
        bams = expand(os.path.join(RMD_BAM_DIR, "{sample}.nobl.bam"), sample=SAMPLES)
    output:
        plot = os.path.join(DEEPTOOLS_DIR, "ATACseq_fingerprint.png"),
        table = os.path.join(DEEPTOOLS_DIR, "ATACseq_fingerprint.tab"),
        metrics = os.path.join(DEEPTOOLS_DIR, "fingerprint_quality_metrics.tab")
    params:
        # JS distance / % enriched need one reference; use the first IgG if present
        jsd = (f"--JSDsample {BLACKLIST_FILTERED_DIR}/{CONTROL_SAMPLES[0]}.nobl.bam"
               if CONTROL_SAMPLES else "")
    threads: 12
    conda:
        "../envs/deeptools.yaml"
    log:
        "logs/deeptools_plotfingerprint/fingerprint.log"
    shell:
        r"""
        mkdir -p {DEEPTOOLS_DIR} logs/deeptools_plotfingerprint
        plotFingerprint -p {threads} \
            -b {input.bams} \
            --ignoreDuplicates \
            -T "Fingerprints" \
            --skipZeros \
            --plotFileFormat png \
            -plot {output.plot} \
            --outRawCounts {output.table} \
            --outQualityMetrics {output.metrics} \
            {params.jsd} 2> {log}
        """
```

- [ ] **Step 2: Extend `qc_all` with the metrics table**

In `rule qc_all` input list, add:
```python
        os.path.join(DEEPTOOLS_DIR, "fingerprint_quality_metrics.tab"),
```

- [ ] **Step 3: Dry-run**

Run: `snakemake -s workflow/Snakefile -d .test -n results/deeptools/fingerprint_quality_metrics.tab`
Expected: DAG builds; `deeptools_plotfingerprint` job; no exceptions. (Optional: `-np` and confirm `--JSDsample …/test_igg.nobl.bam` and `--outQualityMetrics` appear in the command.)

- [ ] **Step 4: Commit**

```bash
git add -A && git commit -m "feat: deepTools fingerprint quality metrics (JS distance vs IgG, % enriched)"
```

---

## Task 3: QC report — NSC/RSC + fingerprint-quality sections

**Files:**
- Modify: `workflow/scripts/build_qc_report.py`
- Modify: `workflow/rules/qc.smk` (`qc_report` rule inputs)

**Interfaces:**
- Consumes: `results/xcor/xcor_summary.tsv`, `results/deeptools/fingerprint_quality_metrics.tab`.

- [ ] **Step 1: Add a NSC/RSC threshold to `build_qc_report.py`**

In the `THRESHOLDS` dict, add (higher-is-better; good, warn):
```python
    "NSC":             (1.05, 1.0, True),
    "RSC":             (0.8, 0.5, True),
```

- [ ] **Step 2: Add the two sections in `build_data`**

After the existing `sections["annotation"] = …` block in `build_data`, add:
```python
    # Cross-correlation (ENCODE NSC/RSC)
    xc = read_tsv(f"{R}/xcor/xcor_summary.tsv") if os.path.exists(
        f"{R}/xcor/xcor_summary.tsv") else []
    for row in xc:
        for k in ("NSC", "RSC"):
            try:
                summary[row["sample"]][k] = flag(float(row[k]), *THRESHOLDS[k])
            except (KeyError, ValueError):
                pass
    sections["xcor"] = {"rows": xc, "flagkey": "NSC", "flagspec": "NSC"}

    # deepTools fingerprint quality metrics (informational)
    fpq = read_tsv(f"{R}/deeptools/fingerprint_quality_metrics.tab") if os.path.exists(
        f"{R}/deeptools/fingerprint_quality_metrics.tab") else []
    sections["fingerprint_quality"] = {"rows": fpq, "flagkey": None, "flagspec": None}
```

- [ ] **Step 3: Register the two sections in the UNITS map + nav order**

In the `UNITS` dict (build_data), add:
```python
        "xcor":       ("ratio · unit-free", "ratio"),
        "fingerprint_quality": ("reads · fragments", "frag"),
```
In the JS `order` array (the numeric-sections list, near the retired `spikein` slot), add after `['complexity', …]`:
```javascript
    ['xcor','Cross-correlation (NSC/RSC)'],['fingerprint_quality','Fingerprint quality (JS distance / % enriched)'],
```

- [ ] **Step 4: Add the report inputs to the `qc_report` rule**

In `qc.smk`'s `rule qc_report` input list, add:
```python
        os.path.join(XCOR_DIR, "xcor_summary.tsv"),
        os.path.join(DEEPTOOLS_DIR, "fingerprint_quality_metrics.tab"),
```

- [ ] **Step 5: Smoke-test the report builder with synthetic xcor input**

```bash
export PATH=/users/jiwang1/.conda/envs/crun_smk/bin:$PATH
SB=$(mktemp -d)/results; mkdir -p "$SB/xcor" "$SB/deeptools" "$SB/qc"
printf "sample\test_frag_len\tNSC\tRSC\tquality_tag\na\t180\t1.20\t0.95\t2\nb\t150\t1.01\t0.40\t-1\n" > "$SB/xcor/xcor_summary.tsv"
printf "Sample\tAUC\tSynthetic AUC\tX-intercept\tElbow Point\tJS Distance\t%% genome enriched\na\t0.7\t0.5\t0.1\t0.8\t0.35\t12\n" > "$SB/deeptools/fingerprint_quality_metrics.tab"
python workflow/scripts/build_qc_report.py --results-dir "$SB" --out "$SB/qc/r.html" --samples a,b --generated x
grep -c "Cross-correlation (NSC/RSC)" "$SB/qc/r.html"
```
Expected: exit 0; grep prints `1`.

- [ ] **Step 6: Dry-run full QC + commit**

Run: `snakemake -s workflow/Snakefile -d .test -n qc_all` → no exceptions.
```bash
git add -A && git commit -m "feat: QC report NSC/RSC + fingerprint-quality sections"
```

---

# PART 2 — DIFFOPEN DIFFERENTIAL-BINDING STAGE

## Task 4: Scaffold — copy verbatim scripts/env/anchor BED; retire notebook

**Files:**
- Copy: `diffopen_annotate.R`, `diffopen_enrich.R`, `diffopen_tracks.R`, `build_constitutive_ctcf.py`, `r-diffopen.yaml`, `ref/constitutive_ctcf_hg38.bed`
- Delete: `cutandrun_Dx.ipynb`, `workflow/scripts/build_diffbind_notebook.py`, `workflow/scripts/diffbind_helpers.R`
- Create: `.test/ref/constitutive_ctcf_hg38.bed` placeholder

- [ ] **Step 1: Copy verbatim assets**

```bash
A=/easley/scratch/projects/amitra/amitra2016502/snakemake_ATACseq_spikein
cp $A/workflow/scripts/{diffopen_annotate.R,diffopen_enrich.R,diffopen_tracks.R,build_constitutive_ctcf.py} workflow/scripts/
cp $A/workflow/envs/r-diffopen.yaml workflow/envs/
cp $A/ref/constitutive_ctcf_hg38.bed ref/
: > .test/ref/constitutive_ctcf_hg38.bed
```

- [ ] **Step 2: Retire the notebook-era files**

```bash
git rm -q cutandrun_Dx.ipynb workflow/scripts/build_diffbind_notebook.py workflow/scripts/diffbind_helpers.R
rm -f tests/test_diffbind*.py 2>/dev/null || true
```

- [ ] **Step 3: Ship the anchor BED (gitignore exception)**

In `.gitignore`, add under the existing `ref/` exceptions:
```gitignore
!ref/constitutive_ctcf_hg38.bed
```

- [ ] **Step 4: Verify scripts parse**

```bash
python3 -c "import ast; ast.parse(open('workflow/scripts/build_constitutive_ctcf.py').read()); print('py ok')"
Rscript -e "invisible(lapply(c('diffopen_annotate.R','diffopen_enrich.R','diffopen_tracks.R'), function(f) parse(file.path('workflow/scripts',f)))); cat('R parse ok\n')"
```
Expected: `py ok` and `R parse ok`.

- [ ] **Step 5: Commit**

```bash
git add -A && git commit -m "chore: diffopen scaffold (R scripts, env, anchor BED); retire notebook"
```

---

## Task 5: Config schema + config.yaml + .test config (remove differential_counts, add diffopen keys)

**Files:**
- Modify: `workflow/schemas/config.schema.yaml`, `config/config.yaml`, `.test/config/config.yaml`, `config/README.md`
- Test: `tests/test_config_schema.py` (unchanged; must stay green)

- [ ] **Step 1: Schema — remove `differential_counts`, add diffopen keys**

In `workflow/schemas/config.schema.yaml`: delete the `differential_counts:` property block and its entry in `required:`. Add these properties (none are `required` — all have defaults / are opt-in):
```yaml
  diffopen_callers:
    type: array
    items: {type: string, enum: [macs2, seacr]}
    description: Which consensus count matrices differential binding runs on.
    default: [macs2, seacr]
  diffopen_modes:
    type: array
    items: {type: string, enum: [none, anchor, rnastable]}
    description: Normalizations for diffopen_all (rnastable needs diffopen_rna_table).
    default: [none, anchor]
  diffopen_ref_label:
    type: string
    description: Reference level of the condition column for the differential test.
    default: Control
  anchor_bed:
    type: string
    description: Invariant reference regions for the `anchor` normalization mode.
    default: ref/constitutive_ctcf_hg38.bed
  anchor_trim_k:
    type: number
    description: MAD multiplier for trimming anchors that move between conditions.
    default: 2.5
  anchor_trim_iter:
    type: integer
    minimum: 0
    description: Trim/re-estimate iterations for the anchor mode.
    default: 2
  anchor_min_anchors:
    type: integer
    minimum: 1
    description: Refuse to normalize on fewer than this many anchors.
    default: 100
  diffopen_min_genes:
    type: integer
    minimum: 0
    description: GO-enrichment / tracks gate; gene sets at or below this are skipped.
    default: 10
  diffopen_go_ont:
    type: string
    enum: [BP, MF, CC]
    description: GO ontology for clusterProfiler enrichment.
    default: BP
  diffopen_track_tier:
    type: string
    description: Which significance tier to draw Gviz tracks for (e.g. p01).
    default: p01
  diffopen_track_top:
    type: integer
    minimum: 1
    description: Top N up and N down regions per class for browser tracks.
    default: 5
  diffopen_rna_table:
    type: string
    description: Path to an RNA-seq DE results table (required only for rnastable mode).
```

- [ ] **Step 2: `config/config.yaml` — swap `differential_counts` for the diffopen block**

Remove the `differential_counts:` line and its comment. Replace with:
```yaml
# ── Differential binding (opt-in: `snakemake diffopen_all`) ──────────────
diffopen_callers: [macs2, seacr]    # consensus matrices to test (macs2 / seacr)
diffopen_modes: [none, anchor]      # add "rnastable" (needs diffopen_rna_table) to opt in
diffopen_ref_label: "cJUN_Ctrl"     # reference level of the `condition` column (set to yours)
anchor_bed: "ref/constitutive_ctcf_hg38.bed"   # invariant reference for the `anchor` mode
anchor_trim_k: 2.5
anchor_trim_iter: 2
anchor_min_anchors: 100
diffopen_min_genes: 10
diffopen_go_ont: "BP"
diffopen_track_tier: "p01"
diffopen_track_top: 5
# rnastable mode (opt-in): point at an RNA-seq DESeq2/edgeR table and add "rnastable" above
# diffopen_rna_table: "ref/rna_de_results.tsv"
```

- [ ] **Step 3: `.test/config/config.yaml` — same swap, `.test`-relative**

Remove `differential_counts:`. Append:
```yaml
diffopen_callers: [macs2, seacr]
diffopen_modes: [none, anchor]
diffopen_ref_label: "cJUN"          # a treatment condition in .test/config/samples.csv
anchor_bed: "ref/constitutive_ctcf_hg38.bed"
anchor_trim_k: 2.5
anchor_trim_iter: 2
anchor_min_anchors: 1               # tiny .test anchor placeholder
diffopen_min_genes: 10
diffopen_go_ont: "BP"
diffopen_track_tier: "p01"
diffopen_track_top: 5
```

- [ ] **Step 4: Update `config/README.md`**

Replace the `differential_counts` mention (line ~69) with a note that differential binding is the opt-in `diffopen_all` target (see README), keyed off `diffopen_callers` / `diffopen_modes` / `diffopen_ref_label`.

- [ ] **Step 5: Run schema tests**

```bash
export PATH=/users/jiwang1/.conda/envs/crun_smk/bin:$PATH
python -m pytest tests/test_config_schema.py -q
```
Expected: 2 passed. (The default + .test configs still validate.)

- [ ] **Step 6: Commit**

```bash
git add -A && git commit -m "feat: config — retire differential_counts, add diffopen keys"
```

---

## Task 6: Adapt `diffopen.R` (read_design, drop spikein, ctcf→anchor)

**Files:**
- Create (copy + edit): `workflow/scripts/diffopen.R`

**Interfaces:**
- Produces (per caller×mode, unchanged filenames): `differential_openness.tsv`, `diffopen_promoter.tsv`, `diffopen_enhancer.tsv`, `*_nominal_p05.tsv`, `*_nominal_p01.tsv`, `size_factors.tsv`, `run_summary.txt`, `MA_plot.png`. CLI: `--mode {none|anchor|rnastable} --counts --samples --outdir --ref-label --trim-k --trim-iter --min-anchors --promoter-bed --enhancer-bed [--anchor <bed>] [--rna-table … --models …]`.

- [ ] **Step 1: Copy the ATAC script**

```bash
cp /easley/scratch/projects/amitra/amitra2016502/snakemake_ATACseq_spikein/workflow/scripts/diffopen.R workflow/scripts/diffopen.R
```

- [ ] **Step 2: Rewrite `read_design` for the CUT&RUN sheet**

Replace the whole `read_design <- function(path, ref_label) { … }` (uses `s$type` + regex pairing) with:
```r
read_design <- function(path, ref_label) {
  s <- read.csv(path, stringsAsFactors = FALSE)
  s$peak_mode <- trimws(ifelse(is.na(s$peak_mode), "", as.character(s$peak_mode)))
  s <- s[s$peak_mode != "", ]                        # treatment samples only
  cond_raw <- s$condition
  lv <- c(ref_label, setdiff(unique(cond_raw), ref_label))     # reference level first
  condition <- factor(cond_raw, levels = lv)
  pair <- factor(as.character(s$replicate))
  data.frame(sample_id = s$sample_id, condition = condition, pair = pair,
             stringsAsFactors = FALSE)
}
```

- [ ] **Step 3: Global identifier renames (safe — code tokens only)**

```bash
sed -i 's/ctcf_overlap/anchor_overlap/g; s/size_factors_ctcf/size_factors_anchor/g' workflow/scripts/diffopen.R
```

- [ ] **Step 4: Remove the spikein IO helper**

Delete the `read_spikein_reads <- function(path) { … }` function (3–4 lines around the original line 79).

- [ ] **Step 5: Remove the spikein size-factor helper**

Delete the `size_factors_spikein <- function(spike) { … }` function (around original line 187).

- [ ] **Step 6: Fix the mode whitelist + validations in `main()`**

Replace:
```r
  if (!mode %in% c("none", "spikein", "ctcf", "rnastable"))
    stop("--mode must be one of: none, spikein, ctcf, rnastable (got '", mode, "')")
```
with:
```r
  if (!mode %in% c("none", "anchor", "rnastable"))
    stop("--mode must be one of: none, anchor, rnastable (got '", mode, "')")
```
Replace the two validation lines:
```r
  if (mode == "spikein" && is.null(a$spikein)) stop("--mode spikein requires --spikein")
  if (mode == "ctcf"    && is.null(a$ctcf))    stop("--mode ctcf requires --ctcf")
```
with:
```r
  if (mode == "anchor"  && is.null(a$anchor))  stop("--mode anchor requires --anchor")
```
And change the ">= 2 conditions" message `'type' column` → `'condition' column`.

- [ ] **Step 7: Replace the spikein + ctcf branches in the size-factor dispatch**

Replace:
```r
  } else if (mode == "spikein") {
    spike <- read_spikein_reads(a$spikein)
    stopifnot(all(samp %in% names(spike)))
    sf <- size_factors_spikein(spike[samp])
  } else if (mode == "ctcf") {
    idx <- which(anchor_overlap(fc$coords, a$ctcf))
    n_anchors <- length(idx)
    message(sprintf("CTCF-overlapping consensus peaks: %d / %d", n_anchors, nrow(counts)))
```
with:
```r
  } else if (mode == "anchor") {
    idx <- which(anchor_overlap(fc$coords, a$anchor))
    n_anchors <- length(idx)
    message(sprintf("anchor-overlapping consensus peaks: %d / %d", n_anchors, nrow(counts)))
```
(The following lines — the `if (n_anchors < …) stop(…)` and `fit <- size_factors_anchor(…)` — stay as-is after the sed rename.)

- [ ] **Step 8: Fix the run_summary label**

Replace:
```r
    if (mode == "ctcf")
      sprintf("CTCF anchors              : %d overlapping -> %d kept after invariance trim (%.1f%% dropped)",
```
with:
```r
    if (mode == "anchor")
      sprintf("anchors                   : %d overlapping -> %d kept after invariance trim (%.1f%% dropped)",
```
(Cosmetic prose elsewhere — the docstring header and the `none`/`ctcf` "should broadly agree" line — may be left; they don't affect execution. Optionally s/none` and `ctcf`/none` and `anchor`/ in the interpretation text.)

- [ ] **Step 9: Verify no stray spikein/ctcf mode tokens + R parses**

```bash
grep -nE "\"spikein\"|\"ctcf\"|a\$spikein|a\$ctcf|read_spikein_reads|size_factors_spikein|mode == .ctcf|mode == .spikein" workflow/scripts/diffopen.R || echo "clean"
Rscript -e "invisible(parse('workflow/scripts/diffopen.R')); cat('diffopen.R parse ok\n')"
```
Expected: `clean` and `diffopen.R parse ok`.

- [ ] **Step 10: Commit**

```bash
git add -A && git commit -m "feat: adapt diffopen.R for CUT&RUN (condition design, drop spikein, ctcf->anchor)"
```

---

## Task 7: Adapt `build_diffopen_report.py` (per-caller, mode set)

**Files:**
- Create (copy + edit): `workflow/scripts/build_diffopen_report.py`

- [ ] **Step 1: Copy + edit**

```bash
cp /easley/scratch/projects/amitra/amitra2016502/snakemake_ATACseq_spikein/workflow/scripts/build_diffopen_report.py workflow/scripts/build_diffopen_report.py
```
Change the `MODES` tuple from `("none", "spikein", "ctcf", "rnastable", "anchor_shape")` to:
```python
MODES = ("none", "anchor", "rnastable")
```
The report already auto-detects present modes via `os.path.isdir(join(root, m))` and takes `--diffopen-dir`, so pointing it at `results/diffopen/<caller>` (done by the rule) yields a per-caller report. No other change needed.

- [ ] **Step 2: Verify parse + smoke-test on a synthetic caller dir**

```bash
export PATH=/users/jiwang1/.conda/envs/crun_smk/bin:$PATH
python3 -c "import ast; ast.parse(open('workflow/scripts/build_diffopen_report.py').read()); print('py ok')"
D=$(mktemp -d)/macs2; mkdir -p "$D/none"
printf "normalization mode        : none\ncontrast                  : cJUN vs Control\nconsensus peaks           : 100\ndifferential (padj<0.05)  : 0\nnominal (pvalue<0.05)     : 4\n" > "$D/none/run_summary.txt"
python workflow/scripts/build_diffopen_report.py --diffopen-dir "$D" --out "$D/diffopen_report.html" && echo "report ok"
```
Expected: `py ok` and `report ok` (writes the HTML).

- [ ] **Step 3: Commit**

```bash
git add -A && git commit -m "feat: adapt build_diffopen_report.py (per-caller, none/anchor/rnastable)"
```

---

## Task 8: `common.smk` diffopen constants + helpers

**Files:**
- Modify: `workflow/rules/common.smk`

**Interfaces:**
- Produces: `DIFFOPEN_DIR`, `DIFFOPEN_CALLERS`, `DIFFOPEN_MODES`, `diffopen_counts(caller)`, `_diffopen_extra_input`, `_diffopen_track_bigwigs`, `_diffopen_track_bwdir`.

- [ ] **Step 1: Append the diffopen block to `common.smk`**

After the QC dir constants, add:
```python
# ── Differential binding (opt-in stage; see rules/diffopen.smk) ──────────
DIFFOPEN_DIR     = f"{RESULT_DIR}/diffopen"
DIFFOPEN_CALLERS = config.get("diffopen_callers", ["macs2", "seacr"])
DIFFOPEN_MODES   = config.get("diffopen_modes", ["none", "anchor"])
DIFFOPEN_COUNTS  = {
    "macs2": f"{CONSENSUS_DIR}/consensus_counts.txt",
    "seacr": f"{SEACR_CONSENSUS_DIR}/consensus_counts.txt",
}

# rnastable needs an RNA-seq DE table; fail fast at DAG-build time.
if "rnastable" in DIFFOPEN_MODES and not config.get("diffopen_rna_table"):
    raise ValueError(
        "diffopen_modes includes 'rnastable' but 'diffopen_rna_table' is unset. "
        "Point diffopen_rna_table at your RNA-seq DESeq2/edgeR results table."
    )

def diffopen_counts(caller):
    return DIFFOPEN_COUNTS[caller]

def _diffopen_extra_input(wildcards):
    """Mode-specific extra input for the `diffopen` rule."""
    if wildcards.mode == "anchor":
        return {"anchor": config.get("anchor_bed", "ref/constitutive_ctcf_hg38.bed")}
    if wildcards.mode == "rnastable":
        rna_table = config.get("diffopen_rna_table")
        if not rna_table:
            raise ValueError(
                "diffopen mode 'rnastable' requires 'diffopen_rna_table' to be set."
            )
        return {"rna_table": rna_table, "models": f"{DIFFOPEN_DIR}/gene_models.rds"}
    return {}

def _diffopen_track_bigwigs(wildcards):
    """Per-mode size-factor-scaled bigWigs for the Gviz tracks (treatment samples)."""
    return expand(
        f"{DIFFOPEN_DIR}/{wildcards.caller}/{wildcards.mode}/bigwig/{{s}}.bw",
        s=TREATMENT_SAMPLES,
    )

def _diffopen_track_bwdir(wildcards):
    return f"{DIFFOPEN_DIR}/{wildcards.caller}/{wildcards.mode}/bigwig"
```

- [ ] **Step 2: Verify common.smk still loads (skeleton dry-run)**

Run: `snakemake -s workflow/Snakefile -d .test -n 2>&1 | tail -3`
Expected: no exceptions (diffopen.smk not yet included, so still just cutandrun/qc; this confirms the new helpers parse).

- [ ] **Step 3: Commit**

```bash
git add -A && git commit -m "feat: common.smk diffopen constants + helpers (caller x mode)"
```

---

## Task 9: `diffopen.smk` rules + Snakefile include

**Files:**
- Create: `workflow/rules/diffopen.smk`
- Modify: `workflow/Snakefile` (include diffopen.smk; keep it OUT of `rule all`)

**Interfaces:**
- Produces: `diffopen_all` target + all per caller×mode outputs.

- [ ] **Step 1: Write `workflow/rules/diffopen.smk`**

```python
# ── Differential binding (OPT-IN) ────────────────────────────────────────
# Not part of `rule all`: a differential test needs >=2 conditions. Request it:
#   snakemake --use-conda --cores N diffopen_all
# Runs each configured normalization on each configured caller's consensus
# matrix. Shared config/helpers live in common.smk.

rule diffopen:
    wildcard_constraints:
        caller="macs2|seacr",
        mode="none|anchor|rnastable",
    input:
        unpack(_diffopen_extra_input),
        script="workflow/scripts/diffopen.R",
        counts=lambda w: diffopen_counts(w.caller),
        samples=config["samples_table"],
        promoter=config["promoter_bed"],
        enhancer=config["enhancer_bed"],
    output:
        table=f"{DIFFOPEN_DIR}/{{caller}}/{{mode}}/differential_openness.tsv",
        promoter=f"{DIFFOPEN_DIR}/{{caller}}/{{mode}}/diffopen_promoter.tsv",
        enhancer=f"{DIFFOPEN_DIR}/{{caller}}/{{mode}}/diffopen_enhancer.tsv",
        all_p05=f"{DIFFOPEN_DIR}/{{caller}}/{{mode}}/differential_openness_nominal_p05.tsv",
        all_p01=f"{DIFFOPEN_DIR}/{{caller}}/{{mode}}/differential_openness_nominal_p01.tsv",
        prom_p05=f"{DIFFOPEN_DIR}/{{caller}}/{{mode}}/diffopen_promoter_nominal_p05.tsv",
        prom_p01=f"{DIFFOPEN_DIR}/{{caller}}/{{mode}}/diffopen_promoter_nominal_p01.tsv",
        enh_p05=f"{DIFFOPEN_DIR}/{{caller}}/{{mode}}/diffopen_enhancer_nominal_p05.tsv",
        enh_p01=f"{DIFFOPEN_DIR}/{{caller}}/{{mode}}/diffopen_enhancer_nominal_p01.tsv",
        factors=f"{DIFFOPEN_DIR}/{{caller}}/{{mode}}/size_factors.tsv",
        summary=f"{DIFFOPEN_DIR}/{{caller}}/{{mode}}/run_summary.txt",
        ma=f"{DIFFOPEN_DIR}/{{caller}}/{{mode}}/MA_plot.png",
    params:
        outdir=lambda w: f"{DIFFOPEN_DIR}/{w.caller}/{w.mode}",
        ref_label=config.get("diffopen_ref_label", "Control"),
        trim_k=lambda w: (config.get("rnastable_trim_k", 2.5)
                          if w.mode == "rnastable" else config.get("anchor_trim_k", 2.5)),
        trim_iter=lambda w: (config.get("rnastable_trim_iter", 2)
                             if w.mode == "rnastable" else config.get("anchor_trim_iter", 2)),
        min_anchors=config.get("anchor_min_anchors", 100),
        extra=lambda w, input: (
            f"--anchor {input.anchor}" if w.mode == "anchor"
            else (
                (f"--rna-table {input.rna_table} --models {input.models} "
                 f"--tss-window {config.get('diffopen_rna_tss_window', 2000)} "
                 f"--rna-gene-col {config.get('diffopen_rna_gene_col', 'gene')} "
                 f"--rna-lfc-col {config.get('diffopen_rna_lfc_col', 'log2FoldChange')} "
                 f"--rna-padj-col {config.get('diffopen_rna_padj_col', 'padj')} "
                 f"--rna-basemean-col {config.get('diffopen_rna_basemean_col', 'baseMean')} "
                 f"--rna-basemean-min {config.get('diffopen_rna_basemean_min', 10)} "
                 f"--rna-padj-min {config.get('diffopen_rna_padj_min', 0.5)} "
                 f"--rna-lfc-max {config.get('diffopen_rna_lfc_max', 0.5)} "
                 f"--promoter-class-required "
                 f"{str(config.get('diffopen_rna_promoter_class_required', True)).lower()}")
                if w.mode == "rnastable" else ""
            )
        ),
    conda:
        "../envs/r-diffopen.yaml"
    log:
        "logs/diffopen/{caller}_{mode}.log",
    shell:
        r"""
        mkdir -p {params.outdir} logs/diffopen
        Rscript workflow/scripts/diffopen.R \
            --mode {wildcards.mode} \
            --counts {input.counts} \
            --samples {input.samples} \
            --outdir {params.outdir} \
            --ref-label '{params.ref_label}' \
            --trim-k {params.trim_k} --trim-iter {params.trim_iter} \
            --min-anchors {params.min_anchors} \
            --promoter-bed {input.promoter} --enhancer-bed {input.enhancer} \
            {params.extra} > {log} 2>&1
        """


# Parse the GTF once into a compact RDS shared by every caller×mode.
rule diffopen_gene_models:
    input:
        script="workflow/scripts/diffopen_annotate.R",
        gtf=config["gtf"],
    output:
        models=f"{DIFFOPEN_DIR}/gene_models.rds",
    conda:
        "../envs/r-diffopen.yaml"
    log:
        "logs/diffopen/gene_models.log",
    shell:
        r"""
        mkdir -p logs/diffopen
        Rscript workflow/scripts/diffopen_annotate.R --models-only \
            --gtf {input.gtf} --models {output.models} > {log} 2>&1
        """


rule diffopen_annotate:
    wildcard_constraints:
        caller="macs2|seacr",
        mode="none|anchor|rnastable",
    input:
        script="workflow/scripts/diffopen_annotate.R",
        promoter=f"{DIFFOPEN_DIR}/{{caller}}/{{mode}}/diffopen_promoter.tsv",
        enhancer=f"{DIFFOPEN_DIR}/{{caller}}/{{mode}}/diffopen_enhancer.tsv",
        gtf=config["gtf"],
        models=f"{DIFFOPEN_DIR}/gene_models.rds",
    output:
        summary=f"{DIFFOPEN_DIR}/{{caller}}/{{mode}}/genes/annotation_summary.tsv",
        universe=f"{DIFFOPEN_DIR}/{{caller}}/{{mode}}/genes/universe_genes.txt",
    params:
        indir=lambda w: f"{DIFFOPEN_DIR}/{w.caller}/{w.mode}",
        outdir=lambda w: f"{DIFFOPEN_DIR}/{w.caller}/{w.mode}/genes",
    conda:
        "../envs/r-diffopen.yaml"
    log:
        "logs/diffopen/annotate_{caller}_{mode}.log",
    shell:
        r"""
        mkdir -p logs/diffopen
        Rscript workflow/scripts/diffopen_annotate.R \
            --indir {params.indir} --gtf {input.gtf} \
            --outdir {params.outdir} --models {input.models} > {log} 2>&1
        """


rule diffopen_enrich:
    wildcard_constraints:
        caller="macs2|seacr",
        mode="none|anchor|rnastable",
    input:
        script="workflow/scripts/diffopen_enrich.R",
        universe=f"{DIFFOPEN_DIR}/{{caller}}/{{mode}}/genes/universe_genes.txt",
    output:
        summary=f"{DIFFOPEN_DIR}/{{caller}}/{{mode}}/enrichment/enrichment_summary.tsv",
    params:
        genedir=lambda w: f"{DIFFOPEN_DIR}/{w.caller}/{w.mode}/genes",
        outdir=lambda w: f"{DIFFOPEN_DIR}/{w.caller}/{w.mode}/enrichment",
        min_genes=config.get("diffopen_min_genes", 10),
        ont=config.get("diffopen_go_ont", "BP"),
    conda:
        "../envs/r-diffopen.yaml"
    log:
        "logs/diffopen/enrich_{caller}_{mode}.log",
    shell:
        r"""
        mkdir -p logs/diffopen
        Rscript workflow/scripts/diffopen_enrich.R \
            --genedir {params.genedir} --outdir {params.outdir} \
            --min-genes {params.min_genes} --ont {params.ont} > {log} 2>&1
        """


rule diffopen_bigwig:
    wildcard_constraints:
        caller="macs2|seacr",
        mode="none|anchor|rnastable",
    input:
        bam=f"{BLACKLIST_FILTERED_DIR}/{{sample}}.nobl.bam",
        bai=f"{BLACKLIST_FILTERED_DIR}/{{sample}}.nobl.bam.bai",
        factors=f"{DIFFOPEN_DIR}/{{caller}}/{{mode}}/size_factors.tsv",
    output:
        bw=f"{DIFFOPEN_DIR}/{{caller}}/{{mode}}/bigwig/{{sample}}.bw",
    params:
        bin_size=config["bin_size"],
        blacklist=config["blacklist"],
    threads: 8
    conda:
        "../envs/deeptools.yaml"
    log:
        "logs/diffopen/bigwig_{caller}_{mode}_{sample}.log",
    shell:
        r"""
        mkdir -p $(dirname {output.bw}) logs/diffopen
        SF=$(awk -F'\t' -v s="{wildcards.sample}" \
               'NR>1 && $1==s && $2+0>0 {{printf "%.10f", 1/$2}}' {input.factors})
        if [ -z "$SF" ]; then
            echo "no usable size factor for {wildcards.sample} in {input.factors}" >&2
            exit 1
        fi
        bamCoverage --bam {input.bam} \
            --scaleFactor $SF \
            --binSize {params.bin_size} \
            --numberOfProcessors {threads} \
            --extendReads \
            --blackListFileName {params.blacklist} \
            --outFileName {output.bw} > {log} 2>&1
        """


rule diffopen_tracks:
    wildcard_constraints:
        caller="macs2|seacr",
        mode="none|anchor|rnastable",
    input:
        script="workflow/scripts/diffopen_tracks.R",
        summary=f"{DIFFOPEN_DIR}/{{caller}}/{{mode}}/genes/annotation_summary.tsv",
        models=f"{DIFFOPEN_DIR}/gene_models.rds",
        bigwigs=_diffopen_track_bigwigs,
    output:
        done=touch(f"{DIFFOPEN_DIR}/{{caller}}/{{mode}}/tracks/.tracks_done"),
    params:
        genedir=lambda w: f"{DIFFOPEN_DIR}/{w.caller}/{w.mode}/genes",
        outdir=lambda w: f"{DIFFOPEN_DIR}/{w.caller}/{w.mode}/tracks",
        bwdir=_diffopen_track_bwdir,
        tier=config.get("diffopen_track_tier", "p01"),
        top=config.get("diffopen_track_top", 5),
        min_genes=config.get("diffopen_min_genes", 10),
    conda:
        "../envs/r-diffopen.yaml"
    log:
        "logs/diffopen/tracks_{caller}_{mode}.log",
    shell:
        r"""
        mkdir -p logs/diffopen {params.outdir}
        Rscript workflow/scripts/diffopen_tracks.R \
            --genedir {params.genedir} --bigwigdir {params.bwdir} \
            --models {input.models} --outdir {params.outdir} \
            --tier {params.tier} --top {params.top} \
            --min-genes {params.min_genes} > {log} 2>&1
        """


# Per-caller HTML comparing that caller's normalization modes side by side.
rule diffopen_report:
    wildcard_constraints:
        caller="macs2|seacr",
    input:
        script="workflow/scripts/build_diffopen_report.py",
        summaries=lambda w: expand(
            f"{DIFFOPEN_DIR}/{w.caller}/{{mode}}/run_summary.txt", mode=DIFFOPEN_MODES),
    output:
        html=f"{DIFFOPEN_DIR}/{{caller}}/diffopen_report.html",
    params:
        indir=lambda w: f"{DIFFOPEN_DIR}/{w.caller}",
    conda:
        "../envs/snakemake.yaml"
    log:
        "logs/diffopen/report_{caller}.log",
    shell:
        r"""
        mkdir -p logs/diffopen
        python workflow/scripts/build_diffopen_report.py \
            --diffopen-dir {params.indir} --out {output.html} > {log} 2>&1
        """


# Aggregate opt-in target: every caller × mode + downstream + per-caller report.
rule diffopen_all:
    input:
        expand(f"{DIFFOPEN_DIR}/{{caller}}/{{mode}}/differential_openness.tsv",
               caller=DIFFOPEN_CALLERS, mode=DIFFOPEN_MODES),
        expand(f"{DIFFOPEN_DIR}/{{caller}}/{{mode}}/genes/annotation_summary.tsv",
               caller=DIFFOPEN_CALLERS, mode=DIFFOPEN_MODES),
        expand(f"{DIFFOPEN_DIR}/{{caller}}/{{mode}}/enrichment/enrichment_summary.tsv",
               caller=DIFFOPEN_CALLERS, mode=DIFFOPEN_MODES),
        expand(f"{DIFFOPEN_DIR}/{{caller}}/{{mode}}/tracks/.tracks_done",
               caller=DIFFOPEN_CALLERS, mode=DIFFOPEN_MODES),
        expand(f"{DIFFOPEN_DIR}/{{caller}}/diffopen_report.html", caller=DIFFOPEN_CALLERS),
```

- [ ] **Step 2: Include diffopen.smk in the Snakefile (NOT in `rule all`)**

In `workflow/Snakefile`, add after `include: "rules/qc.smk"`:
```python
include: "rules/diffopen.smk"
```
Leave `rule all` unchanged (diffopen stays opt-in).

- [ ] **Step 3: Dry-run diffopen_all (caller×mode expansion)**

Run: `snakemake -s workflow/Snakefile -d .test -n diffopen_all 2>&1 | sed -n '/Job stats:/,/total/p'`
Expected: DAG builds; jobs include `diffopen` (2 callers × 2 modes = 4), `diffopen_gene_models` (1), `diffopen_annotate`/`diffopen_enrich`/`diffopen_tracks` (4 each), `diffopen_bigwig` (2 callers × 2 modes × 3 treatment samples = 12), `diffopen_report` (2); no exceptions/ambiguity.

- [ ] **Step 4: Confirm `rule all` did NOT pull diffopen in**

Run: `snakemake -s workflow/Snakefile -d .test -n 2>&1 | grep -c diffopen || true`
Expected: `0` (default target still just cutandrun_all + qc_all).

- [ ] **Step 5: Commit**

```bash
git add -A && git commit -m "feat: diffopen.smk stage (caller x mode) wired opt-in into the Snakefile"
```

---

## Task 10: README + docs + final validation

**Files:**
- Modify: `README.md`, `.snakemake-workflow-catalog.yml`

- [ ] **Step 1: README — replace the notebook section with the diffopen stage**

Replace the "Differential binding" subsection (which references `build_diffbind_notebook.py` / `cutandrun_Dx.ipynb`) with a description of the opt-in `diffopen_all` target: the caller (`diffopen_callers`) × mode (`none`/`anchor`/`rnastable`) matrix, `diffopen_ref_label`, the shipped constitutive-CTCF `anchor_bed`, and the outputs under `results/diffopen/<caller>/<mode>/` (+ per-caller `diffopen_report.html`). Add the two new ENCODE QC metrics (NSC/RSC, fingerprint quality) to the QC feature list and the `results/xcor/` line to the output map.

- [ ] **Step 2: Catalog metadata — mention diffopen is opt-in**

In `.snakemake-workflow-catalog.yml`, extend the description to note `diffopen_all` as an opt-in differential-binding target (alongside `cutandrun_all` / `qc_all`).

- [ ] **Step 3: Full validation sweep**

```bash
export PATH=/users/jiwang1/.conda/envs/crun_smk/bin:$PATH
python -m pytest tests/ -q                                   # all unit tests green
snakemake -s workflow/Snakefile -d .test -n                  # default (cutandrun+qc), no diffopen
snakemake -s workflow/Snakefile -d .test -n qc_all           # incl. xcor + fingerprint metrics
snakemake -s workflow/Snakefile -d .test -n diffopen_all     # opt-in stage resolves
snakemake -s workflow/Snakefile -d .test cutandrun_all qc_all diffopen_all --rulegraph dot > /tmp/rg.dot && head -1 /tmp/rg.dot
```
Expected: pytest all pass; each dry-run "Building DAG of jobs…" with no exceptions; rulegraph prints `digraph snakemake_dag {`.

- [ ] **Step 4: Commit**

```bash
git add -A && git commit -m "docs: diffopen + ENCODE QC in README/catalog; full dry-run green"
```

---

## Self-Review

**Spec coverage:**
- §1a cross-correlation → Task 1. §1b fingerprint metrics → Task 2. §1c report integration → Task 3. §1d qc_all → Tasks 1–3.
- §2a caller×mode → Tasks 5 (config), 8 (helpers), 9 (rules). §2b rules → Task 9. §2c diffopen.R adaptations → Task 6. §2d report → Task 7. §2e config → Task 5. §2f retire notebook → Task 4.
- Envs → Tasks 1 (phantompeak), 4 (r-diffopen). Anchor BED shipped → Task 4. Testing → each task's dry-run + Task 10 sweep.

**Placeholder scan:** No TBD/TODO. diffopen.R edits are given as concrete copy + sed + explicit old→new Edit blocks with the actual code; downstream R scripts and the report are "copy verbatim" + one enumerated change. No vague "handle errors" steps.

**Type/name consistency:** `XCOR_DIR` defined in common.smk (Task 1) and used in qc.smk (Tasks 1,3). `CONTROL_SAMPLES` (already in common.smk) used by the fingerprint JSD param (Task 2). diffopen helpers `diffopen_counts`, `_diffopen_extra_input`, `_diffopen_track_bigwigs`, `_diffopen_track_bwdir`, `DIFFOPEN_DIR/CALLERS/MODES` defined in Task 8 and consumed in Task 9. `diffopen.R` mode tokens (`none`/`anchor`/`rnastable`) and `--anchor` arg (Task 6) match the rule's `--mode`/`{params.extra}` (Task 9). Output filenames (`differential_openness.tsv`, `diffopen_promoter.tsv`, …) match between diffopen.R (unchanged names) and the rule outputs. `build_diffopen_report.py` `MODES` (Task 7) matches `DIFFOPEN_MODES` default (Task 5/8). Report input `run_summary.txt` per mode matches the `diffopen` rule output (Task 9).
