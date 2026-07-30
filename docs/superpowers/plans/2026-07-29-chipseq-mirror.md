# CUT&RUN ← ChIP-seq mirror — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Bring four ChIP-seq capabilities into the CUT&RUN workflow: a dual-control sample sheet (`input_control` + `igg_control` + `control_type`), an ENCODE IDR self-consistency/rescue-ratio reproducibility QC, a per-IP-vs-control fingerprint JSD, and an opt-in downstream module (ChIPseeker annotation+GO, HOMER motifs, peak Jaccard/overlap, deepTools peak-heatmap + gene-body metagene).

**Architecture:** Part A generalizes the control model in the pure `samplesheet.py` module + `common.smk` and updates the control-consuming rules. Parts B/C add QC rules to `qc.smk` (+ one copied script + report sections). Part D adds `rules/downstream.smk` (+ two copied scripts + two envs), wired opt-in. ChIP's DESeq2 `differential_binding` is deliberately NOT ported (the existing `diffopen` stage supersedes it).

**Tech Stack:** snakemake, samtools/MACS2/IDR/bedtools/deepTools, R/Bioconductor (ChIPseeker, clusterProfiler, DESeq2), HOMER, phantompeakqualtools, pandas, pytest.

## Global Constraints

- **Source of truth for mirrored code:** `/easley/scratch/projects/amitra/amitra2016502/snakemake_ChIP_seq` (`$CHIP`). "Copy verbatim" = byte-for-byte; "copy + rename" applies only the listed token substitutions.
- **Rename map (ChIP → CUT&RUN)** applied wherever a ChIP rule/script is adapted:
  `IP_SAMPLES`→`TREATMENT_SAMPLES`, `RATIO_SAMPLES`→`CONTROLLED_SAMPLES`, `ratio_input_bam`→`control_bam`, `config["human_fasta"]`→`config["genome_fasta"]`, `MACS2_GENOME`→`config["macs2_genome"]`, `BROAD_CUTOFF`→`config["macs2_broad_cutoff"]`, `MACS2_QVALUE`→`config["macs2_qvalue"]`, `PEAK_MODE[x]`→`SS.peak_mode(x)`, `peak_file`→`macs2_peak` (aliased). `IDR_SAMPLES`/`IDR_GROUPS`/`GROUPS` already exist and keep their names.
- **`control_type` default is `igg`** (CUT&RUN norm), unlike ChIP's `input`.
- **MACS2 is the backbone:** downstream + reproducibility operate on MACS2 per-sample peaks (`macs2_peak`) and the MACS2 consensus (`results/consensus/`), not SEACR.
- **Opt-in stages:** `downstream_all` is included in the Snakefile but NOT in `rule all` (like `diffopen_all`). `diffopen` is unaffected by this work.
- **Validation env:** `/users/jiwang1/.conda/envs/crun_smk/bin` on PATH (snakemake 9.23, pandas, jsonschema, pytest). `Rscript` available for R parse checks.
- **Dry-run harness:** `snakemake -s workflow/Snakefile -d .test -n <target>`; expect "Building DAG of jobs…" with no exceptions. `.test/config/samples.csv` after Task 2: 2-rep `cJUN` (narrow) + 1-rep `H3K27me3` (broad) + `test_igg` control, with the IgG moved into `igg_control`.
- Commit after each task; repo on `main`, no remote.

---

## File Structure

```
workflow/
├── envs/
│   ├── chipseeker.yaml           # Task 6 — copy verbatim from $CHIP
│   ├── homer.yaml                # Task 6 — copy verbatim
│   └── phantompeak.yaml          # Task 6 — update (add r-catools r-snow r-snowfall r-bitops gawk)
├── rules/
│   ├── common.smk                # Tasks 3,4,5,6 — control helpers, JSD/REPRO/downstream dirs+helpers
│   ├── cutandrun.smk             # Task 3 — igg_bam→control_bam; SEACR control = resolved
│   ├── qc.smk                    # Tasks 4,5 — fingerprint_jsd + repro_* rules
│   └── downstream.smk            # Task 7 — new (opt-in downstream_all)
├── schemas/config.schema.yaml    # Task 2 — add control_type
├── scripts/
│   ├── samplesheet.py            # Task 1 — dual-control resolution
│   ├── idr_reproducibility.py    # Task 5 — copy verbatim
│   ├── peak_annotation.R         # Task 6 — copy verbatim
│   ├── peak_overlap.py           # Task 6 — copy verbatim
│   └── build_qc_report.py        # Tasks 4,5 — fingerprint repoint + idr_reproducibility section
├── Snakefile                     # Task 7 — include downstream.smk (not in rule all)
config/{config.yaml,samples.csv,README.md}   # Task 2
.test/config/{config.yaml,samples.csv}        # Task 2
tests/{test_samplesheet.py,test_config_schema.py}   # Tasks 1,2
```

---

# PART A — DUAL-CONTROL SAMPLE SHEET

## Task 1: `samplesheet.py` dual-control resolution + tests

**Files:** Modify `workflow/scripts/samplesheet.py`; Test `tests/test_samplesheet.py`.

**Interfaces produced:** `SampleSheet.igg_control(s)`, `SampleSheet.resolved_control(s, control_type)`; `igg_control` added to `REQUIRED_COLUMNS`; `validate()` checks both control columns.

- [ ] **Step 1: Add the failing tests** (append to `tests/test_samplesheet.py`)

```python
DUAL = ("t1,cJUN,1,in1,igg1,narrow,a\n"
        "t2,cJUN,2,,igg1,narrow,b\n"
        "in1,Input,1,,,,c\n"
        "igg1,IgG,1,,,,d\n")

def _dual(tmp_path, rows):
    p = tmp_path / "s.csv"
    p.write_text("sample_id,condition,replicate,input_control,igg_control,peak_mode,notes\n" + rows)
    return ss.SampleSheet(str(p))

def test_resolved_control_igg_primary(tmp_path):
    s = _dual(tmp_path, DUAL)
    assert s.resolved_control("t1", "igg") == "igg1"   # igg primary
    assert s.resolved_control("t1", "input") == "in1"  # input primary

def test_resolved_control_fallback(tmp_path):
    s = _dual(tmp_path, DUAL)
    # t2 has only igg1; input requested -> falls back to igg1
    assert s.resolved_control("t2", "input") == "igg1"
    assert s.resolved_control("t2", "igg") == "igg1"

def test_resolved_control_none(tmp_path):
    s = _dual(tmp_path, "t1,cJUN,1,,,narrow,a\nc1,Ctl,1,,,,x\n")
    assert s.resolved_control("t1", "igg") == ""

def test_validate_bad_igg_ref(tmp_path):
    s = _dual(tmp_path, "t1,cJUN,1,,missing,narrow,a\nigg1,IgG,1,,,,d\n")
    with pytest.raises(ValueError, match="igg_control|control"):
        s.validate()
```

Also update the existing `_sheet` helper's rows (the old single-control header) — the old tests use a 6-column sheet without `igg_control`; add an empty `igg_control` so `load_samples` still finds the column. Change the old `_sheet` header to include it and insert an empty field in each row:
```python
def _sheet(tmp_path, rows):
    p = tmp_path / "samples.csv"
    p.write_text("sample_id,condition,replicate,input_control,igg_control,peak_mode,notes\n" + rows)
    return ss.SampleSheet(str(p))
```
and update `REF`/other inline rows to add the extra empty `igg_control` column (one extra comma before `peak_mode`), e.g. `"t1,cJUN,1,igg,,narrow,a\n"`.

- [ ] **Step 2: Run tests — expect failures** (`resolved_control`/`igg_control` missing, column count)

Run: `export PATH=/users/jiwang1/.conda/envs/crun_smk/bin:$PATH && python -m pytest tests/test_samplesheet.py -q`
Expected: FAIL (AttributeError `resolved_control`, or missing-column error).

- [ ] **Step 3: Implement dual control in `samplesheet.py`**

Change `REQUIRED_COLUMNS`:
```python
REQUIRED_COLUMNS = ["sample_id", "condition", "replicate", "input_control",
                    "igg_control", "peak_mode", "notes"]
```
In `load_samples`, strip `igg_control` too:
```python
    for c in ("sample_id", "condition", "input_control", "igg_control", "peak_mode"):
        df[c] = df[c].astype(str).str.strip()
```
In `SampleSheet.__init__`, add the igg map (next to `self._input`):
```python
        self._igg = dict(zip(d["sample_id"], d["igg_control"]))
```
Add methods (next to `input_control`):
```python
    def igg_control(self, sample):
        return self._igg.get(sample, "")

    def resolved_control(self, sample, control_type):
        """Effective MACS2/SEACR/track control for an IP sample: the column named
        by control_type (input|igg), else the other column, else '' (none)."""
        inp, igg = self._input.get(sample, ""), self._igg.get(sample, "")
        primary, secondary = (inp, igg) if control_type == "input" else (igg, inp)
        return primary or secondary or ""
```
In `validate()`, replace the single-control check with one covering both columns:
```python
        controls = set(self.control_samples)
        for s in self.treatment_samples:
            for col, val in (("input_control", self._input[s]), ("igg_control", self._igg[s])):
                if val and val not in controls:
                    raise ValueError(
                        f"{col} '{val}' for sample '{s}' is not an existing "
                        f"control (empty peak_mode) sample")
```
(Keep the mixed-peak_mode and valid-peak_mode checks unchanged.)

- [ ] **Step 4: Run tests — expect pass**

Run: `python -m pytest tests/test_samplesheet.py -q`
Expected: all pass (old + 4 new).

- [ ] **Step 5: Commit** — `git add -A && git commit -m "feat: samplesheet.py dual-control (input+igg) resolution + tests"`

---

## Task 2: Config `control_type` + sheet migration

**Files:** `workflow/schemas/config.schema.yaml`, `config/config.yaml`, `config/samples.csv`, `config/README.md`, `.test/config/config.yaml`, `.test/config/samples.csv`.

- [ ] **Step 1: Schema — add `control_type`** (property + required)

Add property (near `macs2_genome`):
```yaml
  control_type:
    type: string
    enum: [input, igg]
    description: Which control column drives peak calling (falls back to the other). CUT&RUN default igg.
    default: igg
```
Add `- control_type` to the `required:` list.

- [ ] **Step 2: `config/config.yaml` — add the key**

After the MACS2 block, add:
```yaml
control_type: "igg"                 # which control drives peak calling: input | igg (fallback to the other)
```

- [ ] **Step 3: Migrate `config/samples.csv` to dual control**

```csv
sample_id,condition,replicate,input_control,igg_control,peak_mode,notes
GSF2801-ChIPseq-OVCAR3-3D-IP-cJun_S4,cJUN_3D,1,,GSF2801-ChIPseq-OVCAR3-3D-IP-IgG_S5,narrow,3D-cJUN
GSF2801-ChIPseq-OVCAR3-3D-IP-IgG_S5,IgG_3D,1,,,,3D-Igg
GSF2801-ChIPseq-OVCAR3-Control-IP-cJun_S1,cJUN_Ctrl,1,,GSF2801-ChIPseq-OVCAR3-Control-IP-IgG_S2,narrow,Ctrl-cJUN
GSF2801-ChIPseq-OVCAR3-Control-IP-IgG_S2,IgG_Ctrl,1,,,,Ctrl-Igg
```

- [ ] **Step 4: `.test/config/config.yaml` + `.test/config/samples.csv`**

Add `control_type: "igg"` to `.test/config/config.yaml`. Rewrite `.test/config/samples.csv`:
```csv
sample_id,condition,replicate,input_control,igg_control,peak_mode,notes
test_cjun_r1,cJUN,1,,test_igg,narrow,narrow-rep1
test_cjun_r2,cJUN,2,,test_igg,narrow,narrow-rep2
test_k27_r1,H3K27me3,1,,test_igg,broad,broad-rep1
test_igg,IgG,1,,,,control
```

- [ ] **Step 5: `config/README.md`** — document the new column + `control_type`

In the sample-sheet table add an `igg_control` row ("`sample_id` of the matched IgG control, or empty") and note: the effective control per IP is chosen by `control_type` (input|igg, default igg) with fallback to the other column.

- [ ] **Step 6: Config schema test stays green**

Run: `export PATH=/users/jiwang1/.conda/envs/crun_smk/bin:$PATH && python -m pytest tests/test_config_schema.py -q`
Expected: 2 passed.

- [ ] **Step 7: Commit** — `git add -A && git commit -m "feat: config control_type + dual-control sample sheets"`

---

## Task 3: `common.smk` control helpers + `cutandrun.smk` consumers

**Files:** `workflow/rules/common.smk`, `workflow/rules/cutandrun.smk`.

- [ ] **Step 1: `common.smk` — replace `igg_bam` with resolved `control_bam`**

Add near the sample-set block:
```python
CONTROL_TYPE = config["control_type"]

def resolved_control(sample):
    return SS.resolved_control(sample, CONTROL_TYPE)
```
Replace the existing `igg_bam` function with:
```python
def control_bam(sample):
    """Blacklist-filtered BAM of a sample's resolved control ('' if none)."""
    c = resolved_control(sample)
    return f"{BLACKLIST_FILTERED_DIR}/{c}.nobl.bam" if c else ""
```
Change `CONTROLLED_SAMPLES` to use the resolved control:
```python
CONTROLLED_SAMPLES = [s for s in TREATMENT_SAMPLES if resolved_control(s)]
```
(The `SEACR_STRINGENT_SAMPLES`/`SEACR_RELAXED_SAMPLES` derive from `CONTROLLED_SAMPLES` and stay as-is.)

- [ ] **Step 2: `cutandrun.smk` — rename `igg_bam`→`control_bam`, SEACR control→resolved**

```bash
cd /easley/scratch/projects/amitra/amitra2016502/snakemake_CutandRun_seq
sed -i 's/igg_bam(/control_bam(/g' workflow/rules/cutandrun.smk
# SEACR control lookup: input_control -> resolved_control
sed -i 's/SS\.input_control(w\.sample)/resolved_control(w.sample)/g' workflow/rules/cutandrun.smk
grep -nE "igg_bam|SS\.input_control|control_bam|resolved_control" workflow/rules/cutandrun.smk
```
Expected: no `igg_bam` / `SS.input_control` left; `control_bam`/`resolved_control` present in the MACS2, log2, and SEACR rules.

- [ ] **Step 3: Dry-run the primary pipeline (control now resolved from igg_control)**

Run: `snakemake -s workflow/Snakefile -d .test -np results/peaks/test_cjun_r1_peaks.narrowPeak results/seacr/test_cjun_r1.stringent.bed results/log2ratio_bigwig/test_cjun_r1.log2ratio.bw 2>&1 | grep -E "macs2 callpeak|SEACR_1.3.sh|bamCompare" | head`
Expected: MACS2 shows `-c results/blacklist_filtered/test_igg.nobl.bam`; SEACR consumes `test_igg` bedgraph; bamCompare `-b2 …/test_igg.nobl.bam`. (The IgG now comes from `igg_control`.)

- [ ] **Step 4: Full primary dry-run** — `snakemake -s workflow/Snakefile -d .test -n cutandrun_all 2>&1 | grep -iE "error|exception" || echo ok`
Expected: `ok`.

- [ ] **Step 5: Commit** — `git add -A && git commit -m "feat: resolve MACS2/SEACR/track control via control_type (input+igg)"`

---

# PART C — PER-IP-vs-CONTROL FINGERPRINT JSD

## Task 4: fingerprint JSD rules + report repoint

**Files:** `workflow/rules/common.smk`, `workflow/rules/qc.smk`, `workflow/scripts/build_qc_report.py`.

- [ ] **Step 1: `common.smk` — add `JSD_DIR`**

Next to the QC dir constants: `JSD_DIR = f"{RESULT_DIR}/qc_fingerprint"`.

- [ ] **Step 2: Revert `deeptools_plotfingerprint` to plot + table only**

In `qc.smk`, remove the `metrics` output line, the `params: jsd=…` block, and the
`--outQualityMetrics {output.metrics}` + `{params.jsd}` lines added earlier — restore the rule to
its plot + `--outRawCounts` form (title comment back to "# 3. Fingerprint (signal-to-noise)").

- [ ] **Step 3: Add `fingerprint_jsd` + `fingerprint_jsd_summary`** (append to `qc.smk`)

```python
# Per-IP-vs-control fingerprint JSD (ENCODE): JS distance + % genome enriched
rule fingerprint_jsd:
    wildcard_constraints:
        sample = _alt(CONTROLLED_SAMPLES)
    input:
        ip = os.path.join(RMD_BAM_DIR, "{sample}.nobl.bam"),
        ip_bai = os.path.join(RMD_BAM_DIR, "{sample}.nobl.bam.bai"),
        control = lambda w: control_bam(w.sample),
        control_bai = lambda w: control_bam(w.sample) + ".bai"
    output:
        metrics = os.path.join(JSD_DIR, "{sample}.jsd.txt"),
        plot = os.path.join(JSD_DIR, "{sample}.fingerprint.png")
    threads: 8
    conda:
        "../envs/deeptools.yaml"
    log:
        "logs/fingerprint_jsd/{sample}.log"
    shell:
        r"""
        mkdir -p {JSD_DIR} logs/fingerprint_jsd
        plotFingerprint -b {input.ip} {input.control} \
            --labels {wildcards.sample} control \
            --JSDsample {input.control} \
            --ignoreDuplicates --skipZeros -p {threads} \
            --outQualityMetrics {output.metrics} \
            --plotFile {output.plot} > {log} 2>&1
        """


rule fingerprint_jsd_summary:
    input:
        expand(os.path.join(JSD_DIR, "{s}.jsd.txt"), s=CONTROLLED_SAMPLES)
    output:
        tsv = os.path.join(QC_DIR, "fingerprint_jsd.tsv"),
        mqc = os.path.join(QC_DIR, "fingerprint_jsd_mqc.txt")
    params:
        samples = CONTROLLED_SAMPLES,
        jsddir = JSD_DIR
    conda:
        "../envs/snakemake.yaml"
    log:
        "logs/fingerprint_jsd/summary.log"
    shell:
        r"""
        mkdir -p {QC_DIR} logs/fingerprint_jsd
        printf "# id: fingerprint_jsd\n# section_name: 'Fingerprint JSD (ENCODE)'\n# description: 'deepTools plotFingerprint IP-vs-control: JS distance and %% genome enriched.'\n# plot_type: 'table'\nSample\tJS distance\tPct genome enriched\tAUC\n" > {output.mqc}
        echo -e "sample\tjs_distance\tpct_genome_enriched\tauc" > {output.tsv}
        for s in {params.samples}; do
            awk -v s="$s" 'BEGIN{{FS="\t";OFS="\t"}}
                NR==1{{for(i=1;i<=NF;i++){{h=$i; gsub(/^ +| +$/,"",h); col[h]=i}} next}}
                $1==s{{print s, $(col["JS Distance"]), $(col["% genome enriched"]), $(col["AUC"])}}' \
                {params.jsddir}/$s.jsd.txt | tee -a {output.tsv} >> {output.mqc}
        done 2> {log}
        """
```

- [ ] **Step 4: `qc_all` + `qc_report` — swap the fingerprint source**

In `qc_all` input list, replace `os.path.join(DEEPTOOLS_DIR, "fingerprint_quality_metrics.tab")` with `os.path.join(QC_DIR, "fingerprint_jsd.tsv")`. In the `qc_report` rule inputs, make the same replacement.

- [ ] **Step 5: `build_qc_report.py` — repoint the fingerprint section**

Change the `fingerprint_quality` loader from
`read_tsv(f"{R}/deeptools/fingerprint_quality_metrics.tab")` (path check + read) to
`read_tsv(f"{R}/qc/fingerprint_jsd.tsv")` (both the `os.path.exists` guard and the read). Leave the
section key/label/flag (informational) unchanged.

- [ ] **Step 6: Smoke-test the report + dry-run**

```bash
export PATH=/users/jiwang1/.conda/envs/crun_smk/bin:$PATH
SB=$(mktemp -d)/results; mkdir -p "$SB/qc"
printf "sample\tjs_distance\tpct_genome_enriched\tauc\na\t0.35\t12\t0.7\n" > "$SB/qc/fingerprint_jsd.tsv"
python workflow/scripts/build_qc_report.py --results-dir "$SB" --out "$SB/qc/r.html" --samples a --generated x && echo "report ok"
snakemake -s workflow/Snakefile -d .test -n results/qc/fingerprint_jsd.tsv 2>&1 | grep -iE "error|exception" || echo "dryrun ok"
```
Expected: `report ok` and `dryrun ok`; `fingerprint_jsd` runs for the 3 controlled samples.

- [ ] **Step 7: Commit** — `git add -A && git commit -m "feat: per-IP-vs-control fingerprint JSD (replaces single-reference metric)"`

---

# PART B — ENCODE IDR SELF-CONSISTENCY REPRODUCIBILITY

## Task 5: reproducibility rules + summary + report section

**Files:** `workflow/rules/common.smk`, `workflow/rules/qc.smk`, `workflow/scripts/idr_reproducibility.py` (copy), `workflow/scripts/build_qc_report.py`.

- [ ] **Step 1: Copy the summary script** — `cp $CHIP/workflow/scripts/idr_reproducibility.py workflow/scripts/` (verbatim; it is a `script:`-directive consumer of `sm.params.groups/members/true_idr/repro_idr_dir`).

- [ ] **Step 2: `common.smk` — reproducibility constants + helpers**

Add:
```python
REPRO_DIR = f"{RESULT_DIR}/idr_reproducibility"

# sample -> condition (for pseudo-rep unit -> group resolution)
SAMPLE_CONDITION = {s: g for g, members in GROUPS.items() for s in members}

def idr_peak_file(group):
    """Module-B IDR consensus peaks for a 2-replicate group (already produced)."""
    ext = "broadPeak" if SS.peak_mode(GROUPS[group][0]) == "broad" else "narrowPeak"
    return f"{CONSENSUS_DIR}/idr/{group}.idr_peaks.{ext}"

# A pseudo-rep "unit" is split into two halves: a single replicate (self) or a pooled condition.
REPRO_UNITS = [f"self__{s}" for s in IDR_SAMPLES] + [f"pool__{g}" for g in IDR_GROUPS]

def _unit_group(unit):
    return SAMPLE_CONDITION[unit[len("self__"):]] if unit.startswith("self__") else unit[len("pool__"):]

def _unit_mode(unit):
    return SS.peak_mode(GROUPS[_unit_group(unit)][0])

def _unit_rep_sample(unit):
    return unit[len("self__"):] if unit.startswith("self__") else GROUPS[_unit_group(unit)][0]

REPRO_NARROW_UNITS = [u for u in REPRO_UNITS if _unit_mode(u) == "narrow"]
REPRO_BROAD_UNITS  = [u for u in REPRO_UNITS if _unit_mode(u) == "broad"]

def pseudo_source_bam(wildcards):
    u = wildcards.unit
    if u.startswith("self__"):
        return f"{BLACKLIST_FILTERED_DIR}/{u[len('self__'):]}.nobl.bam"
    return f"{REPRO_DIR}/pool/{u[len('pool__'):]}.bam"

def unit_control_bam(wildcards):
    c = resolved_control(_unit_rep_sample(wildcards.unit))
    return f"{BLACKLIST_FILTERED_DIR}/{c}.nobl.bam" if c else []

def unit_control_arg(wildcards):
    c = resolved_control(_unit_rep_sample(wildcards.unit))
    return f"-c {BLACKLIST_FILTERED_DIR}/{c}.nobl.bam" if c else ""
```

- [ ] **Step 3: `qc.smk` — reproducibility rules**

Copy the six rules `repro_pool`, `repro_split`, `repro_relaxed_narrow`, `repro_relaxed_broad`,
`repro_idr_narrow`, `repro_idr_broad`, and `idr_reproducibility_summary` from
`$CHIP/workflow/rules/qc.smk` (they are self-contained). Apply the **rename map**: `MACS2_GENOME`→
`config["macs2_genome"]`, `BROAD_CUTOFF`→`config["macs2_broad_cutoff"]` (in the `relaxed` rules'
`params.genome`/`params.broad_cutoff`). Everything else (`REPRO_DIR`, `REPRO_*_UNITS`,
`pseudo_source_bam`, `unit_control_bam`, `unit_control_arg`, `idr_peak_file`, `GROUPS`,
`IDR_GROUPS`, `config["idr_*"]`) already exists after Step 2. The `idr_reproducibility_summary`
rule's `params` are `groups=IDR_GROUPS`, `members={g: GROUPS[g] for g in IDR_GROUPS}`,
`true_idr={g: idr_peak_file(g) for g in IDR_GROUPS}`, `repro_idr_dir=os.path.join(REPRO_DIR,"idr")`,
`input.true_idr=[idr_peak_file(g) for g in IDR_GROUPS]`,
`input.counts=[os.path.join(REPRO_DIR,"idr",f"{u}.n.txt") for u in REPRO_UNITS]`,
`output.tsv=os.path.join(QC_DIR,"idr_reproducibility.tsv")`, `output.mqc=…_mqc.txt`,
`script="../scripts/idr_reproducibility.py"`.

- [ ] **Step 4: `qc_all` — add the reproducibility target (guarded)**

In `rule qc_all` input, add:
```python
        ([os.path.join(QC_DIR, "idr_reproducibility.tsv")] if IDR_GROUPS else []),
```

- [ ] **Step 5: `build_qc_report.py` — add the `idr_reproducibility` section**

In `build_data`, after the fingerprint section, add:
```python
    # ENCODE IDR reproducibility (self-consistency + rescue ratios)
    repro = read_tsv(f"{R}/qc/idr_reproducibility.tsv") if os.path.exists(
        f"{R}/qc/idr_reproducibility.tsv") else []
    for row in repro:
        st = row.get("status", "")
        summary.setdefault(row.get("condition", ""), {})
    sections["idr_reproducibility"] = {"rows": repro, "flagkey": "status", "flagspec": None}
```
Add to the `UNITS` map: `"idr_reproducibility": ("ratio · unit-free", "ratio"),`. Add to the JS
`order` array (after `xcor`/`fingerprint_quality`):
`['idr_reproducibility','IDR reproducibility (ENCODE)'],`.

- [ ] **Step 6: `qc_report` inputs — add the guarded target**

In the `qc_report` rule input, add:
```python
        ([os.path.join(QC_DIR, "idr_reproducibility.tsv")] if IDR_GROUPS else []),
```

- [ ] **Step 7: Smoke-test report + dry-run qc_all**

```bash
export PATH=/users/jiwang1/.conda/envs/crun_smk/bin:$PATH
python3 -c "import ast; ast.parse(open('workflow/scripts/idr_reproducibility.py').read()); print('py ok')"
SB=$(mktemp -d)/results; mkdir -p "$SB/qc"
printf "condition\tNt\tN1\tN2\tNp\tself_consistency_ratio\trescue_ratio\tstatus\ncJUN\t120\t110\t100\t130\t1.10\t1.08\tpass\n" > "$SB/qc/idr_reproducibility.tsv"
python workflow/scripts/build_qc_report.py --results-dir "$SB" --out "$SB/qc/r.html" --samples a --generated x >/dev/null && grep -c "IDR reproducibility (ENCODE)" "$SB/qc/r.html"
snakemake -s workflow/Snakefile -d .test -n qc_all 2>&1 | grep -iE "error|exception|ambiguous" || echo "qc_all dryrun ok"
snakemake -s workflow/Snakefile -d .test -n results/qc/idr_reproducibility.tsv 2>&1 | grep -E "repro_pool|repro_split|repro_idr_narrow|idr_reproducibility_summary" | head
```
Expected: `py ok`; grep prints `1`; `qc_all dryrun ok`; repro rules present (self__ + pool__ units for the 2-rep cJUN condition).

- [ ] **Step 8: Commit** — `git add -A && git commit -m "feat: ENCODE IDR self-consistency/rescue reproducibility QC"`

---

# PART D — DOWNSTREAM MODULE (opt-in)

## Task 6: Downstream scaffold — envs, scripts, common.smk dirs/helpers

**Files:** copy `chipseeker.yaml`, `homer.yaml`, `peak_annotation.R`, `peak_overlap.py`; edit `phantompeak.yaml`, `common.smk`.

- [ ] **Step 1: Copy envs + scripts**

```bash
cd /easley/scratch/projects/amitra/amitra2016502/snakemake_CutandRun_seq
C=/easley/scratch/projects/amitra/amitra2016502/snakemake_ChIP_seq
cp $C/workflow/envs/{chipseeker.yaml,homer.yaml} workflow/envs/
cp $C/workflow/scripts/{peak_annotation.R,peak_overlap.py} workflow/scripts/
```

- [ ] **Step 2: Update `phantompeak.yaml`** — add the run_spp.R runtime deps

Append to `workflow/envs/phantompeak.yaml` dependencies: `r-catools`, `r-snow`, `r-snowfall`,
`r-bitops`, `gawk`.

- [ ] **Step 3: `common.smk` — downstream dirs + peak-file aliases**

Add:
```python
# Downstream-analysis directories (workflow/rules/downstream.smk)
ANNOTATION_DIR = f"{RESULT_DIR}/annotation"   # ChIPseeker peak annotation + GO
MOTIF_DIR      = f"{RESULT_DIR}/motifs"        # HOMER motif enrichment
OVERLAP_DIR    = f"{RESULT_DIR}/peak_overlap"  # peak-set Jaccard / overlap

def peak_file(sample):
    """Alias used by the downstream rules/scripts: the MACS2 peak for an IP sample."""
    return macs2_peak(sample)

def all_peak_files():
    return [macs2_peak(s) for s in TREATMENT_SAMPLES]
```

- [ ] **Step 4: Verify parses**

```bash
export PATH=/users/jiwang1/.conda/envs/crun_smk/bin:$PATH
python3 -c "import ast; ast.parse(open('workflow/scripts/peak_overlap.py').read()); print('py ok')"
Rscript -e "invisible(parse('workflow/scripts/peak_annotation.R')); cat('R ok\n')" 2>&1 | tail -1
```
Expected: `py ok` and `R ok`.

- [ ] **Step 5: Commit** — `git add -A && git commit -m "chore: downstream scaffold (chipseeker/homer envs, peak scripts, dirs)"`

---

## Task 7: `downstream.smk` rules + Snakefile include

**Files:** Create `workflow/rules/downstream.smk`; modify `workflow/Snakefile`.

- [ ] **Step 1: Write `workflow/rules/downstream.smk`**

Copy `$CHIP/workflow/rules/downstream.smk` and apply: (a) delete the `differential_binding` rule
and its line in `downstream_all`; (b) apply the rename map (`IP_SAMPLES`→`TREATMENT_SAMPLES`,
`config["human_fasta"]`→`config["genome_fasta"]`, `PEAK_MODE[w.sample]`→`SS.peak_mode(w.sample)`);
(c) keep `peak_file`/`all_peak_files`/`CONSENSUS_DIR`/`BIGWIG_DIR`/`GTF_FILE`/`ANNOTATION_DIR`/
`MOTIF_DIR`/`OVERLAP_DIR`/`DEEPTOOLS_DIR`/`PEAKS_DIR`/`TMP_DIR` (all now defined). Resulting rules:
`annotate_peaks`, `motif_enrichment`, `peak_jaccard`, `peak_overlap_matrix`,
`deeptools_peak_heatmap`, `deeptools_metagene`, and `downstream_all` (without the diff-binding
target). The `downstream_all` header:
```python
def _overlap_targets():
    if len(TREATMENT_SAMPLES) >= 2:
        return [os.path.join(OVERLAP_DIR, "jaccard_matrix.tsv"),
                os.path.join(OVERLAP_DIR, "jaccard_heatmap.png")]
    return []

rule downstream_all:
    input:
        expand(os.path.join(ANNOTATION_DIR, "{sample}", "{sample}.annotation.tsv"), sample=TREATMENT_SAMPLES),
        expand(os.path.join(ANNOTATION_DIR, "{sample}", "{sample}.GO.tsv"), sample=TREATMENT_SAMPLES),
        expand(os.path.join(MOTIF_DIR, "{sample}", "knownResults.txt"), sample=TREATMENT_SAMPLES),
        _overlap_targets(),
        os.path.join(DEEPTOOLS_DIR, "Heatmap_peaks.png"),
        os.path.join(DEEPTOOLS_DIR, "Profile_peaks.png"),
        os.path.join(DEEPTOOLS_DIR, "Heatmap_genebody.png"),
        os.path.join(DEEPTOOLS_DIR, "Profile_genebody.png")
```
(`annotate_peaks` and `motif_enrichment` keep `wildcard_constraints: sample = _alt(TREATMENT_SAMPLES)`.)

- [ ] **Step 2: Snakefile — include downstream (NOT in `rule all`)**

Add after the diffopen include:
```python
include: "rules/downstream.smk"   # opt-in (downstream_all); NOT wired into rule all
```
Leave `rule all` unchanged.

- [ ] **Step 3: Dry-run `downstream_all` + confirm default excludes it**

```bash
export PATH=/users/jiwang1/.conda/envs/crun_smk/bin:$PATH
snakemake -s workflow/Snakefile -d .test -n downstream_all 2>&1 | sed -n '/Job stats:/,/total/p'
snakemake -s workflow/Snakefile -d .test -n 2>&1 | grep -c "annotate_peaks\|motif_enrichment\|downstream" || echo "0 (downstream not in default)"
```
Expected: DAG builds with `annotate_peaks`/`motif_enrichment` (×3 IP samples), `peak_jaccard`+`peak_overlap_matrix` (≥2 IP), `deeptools_peak_heatmap`, `deeptools_metagene`; no exceptions. Default-target grep prints `0`.

- [ ] **Step 4: Commit** — `git add -A && git commit -m "feat: downstream.smk (annotation/motifs/overlap/heatmap/metagene), opt-in"`

---

## Task 8: README/catalog + full validation sweep

**Files:** `README.md`, `.snakemake-workflow-catalog.yml`.

- [ ] **Step 1: README** — document the additions

- In the QC feature list add: ENCODE IDR self-consistency/rescue reproducibility and per-IP-vs-control fingerprint JSD.
- Add a **Downstream analysis** section: opt-in `downstream_all` (ChIPseeker annotation + GO, HOMER motifs, peak Jaccard/overlap, deepTools peak-heatmap + gene-body metagene), on the MACS2 backbone.
- Sample-sheet section: add `igg_control` column + `control_type` (input|igg, default igg, fallback).
- Output map: add `results/{annotation,motifs,peak_overlap,idr_reproducibility,qc_fingerprint}/`.
- Citations: add HOMER (Heinz et al. 2010) and ChIPseeker (already cited).

- [ ] **Step 2: Catalog metadata** — note `downstream_all` as an opt-in target alongside `cutandrun_all`/`qc_all`/`diffopen_all`.

- [ ] **Step 3: Full validation sweep**

```bash
export PATH=/users/jiwang1/.conda/envs/crun_smk/bin:$PATH
python -m pytest tests/ -q                                    # all unit tests
snakemake -s workflow/Snakefile -d .test -n                   # default (no diffopen/downstream)
snakemake -s workflow/Snakefile -d .test -n qc_all            # incl. fingerprint_jsd + repro
snakemake -s workflow/Snakefile -d .test -n downstream_all    # opt-in downstream
snakemake -s workflow/Snakefile -d .test -n diffopen_all      # still resolves (unaffected)
snakemake -s workflow/Snakefile -d .test cutandrun_all qc_all downstream_all diffopen_all --rulegraph dot > /tmp/rg.dot && head -1 /tmp/rg.dot
```
Expected: pytest all pass; every dry-run green (no exceptions); rulegraph prints `digraph snakemake_dag {`.

- [ ] **Step 4: Commit** — `git add -A && git commit -m "docs: ChIP-seq-mirror additions in README/catalog; full dry-run green"`

---

## Self-Review

**Spec coverage:** Part A → Tasks 1–3 (samplesheet, config, common/cutandrun consumers). Part B → Task 5 (+ common.smk helpers). Part C → Task 4. Part D → Tasks 6–7. Envs → Tasks 4/6. Report sections → Tasks 4/5. README/catalog + validation → Task 8.

**Placeholder scan:** No TBD/TODO. The two ChIP "copy + rename" tasks (5 repro rules, 7 downstream rules) name the exact source file and an explicit rename map + list the resulting rules; novel/adapted code (samplesheet dual control, fingerprint_jsd, common.smk helpers, downstream_all header, report sections, config) is given in full.

**Type/name consistency:** `resolved_control`/`control_bam`/`CONTROL_TYPE` defined in Tasks 1/3 and consumed in Tasks 3/4/5. `CONTROLLED_SAMPLES` redefined via resolved control (Task 3) and used by fingerprint_jsd (Task 4) + SEACR (existing). `JSD_DIR` (Task 4), `REPRO_DIR`/`REPRO_*_UNITS`/`idr_peak_file`/`pseudo_source_bam`/`unit_control_*` (Task 5) match the qc.smk rules that use them. `peak_file`/`all_peak_files`/`ANNOTATION_DIR`/`MOTIF_DIR`/`OVERLAP_DIR` (Task 6) match downstream.smk (Task 7). Report section keys `fingerprint_quality`(repointed)/`idr_reproducibility` match the TSV paths written by the rules. `control_type` schema key (Task 2) matches `config["control_type"]` (Task 3).
