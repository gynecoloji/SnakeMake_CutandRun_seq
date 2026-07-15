# CUT&RUN Snakemake Workflow Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a paired-end CUT&RUN Snakemake workflow that mirrors the sibling ATAC-seq spike-in pipeline, but without spike-in, with per-sample narrow/broad peak mode and per-sample IgG/Input controls, using both MACS2 and SEACR (each producing a consensus count matrix), IgG-subtracted signal tracks, and a generic DESeq2 differential notebook.

**Architecture:** Standard Snakemake Workflow Catalog layout (`workflow/Snakefile` + `rules/common.smk`, `rules/cutandrun.smk`, `rules/qc.smk`; params validated against `workflow/schemas/config.schema.yaml`; per-rule conda envs). Sample-sheet parsing/validation is factored into a pure, unit-tested Python module `workflow/scripts/samplesheet.py` that `common.smk` imports. Most rules are adapted verbatim from `../snakemake_ATACseq_spikein`; novel logic (sample sets, MACS2 control routing, SEACR, log2 tracks, SEACR consensus) is written fresh.

**Tech Stack:** Snakemake ≥8, Bowtie2, samtools, Picard, fastp/FastQC, MACS2, SEACR 1.3, deepTools, bedtools, IDR, featureCounts (Subread), R/Bioconductor (DESeq2, ChIPseeker), pandas, pytest.

## Global Constraints

- **Source of truth for mirrored code:** `/easley/scratch/projects/amitra/amitra2016502/snakemake_ATACseq_spikein` (referred to below as `$ATAC`). "Copy verbatim" means byte-for-byte unless a change is specified.
- **No spike-in anywhere.** Do not carry over `spikein_*` config keys, the combined-genome index, spike-in read splitting, `compute_spikein_factors.py`, `create_spikein_bigwig`, or `spikein_qc`.
- **Config is validated** against `workflow/schemas/config.schema.yaml` at parse time (the schema is the single source of truth; the catalog renders it). Every config key used by a rule MUST exist in the schema with a default.
- **Sample sheet columns (exact):** `sample_id,condition,replicate,input_control,peak_mode,notes`.
- **Control vs treatment:** a row is a CONTROL (IgG/Input) iff `peak_mode` is empty; otherwise TREATMENT. `peak_mode` ∈ {`narrow`,`broad`,empty}.
- **`condition` is the reproducibility group** for treatment samples; reproducibility method by replicate count: ≥3 → majority, ==2 → IDR, ==1 → single (mirrors ATAC).
- **Read naming:** `data/<sample_id>_R1_001.fastq.gz` / `_R2_001.fastq.gz`.
- **Genome default:** hg38, chr-prefixed UCSC. Config key is `genome_fasta` (not `human_fasta`); index prefix is `genome_index` (not `combined_index`).
- **Git is NOT initialized** in the working dir. "Commit" steps below run `git add`/`git commit`; if `git init` has not been done, the executor should either `git init` first (once, in Task 1) or skip commits. Task 1 Step 0 initializes git.
- **Dry-run test harness:** `.test/` holds placeholder config/samples/ref/data so `snakemake -n -d .test <target>` validates the DAG without real genomes. Placeholders may be empty files; dry-run does not read their contents. The default validation command for a target `T` is: `snakemake -s workflow/Snakefile -d .test -n <T>` → expect "Building DAG of jobs..." and a job listing with **no** exceptions.
- **pytest** runs from repo root: `pytest tests/ -v`.

**`.test/config/samples.csv` used by every dry-run:**
```csv
sample_id,condition,replicate,input_control,peak_mode,notes
test_cjun_r1,cJUN,1,test_igg,narrow,narrow-rep1
test_cjun_r2,cJUN,2,test_igg,narrow,narrow-rep2
test_k27_r1,H3K27me3,1,test_igg,broad,broad-rep1
test_igg,IgG,1,,,control
```
This exercises: narrow (2-rep → IDR), broad (1-rep → single), and a shared IgG control.

---

## File Structure

```
snakemake_CutandRun_seq/
├── config/
│   ├── config.yaml              # Task 2 — default params (schema-valid)
│   ├── samples.csv              # Task 2 — example sheet (reference format)
│   └── README.md                # Task 2 — sample-sheet + config reference
├── workflow/
│   ├── Snakefile                # Task 4 — entry point (targets: cutandrun_all, qc_all)
│   ├── rules/
│   │   ├── common.smk           # Task 4 — imports samplesheet.py, dirs, GROUPS, helpers
│   │   ├── cutandrun.smk        # Tasks 5–11 — primary pipeline
│   │   └── qc.smk               # Tasks 12–13 — QC pipeline
│   ├── schemas/
│   │   └── config.schema.yaml   # Task 2 — parameter schema (single source of truth)
│   ├── scripts/
│   │   ├── samplesheet.py       # Task 3 — pure sample-sheet logic (unit-tested)
│   │   ├── process_sam.py       # Task 1 — copy verbatim from $ATAC
│   │   ├── blacklist-stats-script.py   # Task 1 — copy verbatim
│   │   ├── tss_score.py         # Task 1 — copy verbatim
│   │   ├── downsample_tss_matrix.py    # Task 1 — copy verbatim
│   │   ├── build_promoter_beds.py      # Task 1 — copy verbatim
│   │   ├── consensus_peaks.py   # Task 9 — adapt (narrow+broad loading)
│   │   ├── seacr_consensus.py   # Task 10 — new
│   │   ├── build_qc_report.py   # Task 13 — adapt (drop spike-in, add SEACR)
│   │   ├── build_diffbind_notebook.py  # Task 14 — adapt (generic DESeq2)
│   │   └── diffbind_helpers.R   # Task 14 — adapt (drop spike-in)
│   └── envs/
│       ├── snakemake.yaml       # Task 1 — copy verbatim
│       ├── macs2.yaml           # Task 1 — copy verbatim
│       ├── deeptools.yaml       # Task 1 — copy verbatim
│       ├── bedtools.yaml        # Task 1 — copy verbatim
│       ├── idr.yaml             # Task 1 — copy verbatim
│       ├── seacr.yaml           # Task 1 — new
│       └── cutandrun_Dx.yaml    # Task 1 — rename of ATACSeq_Dx.yaml
├── ref/
│   ├── hg38_blacklist_regions.bed      # Task 1 — copy from $ATAC
│   ├── promoter_chr1-22X.bed           # Task 1 — copy from $ATAC
│   ├── enhancer_chr1-22X.bed           # Task 1 — copy from $ATAC
│   └── .gitkeep
├── data/.gitkeep                # Task 1
├── tests/
│   ├── test_samplesheet.py      # Task 3
│   ├── test_consensus_peaks.py  # Task 9
│   └── test_seacr_consensus.py  # Task 10
├── .test/                       # Task 1 — dry-run harness
│   ├── config/{config.yaml,samples.csv}
│   ├── ref/{genome.fa,hg38_blacklist_regions.bed,gencode.gtf,hg38.2bit,promoter_chr1-22X.bed,enhancer_chr1-22X.bed}
│   └── data/*.fastq.gz          # 8 empty placeholders
├── .snakemake-workflow-catalog.yml     # Task 15
├── .gitignore                  # Task 1
├── run_pipeline.sh             # Task 15
└── README.md                   # Task 15
```

---

## Task 1: Project scaffold, envs, shipped scripts/refs, and `.test/` harness

**Files:**
- Create dirs: `config/ workflow/{rules,schemas,scripts,envs} ref/ data/ tests/ .test/{config,ref,data}`
- Copy: five envs + five scripts + three ref BEDs from `$ATAC`
- Create: `workflow/envs/seacr.yaml`, `workflow/envs/cutandrun_Dx.yaml`, `.gitignore`, `.gitkeep`s, `.test/` placeholders

**Interfaces:**
- Produces: the directory layout and all verbatim-copied assets that later tasks import/reference.

- [ ] **Step 0: Init git**

```bash
cd /easley/scratch/projects/amitra/amitra2016502/snakemake_CutandRun_seq
git init -q
```

- [ ] **Step 1: Make directories**

```bash
mkdir -p config workflow/rules workflow/schemas workflow/scripts workflow/envs \
        ref data tests .test/config .test/ref .test/data \
        docs/superpowers/specs docs/superpowers/plans
```

- [ ] **Step 2: Copy verbatim assets from the ATAC repo**

```bash
A=/easley/scratch/projects/amitra/amitra2016502/snakemake_ATACseq_spikein
# envs reused unchanged
cp $A/workflow/envs/{snakemake.yaml,macs2.yaml,deeptools.yaml,bedtools.yaml,idr.yaml} workflow/envs/
# Dx env renamed
cp $A/workflow/envs/ATACSeq_Dx.yaml workflow/envs/cutandrun_Dx.yaml
# scripts reused unchanged
cp $A/workflow/scripts/{process_sam.py,blacklist-stats-script.py,tss_score.py,downsample_tss_matrix.py,build_promoter_beds.py} workflow/scripts/
# shipped reference BEDs
cp $A/ref/{hg38_blacklist_regions.bed,promoter_chr1-22X.bed,enhancer_chr1-22X.bed} ref/
touch ref/.gitkeep data/.gitkeep
```

- [ ] **Step 3: Fix the `name:` in the renamed Dx env**

Edit `workflow/envs/cutandrun_Dx.yaml` line 1: change `name: ATACseq_Dx` → `name: cutandrun_Dx`. Leave the `prefix:` line and all dependencies unchanged.

- [ ] **Step 4: Create `workflow/envs/seacr.yaml`**

```yaml
name: seacr
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - seacr=1.3
  - bedtools=2.31.1
  - samtools=1.21
  - r-base>=4.0
```

- [ ] **Step 5: Create `.gitignore`**

```gitignore
results/
logs/
.snakemake/
ref/*
!ref/.gitkeep
!ref/hg38_blacklist_regions.bed
!ref/promoter_chr1-22X.bed
!ref/enhancer_chr1-22X.bed
data/*
!data/.gitkeep
__pycache__/
*.pyc
.pytest_cache/
```

- [ ] **Step 6: Create the `.test/` dry-run harness**

Placeholder reference + data files (empty is fine — dry-run never reads them):
```bash
touch .test/ref/genome.fa .test/ref/hg38_blacklist_regions.bed \
      .test/ref/gencode.gtf .test/ref/hg38.2bit \
      .test/ref/promoter_chr1-22X.bed .test/ref/enhancer_chr1-22X.bed
for s in test_cjun_r1 test_cjun_r2 test_k27_r1 test_igg; do
  : > .test/data/${s}_R1_001.fastq.gz
  : > .test/data/${s}_R2_001.fastq.gz
done
```

Create `.test/config/samples.csv` with the exact content from Global Constraints.

`.test/config/config.yaml` is created in Task 2 Step 4 (needs the schema keys).

- [ ] **Step 7: Verify copies imported cleanly**

Run: `python -c "import ast; [ast.parse(open('workflow/scripts/'+f).read()) for f in ['process_sam.py','tss_score.py','downsample_tss_matrix.py','build_promoter_beds.py']]; print('ok')"`
Expected: `ok`

- [ ] **Step 8: Commit**

```bash
git add -A
git commit -m "chore: scaffold CUT&RUN workflow (envs, shipped scripts/refs, .test harness)"
```

---

## Task 2: Config schema, default config, example sample sheet, config README

**Files:**
- Create: `workflow/schemas/config.schema.yaml`, `config/config.yaml`, `config/samples.csv`, `config/README.md`, `.test/config/config.yaml`
- Test: `tests/test_config_schema.py`

**Interfaces:**
- Produces: config keys consumed everywhere. Exact keys (with defaults): `samples_table`, `adapter_r1`, `adapter_r2`, `genome_fasta`, `genome_index`, `align_chroms`, `keep_chroms`, `blacklist`, `effective_genome_size`, `bin_size`, `max_fragment_length`, `remove_duplicates`, `macs2_genome`, `macs2_qvalue`, `macs2_broad_cutoff`, `seacr_norm`, `seacr_narrow_stringency`, `seacr_broad_stringency`, `consensus_window`, `consensus_min_replicates`, `idr_threshold`, `idr_relaxed_pvalue`, `idr_top_n_peaks`, `keep_chroms_regex`, `differential_counts`, `gtf`, `promoter_bed`, `enhancer_bed`.

- [ ] **Step 1: Write `workflow/schemas/config.schema.yaml`**

Start from `$ATAC/workflow/schemas/config.schema.yaml`. Remove `spikein_fasta`, `spikein_prefix`, `combined_index`, `spikein_pct_min`, `spikein_pct_max`, and the `peak_types` block. Rename `human_fasta`→`genome_fasta`, `combined_index`→`genome_index`. Add the new keys. Full file:

```yaml
"$schema": "http://json-schema.org/draft-07/schema#"
title: CUT&RUN peak-calling + consensus + QC — configuration
description: >-
  Parameters for the CUT&RUN pipeline (no spike-in). File paths are relative to
  the working directory (the repo root).
type: object

properties:

  samples_table:
    type: string
    description: >-
      Path to the sample sheet CSV. Columns: sample_id, condition, replicate,
      input_control, peak_mode, notes.
    default: config/samples.csv

  adapter_r1:
    type: string
    description: Optional explicit R1 adapter that OVERRIDES fastp auto-detection.
  adapter_r2:
    type: string
    description: Optional explicit R2 adapter (used with adapter_r1).

  genome_fasta:
    type: string
    description: Genome FASTA, chr-prefixed UCSC (hg38) to match the blacklist.
    default: ref/genome.fa
  genome_index:
    type: string
    description: Bowtie2 index prefix, built automatically by build_genome_index.
    default: ref/genome/genome
  align_chroms:
    type: array
    items: {type: string}
    description: Chromosomes kept when building the index ([] = keep all).
  keep_chroms:
    type: array
    items: {type: string}
    description: Analysis keep-set for the final BAM (mito-% QC recorded first).
  blacklist:
    type: string
    description: ENCODE-style blacklist BED (chr-prefixed).
    default: ref/hg38_blacklist_regions.bed

  effective_genome_size:
    type: integer
    minimum: 1
    description: Effective genome size for deepTools RPGC normalization (hg38).
    default: 2913022398
  bin_size:
    type: integer
    minimum: 1
    description: bigWig bin size in bp.
    default: 25
  max_fragment_length:
    type: integer
    minimum: 1
    description: Bowtie2 -X maximum fragment length (CUT&RUN default 700).
    default: 700
  remove_duplicates:
    type: boolean
    description: Picard REMOVE_DUPLICATES. false keeps (marks) dups for low-input CUT&RUN.
    default: true

  macs2_genome:
    type: string
    description: MACS2 -g effective genome preset (hs, mm, ce, dm).
    default: hs
  macs2_qvalue:
    type: number
    minimum: 0
    maximum: 1
    description: MACS2 -q cutoff for narrow/broad peak calls.
    default: 0.05
  macs2_broad_cutoff:
    type: number
    minimum: 0
    maximum: 1
    description: MACS2 --broad-cutoff for broad peak calls.
    default: 0.1

  seacr_norm:
    type: string
    enum: [norm, non]
    description: SEACR normalization mode (norm when using an IgG control).
    default: norm
  seacr_narrow_stringency:
    type: string
    enum: [stringent, relaxed]
    description: SEACR stringency for peak_mode=narrow samples.
    default: stringent
  seacr_broad_stringency:
    type: string
    enum: [stringent, relaxed]
    description: SEACR stringency for peak_mode=broad samples.
    default: relaxed

  consensus_window:
    type: integer
    minimum: 1
    description: Fixed MACS2-consensus peak width around each summit, bp.
    default: 500
  consensus_min_replicates:
    type: integer
    minimum: 1
    description: Majority-vote threshold for conditions with >=3 replicates.
    default: 2
  idr_threshold:
    type: number
    minimum: 0
    maximum: 1
    description: IDR threshold for conditions with exactly 2 replicates.
    default: 0.05
  idr_relaxed_pvalue:
    type: number
    minimum: 0
    maximum: 1
    description: MACS2 -p for the relaxed peak calls used as IDR input.
    default: 0.1
  idr_top_n_peaks:
    type: integer
    minimum: 1
    description: Top relaxed peaks retained per replicate for IDR.
    default: 150000
  keep_chroms_regex:
    type: string
    description: Regex used by the consensus step to filter chromosomes.
    default: "^chr([1-9]|1[0-9]|2[0-2]|X)$"

  differential_counts:
    type: string
    enum: [macs2, seacr, both]
    description: Which consensus count matrix the differential notebook analyzes.
    default: macs2

  gtf:
    type: string
    description: GENCODE GTF (chr-prefixed) used for TSS-enrichment QC.
    default: ref/gencode.v36.annotation.gtf
  promoter_bed:
    type: string
    description: Promoter BED for reads-in-annotation QC and promoter/distal split.
    default: ref/promoter_chr1-22X.bed
  enhancer_bed:
    type: string
    description: Enhancer BED for reads-in-annotation QC.
    default: ref/enhancer_chr1-22X.bed

required:
  - samples_table
  - genome_fasta
  - genome_index
  - align_chroms
  - keep_chroms
  - blacklist
  - effective_genome_size
  - bin_size
  - max_fragment_length
  - remove_duplicates
  - macs2_genome
  - macs2_qvalue
  - macs2_broad_cutoff
  - seacr_norm
  - seacr_narrow_stringency
  - seacr_broad_stringency
  - consensus_window
  - consensus_min_replicates
  - idr_threshold
  - idr_relaxed_pvalue
  - idr_top_n_peaks
  - keep_chroms_regex
  - differential_counts
  - gtf
  - promoter_bed
  - enhancer_bed
```

- [ ] **Step 2: Write `config/config.yaml`**

```yaml
# CUT&RUN analysis configuration (no spike-in)

samples_table: "config/samples.csv"

# fastp auto-detects adapters for paired-end reads; uncomment to override.
# adapter_r1: "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
# adapter_r2: "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

# ── Alignment (human-only Bowtie2 index) ────────────────────────────────
genome_fasta: "ref/genome.fa"
genome_index: "ref/genome/genome"
align_chroms: [chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, chr20, chr21, chr22, chrX, chrM]
keep_chroms: [chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, chr20, chr21, chr22, chrX]
blacklist: "ref/hg38_blacklist_regions.bed"

max_fragment_length: 700      # Bowtie2 -X (CUT&RUN fragments are short)
remove_duplicates: true       # false keeps duplicates (low-input CUT&RUN)

# ── deepTools tracks ────────────────────────────────────────────────────
effective_genome_size: 2913022398
bin_size: 25

# ── MACS2 ───────────────────────────────────────────────────────────────
macs2_genome: "hs"
macs2_qvalue: 0.05
macs2_broad_cutoff: 0.1

# ── SEACR ───────────────────────────────────────────────────────────────
seacr_norm: "norm"                 # norm when an IgG control is provided
seacr_narrow_stringency: "stringent"
seacr_broad_stringency: "relaxed"

# ── MACS2 consensus + reproducibility ───────────────────────────────────
consensus_window: 500
consensus_min_replicates: 2
idr_threshold: 0.05
idr_relaxed_pvalue: 0.1
idr_top_n_peaks: 150000
keep_chroms_regex: "^chr([1-9]|1[0-9]|2[0-2]|X)$"

# ── Differential binding ────────────────────────────────────────────────
differential_counts: "macs2"       # macs2 | seacr | both

# ── QC annotation ───────────────────────────────────────────────────────
gtf: "ref/gencode.v36.annotation.gtf"
promoter_bed: "ref/promoter_chr1-22X.bed"
enhancer_bed: "ref/enhancer_chr1-22X.bed"
```

- [ ] **Step 3: Write `config/samples.csv` (reference example)**

```csv
sample_id,condition,replicate,input_control,peak_mode,notes
GSF2801-ChIPseq-OVCAR3-3D-IP-cJun_S4,cJUN_3D,1,GSF2801-ChIPseq-OVCAR3-3D-IP-IgG_S5,narrow,3D-cJUN
GSF2801-ChIPseq-OVCAR3-3D-IP-IgG_S5,IgG_3D,1,,,3D-Igg
GSF2801-ChIPseq-OVCAR3-Control-IP-cJun_S1,cJUN_Ctrl,1,GSF2801-ChIPseq-OVCAR3-Control-IP-IgG_S2,narrow,Ctrl-cJUN
GSF2801-ChIPseq-OVCAR3-Control-IP-IgG_S2,IgG_Ctrl,1,,,Ctrl-Igg
```
(Note: the two cJUN rows are given distinct `condition` labels `cJUN_3D`/`cJUN_Ctrl` so they are not mistaken for replicates — see config/README.md.)

- [ ] **Step 4: Write `.test/config/config.yaml`**

Copy `config/config.yaml` and change only these paths so dry-run finds the `.test/ref` placeholders:
```yaml
genome_fasta: "ref/genome.fa"
genome_index: "ref/genome/genome"
blacklist: "ref/hg38_blacklist_regions.bed"
gtf: "ref/gencode.gtf"
promoter_bed: "ref/promoter_chr1-22X.bed"
enhancer_bed: "ref/enhancer_chr1-22X.bed"
```
(These are relative to `.test/` because dry-runs pass `-d .test`.) Keep every other key identical to `config/config.yaml`.

- [ ] **Step 5: Write `tests/test_config_schema.py`**

```python
import yaml, json
from pathlib import Path
from jsonschema import validate

def _load(p):
    return yaml.safe_load(Path(p).read_text())

def test_default_config_matches_schema():
    schema = _load("workflow/schemas/config.schema.yaml")
    config = _load("config/config.yaml")
    # Fill required keys that are commented-out optionals
    validate(instance=config, schema=schema)

def test_test_config_matches_schema():
    schema = _load("workflow/schemas/config.schema.yaml")
    config = _load(".test/config/config.yaml")
    validate(instance=config, schema=schema)
```

- [ ] **Step 6: Run the schema tests**

Run: `pytest tests/test_config_schema.py -v`
Expected: 2 passed. (If `jsonschema` is missing: `pip install jsonschema pyyaml`.)

- [ ] **Step 7: Write `config/README.md`**

Adapt `$ATAC/config/README.md`: replace the sample-sheet table with the 6-column CUT&RUN schema; document control-vs-treatment (empty `peak_mode` = control), the per-condition same-`peak_mode` rule, and the "give biologically distinct groups distinct `condition` labels" note. Point reference-data section at `genome_fasta`, `blacklist`, `gtf`, `promoter_bed`, `enhancer_bed` (drop dm6/spike-in). State the `genome_index` and `genome.chrom.sizes` are built automatically.

- [ ] **Step 8: Commit**

```bash
git add -A && git commit -m "feat: config schema, default config, sample sheet, config README"
```

---

## Task 3: `samplesheet.py` — sample-sheet parsing, derived sets, validation, helpers

**Files:**
- Create: `workflow/scripts/samplesheet.py`
- Test: `tests/test_samplesheet.py`

**Interfaces:**
- Produces (imported by `common.smk` in Task 4):
  - `load_samples(path) -> pandas.DataFrame` (columns coerced: `sample_id,condition,replicate:int,input_control:str,peak_mode:str,notes`; NaN→"").
  - `SampleSheet` dataclass-like object with attributes/methods:
    - `.all_samples -> list[str]`
    - `.treatment_samples -> list[str]` (peak_mode non-empty)
    - `.control_samples -> list[str]` (peak_mode empty)
    - `.narrow_samples -> list[str]`, `.broad_samples -> list[str]`
    - `.groups -> dict[str, list[str]]` (condition → treatment sample_ids)
    - `.group_method -> dict[str, str]` ("majority"|"idr"|"single")
    - `.idr_groups -> list[str]`, `.idr_samples -> list[str]`
    - `.idr_pairs -> list[tuple[str,str,str]]` (group, rep_a, rep_b)
    - `.peak_mode(sample) -> str`
    - `.macs2_ext(sample) -> str` ("narrowPeak"|"broadPeak")
    - `.input_control(sample) -> str` ("" if none)
    - `.seacr_stringency(sample, cfg) -> str` (from cfg["seacr_narrow_stringency"]/broad)
    - `.validate()` raises `ValueError` on the three validation rules.
  - Module constant `REQUIRED_COLUMNS = ["sample_id","condition","replicate","input_control","peak_mode","notes"]`.

- [ ] **Step 1: Write failing tests `tests/test_samplesheet.py`**

```python
import textwrap
import pytest
from pathlib import Path
import sys
sys.path.insert(0, "workflow/scripts")
import samplesheet as ss

def _sheet(tmp_path, rows):
    p = tmp_path / "samples.csv"
    p.write_text("sample_id,condition,replicate,input_control,peak_mode,notes\n" + rows)
    return ss.SampleSheet(str(p))

REF = ("t1,cJUN,1,igg,narrow,a\n"
       "t2,cJUN,2,igg,narrow,b\n"
       "t3,K27,1,igg,broad,c\n"
       "igg,IgG,1,,,ctrl\n")

def test_sets(tmp_path):
    s = _sheet(tmp_path, REF)
    assert s.all_samples == ["t1","t2","t3","igg"]
    assert s.treatment_samples == ["t1","t2","t3"]
    assert s.control_samples == ["igg"]
    assert s.narrow_samples == ["t1","t2"]
    assert s.broad_samples == ["t3"]

def test_groups_and_method(tmp_path):
    s = _sheet(tmp_path, REF)
    assert s.groups == {"cJUN": ["t1","t2"], "K27": ["t3"]}
    assert s.group_method == {"cJUN": "idr", "K27": "single"}
    assert s.idr_groups == ["cJUN"]
    assert s.idr_samples == ["t1","t2"]
    assert s.idr_pairs == [("cJUN","t1","t2")]

def test_helpers(tmp_path):
    s = _sheet(tmp_path, REF)
    assert s.macs2_ext("t1") == "narrowPeak"
    assert s.macs2_ext("t3") == "broadPeak"
    assert s.input_control("t1") == "igg"
    assert s.input_control("igg") == ""
    cfg = {"seacr_narrow_stringency": "stringent", "seacr_broad_stringency": "relaxed"}
    assert s.seacr_stringency("t1", cfg) == "stringent"
    assert s.seacr_stringency("t3", cfg) == "relaxed"

def test_majority_method(tmp_path):
    s = _sheet(tmp_path, "a,X,1,igg,narrow,\nb,X,2,igg,narrow,\nc,X,3,igg,narrow,\nigg,IgG,1,,,\n")
    assert s.group_method["X"] == "majority"

def test_validate_bad_control_ref(tmp_path):
    s = _sheet(tmp_path, "t1,cJUN,1,missing,narrow,\nigg,IgG,1,,,\n")
    with pytest.raises(ValueError, match="input_control"):
        s.validate()

def test_validate_mixed_peakmode(tmp_path):
    s = _sheet(tmp_path, "t1,cJUN,1,igg,narrow,\nt2,cJUN,2,igg,broad,\nigg,IgG,1,,,\n")
    with pytest.raises(ValueError, match="peak_mode"):
        s.validate()

def test_validate_bad_peakmode_value(tmp_path):
    s = _sheet(tmp_path, "t1,cJUN,1,igg,wide,\nigg,IgG,1,,,\n")
    with pytest.raises(ValueError, match="narrow|broad"):
        s.validate()

def test_validate_ok(tmp_path):
    _sheet(tmp_path, REF).validate()  # no raise
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_samplesheet.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'samplesheet'` (or AttributeError).

- [ ] **Step 3: Implement `workflow/scripts/samplesheet.py`**

```python
#!/usr/bin/env python3
"""Sample-sheet parsing, derived sample sets, validation and helpers for the
CUT&RUN workflow. Pure module (no snakemake import) so it is unit-testable and
importable from common.smk."""
from itertools import combinations
import pandas as pd

REQUIRED_COLUMNS = ["sample_id", "condition", "replicate", "input_control",
                    "peak_mode", "notes"]
VALID_PEAK_MODES = {"narrow", "broad"}


def load_samples(path):
    df = pd.read_csv(path, dtype=str).fillna("")
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"samples sheet {path} missing columns: {missing}")
    for c in ("sample_id", "condition", "input_control", "peak_mode"):
        df[c] = df[c].astype(str).str.strip()
    df["replicate"] = df["replicate"].astype(str).str.strip()
    return df


def _repro_method(n):
    if n >= 3:
        return "majority"
    if n == 2:
        return "idr"
    return "single"


class SampleSheet:
    def __init__(self, path):
        self.df = load_samples(path)
        d = self.df
        self.all_samples = d["sample_id"].tolist()
        self._peak_mode = dict(zip(d["sample_id"], d["peak_mode"]))
        self._input = dict(zip(d["sample_id"], d["input_control"]))
        self.treatment_samples = [s for s in self.all_samples if self._peak_mode[s]]
        self.control_samples = [s for s in self.all_samples if not self._peak_mode[s]]
        self.narrow_samples = [s for s in self.treatment_samples if self._peak_mode[s] == "narrow"]
        self.broad_samples = [s for s in self.treatment_samples if self._peak_mode[s] == "broad"]
        trt = d[d["sample_id"].isin(self.treatment_samples)]
        self.groups = {g: m["sample_id"].tolist()
                       for g, m in trt.groupby("condition", sort=False)}
        self.group_method = {g: _repro_method(len(m)) for g, m in self.groups.items()}
        self.idr_groups = [g for g in self.groups if self.group_method[g] == "idr"]
        self.idr_samples = [s for g in self.idr_groups for s in self.groups[g]]
        self.idr_pairs = []
        for g, members in self.groups.items():
            for a, b in combinations(members, 2):
                self.idr_pairs.append((g, a, b))

    def peak_mode(self, sample):
        return self._peak_mode.get(sample, "")

    def macs2_ext(self, sample):
        return "broadPeak" if self._peak_mode.get(sample) == "broad" else "narrowPeak"

    def input_control(self, sample):
        return self._input.get(sample, "")

    def seacr_stringency(self, sample, cfg):
        if self._peak_mode.get(sample) == "broad":
            return cfg["seacr_broad_stringency"]
        return cfg["seacr_narrow_stringency"]

    def validate(self):
        controls = set(self.control_samples)
        for s in self.treatment_samples:
            ic = self._input[s]
            if ic and ic not in controls:
                raise ValueError(
                    f"input_control '{ic}' for sample '{s}' is not an existing "
                    f"control (empty peak_mode) sample")
        for g, members in self.groups.items():
            modes = {self._peak_mode[s] for s in members}
            if len(modes) > 1:
                raise ValueError(
                    f"condition '{g}' mixes peak_mode values {modes}; all "
                    f"replicates of a condition must share narrow or broad")
        for s in self.all_samples:
            pm = self._peak_mode[s]
            if pm and pm not in VALID_PEAK_MODES:
                raise ValueError(
                    f"sample '{s}' has peak_mode '{pm}'; must be narrow, broad or empty")
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_samplesheet.py -v`
Expected: 8 passed.

- [ ] **Step 5: Commit**

```bash
git add -A && git commit -m "feat: samplesheet.py sample-set parsing + validation with tests"
```

---

## Task 4: `Snakefile`, `common.smk`, and loadable rule-file skeletons

**Files:**
- Create: `workflow/Snakefile`, `workflow/rules/common.smk`, `workflow/rules/cutandrun.smk` (skeleton), `workflow/rules/qc.smk` (skeleton)

**Interfaces:**
- Consumes: `samplesheet.py` (Task 3), config schema (Task 2).
- Produces: global names for all rules — directory constants (`RESULT_DIR, FASTQC_DIR, FASTP_DIR, ALIGN_DIR, TMP_DIR, FILTERED_DIR, DEDUP_DIR, BLACKLIST_FILTERED_DIR, PEAKS_DIR, SEACR_DIR, SEACR_BEDGRAPH_DIR, BIGWIG_DIR, LOG2_BIGWIG_DIR, CONSENSUS_DIR, SEACR_CONSENSUS_DIR, QC_DIR, RELAXED_PEAKS_DIR` and the QC dirs mirrored from ATAC), the `SS` (SampleSheet) object and its exported lists (`SAMPLES, TREATMENT_SAMPLES, CONTROL_SAMPLES, NARROW_SAMPLES, BROAD_SAMPLES, GROUPS, GROUP_METHOD, IDR_GROUPS, IDR_SAMPLES, IDR_PAIRS`), reference constants (`GENOME_FASTA, GENOME_INDEX, CHROM_SIZES, GENOME_2BIT, GTF_FILE, PROMOTER_BED, ENHANCER_BED, EGS`), and helpers `_alt`, `FASTP_ADAPTER_ARGS`, `igg_bam(sample)`, `macs2_peak(sample)`, `_group_relaxed_inputs`.

- [ ] **Step 1: Write `workflow/Snakefile`**

```python
# CUT&RUN pipeline (no spike-in): peak calling (MACS2 + SEACR) + consensus + QC.
from snakemake.utils import min_version

min_version("8.0")

configfile: "config/config.yaml"

include: "rules/common.smk"
include: "rules/cutandrun.smk"
include: "rules/qc.smk"

rule all:
    input:
        rules.cutandrun_all.input,
        rules.qc_all.input,
    default_target: True
```

- [ ] **Step 2: Write `workflow/rules/common.smk`**

```python
# Shared setup for cutandrun.smk and qc.smk: config validation, sample sheet
# (via workflow/scripts/samplesheet.py), directory constants, and helpers.
import os
import re
import sys
from snakemake.utils import validate

validate(config, "../schemas/config.schema.yaml")

sys.path.insert(0, os.path.join(os.path.dirname(workflow.snakefile), "scripts"))
import samplesheet

# ── Samples ─────────────────────────────────────────────────────────────
SS = samplesheet.SampleSheet(config["samples_table"])
SS.validate()
SAMPLES           = SS.all_samples
TREATMENT_SAMPLES = SS.treatment_samples
CONTROL_SAMPLES   = SS.control_samples
NARROW_SAMPLES    = SS.narrow_samples
BROAD_SAMPLES     = SS.broad_samples
GROUPS            = SS.groups
GROUP_METHOD      = SS.group_method
IDR_GROUPS        = SS.idr_groups
IDR_SAMPLES       = SS.idr_samples
IDR_PAIRS         = SS.idr_pairs
# Treatment samples that have an IgG/Input control (for log2 tracks + SEACR)
CONTROLLED_SAMPLES = [s for s in TREATMENT_SAMPLES if SS.input_control(s)]

# ── Output directories ──────────────────────────────────────────────────
RESULT_DIR             = "results"
FASTQC_DIR             = f"{RESULT_DIR}/fastqc"
FASTP_DIR              = f"{RESULT_DIR}/fastp"
ALIGN_DIR              = f"{RESULT_DIR}/aligned"
TMP_DIR                = f"{RESULT_DIR}/tmp"
FILTERED_DIR           = f"{RESULT_DIR}/filtered"
DEDUP_DIR              = f"{RESULT_DIR}/dedup"
BLACKLIST_FILTERED_DIR = f"{RESULT_DIR}/blacklist_filtered"
PEAKS_DIR              = f"{RESULT_DIR}/peaks"
BIGWIG_DIR             = f"{RESULT_DIR}/bigwig"
LOG2_BIGWIG_DIR        = f"{RESULT_DIR}/log2ratio_bigwig"
SEACR_DIR              = f"{RESULT_DIR}/seacr"
SEACR_BEDGRAPH_DIR     = f"{SEACR_DIR}/bedgraph"
RELAXED_PEAKS_DIR      = f"{RESULT_DIR}/peaks_relaxed"
CONSENSUS_DIR          = f"{RESULT_DIR}/consensus"
SEACR_CONSENSUS_DIR    = f"{RESULT_DIR}/consensus_seacr"
QC_DIR                 = f"{RESULT_DIR}/qc"

# QC aliases + QC-only dirs (mirror ATAC)
RMD_BAM_DIR    = BLACKLIST_FILTERED_DIR
PEAK_DIR       = PEAKS_DIR
BEDGRAPH_DIR   = f"{RESULT_DIR}/bedgraph"
DEEPTOOLS_DIR  = f"{RESULT_DIR}/deeptools"
FRIP_DIR       = f"{RESULT_DIR}/FRiP"
IDR_DIR        = f"{RESULT_DIR}/idr"
RELAXED_DIR    = f"{RESULT_DIR}/qc_relaxed_peaks"
COMPLEXITY_DIR = f"{RESULT_DIR}/library_complexity"
ANNOT_DIR      = f"{RESULT_DIR}/peak_annotation"

# ── Reference data / config ─────────────────────────────────────────────
GENOME_FASTA = config["genome_fasta"]
GENOME_INDEX = config["genome_index"]
CHROM_SIZES  = os.path.join("ref", "genome.chrom.sizes")
GENOME_2BIT  = os.path.join("ref", "hg38.2bit")
GTF_FILE     = config["gtf"]
PROMOTER_BED = config["promoter_bed"]
ENHANCER_BED = config["enhancer_bed"]
EGS          = config["effective_genome_size"]

def _alt(names):
    return "|".join(re.escape(n) for n in names) if names else "a^"

def _fastp_adapter_args():
    r1 = str(config.get("adapter_r1") or "").strip()
    r2 = str(config.get("adapter_r2") or "").strip()
    if r1:
        args = f"--adapter_sequence {r1}"
        if r2:
            args += f" --adapter_sequence_r2 {r2}"
        return args
    return "--detect_adapter_for_pe"

FASTP_ADAPTER_ARGS = _fastp_adapter_args()

def igg_bam(sample):
    """Blacklist-filtered BAM of a sample's IgG control ('' if none)."""
    ic = SS.input_control(sample)
    return f"{BLACKLIST_FILTERED_DIR}/{ic}.nobl.bam" if ic else ""

def macs2_peak(sample):
    """Per-sample MACS2 peak path with the correct narrow/broad extension."""
    return f"{PEAKS_DIR}/{sample}_peaks.{SS.macs2_ext(sample)}"

def seacr_peak(sample):
    """Per-sample SEACR output BED path."""
    return f"{SEACR_DIR}/{sample}.{SS.seacr_stringency(sample, config)}.bed"

def _group_relaxed_inputs(wildcards):
    ext = "broadPeak" if SS.peak_mode(GROUPS[wildcards.group][0]) == "broad" else "narrowPeak"
    return [f"{RELAXED_PEAKS_DIR}/{s}_relaxed.{ext}" for s in GROUPS[wildcards.group]]
```

- [ ] **Step 3: Write skeleton `workflow/rules/cutandrun.smk`**

```python
# Primary CUT&RUN pipeline. Rules added task-by-task; cutandrun_all aggregates.
rule cutandrun_all:
    input:
        []
```

- [ ] **Step 4: Write skeleton `workflow/rules/qc.smk`**

```python
# CUT&RUN QC pipeline. Rules added task-by-task; qc_all aggregates.
rule qc_all:
    input:
        []
```

- [ ] **Step 5: Verify the skeleton loads and validates config + sample sheet**

Run: `snakemake -s workflow/Snakefile -d .test -n`
Expected: "Building DAG of jobs..." then "Nothing to be done" (0 jobs) — no schema or sample-sheet exceptions. This confirms `common.smk` parses, the schema validates `.test/config/config.yaml`, and `SS.validate()` passes.

- [ ] **Step 6: Commit**

```bash
git add -A && git commit -m "feat: Snakefile + common.smk (sample sets, dirs, helpers) + rule skeletons"
```

---

## Task 5: Preprocessing rules — index, chrom sizes, FastQC, fastp, align, filter, dedup, blacklist

**Files:**
- Modify: `workflow/rules/cutandrun.smk` (append rules; extend `cutandrun_all`)

**Interfaces:**
- Consumes: `common.smk` names.
- Produces: `{BLACKLIST_FILTERED_DIR}/{sample}.nobl.bam(.bai)`, `{DEDUP_DIR}/{sample}.dedup.bam` + metrics, `{FILTERED_DIR}/{sample}.sorted.filtered.bam` + idxstats, `{ALIGN_DIR}/{sample}.bowtie2.log`, `{CHROM_SIZES}`, `{GENOME_INDEX}.build.done`.

- [ ] **Step 1: Append the preprocessing rules**

Adapt from `$ATAC/workflow/rules/atacseq.smk`. Use these rules (copy the ATAC rule bodies except where noted):

`genome_chrom_sizes` (NEW):
```python
rule genome_chrom_sizes:
    input: fasta = GENOME_FASTA
    output: sizes = CHROM_SIZES
    conda: "../envs/snakemake.yaml"
    log: "logs/genome_chrom_sizes/sizes.log"
    shell:
        r"""
        mkdir -p $(dirname {output.sizes}) logs/genome_chrom_sizes
        samtools faidx {input.fasta} 2> {log}
        cut -f1,2 {input.fasta}.fai > {output.sizes} 2>> {log}
        """
```

`build_genome_index` (REPLACES build_combined_genome — human only, no spike-in):
```python
rule build_genome_index:
    input:
        fasta = GENOME_FASTA
    output:
        done = touch(f"{GENOME_INDEX}.build.done")
    params:
        index  = GENOME_INDEX,
        chroms = config["align_chroms"],
        subset = f"{GENOME_INDEX}.subset.fa"
    threads: 8
    conda: "../envs/snakemake.yaml"
    log: "logs/build_genome_index/build.log"
    shell:
        r"""
        mkdir -p $(dirname {params.index}) logs/build_genome_index
        if [ -n "{params.chroms}" ]; then
            samtools faidx {input.fasta} {params.chroms} > {params.subset} 2>> {log}
            GENOME={params.subset}
        else
            GENOME={input.fasta}
        fi
        bowtie2-build --threads {threads} $GENOME {params.index} >> {log} 2>&1
        if [ -n "{params.chroms}" ]; then rm -f {params.subset}; fi
        """
```

`fastqc`, `fastp` — copy verbatim from `$ATAC` (`atacseq.smk` lines 53–101).

`bowtie2_align` (CUT&RUN params, human index):
```python
rule bowtie2_align:
    input:
        r1   = f"{FASTP_DIR}/{{sample}}_R1.trimmed.fastq.gz",
        r2   = f"{FASTP_DIR}/{{sample}}_R2.trimmed.fastq.gz",
        done = f"{GENOME_INDEX}.build.done"
    output:
        bam = f"{ALIGN_DIR}/{{sample}}.bam",
        bai = f"{ALIGN_DIR}/{{sample}}.bam.bai",
        summary = f"{ALIGN_DIR}/{{sample}}.bowtie2.log"
    params:
        index = GENOME_INDEX,
        maxfrag = config["max_fragment_length"]
    threads: 20
    conda: "../envs/snakemake.yaml"
    shell:
        r"""
        mkdir -p {ALIGN_DIR} {TMP_DIR}
        bowtie2 -x {params.index} -1 {input.r1} -2 {input.r2} \
               -p {threads} \
               --end-to-end --very-sensitive --no-mixed --no-discordant --no-unal \
               -I 10 -X {params.maxfrag} \
               2> {output.summary} \
            | samtools sort -@ 4 -m 2G -T {TMP_DIR}/{wildcards.sample}.rawsort -o {output.bam} -
        samtools index -@ {threads} {output.bam}
        """
```

`samtools_sort_filter_index` — copy the ATAC rule (`atacseq.smk` 136–189) but REMOVE the spike-in prefix logic: drop `params.spikein_prefix` and change the unique-read `awk` from `awk -v p="^{params.spikein_prefix}" '/^@/ || $3 !~ p'` to a plain header/keep pass. Replace the filter `samtools view` pipeline body with:
```bash
        samtools view -@ {threads} -hS -f 2 -F 2316 {input} | grep -v "XS:i:" \
            > {TMP_DIR}/temp_{wildcards.sample}.unsorted.sam 2>> {log}
```
Keep everything else (name-sort → `process_sam.py` orphan removal → coordinate-sort prekeep → idxstats mito-% → `keep_chroms` restriction → flagstat → cleanup) identical, with `params.keep_chroms = config["keep_chroms"]`.

`remove_duplicates` (toggle-aware):
```python
rule remove_duplicates:
    input:
        filtered_bam = f"{FILTERED_DIR}/{{sample}}.sorted.filtered.bam"
    output:
        dedup_bam = f"{DEDUP_DIR}/{{sample}}.dedup.bam",
        metrics = f"{DEDUP_DIR}/{{sample}}.dedup.metrics.txt"
    params:
        remove = "true" if config["remove_duplicates"] else "false"
    threads: 4
    conda: "../envs/snakemake.yaml"
    log: "logs/dedup/{sample}.log"
    shell:
        r"""
        mkdir -p {DEDUP_DIR}
        java -jar ref/picard.jar MarkDuplicates \
               INPUT={input.filtered_bam} OUTPUT={output.dedup_bam} \
               METRICS_FILE={output.metrics} REMOVE_DUPLICATES={params.remove} \
               ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=tmp 2> {log}
        samtools index {output.dedup_bam}
        """
```

`filter_blacklist` — copy verbatim from `$ATAC` (`atacseq.smk` 218–292).

- [ ] **Step 2: Extend `cutandrun_all` input**

Replace the `cutandrun_all` input list with:
```python
rule cutandrun_all:
    input:
        expand(f"{FASTQC_DIR}/{{s}}_R1_001_fastqc.html", s=SAMPLES),
        expand(f"{FASTQC_DIR}/{{s}}_R2_001_fastqc.html", s=SAMPLES),
        expand(f"{FASTP_DIR}/{{s}}.html", s=SAMPLES),
        expand(f"{FILTERED_DIR}/{{s}}.sorted.filtered.bam", s=SAMPLES),
        expand(f"{DEDUP_DIR}/{{s}}.dedup.bam", s=SAMPLES),
        expand(f"{BLACKLIST_FILTERED_DIR}/{{s}}.nobl.bam", s=SAMPLES),
        expand(f"{BLACKLIST_FILTERED_DIR}/{{s}}.nobl.bam.bai", s=SAMPLES),
```
(Leave a trailing comma; later tasks append more targets.)

- [ ] **Step 3: Dry-run the preprocessing chain**

Run: `snakemake -s workflow/Snakefile -d .test -n results/blacklist_filtered/test_igg.nobl.bam`
Expected: DAG builds; job list includes `build_genome_index, fastqc, fastp, bowtie2_align, samtools_sort_filter_index, remove_duplicates, filter_blacklist`; no exceptions.

- [ ] **Step 4: Commit**

```bash
git add -A && git commit -m "feat: preprocessing rules (index, align, filter, dedup, blacklist)"
```

---

## Task 6: Signal tracks — RPGC bigWig + IgG-subtracted log2 bigWig

**Files:**
- Modify: `workflow/rules/cutandrun.smk`

**Interfaces:**
- Produces: `{BIGWIG_DIR}/{sample}.bw` (all samples), `{LOG2_BIGWIG_DIR}/{sample}.log2ratio.bw` (CONTROLLED_SAMPLES).

- [ ] **Step 1: Append `create_bigwig`**

Copy the ATAC `create_bigwig` rule (`atacseq.smk` 427–453) verbatim (RPGC, blacklist, extendReads) — it already uses `BLACKLIST_FILTERED_DIR`, `EGS`, `bin_size`.

- [ ] **Step 2: Append `create_log2ratio_bigwig` (NEW)**

```python
rule create_log2ratio_bigwig:
    wildcard_constraints:
        sample = _alt(CONTROLLED_SAMPLES)
    input:
        treat = f"{BLACKLIST_FILTERED_DIR}/{{sample}}.nobl.bam",
        treat_bai = f"{BLACKLIST_FILTERED_DIR}/{{sample}}.nobl.bam.bai",
        ctrl = lambda w: igg_bam(w.sample),
        ctrl_bai = lambda w: igg_bam(w.sample) + ".bai"
    output:
        bw = f"{LOG2_BIGWIG_DIR}/{{sample}}.log2ratio.bw"
    params:
        bin_size = config["bin_size"],
        blacklist = config["blacklist"]
    threads: 8
    conda: "../envs/deeptools.yaml"
    log: "logs/log2ratio_bigwig/{sample}.log"
    shell:
        r"""
        mkdir -p {LOG2_BIGWIG_DIR} logs/log2ratio_bigwig
        bamCompare -b1 {input.treat} -b2 {input.ctrl} \
            --operation log2 --normalizeUsing CPM \
            --binSize {params.bin_size} --numberOfProcessors {threads} \
            --extendReads --blackListFileName {params.blacklist} \
            --outFileName {output.bw} > {log} 2>&1
        """
```

- [ ] **Step 3: Extend `cutandrun_all` with track targets**

Append:
```python
        expand(f"{BIGWIG_DIR}/{{s}}.bw", s=SAMPLES),
        expand(f"{LOG2_BIGWIG_DIR}/{{s}}.log2ratio.bw", s=CONTROLLED_SAMPLES),
```

- [ ] **Step 4: Dry-run tracks**

Run: `snakemake -s workflow/Snakefile -d .test -n results/log2ratio_bigwig/test_cjun_r1.log2ratio.bw results/bigwig/test_igg.bw`
Expected: DAG builds; `create_log2ratio_bigwig` job shows `-b1 ...test_cjun_r1.nobl.bam -b2 ...test_igg.nobl.bam`; no exceptions.

- [ ] **Step 5: Commit**

```bash
git add -A && git commit -m "feat: RPGC bigWig + IgG-subtracted log2 bamCompare tracks"
```

---

## Task 7: MACS2 peak calling (narrow + broad, IgG-control aware)

**Files:**
- Modify: `workflow/rules/cutandrun.smk`

**Interfaces:**
- Produces: `{PEAKS_DIR}/{sample}_peaks.narrowPeak` (NARROW_SAMPLES), `{PEAKS_DIR}/{sample}_peaks.broadPeak` (BROAD_SAMPLES).

- [ ] **Step 1: Append `call_peaks_macs2_narrow`**

```python
rule call_peaks_macs2_narrow:
    wildcard_constraints:
        sample = _alt(NARROW_SAMPLES)
    input:
        treatment = f"{BLACKLIST_FILTERED_DIR}/{{sample}}.nobl.bam",
        control = lambda w: [igg_bam(w.sample)] if igg_bam(w.sample) else []
    output:
        peaks = f"{PEAKS_DIR}/{{sample}}_peaks.narrowPeak"
    params:
        outdir = PEAKS_DIR, name = "{sample}",
        genome = config["macs2_genome"], q = config["macs2_qvalue"],
        control_arg = lambda w: f"-c {igg_bam(w.sample)}" if igg_bam(w.sample) else ""
    conda: "../envs/macs2.yaml"
    log: "logs/macs2/{sample}.log"
    shell:
        r"""
        mkdir -p {params.outdir} logs/macs2
        macs2 callpeak -t {input.treatment} {params.control_arg} \
              -f BAMPE -g {params.genome} --outdir {params.outdir} \
              -n {params.name} --nomodel -q {params.q} > {log} 2>&1
        """
```

- [ ] **Step 2: Append `call_peaks_macs2_broad`**

Same as narrow but `wildcard_constraints: sample = _alt(BROAD_SAMPLES)`, output `..._peaks.broadPeak`, and the shell adds `--broad --broad-cutoff {params.broad_cutoff}` (add `params.broad_cutoff = config["macs2_broad_cutoff"]`).

- [ ] **Step 3: Extend `cutandrun_all` with per-sample MACS2 peaks**

Append (uses the `macs2_peak` helper for correct extension):
```python
        [macs2_peak(s) for s in TREATMENT_SAMPLES],
```

- [ ] **Step 4: Dry-run both peak modes**

Run: `snakemake -s workflow/Snakefile -d .test -n results/peaks/test_cjun_r1_peaks.narrowPeak results/peaks/test_k27_r1_peaks.broadPeak`
Expected: DAG builds; narrow job has `-c ...test_igg.nobl.bam` and no `--broad`; broad job has `--broad --broad-cutoff 0.1`; no exceptions.

- [ ] **Step 5: Commit**

```bash
git add -A && git commit -m "feat: MACS2 narrow/broad peak calling with per-sample IgG control"
```

---

## Task 8: SEACR peak calling (bedgraph + control-aware SEACR)

**Files:**
- Modify: `workflow/rules/cutandrun.smk`

**Interfaces:**
- Consumes: `{CHROM_SIZES}` (Task 5).
- Produces: `{SEACR_BEDGRAPH_DIR}/{sample}.fragments.bedgraph` (treatment + referenced controls), `{SEACR_DIR}/{sample}.{stringency}.bed` (CONTROLLED_SAMPLES).

- [ ] **Step 1: Append `seacr_bedgraph`**

Needed for each treatment sample AND its control. Build for `set(CONTROLLED_SAMPLES) | {SS.input_control(s) for s in CONTROLLED_SAMPLES}`. Rule keyed on sample wildcard:
```python
rule seacr_bedgraph:
    input:
        bam = f"{BLACKLIST_FILTERED_DIR}/{{sample}}.nobl.bam",
        sizes = CHROM_SIZES
    output:
        bg = f"{SEACR_BEDGRAPH_DIR}/{{sample}}.fragments.bedgraph"
    params:
        tmp_bed = f"{TMP_DIR}/{{sample}}.seacr.frag.bed",
        tmp_ns  = f"{TMP_DIR}/{{sample}}.seacr.namesorted.bam"
    threads: 8
    conda: "../envs/bedtools.yaml"
    log: "logs/seacr_bedgraph/{sample}.log"
    shell:
        r"""
        mkdir -p {SEACR_BEDGRAPH_DIR} {TMP_DIR} logs/seacr_bedgraph
        samtools sort -n -@ {threads} -o {params.tmp_ns} {input.bam} 2> {log}
        bedtools bamtobed -bedpe -i {params.tmp_ns} 2>> {log} \
          | awk 'BEGIN{{OFS="\t"}} $1==$4 {{print $1, ($2<$5?$2:$5), ($3>$6?$3:$6)}}' \
          | sort -k1,1 -k2,2n > {params.tmp_bed} 2>> {log}
        bedtools genomecov -bg -i {params.tmp_bed} -g {input.sizes} > {output.bg} 2>> {log}
        rm -f {params.tmp_bed} {params.tmp_ns}
        """
```

- [ ] **Step 2: Append `call_peaks_seacr`**

```python
rule call_peaks_seacr:
    wildcard_constraints:
        sample = _alt(CONTROLLED_SAMPLES)
    input:
        treat = f"{SEACR_BEDGRAPH_DIR}/{{sample}}.fragments.bedgraph",
        ctrl = lambda w: f"{SEACR_BEDGRAPH_DIR}/{SS.input_control(w.sample)}.fragments.bedgraph"
    output:
        bed = lambda w: f"{SEACR_DIR}/{w.sample}.{SS.seacr_stringency(w.sample, config)}.bed"
    params:
        norm = config["seacr_norm"],
        stringency = lambda w: SS.seacr_stringency(w.sample, config),
        prefix = lambda w: f"{SEACR_DIR}/{w.sample}"
    conda: "../envs/seacr.yaml"
    log: "logs/seacr/{sample}.log"
    shell:
        r"""
        mkdir -p {SEACR_DIR} logs/seacr
        SEACR_1.3.sh {input.treat} {input.ctrl} {params.norm} {params.stringency} {params.prefix} > {log} 2>&1
        """
```
Note: SEACR writes `<prefix>.<stringency>.bed`; the output lambda names exactly that file.

- [ ] **Step 3: Extend `cutandrun_all` with SEACR peaks**

Append:
```python
        [seacr_peak(s) for s in CONTROLLED_SAMPLES],
```

- [ ] **Step 4: Dry-run SEACR**

Run: `snakemake -s workflow/Snakefile -d .test -n results/seacr/test_cjun_r1.stringent.bed results/seacr/test_k27_r1.relaxed.bed`
Expected: DAG builds; `call_peaks_seacr` for `test_cjun_r1` consumes its own + `test_igg` bedgraph; narrow→stringent, broad→relaxed; no exceptions.

- [ ] **Step 5: Commit**

```bash
git add -A && git commit -m "feat: SEACR peak calling (fragment bedgraph + IgG control)"
```

---

## Task 9: MACS2 reproducibility + fixed-width consensus + counts (narrow/broad aware)

**Files:**
- Modify: `workflow/scripts/consensus_peaks.py` (narrow+broad loading), `workflow/rules/cutandrun.smk`
- Test: `tests/test_consensus_peaks.py`

**Interfaces:**
- Consumes: `macs2_peak(sample)`, IDR peaks.
- Produces: `{RELAXED_PEAKS_DIR}/{sample}_relaxed.{ext}`, `{CONSENSUS_DIR}/idr/{group}.idr_peaks.{ext}`, `{CONSENSUS_DIR}/consensus_peaks.{bed,saf}`, `{CONSENSUS_DIR}/consensus_counts.txt`.

- [ ] **Step 1: Copy `consensus_peaks.py` then write the failing broad test**

```bash
cp /easley/scratch/projects/amitra/amitra2016502/snakemake_ATACseq_spikein/workflow/scripts/consensus_peaks.py workflow/scripts/consensus_peaks.py
```
`tests/test_consensus_peaks.py`:
```python
import sys
sys.path.insert(0, "workflow/scripts")
import consensus_peaks as cp

def _write(p, lines):
    p.write_text("\n".join("\t".join(map(str, r)) for r in lines) + "\n")

def test_load_narrowpeak_uses_summit_offset(tmp_path):
    f = tmp_path / "a_peaks.narrowPeak"
    # chrom start end name score strand signal p q summitoffset
    _write(f, [["chr1", 100, 300, "p1", 50, ".", 5.0, 9.0, 7.5, 40]])
    peaks = cp.load_peaks(str(f), "a")
    assert peaks[0].summit == 140          # start + offset
    assert peaks[0].score == 7.5           # col 9 (-log10 q)

def test_load_broadpeak_uses_midpoint(tmp_path):
    f = tmp_path / "b_peaks.broadPeak"
    # broadPeak has 9 columns, no summit offset
    _write(f, [["chr1", 100, 300, "p1", 50, ".", 5.0, 9.0, 7.5]])
    peaks = cp.load_peaks(str(f), "b")
    assert peaks[0].summit == 200          # midpoint (100+300)//2
    assert peaks[0].score == 7.5
```

- [ ] **Step 2: Run to verify failure**

Run: `pytest tests/test_consensus_peaks.py -v`
Expected: FAIL — `AttributeError: module has no attribute 'load_peaks'`.

- [ ] **Step 3: Adapt `consensus_peaks.py`**

Rename `load_narrowpeak` → `load_peaks` and make it column-count aware. Replace the function body with:
```python
def load_peaks(path, sample):
    """Parse a MACS2 narrowPeak (10-col) or broadPeak (9-col) file into Peak
    objects tagged with `sample`. narrowPeak: summit = start + col10 offset;
    broadPeak: summit = peak midpoint. score = col9 (-log10 q) for both."""
    peaks = []
    for line in Path(path).read_text().splitlines():
        if not line.strip():
            continue
        f = line.split("\t")
        start, end = int(f[1]), int(f[2])
        score = float(f[8])
        summit = start + int(f[9]) if len(f) >= 10 else (start + end) // 2
        peaks.append(Peak(f[0], start, end, score, summit, sample))
    return peaks
```
Update the three internal callers (`majority_keep` inputs come via `build_consensus`; find every `load_narrowpeak(` and rename to `load_peaks(`). In the `if "snakemake" in globals()` block, change `_npaths` to use the per-sample extension: the rule will pass a dict, so replace the hard-coded `_peaks.narrowPeak` construction with values passed via `sm.params` (see Step 4). Concretely, change the block to read `_npaths = dict(sm.params.narrowpeak_paths)` and `_ipaths = dict(sm.params.idr_paths)`.

- [ ] **Step 4: Run tests to verify pass**

Run: `pytest tests/test_consensus_peaks.py -v`
Expected: 2 passed.

- [ ] **Step 5: Append reproducibility + consensus rules to `cutandrun.smk`**

`relaxed_peaks` (narrow/broad aware) — start from ATAC `relaxed_peaks` (`atacseq.smk` 489–522). Change: `wildcard_constraints.sample = _alt(IDR_SAMPLES)`; compute per-sample ext; output `f"{RELAXED_PEAKS_DIR}/{{sample}}_relaxed.{{ext}}"` is awkward with two wildcards — instead split into two rules `relaxed_peaks_narrow` (constraint `_alt([s for s in IDR_SAMPLES if s in NARROW_SAMPLES])`, output `_relaxed.narrowPeak`, no `--broad`) and `relaxed_peaks_broad` (broad set, output `_relaxed.broadPeak`, add `--broad --broad-cutoff {broad_cutoff}`). Keep the sort `-k8,8gr` + `head -n {top_n}` top-N trimming from ATAC.

`reproducible_idr` — start from ATAC (`atacseq.smk` 526–553). Change output to `f"{CONSENSUS_DIR}/idr/{{group}}.idr_peaks.{{ext}}"`? Two wildcards again — instead compute ext inside via a params lambda and keep output `f"{CONSENSUS_DIR}/idr/{{group}}.idr_peaks.narrowPeak"` for narrow IDR groups and `.broadPeak` for broad groups by splitting into `reproducible_idr_narrow`/`reproducible_idr_broad` with group-set constraints:
```python
NARROW_IDR_GROUPS = [g for g in IDR_GROUPS if SS.peak_mode(GROUPS[g][0]) == "narrow"]
BROAD_IDR_GROUPS  = [g for g in IDR_GROUPS if SS.peak_mode(GROUPS[g][0]) == "broad"]
```
(add these to `common.smk` Step 2 alongside IDR_GROUPS, or compute inline). Each variant sets `--input-file-type narrowPeak|broadPeak`. Body otherwise identical to ATAC.

`consensus_peaks` rule — start from ATAC (`atacseq.smk` 556–577). Change `input.narrowpeaks` → per-sample MACS2 peaks with correct ext: `input.peaks = [macs2_peak(s) for s in TREATMENT_SAMPLES]`, `input.idr = [f"{CONSENSUS_DIR}/idr/{g}.idr_peaks.{'broadPeak' if g in BROAD_IDR_GROUPS else 'narrowPeak'}" for g in IDR_GROUPS]`. Pass `params.narrowpeak_paths = {s: macs2_peak(s) for s in TREATMENT_SAMPLES}` and `params.idr_paths = {g: (idr path) for g in IDR_GROUPS}` (matching the script's Step-3 change). Keep `min_reps`, `window`, `keep_regex` params.

`count_fragments_consensus` — copy ATAC (`atacseq.smk` 580–601) but `bams = expand(f"{BLACKLIST_FILTERED_DIR}/{{s}}.nobl.bam", s=TREATMENT_SAMPLES)` (treatment samples only).

- [ ] **Step 6: Extend `cutandrun_all`**

Append:
```python
        f"{CONSENSUS_DIR}/consensus_peaks.bed",
        f"{CONSENSUS_DIR}/consensus_counts.txt",
```

- [ ] **Step 7: Dry-run consensus (narrow IDR + broad single)**

Run: `snakemake -s workflow/Snakefile -d .test -n results/consensus/consensus_counts.txt`
Expected: DAG builds; includes `relaxed_peaks_narrow` + `reproducible_idr_narrow` for the cJUN pair, `call_peaks_macs2_broad` for K27 (single, no IDR), `consensus_peaks`, `count_fragments_consensus` over 3 treatment BAMs; no exceptions.

- [ ] **Step 8: Commit**

```bash
git add -A && git commit -m "feat: MACS2 reproducibility + fixed-width consensus + counts (narrow/broad)"
```

---

## Task 10: SEACR consensus + counts (overlap-based variable-width union)

**Files:**
- Create: `workflow/scripts/seacr_consensus.py`
- Modify: `workflow/rules/cutandrun.smk`
- Test: `tests/test_seacr_consensus.py`

**Interfaces:**
- Consumes: per-sample SEACR beds, `GROUPS`, `GROUP_METHOD`, `consensus_min_replicates`, blacklist, `keep_chroms_regex`.
- Produces: `{SEACR_CONSENSUS_DIR}/consensus_peaks.{bed,saf}`, `{SEACR_CONSENSUS_DIR}/consensus_counts.txt`.

- [ ] **Step 1: Write failing tests `tests/test_seacr_consensus.py`**

```python
import sys
sys.path.insert(0, "workflow/scripts")
import seacr_consensus as sc

def test_seacr_bed_load(tmp_path):
    f = tmp_path / "s.bed"
    f.write_text("chr1\t100\t300\t12.5\t9.0\tchr1:180-200\n")
    ivs = sc.load_seacr_bed(str(f))
    assert ivs == [("chr1", 100, 300)]

def test_reproducible_union_majority():
    # condition X, 3 reps, min_reps=2: peak present in >=2 reps kept
    groups = {"X": ["a", "b", "c"]}
    method = {"X": "majority"}
    peaks = {
        "a": [("chr1", 100, 200)],
        "b": [("chr1", 150, 250)],           # overlaps a
        "c": [("chr2", 500, 600)],           # alone
    }
    out = sc.reproducible_union(groups, method, peaks, min_reps=2)
    # chr1 region reproducible (a,b overlap); chr2 not (only c)
    assert any(c == "chr1" for c, s, e in out)
    assert not any(c == "chr2" for c, s, e in out)

def test_reproducible_union_single_passthrough():
    groups = {"Y": ["a"]}
    method = {"Y": "single"}
    peaks = {"a": [("chr3", 10, 20)]}
    out = sc.reproducible_union(groups, method, peaks, min_reps=2)
    assert ("chr3", 10, 20) in out

def test_merge_and_filter():
    ivs = [("chr1", 100, 200), ("chr1", 150, 300), ("chrM", 0, 50)]
    merged = sc.merge_and_filter(ivs, keep_regex=r"^chr([1-9]|1[0-9]|2[0-2]|X)$",
                                 blacklist={})
    assert ("chr1", 100, 300) in merged
    assert not any(c == "chrM" for c, s, e in merged)
```

- [ ] **Step 2: Run to verify failure**

Run: `pytest tests/test_seacr_consensus.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'seacr_consensus'`.

- [ ] **Step 3: Implement `workflow/scripts/seacr_consensus.py`**

```python
#!/usr/bin/env python3
"""Overlap-based reproducible variable-width union of SEACR peaks → BED + SAF.

SEACR emits enrichment domains (no summit / p-value), so reproducibility is
overlap-based, not IDR/fixed-width. Per condition: keep a peak if it overlaps
peaks in >=K replicates (K = min_reps for >=3 reps; 2-of-2 for 2 reps;
passthrough for 1 rep). Union across conditions, merge overlapping, drop
blacklist / off-target chroms.
"""
import re
from pathlib import Path
from collections import defaultdict


def load_seacr_bed(path):
    """SEACR .bed: chrom start end total_signal max_signal max_region. Return
    [(chrom, start, end)]."""
    out = []
    for line in Path(path).read_text().splitlines():
        if not line.strip():
            continue
        f = line.split("\t")
        out.append((f[0], int(f[1]), int(f[2])))
    return out


def _overlaps(a, intervals):
    c, s, e = a
    for c2, s2, e2 in intervals:
        if c2 == c and s2 < e and s < e2:
            return True
    return False


def reproducible_union(groups, method, peaks_by_sample, min_reps):
    """Return list of (chrom,start,end) reproducible peaks unioned across
    conditions (not yet merged)."""
    kept = []
    for g, members in groups.items():
        m = method[g]
        rep_lists = [peaks_by_sample.get(s, []) for s in members]
        if m == "single":
            for s in members:
                kept.extend(peaks_by_sample.get(s, []))
            continue
        k = min_reps if m == "majority" else 2
        for i, s in enumerate(members):
            others = [rep_lists[j] for j in range(len(rep_lists)) if j != i]
            for p in peaks_by_sample.get(s, []):
                cover = 1 + sum(1 for o in others if _overlaps(p, o))
                if cover >= k:
                    kept.append(p)
    return kept


def merge_and_filter(intervals, keep_regex, blacklist):
    """Sort, chrom/blacklist filter, then merge overlapping into non-overlapping
    (chrom,start,end)."""
    keep_re = re.compile(keep_regex)
    bl = blacklist or {}

    def _bl_hit(c, s, e):
        for bs, be in bl.get(c, []):
            if bs < e and s < be:
                return True
        return False

    ivs = sorted((c, s, e) for (c, s, e) in intervals
                 if keep_re.fullmatch(c) and not _bl_hit(c, s, e))
    merged = []
    for c, s, e in ivs:
        if merged and merged[-1][0] == c and s <= merged[-1][2]:
            pc, ps, pe = merged[-1]
            merged[-1] = (pc, ps, max(pe, e))
        else:
            merged.append((c, s, e))
    return merged


def _load_blacklist(path):
    by = defaultdict(list)
    if not path:
        return by
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        f = line.split("\t")
        by[f[0]].append((int(f[1]), int(f[2])))
    return by


def write_bed(consensus, path):
    lines = [f"{c}\t{s}\t{e}\tseacr_consensus_peak_{i}\t.\t."
             for i, (c, s, e) in enumerate(consensus, 1)]
    Path(path).write_text("\n".join(lines) + ("\n" if lines else ""))


def write_saf(consensus, path):
    lines = ["GeneID\tChr\tStart\tEnd\tStrand"]
    for i, (c, s, e) in enumerate(consensus, 1):
        lines.append(f"seacr_consensus_peak_{i}\t{c}\t{s + 1}\t{e}\t.")
    Path(path).write_text("\n".join(lines) + "\n")


if "snakemake" in globals():  # pragma: no cover
    sm = snakemake  # noqa: F821
    _groups = dict(sm.params.groups)
    _method = dict(sm.params.group_method)
    _paths = dict(sm.params.seacr_paths)  # sample -> bed path
    _peaks = {s: load_seacr_bed(p) for s, p in _paths.items()}
    _union = reproducible_union(_groups, _method, _peaks, int(sm.params.min_reps))
    _blacklist = _load_blacklist(str(sm.input.blacklist))
    _consensus = merge_and_filter(_union, sm.params.keep_regex, _blacklist)
    write_bed(_consensus, str(sm.output.bed))
    write_saf(_consensus, str(sm.output.saf))
```

- [ ] **Step 4: Run tests to verify pass**

Run: `pytest tests/test_seacr_consensus.py -v`
Expected: 4 passed.

- [ ] **Step 5: Append `seacr_consensus_peaks` + `count_fragments_seacr_consensus` rules**

```python
rule seacr_consensus_peaks:
    input:
        seacr = [seacr_peak(s) for s in CONTROLLED_SAMPLES],
        blacklist = config["blacklist"]
    output:
        bed = f"{SEACR_CONSENSUS_DIR}/consensus_peaks.bed",
        saf = f"{SEACR_CONSENSUS_DIR}/consensus_peaks.saf"
    params:
        groups = {g: [s for s in m if s in CONTROLLED_SAMPLES] for g, m in GROUPS.items()},
        group_method = GROUP_METHOD,
        seacr_paths = {s: seacr_peak(s) for s in CONTROLLED_SAMPLES},
        min_reps = config["consensus_min_replicates"],
        keep_regex = config["keep_chroms_regex"]
    conda: "../envs/snakemake.yaml"
    log: "logs/seacr_consensus/consensus.log"
    script: "../scripts/seacr_consensus.py"

rule count_fragments_seacr_consensus:
    input:
        saf = f"{SEACR_CONSENSUS_DIR}/consensus_peaks.saf",
        bams = expand(f"{BLACKLIST_FILTERED_DIR}/{{s}}.nobl.bam", s=TREATMENT_SAMPLES),
        bais = expand(f"{BLACKLIST_FILTERED_DIR}/{{s}}.nobl.bam.bai", s=TREATMENT_SAMPLES)
    output:
        counts = f"{SEACR_CONSENSUS_DIR}/consensus_counts.txt",
        summary = f"{SEACR_CONSENSUS_DIR}/consensus_counts.txt.summary"
    threads: 8
    conda: "../envs/snakemake.yaml"
    log: "logs/seacr_consensus_counts/featurecounts.log"
    shell:
        r"""
        mkdir -p {SEACR_CONSENSUS_DIR} logs/seacr_consensus_counts
        featureCounts -F SAF -a {input.saf} -p --countReadPairs \
            -T {threads} -o {output.counts} {input.bams} > {log} 2>&1
        """
```
Note: `params.groups` filters to CONTROLLED_SAMPLES because only treatments with an IgG get SEACR peaks; a condition whose members lack SEACR peaks is dropped from the SEACR consensus.

- [ ] **Step 6: Extend `cutandrun_all`**

Append:
```python
        f"{SEACR_CONSENSUS_DIR}/consensus_peaks.bed",
        f"{SEACR_CONSENSUS_DIR}/consensus_counts.txt",
```

- [ ] **Step 7: Dry-run SEACR consensus**

Run: `snakemake -s workflow/Snakefile -d .test -n results/consensus_seacr/consensus_counts.txt`
Expected: DAG builds; `seacr_consensus_peaks` consumes the SEACR beds; `count_fragments_seacr_consensus` over 3 treatment BAMs; no exceptions.

- [ ] **Step 8: Commit**

```bash
git add -A && git commit -m "feat: SEACR overlap-based consensus + featureCounts matrix"
```

---

## Task 11: Blacklist stats + finalize `cutandrun_all`

**Files:**
- Modify: `workflow/rules/cutandrun.smk`

- [ ] **Step 1: Append `blacklist_stats`**

Copy ATAC `blacklist_stats` (`atacseq.smk` 604–619) verbatim; it already uses `DEDUP_DIR`/`BLACKLIST_FILTERED_DIR` and `params.samples = SAMPLES`.

- [ ] **Step 2: Add its target to `cutandrun_all`**

Append: `f"{QC_DIR}/blacklist_filtering_stats.txt",`

- [ ] **Step 3: Dry-run the whole primary pipeline**

Run: `snakemake -s workflow/Snakefile -d .test -n cutandrun_all`
Expected: DAG builds end-to-end for all 4 samples across preprocessing, tracks, MACS2 (narrow+broad), SEACR, both consensus matrices, blacklist stats; no exceptions and no missing-rule errors.

- [ ] **Step 4: Commit**

```bash
git add -A && git commit -m "feat: blacklist stats; finalize cutandrun_all target"
```

---

## Task 12: QC pipeline rules (deepTools, FRiP, IDR, complexity, annotations, peak summary, MultiQC)

**Files:**
- Modify: `workflow/rules/qc.smk`

**Interfaces:**
- Produces: all QC outputs listed in `qc_all` except the interactive report (Task 13).

- [ ] **Step 1: Copy the ATAC QC rules that are unchanged**

From `$ATAC/workflow/rules/qc.smk`, copy verbatim into `qc.smk` (above the `qc_all` rule): `deeptools_bedgraph`, `deeptools_fragmentsize`, `deeptools_plotfingerprint`, `deeptools_cor_multibam`, `deeptools_cor_scatterplot`, `deeptools_cor_heatmap`, `deeptools_cor_pca`, `deeptools_gc_bias`, `deeptools_tss`, `deeptools_tss_heatmap_downsample`, `tss_enrichment_score`, `qc_relaxed_peaks`, `idr`, `calculate_library_complexity`, `reads_in_annotations`, `multiqc_fastqc`, and `localrules: multiqc_fastqc`. These all reference `RMD_BAM_DIR`, `DEEPTOOLS_DIR`, `EGS`, `SAMPLES`, `IDR_PAIRS`, `PROMOTER_BED`, `ENHANCER_BED`, `GENOME_2BIT`, `GTF_FILE` — all defined in `common.smk`.

- [ ] **Step 2: Add the FRiP rule (MACS2 + SEACR, per-sample ext)**

Replace ATAC's wildcard `{condition}` FRiP with an explicit per-treatment approach. Add:
```python
rule frip_macs2:
    input:
        bamfile = os.path.join(RMD_BAM_DIR, "{sample}.nobl.bam"),
        peakfile = lambda w: macs2_peak(w.sample)
    output:
        fripfile = os.path.join(FRIP_DIR, "{sample}.macs2.frip.txt")
    conda: "../envs/bedtools.yaml"
    log: "logs/FRiP/{sample}.macs2.log"
    shell:
        r"""
        mkdir -p {FRIP_DIR} logs/FRiP
        total=$(samtools view -c {input.bamfile})
        in_peaks=$(bedtools intersect -u -abam {input.bamfile} -b {input.peakfile} | samtools view -c) 2> {log}
        frip=$(echo "scale=4; $in_peaks / $total" | bc)
        echo -e "{wildcards.sample}\t$in_peaks\t$total\t$frip" > {output.fripfile}
        """

rule frip_seacr:
    wildcard_constraints:
        sample = _alt(CONTROLLED_SAMPLES)
    input:
        bamfile = os.path.join(RMD_BAM_DIR, "{sample}.nobl.bam"),
        peakfile = lambda w: seacr_peak(w.sample)
    output:
        fripfile = os.path.join(FRIP_DIR, "{sample}.seacr.frip.txt")
    conda: "../envs/bedtools.yaml"
    log: "logs/FRiP/{sample}.seacr.log"
    shell:
        r"""
        mkdir -p {FRIP_DIR} logs/FRiP
        total=$(samtools view -c {input.bamfile})
        in_peaks=$(bedtools intersect -u -abam {input.bamfile} -b {input.peakfile} | samtools view -c) 2> {log}
        frip=$(echo "scale=4; $in_peaks / $total" | bc)
        echo -e "{wildcards.sample}\t$in_peaks\t$total\t$frip" > {output.fripfile}
        """
```

- [ ] **Step 3: Add the peak-summary rule (MACS2 + SEACR)**

Adapt ATAC `peak_summary` (`qc.smk` 568–598) to write two tables — one for MACS2 (`peak_summary_macs2.tsv/_mqc.txt`, iterating `TREATMENT_SAMPLES`, using `macs2_peak(s)` and `{s}.macs2.frip.txt`) and one for SEACR (`peak_summary_seacr.tsv/_mqc.txt`, iterating `CONTROLLED_SAMPLES`, using `seacr_peak(s)` and `{s}.seacr.frip.txt`; SEACR bed width is `$3-$2`, and FRiP from the seacr frip file). Provide both as two rules `peak_summary_macs2` and `peak_summary_seacr` mirroring the ATAC shell (the awk width/FRiP logic is identical; only the input paths and column source differ).

- [ ] **Step 4: Write the `qc_all` aggregate (minus the report)**

```python
rule qc_all:
    input:
        expand(os.path.join(BEDGRAPH_DIR, "{s}.nobl.RPGC.bedgraph"), s=SAMPLES),
        os.path.join(DEEPTOOLS_DIR, "fragmentSize.png"),
        os.path.join(DEEPTOOLS_DIR, "ATACseq_fingerprint.png"),
        os.path.join(DEEPTOOLS_DIR, "deeptools_scatterplot.png"),
        os.path.join(DEEPTOOLS_DIR, "deeptools_heatmap.png"),
        os.path.join(DEEPTOOLS_DIR, "deeptools_PCA.png"),
        expand(os.path.join(DEEPTOOLS_DIR, "{s}.gc_content.png"), s=SAMPLES),
        os.path.join(DEEPTOOLS_DIR, "Heatmap_TSS.png"),
        os.path.join(DEEPTOOLS_DIR, "Profile_TSS.png"),
        os.path.join(DEEPTOOLS_DIR, "tss_heatmap_downsampled.json"),
        os.path.join(QC_DIR, "tss_enrichment_mqc.txt"),
        expand(os.path.join(FRIP_DIR, "{s}.macs2.frip.txt"), s=TREATMENT_SAMPLES),
        expand(os.path.join(FRIP_DIR, "{s}.seacr.frip.txt"), s=CONTROLLED_SAMPLES),
        [os.path.join(IDR_DIR, f"{g}--{r1}--{r2}--idr_peaks.{('broadPeak' if g in BROAD_IDR_GROUPS else 'narrowPeak')}.txt")
         for g, r1, r2 in IDR_PAIRS],
        expand(os.path.join(COMPLEXITY_DIR, "{s}_complexity.txt"), s=SAMPLES),
        os.path.join(ANNOT_DIR, "reads_in_annotations_mqc.txt"),
        os.path.join(QC_DIR, "peak_summary_macs2_mqc.txt"),
        os.path.join(QC_DIR, "peak_summary_seacr_mqc.txt"),
        os.path.join(QC_DIR, "multiqc_fastqc.html"),
        os.path.join(QC_DIR, "cutandrun_qc_report.html"),
```
The `idr` rule's `condition` wildcard: the ATAC `idr` rule takes a `{condition}` wildcard (narrowPeak/broadPeak) — the target list above supplies the right ext per group, so the rule works unchanged.

- [ ] **Step 5: Dry-run the QC pipeline (report target will fail until Task 13 — test the rest)**

Run: `snakemake -s workflow/Snakefile -d .test -n results/qc/peak_summary_macs2_mqc.txt results/qc/peak_summary_seacr_mqc.txt results/FRiP/test_cjun_r1.macs2.frip.txt results/FRiP/test_cjun_r1.seacr.frip.txt results/deeptools/Heatmap_TSS.png`
Expected: DAG builds for these targets; no exceptions.

- [ ] **Step 6: Commit**

```bash
git add -A && git commit -m "feat: QC pipeline (deepTools, FRiP macs2+seacr, IDR, complexity, peak summary)"
```

---

## Task 13: Interactive QC report (adapt `build_qc_report.py`) + `qc_report` rule

**Files:**
- Create: `workflow/scripts/build_qc_report.py` (adapted copy)
- Modify: `workflow/rules/qc.smk`

**Interfaces:**
- Produces: `{QC_DIR}/cutandrun_qc_report.html`.

- [ ] **Step 1: Copy then read the ATAC report builder**

```bash
cp /easley/scratch/projects/amitra/amitra2016502/snakemake_ATACseq_spikein/workflow/scripts/build_qc_report.py workflow/scripts/build_qc_report.py
```
Read it fully to locate spike-in sections (search `spikein`, `normalization_factors`, `spikein_fraction`, `norm_factor`).

- [ ] **Step 2: Remove spike-in sections from the report builder**

Delete/guard every code path that loads or renders: `results/spikein/normalization_factors.tsv`, `results/spikein_qc/spikein_fraction.tsv`, and any "Spike-in %" / "normalization factor" panel. Where a table row iterates these, remove that row. Keep alignment/mito/dup/blacklist/FRiP/TSS/complexity/fragment-size/GC/correlation/PCA/fingerprint/reads-in-annotation/consensus panels. Make missing files non-fatal (wrap each loader in a `try/except FileNotFoundError` that skips the panel) so the report renders without the removed inputs.

- [ ] **Step 3: Add MACS2-vs-SEACR peak panel**

Add a panel that reads `results/qc/peak_summary_macs2.tsv` and `results/qc/peak_summary_seacr.tsv` and renders both peak-count/width/FRiP tables side by side. Reuse the builder's existing table-render helper (same one used by the ATAC peak_summary panel). Retitle the report "CUT&RUN QC Report" and change the default output name references to `cutandrun_qc_report.html`.

- [ ] **Step 4: Add the `qc_report` rule**

Adapt ATAC `qc_report` (`qc.smk` 623–666): drop the spike-in inputs (`normalization_factors.tsv`, `spikein_fraction.tsv`); change `peak_summary.tsv` → both `peak_summary_macs2.tsv` and `peak_summary_seacr.tsv`; output `os.path.join(QC_DIR, "cutandrun_qc_report.html")`; keep the deepTools PNG/JSON/TSV inputs. Invoke `python workflow/scripts/build_qc_report.py --results-dir {params.results} --out {output.html} --samples {params.samples} --generated ...`.

- [ ] **Step 5: Dry-run the report + full QC**

Run: `snakemake -s workflow/Snakefile -d .test -n qc_all`
Expected: DAG builds end-to-end including `qc_report`; no exceptions.

- [ ] **Step 6: Smoke-test the report builder on synthetic inputs (optional but recommended)**

Create a throwaway `results/` tree with minimal valid TSVs (a couple of rows each) under a temp dir and run `python workflow/scripts/build_qc_report.py --results-dir <tmp> --out /tmp/r.html --samples a,b --generated x`. Expected: exits 0 and writes an HTML file. (Skip if too costly; the dry-run already validates wiring.)

- [ ] **Step 7: Commit**

```bash
git add -A && git commit -m "feat: interactive CUT&RUN QC report (spike-in removed, SEACR added)"
```

---

## Task 14: Generic DESeq2 differential-binding notebook

**Files:**
- Create: `workflow/scripts/build_diffbind_notebook.py` (adapted), `workflow/scripts/diffbind_helpers.R` (adapted)
- Generate: `cutandrun_Dx.ipynb`

**Interfaces:**
- Consumes: `results/consensus/consensus_counts.txt` and/or `results/consensus_seacr/consensus_counts.txt`, `config/samples.csv`, `promoter_bed`.
- Produces: `cutandrun_Dx.ipynb` (executed downstream in the `cutandrun_Dx` env, out of the Snakemake DAG — mirrors ATAC's notebook workflow).

- [ ] **Step 1: Copy the ATAC generators**

```bash
A=/easley/scratch/projects/amitra/amitra2016502/snakemake_ATACseq_spikein
cp $A/workflow/scripts/build_diffbind_notebook.py workflow/scripts/build_diffbind_notebook.py
cp $A/workflow/scripts/diffbind_helpers.R workflow/scripts/diffbind_helpers.R
```
Read both fully.

- [ ] **Step 2: Rewrite `build_diffbind_notebook.py` for generic CUT&RUN DESeq2**

Change the notebook it emits so that it:
- Reads `config/samples.csv` (columns condition/replicate), restricting to treatment samples (non-empty peak_mode); builds a `colData` frame from the featureCounts matrix columns.
- Selects the count matrix by a `differential_counts` argument (`macs2`|`seacr`|`both`) → `results/consensus/consensus_counts.txt` and/or `results/consensus_seacr/consensus_counts.txt`; when `both`, loops the analysis over each.
- Runs DESeq2 with **default median-of-ratios** normalization (NO spike-in size factors); design `~ condition`, upgraded to `~ replicate + condition` when every condition has the same replicate labels (balanced paired design).
- Iterates **all pairwise treatment-condition contrasts** (drop IgG conditions — they have no treatment rows so are already absent from `colData`).
- Splits consensus regions into promoter vs distal by overlap with `promoter_bed` and tests each subset; emits results + significant subset (`padj<0.05 & |log2FC|>1`), MA + volcano + PCA (VST), ChIPseeker nearest-gene annotation → `results/diff_region/`.
- Drops the Notch/HES1/Gviz positive-control section entirely.
- Targets the `cutandrun_Dx` kernel/env; output filename `cutandrun_Dx.ipynb`.

- [ ] **Step 3: Trim `diffbind_helpers.R`**

Remove any spike-in helper (functions referencing `NF`, `sizeFactor = 1/NF`, spike-in normalization). Keep the DESeq2/promoter-split/annotation/plot helpers. Ensure no remaining function references spike-in inputs.

- [ ] **Step 4: Generate and validate the notebook**

Run: `python workflow/scripts/build_diffbind_notebook.py` (or with its expected args) then
`python -c "import nbformat; nbformat.read('cutandrun_Dx.ipynb', as_version=4); print('valid notebook')"`
Expected: `valid notebook`.

- [ ] **Step 5: Commit**

```bash
git add -A && git commit -m "feat: generic DESeq2 differential-binding notebook (no spike-in)"
```

---

## Task 15: Top-level docs, run script, catalog metadata; full-pipeline dry-run

**Files:**
- Create: `README.md`, `run_pipeline.sh`, `.snakemake-workflow-catalog.yml`, `.test/README.md`

- [ ] **Step 1: Write `.snakemake-workflow-catalog.yml`**

Adapt `$ATAC/.snakemake-workflow-catalog.yml`: keep `--use-conda` mandatory-flag; update the description to CUT&RUN (no spike-in), list reference downloads (genome FASTA, ENCODE blacklist, GENCODE GTF, hg38.2bit, picard.jar — drop dm6). Mention targets `cutandrun_all` / `qc_all`.

- [ ] **Step 2: Write `run_pipeline.sh`**

Adapt `$ATAC/run_pipeline.sh`: replace image name and `atacseq_all` references with `cutandrun`/`cutandrun_all`. (If Docker images are not built, this is a convenience wrapper only; note in README it requires building the image — Docker packaging is out of scope for v1.)

- [ ] **Step 3: Write `README.md`**

Write a CUT&RUN-specific README covering: overview (MACS2+SEACR, IgG controls, narrow/broad, IgG-subtracted tracks, dual consensus matrices, generic DESeq2), the 6-column sample sheet, config highlights, reference-data download/generate steps (genome, blacklist, GTF, 2bit, picard.jar; promoter/enhancer BEDs shipped), run commands (`snakemake --use-conda --cores N`, `cutandrun_all`, `qc_all`), output-directory map (`results/{...,seacr,consensus,consensus_seacr,log2ratio_bigwig,...}`), and the differential-notebook usage. Cite MACS2, SEACR (Meers et al. 2019), deepTools, Bowtie2, samtools, IDR, featureCounts, Corces-2018 consensus, bedtools, FastQC, fastp, DESeq2.

- [ ] **Step 4: Write `.test/README.md`**

One paragraph: this directory is a dry-run harness (placeholder ref/data) used by `snakemake -n -d .test <target>` to validate the DAG without real genomes; it is not an execution dataset.

- [ ] **Step 5: Full-pipeline dry-run (primary + QC)**

Run: `snakemake -s workflow/Snakefile -d .test -n`
Expected: "Building DAG of jobs..." then a full job table covering `cutandrun_all` + `qc_all` for all 4 test samples; **no** exceptions, **no** missing-input/rule errors.

- [ ] **Step 6: Rulegraph sanity check**

Run: `snakemake -s workflow/Snakefile -d .test --rulegraph cutandrun_all qc_all > /tmp/rg.dot && head -1 /tmp/rg.dot`
Expected: prints `digraph snakemake_dag {` (valid DOT), confirming the DAG is acyclic and complete.

- [ ] **Step 7: Commit**

```bash
git add -A && git commit -m "docs: README, run script, catalog metadata; full-pipeline dry-run green"
```

---

## Self-Review

**Spec coverage** (each spec section → task):
- §2 sample sheet + derived sets + validation → Tasks 2, 3.
- §3 removed/changed/added → Tasks 2 (config), 5 (index/align/filter/dedup), 6 (log2), 7 (MACS2 control), 8 (SEACR), 9 (consensus narrow/broad), 10 (SEACR consensus).
- §4 primary rules 1–20 → Tasks 5 (1–8 incl. chrom sizes + index + dedup + blacklist), 6 (9–10 tracks), 7 (11–12 MACS2), 8 (13–14 SEACR), 9 (15–17 MACS2 reproducibility/consensus/counts), 10 (18–19 SEACR consensus/counts), 11 (20 blacklist stats).
- §5 QC (spike-in removed, SEACR added) → Tasks 12, 13.
- §6 differential notebook (macs2/seacr/both) → Task 14; `differential_counts` key → Task 2.
- §7 config schema (renames/removes/adds) → Task 2.
- §8 envs (seacr new, Dx rename, reuse) → Task 1.
- §9 reference data (chrom.sizes built, shipped BEDs) → Tasks 1, 5.
- §10 scripts (reuse/adapt/new/drop) → Tasks 1, 9, 10, 13, 14.
- §11 testing (dry-run harness, rulegraph) → Task 1 harness + per-task dry-runs + Task 15.

**Placeholder scan:** No "TBD/TODO". Task 13 (report) and Task 14 (notebook) give targeted edit instructions against concrete copied files rather than full reprints — the exact search terms and required behaviors are specified; these are diffs on 841/244-line files where reprinting verbatim would be error-prone. All novel modules (`samplesheet.py`, `seacr_consensus.py`, consensus `load_peaks`) and all new rules are shown in full.

**Type/name consistency:** `SampleSheet` API (Task 3) matches `common.smk` usage (Task 4): `input_control`, `macs2_ext`, `seacr_stringency(sample, cfg)`, `peak_mode`, `groups`, `group_method`. Helpers `macs2_peak`, `seacr_peak`, `igg_bam` defined in Task 4 and used in Tasks 6–13. `load_peaks` renamed consistently in Task 9 (script + snakemake block + rule params `narrowpeak_paths`/`idr_paths`). SEACR script `seacr_paths`/`groups`/`group_method`/`min_reps`/`keep_regex` params match rule (Task 10). `BROAD_IDR_GROUPS`/`NARROW_IDR_GROUPS` introduced in Task 9 and reused in Task 12 `qc_all` — both tasks must define them in `common.smk` (Task 9 Step 5 adds them to common.smk).
