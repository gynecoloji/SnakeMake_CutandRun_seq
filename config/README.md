# Configuration

This workflow is configured through two files in this directory:

- `config.yaml` — all workflow parameters (see below)
- `samples.csv` — the sample sheet

plus reference data you download into `ref/` (see [Reference data](#reference-data)).

## Sample sheet (`config/samples.csv`)

CSV with one row per sample and these columns:

| column          | description                                                                                     |
|-----------------|-------------------------------------------------------------------------------------------------|
| `sample_id`     | Sample name. Raw reads must be `data/<sample_id>_R1_001.fastq.gz` / `_R2_001.fastq.gz`.          |
| `condition`     | Free-text label. For **treatment** rows this is the replicate group used for reproducibility.    |
| `replicate`     | Replicate index within the condition (integer).                                                 |
| `input_control` | `sample_id` of the matched Input control for this row. Leave empty for control rows.              |
| `igg_control`   | `sample_id` of the matched IgG control for this row. Leave empty for control rows.                |
| `peak_mode`     | `narrow` or `broad`. **Empty marks the row as a control (IgG/Input)** — it is not peak-called.    |
| `notes`         | Free text; ignored by the pipeline.                                                             |

The effective control used for peak calling (MACS2 `-c`, SEACR control track, `bamCompare` `-b2`)
for each treatment row is chosen by `control_type` in `config.yaml` (`input` or `igg`, default
`igg`): the named column is used if non-empty, otherwise the pipeline falls back to the other
column.

Example:

```csv
sample_id,condition,replicate,input_control,igg_control,peak_mode,notes
GSF2801-ChIPseq-OVCAR3-3D-IP-cJun_S4,cJUN_3D,1,,GSF2801-ChIPseq-OVCAR3-3D-IP-IgG_S5,narrow,3D-cJUN
GSF2801-ChIPseq-OVCAR3-3D-IP-IgG_S5,IgG_3D,1,,,,3D-Igg
GSF2801-ChIPseq-OVCAR3-Control-IP-cJun_S1,cJUN_Ctrl,1,,GSF2801-ChIPseq-OVCAR3-Control-IP-IgG_S2,narrow,Ctrl-cJUN
GSF2801-ChIPseq-OVCAR3-Control-IP-IgG_S2,IgG_Ctrl,1,,,,Ctrl-Igg
```

### Treatment vs control

- A row is a **control** (IgG/Input) when `peak_mode` is empty. Controls are aligned,
  filtered, deduplicated, blacklist-filtered and turned into RPGC bigWigs, and are used as the
  MACS2 `-c` control, the SEACR control track, and the `bamCompare` `-b2` for their matched
  treatments — but they are never peak-called themselves.
- Every other row is a **treatment**. Its `peak_mode` (`narrow`/`broad`) selects MACS2 narrow vs
  `--broad` and the SEACR stringency; its `input_control`/`igg_control` name the controls to pair
  with it (see `control_type` above for which one is used).

### Per-condition rules

- **`condition` is the reproducibility group.** Reproducibility handling is derived from the
  number of replicates in each treatment condition:
  - **≥ 3 replicates** → majority vote (kept if a peak recurs in ≥ `consensus_min_replicates`).
  - **exactly 2 replicates** → IDR (`idr_threshold`) for MACS2; 2-of-2 overlap for SEACR.
  - **1 replicate** → the sample's own peaks are used as-is.
- **All replicates of one condition must share the same `peak_mode`** (a condition is either narrow
  or broad; consensus/IDR cannot mix the two). The workflow errors out if they differ.
- **Give biologically distinct groups distinct `condition` labels.** If the same antibody target
  was profiled in two contexts (e.g. cJUN in "3D" and "Control"), label them `cJUN_3D` and
  `cJUN_Ctrl` — otherwise the two single-replicate rows would be treated as two replicates of one
  condition and (incorrectly) run through IDR.
- Each treatment's `input_control`/`igg_control` (whichever is filled in) must reference an
  existing control (`peak_mode`-empty) `sample_id`.

## Parameters (`config/config.yaml`)

Every parameter — with its type, default, and description — is defined once in the config schema,
[`workflow/schemas/config.schema.yaml`](../workflow/schemas/config.schema.yaml). That schema is the
single source of truth: the workflow validates `config.yaml` against it on every run (filling in
defaults for anything you omit).

To configure a run, edit `config.yaml` directly. At minimum, point the reference-file paths
(`genome_fasta`, `blacklist`, `gtf`, `promoter_bed`, `enhancer_bed`) at the files you provide. CUT&RUN
specifics worth reviewing: `max_fragment_length` (Bowtie2 `-X`, default 700), `remove_duplicates`
(set `false` to keep duplicates for low-input libraries), `control_type` (`input`|`igg`, default
`igg` — which sample-sheet control column drives peak calling, falling back to the other), and the
`macs2_*` and `seacr_*` peak-calling knobs. Differential binding is the opt-in `diffopen_all` target (see the top-level
README), keyed off `diffopen_callers` (which consensus matrices), `diffopen_modes`
(`none`/`anchor`/`rnastable`), and `diffopen_ref_label` (the reference `condition`).

## Reference data

Genomes, indexes and large annotations are **not** shipped in the repo. Download / place them under
`ref/` before running, matching the paths in `config.yaml`:

- `ref/genome.fa` — chr-prefixed UCSC human genome (hg38)
- `ref/hg38_blacklist_regions.bed` — ENCODE hg38 blacklist (shipped)
- `ref/gencode.v36.annotation.gtf` — GENCODE annotation (for TSS QC)
- `ref/hg38.2bit` — for `computeGCBias`
- `ref/picard.jar` — Picard (used by MarkDuplicates)
- `ref/promoter_chr1-22X.bed`, `ref/enhancer_chr1-22X.bed` — Ensembl Regulatory Build (shipped)

The Bowtie2 index (`ref/genome/`) and `ref/genome.chrom.sizes` are built automatically by the
`build_genome_index` / `genome_chrom_sizes` rules.

See the top-level `README.md` for full setup and run instructions.
