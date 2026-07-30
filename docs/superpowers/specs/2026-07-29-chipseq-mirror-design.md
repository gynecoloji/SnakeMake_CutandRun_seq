# CUT&RUN update: mirror ChIP-seq (dual control, ENCODE reproducibility, downstream module) — Design

**Date:** 2026-07-29
**Working dir:** `/easley/scratch/projects/amitra/amitra2016502/snakemake_CutandRun_seq`
**Mirrors:** `../snakemake_ChIP_seq` (its `rules/downstream.smk`, the `repro_*`/`fingerprint_jsd`
rules in `rules/qc.smk`, `scripts/{idr_reproducibility.py,peak_annotation.R,peak_overlap.py}`,
`envs/{chipseeker,homer,phantompeakqualtools}.yaml`, and the dual-control sample-sheet logic in
`rules/common.smk`).

The ChIP-seq primary pipeline is already a subset of CUT&RUN's (CUT&RUN adds SEACR), so nothing
in the primary stage changes. This update brings four ChIP-seq capabilities across, all confirmed
with the user:

- **Part A — Dual-control sample sheet:** add an `igg_control` column beside the existing
  `input_control`, plus a `control_type` (input|igg) config with fallback.
- **Part B — ENCODE IDR self-consistency / rescue-ratio reproducibility QC** (pooled & self
  pseudo-replicates per 2-replicate condition).
- **Part C — Per-IP-vs-control fingerprint JSD** (replaces the current single-reference
  fingerprint-quality metric with the more correct per-IP-vs-its-own-control JSD).
- **Part D — Downstream module** (opt-in `downstream_all`): ChIPseeker annotation + GO, HOMER
  motif enrichment, peak Jaccard/overlap matrix, deepTools peak-heatmap + gene-body metagene.
  ChIP's DESeq2 `differential_binding` is **not** ported — the existing `diffopen` stage supersedes it.

Ordering: Part A is the foundation (B/C/D all consume the resolved control), so it lands first.

---

## Part A — Dual-control sample sheet

### Sample sheet
`config/samples.csv` gains an **`igg_control`** column. Final columns:
`sample_id,condition,replicate,input_control,igg_control,peak_mode,notes`.

- `input_control` — sample_id of the matched **Input** control (sonicated chromatin), or empty.
- `igg_control` — sample_id of the matched **IgG** control, or empty.
- A row with empty `peak_mode` is a control-only sample (aligned + tracks, never peak-called),
  as today.

### Control resolution
New config `control_type` (enum `input|igg`, **default `igg`** — IgG is the standard CUT&RUN
control). For an IP sample the effective MACS2/SEACR/track control is: the column named by
`control_type`; if empty, fall back to the other column; else no control. Mirrors ChIP's
`_resolve_control`.

### `samplesheet.py` (extend the pure module)
- Parse `igg_control` (add to `REQUIRED_COLUMNS`; NaN→"").
- `input_control(sample)` / `igg_control(sample)` accessors.
- `resolved_control(sample, control_type)` → the effective control sample_id or "".
- `SampleSheet.validate()` also checks every `input_control` **and** `igg_control` reference
  resolves to an existing sample.
- Backward compatibility: a sheet with only `input_control` filled (today's example) still works —
  with `control_type=igg`, `igg_control` empty → falls back to `input_control`.
- Unit tests (`tests/test_samplesheet.py`): dual-control resolution (input primary, igg primary,
  fallback each way, none), and validation of a bad `igg_control` reference.

### `common.smk`
- Read `CONTROL_TYPE = config["control_type"]`.
- Replace `igg_bam(sample)` with `control_bam(sample)` = blacklist-filtered BAM of
  `SS.resolved_control(sample, CONTROL_TYPE)` ("" if none), and `control_arg(sample)` = `-c <bam>`
  or "".
- `CONTROLLED_SAMPLES` = IP (treatment) samples with a resolved control (drives log2 tracks,
  SEACR, fingerprint JSD, and the pseudo-rep control reuse).

### Consumers to update (rename `igg_bam`→`control_bam`, `input_control`→resolved)
- `cutandrun.smk`: `call_peaks_macs2_narrow`/`broad` (`-c` via `control_arg`),
  `create_log2ratio_bigwig` (`-b2` = `control_bam`), `seacr_bedgraph`/`call_peaks_seacr_*`
  (control bedgraph = resolved control), and the `CONTROLLED_SAMPLES`/`SEACR_*` sets.
- `diffopen.smk` does not use the control directly (it reads the count matrices) — unchanged.

### Config
Add `control_type: "igg"` to `config/config.yaml`, `.test/config/config.yaml`, and the schema
(enum input|igg, default igg, required). Migrate the shipped `config/samples.csv` and
`.test/config/samples.csv` to add the `igg_control` column with the current IgG moved there
(keeping `input_control` empty for IP rows), so the sheets are unambiguous under `control_type=igg`.

---

## Part B — ENCODE IDR self-consistency reproducibility (`qc.smk`)

Copy ChIP's pseudo-replicate reproducibility, adapted to CUT&RUN names/config. Assesses each
2-replicate condition with self-pseudoreplicates (each replicate split in half) and pooled
pseudoreplicates (both replicates pooled then split), plus the true-replicate IDR count, and
reports the ENCODE self-consistency and rescue ratios.

### `common.smk` additions
- `REPRO_DIR = f"{RESULT_DIR}/idr_reproducibility"`.
- `REPRO_UNITS = [f"self__{s}" for s in IDR_SAMPLES] + [f"pool__{g}" for g in IDR_GROUPS]`,
  split into `REPRO_NARROW_UNITS`/`REPRO_BROAD_UNITS` by the unit's condition mode.
- Helpers `pseudo_source_bam`, `unit_control_bam`, `unit_control_arg` (reuse the unit's
  representative sample's resolved control), `idr_peak_file(group)` (Module-B IDR peaks, already
  produced by `reproducible_idr_*`). `IDR_SAMPLES`/`IDR_GROUPS`/`GROUPS`/`GROUP_MODE` already exist.

### `qc.smk` rules (copied from ChIP, adapted)
`repro_pool` (samtools merge the 2 replicate BAMs), `repro_split` (deterministic read-name-hash
split into pr1/pr2), `repro_relaxed_narrow`/`repro_relaxed_broad` (MACS2 relaxed peaks per
pseudo-half, reusing the unit's control; uses `macs2_genome`/`macs2_broad_cutoff`/
`idr_relaxed_pvalue`/`idr_top_n_peaks`), `repro_idr_narrow`/`repro_idr_broad` (IDR between the two
halves → reproducible-peak count), and `idr_reproducibility_summary` (script
`idr_reproducibility.py`) → `results/qc/idr_reproducibility.tsv` + `_mqc.txt` with per-condition
Nt/Np/N1/N2, self-consistency ratio, rescue ratio, and an ENCODE pass/borderline/fail flag.

### Script
`workflow/scripts/idr_reproducibility.py` — copy verbatim (a `script:`-directive consumer of
`snakemake.params` groups/members/true_idr + the per-unit `.n.txt` counts).

### Wiring
Add `results/qc/idr_reproducibility.tsv` to `qc_all` (guarded: only meaningful when
`IDR_GROUPS` is non-empty; when there are no 2-rep conditions the target list is empty and the
stage is skipped).

---

## Part C — Per-IP-vs-control fingerprint JSD (`qc.smk`)

Replace the current bolt-on (the `deeptools_plotfingerprint` rule's `--outQualityMetrics` +
`--JSDsample <first control>`) with ChIP's dedicated per-IP rule.

- **`deeptools_plotfingerprint`** → revert to plot + raw-counts only (drop the `metrics` output,
  `--outQualityMetrics`, and the `jsd` param). It stays the all-sample signal-to-noise plot.
- **`fingerprint_jsd`** (new, per `CONTROLLED_SAMPLES`): `plotFingerprint -b <IP> <control>
  --JSDsample <control> --outQualityMetrics …` → per-IP JS distance vs **its own** control +
  % genome enriched + AUC. `JSD_DIR = f"{RESULT_DIR}/qc_fingerprint"`.
- **`fingerprint_jsd_summary`** → `results/qc/fingerprint_jsd.tsv` + `_mqc.txt`.
- **Report** (`build_qc_report.py`): the `fingerprint_quality` section now reads
  `results/qc/fingerprint_jsd.tsv` (columns sample/js_distance/pct_genome_enriched/auc) instead of
  `results/deeptools/fingerprint_quality_metrics.tab`. Keep it informational (no hard flag).
- **Wiring:** in `qc_all` and the `qc_report` inputs, replace
  `results/deeptools/fingerprint_quality_metrics.tab` with `results/qc/fingerprint_jsd.tsv`.

---

## Part D — Downstream module (`workflow/rules/downstream.smk`, opt-in)

New opt-in stage `downstream_all` (included in the Snakefile but, like `diffopen`, **not** wired
into `rule all` — keeps the default run light; ChIP puts it in `all`, we deliberately don't).
IP samples = `TREATMENT_SAMPLES`; peak files = the MACS2 per-sample peaks (`macs2_peak`); the
consensus is the MACS2 consensus (`results/consensus/`). Rules (copied from ChIP, adapted to
CUT&RUN names — `genome_fasta` for HOMER, `macs2_peak`/`PEAK_MODE` via `SS`):

- `annotate_peaks` (per IP sample): ChIPseeker annotation (feature distribution, dist-to-TSS) + GO
  → `results/annotation/<sample>/`. Script `peak_annotation.R`; env `chipseeker.yaml`.
- `motif_enrichment` (per IP sample): HOMER `findMotifsGenome.pl` on `cut -f1-6` of the peak file,
  `-size 200` for narrow / `-size given` for broad → `results/motifs/<sample>/`. Env `homer.yaml`.
- `peak_jaccard` + `peak_overlap_matrix`: pairwise bedtools Jaccard across IP samples → long TSV →
  `results/peak_overlap/jaccard_matrix.tsv` + heatmap. Script `peak_overlap.py`; guarded to run only
  with ≥2 IP samples.
- `deeptools_peak_heatmap`: `computeMatrix reference-point --referencePoint center` over the MACS2
  consensus peaks with all-sample RPGC bigWigs → `Heatmap_peaks.png`/`Profile_peaks.png`.
- `deeptools_metagene`: `computeMatrix scale-regions --metagene` over the GTF → `Heatmap_genebody.png`/
  `Profile_genebody.png`.
- `downstream_all`: aggregates the above (annotation + GO + motifs per IP sample; jaccard when ≥2
  IP; peak heatmap + metagene). **No** `differential_binding`.

### Scripts / envs / dirs
- Copy verbatim: `workflow/scripts/peak_annotation.R`, `workflow/scripts/peak_overlap.py`.
- New envs: `workflow/envs/chipseeker.yaml`, `workflow/envs/homer.yaml` (copy verbatim).
- Update `workflow/envs/phantompeak.yaml` to match ChIP's fuller `phantompeakqualtools.yaml`
  (add `r-catools r-snow r-snowfall r-bitops gawk`) so `run_spp.R` has all its runtime deps.
- `common.smk` dirs: `ANNOTATION_DIR`, `MOTIF_DIR`, `OVERLAP_DIR` (drop ChIP's `DIFFBIND_DIR` —
  not porting differential_binding).
- `common.smk` helper `all_peak_files()` = `[macs2_peak(s) for s in TREATMENT_SAMPLES]`;
  `peak_file(sample)` alias = `macs2_peak(sample)` (ChIP scripts/rules call `peak_file`).

### Snakefile
Add `include: "rules/downstream.smk"` after the diffopen include; leave `rule all` unchanged
(opt-in via `downstream_all`).

---

## Envs summary

- **New (copy verbatim):** `chipseeker.yaml`, `homer.yaml`.
- **Update:** `phantompeak.yaml` — add `r-catools`, `r-snow`, `r-snowfall`, `r-bitops`, `gawk`.
- Reuse: `deeptools.yaml` (peak heatmap/metagene, fingerprint JSD, jaccard-matrix plot),
  `macs2.yaml` (repro relaxed peaks), `idr.yaml` (repro IDR), `bedtools.yaml` (jaccard),
  `snakemake.yaml` (pool/split/summaries).

## Scripts summary

- Copy verbatim: `peak_annotation.R`, `peak_overlap.py`, `idr_reproducibility.py`.
- Adapt: `samplesheet.py` (dual control), `build_qc_report.py` (fingerprint section source +
  new `idr_reproducibility` section).
- Unchanged: everything else (diffopen.*, consensus_peaks.py, seacr_consensus.py, …).

## Config summary (`config/config.yaml` + schema + `.test`)

Add: `control_type: "igg"` (enum input|igg). Sample sheets gain the `igg_control` column.
No new keys for reproducibility (reuses `idr_*`), fingerprint JSD, or downstream (HOMER/ChIPseeker
use `genome_fasta`/`gtf`, already present). No `contrasts` key (differential_binding not ported).

## QC report (`build_qc_report.py`) integration

- **Fingerprint section:** repoint to `results/qc/fingerprint_jsd.tsv`.
- **New `idr_reproducibility` section:** read `results/qc/idr_reproducibility.tsv` (per-condition
  self-consistency + rescue ratios), flagged pass/borderline/fail per the ENCODE convention (both
  ratios < 2 = pass). Added to the report's section order + UNITS map.

## Testing

- `tests/test_samplesheet.py`: dual-control resolution + validation (extend existing).
- `tests/test_config_schema.py`: default + `.test` configs validate with `control_type` +
  `igg_control` sheets (stays green).
- Dry-runs (`snakemake -n -d .test`): `qc_all` (now incl. `fingerprint_jsd`, `repro_*`,
  `idr_reproducibility_summary`), `downstream_all` (annotation/motif/jaccard/heatmap/metagene
  resolve for the 2 IP conditions), and the default target (unchanged; excludes diffopen +
  downstream). Smoke-test `build_qc_report.py` with synthetic `fingerprint_jsd.tsv` +
  `idr_reproducibility.tsv`. Rscript-parse `peak_annotation.R`; py-parse `peak_overlap.py` /
  `idr_reproducibility.py`.
- `.test` harness: add the `igg_control` column to `.test/config/samples.csv`; the existing 2-rep
  `cJUN` condition already exercises the reproducibility path. No new `.test/ref` placeholders are
  needed — the downstream rules reuse the existing `genome.fa`/`gencode.gtf`/`consensus` placeholders,
  and dry-run never reads their contents.

## Out of scope (v1)

- ChIP's DESeq2 `differential_binding` rule + `contrasts` config (diffopen supersedes it).
- Wiring `downstream_all` into the default `rule all` (kept opt-in, like diffopen).
- Running downstream/reproducibility on the SEACR peak set (uses the MACS2 backbone, as elsewhere).
- Real-data execution (needs the chipseeker/homer/phantompeak/macs2/idr conda envs + real BAMs).

## Notes

- `control_type` defaults to `igg` (CUT&RUN norm), unlike ChIP's `input` default.
- Downstream annotation/motifs run on **MACS2** per-sample peaks and the MACS2 consensus, matching
  the "MACS2 is the backbone" convention already established in this workflow.
- The reproducibility and fingerprint-JSD stages need a resolved control and ≥2-rep conditions
  respectively; both degrade to empty target lists (skipped) when the sample sheet doesn't supply
  them, so single-replicate / control-less runs stay valid.
