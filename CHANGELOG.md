# Changelog

All notable changes to this workflow are documented here. This file is managed
by [release-please](https://github.com/googleapis/release-please) from
[Conventional Commit](https://www.conventionalcommits.org/) messages — a release
PR bumps the version and prepends the sections below on each merge to `main`.

## [0.2.1](https://github.com/gynecoloji/snakemake_CutandRunseq/compare/v0.2.0...v0.2.1) (2026-07-30)


### Documentation

* add snakevision workflow tube-map (images/rulegraph.svg) + README diagram section ([5591c99](https://github.com/gynecoloji/snakemake_CutandRunseq/commit/5591c99c4542d5f38844121e007d2cb8209339ec))
* correct repo slug to snakemake_CutandRunseq in CITATION.cff + release-please package-name ([b1d5368](https://github.com/gynecoloji/snakemake_CutandRunseq/commit/b1d536855f24eb854743aab94bc9f0f17a3ca4de))

## [0.2.0](https://github.com/gynecoloji/snakemake_CutandRunseq/compare/v0.1.0...v0.2.0) (2026-07-30)


### Added

* adapt build_diffopen_report.py (per-caller, none/anchor/rnastable) ([cc3cb21](https://github.com/gynecoloji/snakemake_CutandRunseq/commit/cc3cb219112c06a9f9e9a9aa77ac084616cda379))
* adapt diffopen.R for CUT&RUN (condition design, drop spikein, ctcf-&gt;anchor) ([f11334a](https://github.com/gynecoloji/snakemake_CutandRunseq/commit/f11334a4217d546ae663bffc7e18e649140da1c9))
* blacklist stats; finalize cutandrun_all target (primary pipeline complete) ([7e4092d](https://github.com/gynecoloji/snakemake_CutandRunseq/commit/7e4092d881875f4fe533a59c315d5aa0db10d840))
* common.smk diffopen constants + helpers (caller x mode) ([95eb7d8](https://github.com/gynecoloji/snakemake_CutandRunseq/commit/95eb7d86cb41d7b7ff89e5f37170b7db95ab4d98))
* config — retire differential_counts, add diffopen keys ([8969cb6](https://github.com/gynecoloji/snakemake_CutandRunseq/commit/8969cb66474582c2914dd24dfc8e28e61b7f99e2))
* config control_type + dual-control sample sheets ([cd61010](https://github.com/gynecoloji/snakemake_CutandRunseq/commit/cd6101003699abd45768c3168965f514ef5477bc))
* config schema, default config, sample sheet, config README ([cbb49b3](https://github.com/gynecoloji/snakemake_CutandRunseq/commit/cbb49b32c3cc9212ce7a0148c7f065e2e2b4c7a7))
* deepTools fingerprint quality metrics (JS distance vs IgG, % enriched) ([c28a07b](https://github.com/gynecoloji/snakemake_CutandRunseq/commit/c28a07b6b94faf233be9aa9b3b755d4d6e9b84f7))
* diffopen.smk stage (caller x mode) wired opt-in into the Snakefile ([2d879e7](https://github.com/gynecoloji/snakemake_CutandRunseq/commit/2d879e7407f0776db2a3b5ab7eb860b1c12c245c))
* downstream.smk (annotation/motifs/overlap/heatmap/metagene), opt-in ([c9bb7a5](https://github.com/gynecoloji/snakemake_CutandRunseq/commit/c9bb7a53cf8d53c3a02bfc7cd732b8e24ecc527f))
* ENCODE cross-correlation QC (phantompeakqualtools NSC/RSC) ([7dc7175](https://github.com/gynecoloji/snakemake_CutandRunseq/commit/7dc7175658aeb0936081977c44a4b3d25ec15d72))
* ENCODE IDR self-consistency/rescue reproducibility QC ([2f878eb](https://github.com/gynecoloji/snakemake_CutandRunseq/commit/2f878eb99f700ce2890689f1e279469910872b66))
* generic DESeq2 differential-binding notebook (no spike-in, pairwise contrasts) ([8eb2921](https://github.com/gynecoloji/snakemake_CutandRunseq/commit/8eb292112cc2d84fb7bb7dc0a79fd16ad26b1dec))
* interactive CUT&RUN QC report (spike-in removed, SEACR added) + qc_report rule ([9cbc716](https://github.com/gynecoloji/snakemake_CutandRunseq/commit/9cbc7169a4468ac49e7e14709366641b1eb9190e))
* MACS2 narrow/broad peak calling with per-sample IgG control ([e1bbd48](https://github.com/gynecoloji/snakemake_CutandRunseq/commit/e1bbd48471e08c5d0bac5a974dca10d62ba5c12f))
* MACS2 reproducibility + fixed-width consensus + counts (narrow/broad) ([e751a33](https://github.com/gynecoloji/snakemake_CutandRunseq/commit/e751a338f4ea7c03cb163babe17ba807e4cda273))
* per-IP-vs-control fingerprint JSD (replaces single-reference metric) ([6883392](https://github.com/gynecoloji/snakemake_CutandRunseq/commit/6883392a2b0f0278d2da97b86cee5ad10eae2972))
* preprocessing rules (index, chrom sizes, align, filter, dedup, blacklist) ([316a362](https://github.com/gynecoloji/snakemake_CutandRunseq/commit/316a3628d06d86d0884bd503f8f71f9dd79eaf56))
* QC pipeline (deepTools, FRiP macs2+seacr, IDR, complexity, peak summaries) ([19b48b6](https://github.com/gynecoloji/snakemake_CutandRunseq/commit/19b48b6c27bd64467da6047a0ac92e5ceb112164))
* QC report NSC/RSC + fingerprint-quality sections ([52d19c4](https://github.com/gynecoloji/snakemake_CutandRunseq/commit/52d19c4ff6af32642db07be81b15b5130b4d19d0))
* resolve MACS2/SEACR/track control via control_type (input+igg) ([3612360](https://github.com/gynecoloji/snakemake_CutandRunseq/commit/36123604d3a34213f8a856fbb9f5ad019f495eb1))
* RPGC bigWig + IgG-subtracted log2 bamCompare tracks ([50f555b](https://github.com/gynecoloji/snakemake_CutandRunseq/commit/50f555b5fe3e68c1346e73c5d2ba5878086bbcb5))
* samplesheet.py dual-control (input+igg) resolution + tests ([c63529f](https://github.com/gynecoloji/snakemake_CutandRunseq/commit/c63529fd20469c24c5c606ae46eebe664a5a8843))
* samplesheet.py sample-set parsing + validation with tests ([fcb0c30](https://github.com/gynecoloji/snakemake_CutandRunseq/commit/fcb0c305e0cf931cf8eb377673dde2245abac15a))
* SEACR overlap-based consensus + featureCounts matrix ([f39ef66](https://github.com/gynecoloji/snakemake_CutandRunseq/commit/f39ef666dabc2a43252b6babf279ebce4954150b))
* SEACR peak calling (fragment bedgraph + IgG control, stringent/relaxed) ([208719f](https://github.com/gynecoloji/snakemake_CutandRunseq/commit/208719fccf08a02a41ee8525eb56c1cdee87721a))
* Snakefile + common.smk (sample sets, dirs, helpers) + rule skeletons ([1908f12](https://github.com/gynecoloji/snakemake_CutandRunseq/commit/1908f129dcae1863d302798abb060ac7f74a1a1a))


### Fixed

* control-free pseudo-rep relaxed peaks so IDR rescue ratio is consistent ([64110e3](https://github.com/gynecoloji/snakemake_CutandRunseq/commit/64110e3f99e51a34e88573eeb2d33db05d8b5f51))


### Documentation

* ChIP-seq-mirror additions in README/catalog; full dry-run green ([2fec94c](https://github.com/gynecoloji/snakemake_CutandRunseq/commit/2fec94ca12300997c9cbeb26bc39db23e3270491))
* design spec for ChIP-seq mirror (dual control, ENCODE reproducibility, downstream module) ([f7ff904](https://github.com/gynecoloji/snakemake_CutandRunseq/commit/f7ff9043d6b5dd0c20034d5786e1a3aaf271d49e))
* design spec for ENCODE QC metrics + diffopen differential-binding stage ([0e74387](https://github.com/gynecoloji/snakemake_CutandRunseq/commit/0e74387c62331e64858b24af5022cbb4b7cd3e4a))
* diffopen + ENCODE QC in README/catalog; full dry-run green ([05fc7c6](https://github.com/gynecoloji/snakemake_CutandRunseq/commit/05fc7c688cb2dea22ce9b1956851bf5186a57930))
* fix IDR-reproducibility scope wording (2-rep conditions) ([3350f75](https://github.com/gynecoloji/snakemake_CutandRunseq/commit/3350f7598c95e793a45a760aed9ae7302ed91afd))
* implementation plan for ChIP-seq mirror ([edae9f3](https://github.com/gynecoloji/snakemake_CutandRunseq/commit/edae9f309bccbe439ae5e024ef00a74ca5829f4d))
* implementation plan for ENCODE QC + diffopen stage ([d1dbc6a](https://github.com/gynecoloji/snakemake_CutandRunseq/commit/d1dbc6af3aaccfb1d8145dce9129242674927035))
* README, run script, catalog metadata; full-pipeline dry-run green ([0f51066](https://github.com/gynecoloji/snakemake_CutandRunseq/commit/0f51066337f2f81769e9d1f8e38fe065dad176f2))
* tidy diffopen.R read_design docstring ([c20cc29](https://github.com/gynecoloji/snakemake_CutandRunseq/commit/c20cc2946ce04bf8ad8fc93e1b6b85f331a8b002))

## 0.1.0

Initial release of the CUT&RUN Snakemake workflow: FastQC/fastp → human Bowtie2
alignment → filter/dedup/blacklist → MACS2 & SEACR peak calling (matched
IgG/Input controls) → reproducible consensus peaks + featureCounts matrices →
RPGC and IgG-subtracted tracks; an ENCODE-style QC stage (FRiP, NSC/RSC
cross-correlation, per-IP fingerprint JSD, IDR self-consistency reproducibility,
library complexity, TSS, correlation/PCA, interactive HTML report); an opt-in
DESeq2 differential-binding stage (`diffopen_all`); and an opt-in downstream
module (`downstream_all`: ChIPseeker annotation, HOMER motifs, peak
Jaccard/overlap, deepTools peak/gene-body heatmaps).
