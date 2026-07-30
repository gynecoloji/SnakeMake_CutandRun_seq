# Changelog

All notable changes to this workflow are documented here. This file is managed
by [release-please](https://github.com/googleapis/release-please) from
[Conventional Commit](https://www.conventionalcommits.org/) messages — a release
PR bumps the version and prepends the sections below on each merge to `main`.

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
