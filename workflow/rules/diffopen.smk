# ── Differential binding (OPT-IN) ────────────────────────────────────────
# Not part of `rule all`: a differential test needs >=2 conditions. Request it:
#   snakemake --use-conda --cores N diffopen_all
# Runs each configured normalization (diffopen_modes) on each configured caller's
# consensus matrix (diffopen_callers). Shared config/helpers live in common.smk.

rule diffopen:
    wildcard_constraints:
        caller="macs2|seacr",
        mode="none|anchor|rnastable",
    input:
        unpack(_diffopen_extra_input),
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
