# CUT&RUN QC pipeline: deepTools coverage/correlation/fingerprint/GC/TSS QC,
# FRiP (MACS2 + SEACR), IDR on relaxed peaks, library complexity, reads-in-annotation,
# peak summaries, a FastQC-only MultiQC, and a self-contained interactive HTML report.
#
# Consumes the primary pipeline's results/. Shared config/dirs/helpers live in common.smk.

localrules: multiqc_fastqc

# ── deepTools coverage/QC + TSS (verbatim from ATAC) ─────────────────────
# 1. RPGC bedgraph per sample (coverage QC; main pipeline makes bigWigs, not bedgraphs)
rule deeptools_bedgraph:
    input:
        bam = os.path.join(RMD_BAM_DIR, "{sample}.nobl.bam")
    output:
        bedgraph = os.path.join(BEDGRAPH_DIR, "{sample}.nobl.RPGC.bedgraph")
    threads: 8
    conda:
        "../envs/deeptools.yaml"
    log:
        "logs/deeptools_bedgraph/{sample}.log"
    shell:
        """
        mkdir -p {BEDGRAPH_DIR} logs/deeptools_bedgraph
        bamCoverage -p {threads} \
            --outFileFormat bedgraph \
            --effectiveGenomeSize {EGS} \
            --normalizeUsing RPGC \
            --binSize 10 --extendReads \
            --bam {input.bam} -o {output.bedgraph} 2> {log}
        """


# 2. Fragment size distribution (ATAC nucleosome periodicity)
rule deeptools_fragmentsize:
    input:
        bams = expand(os.path.join(RMD_BAM_DIR, "{sample}.nobl.bam"), sample=SAMPLES)
    output:
        plot = os.path.join(DEEPTOOLS_DIR, "fragmentSize.png"),
        table = os.path.join(DEEPTOOLS_DIR, "fragmentsize.txt"),
        raw = os.path.join(DEEPTOOLS_DIR, "fragment_lengths.txt")
    threads: 12
    conda:
        "../envs/deeptools.yaml"
    log:
        "logs/deeptools_fragmentsize/fragmentsize.log"
    shell:
        """
        mkdir -p {DEEPTOOLS_DIR} logs/deeptools_fragmentsize
        bamPEFragmentSize -p {threads} \
            -hist {output.plot} \
            -T "Fragment size of PE ATACseq data" \
            --maxFragmentLength 1500 \
            -b {input.bams} \
            --outRawFragmentLengths {output.raw} \
            --table {output.table} 2> {log}
        """


# 3. Fingerprint (signal-to-noise)
rule deeptools_plotfingerprint:
    input:
        bams = expand(os.path.join(RMD_BAM_DIR, "{sample}.nobl.bam"), sample=SAMPLES)
    output:
        plot = os.path.join(DEEPTOOLS_DIR, "ATACseq_fingerprint.png"),
        table = os.path.join(DEEPTOOLS_DIR, "ATACseq_fingerprint.tab")
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
            --outRawCounts {output.table} 2> {log}
        """


# 4. Correlation / PCA across samples
rule deeptools_cor_multibam:
    input:
        bams = expand(os.path.join(RMD_BAM_DIR, "{sample}.nobl.bam"), sample=SAMPLES)
    output:
        npz = os.path.join(DEEPTOOLS_DIR, "deeptools_multiBAM.out.npz"),
        counts = os.path.join(DEEPTOOLS_DIR, "deeptools_readCounts.tab")
    threads: 12
    conda:
        "../envs/deeptools.yaml"
    log:
        "logs/deeptools_correlation/multibam.log"
    shell:
        """
        mkdir -p {DEEPTOOLS_DIR} logs/deeptools_correlation
        multiBamSummary bins \
            -bs 5000 \
            --ignoreDuplicates \
            -p {threads} \
            --bamfiles {input.bams} \
            -out {output.npz} \
            --outRawCounts {output.counts} 2> {log}
        """


rule deeptools_cor_scatterplot:
    input:
        npz = os.path.join(DEEPTOOLS_DIR, "deeptools_multiBAM.out.npz")
    output:
        plot = os.path.join(DEEPTOOLS_DIR, "deeptools_scatterplot.png")
    threads: 4
    conda:
        "../envs/deeptools.yaml"
    log:
        "logs/deeptools_correlation/scatterplot.log"
    shell:
        """
        plotCorrelation --corData {input.npz} \
            --whatToPlot scatterplot \
            --skipZero \
            --plotTitle "Scatterplot" \
            --plotFileFormat png \
            --corMethod spearman \
            --log1p \
            --plotFile {output.plot} 2> {log}
        """


rule deeptools_cor_heatmap:
    input:
        npz = os.path.join(DEEPTOOLS_DIR, "deeptools_multiBAM.out.npz")
    output:
        plot = os.path.join(DEEPTOOLS_DIR, "deeptools_heatmap.png"),
        cormat = os.path.join(DEEPTOOLS_DIR, "correlation_matrix.tab")
    threads: 4
    conda:
        "../envs/deeptools.yaml"
    log:
        "logs/deeptools_correlation/heatmap.log"
    shell:
        """
        plotCorrelation --corData {input.npz} \
            --whatToPlot heatmap \
            --skipZero \
            --plotTitle "Heatmap" \
            --plotFileFormat png \
            --corMethod spearman \
            --log1p \
            --outFileCorMatrix {output.cormat} \
            --plotFile {output.plot} 2> {log}
        """


rule deeptools_cor_pca:
    input:
        npz = os.path.join(DEEPTOOLS_DIR, "deeptools_multiBAM.out.npz")
    output:
        plot = os.path.join(DEEPTOOLS_DIR, "deeptools_PCA.png"),
        data = os.path.join(DEEPTOOLS_DIR, "deeptools_PCA.tab")
    threads: 4
    conda:
        "../envs/deeptools.yaml"
    log:
        "logs/deeptools_correlation/pca.log"
    shell:
        """
        plotPCA --corData {input.npz} \
            --plotTitle "PCA" \
            --plotFileFormat png \
            --ntop 1000 \
            --plotFile {output.plot} \
            --outFileNameData {output.data} 2> {log}
        """


# 5. GC bias
rule deeptools_gc_bias:
    input:
        bam = os.path.join(RMD_BAM_DIR, "{sample}.nobl.bam"),
        genome = GENOME_2BIT
    output:
        freq = os.path.join(DEEPTOOLS_DIR, "{sample}.gc_content.txt"),
        plot = os.path.join(DEEPTOOLS_DIR, "{sample}.gc_content.png")
    threads: 8
    conda:
        "../envs/deeptools.yaml"
    log:
        "logs/deeptools_gc_bias/{sample}.log"
    shell:
        """
        mkdir -p {DEEPTOOLS_DIR} logs/deeptools_gc_bias
        computeGCBias -b {input.bam} \
            --effectiveGenomeSize {EGS} \
            -p {threads} \
            --genome {input.genome} \
            -freq {output.freq} \
            --biasPlot {output.plot} \
            --plotFileFormat png 2> {log}
        """


# 6. TSS enrichment: heatmap + profile + profile-data table (reuses the main
#    pipeline's RPGC bigWigs instead of regenerating them).
rule deeptools_tss:
    input:
        bigwigs = expand(os.path.join(BIGWIG_DIR, "{sample}.bw"), sample=SAMPLES),
        gtf = GTF_FILE
    output:
        matrix = os.path.join(DEEPTOOLS_DIR, "matrix.mat.gz"),
        heatmap = os.path.join(DEEPTOOLS_DIR, "Heatmap_TSS.png"),
        profile = os.path.join(DEEPTOOLS_DIR, "Profile_TSS.png"),
        profiledata = os.path.join(DEEPTOOLS_DIR, "Profile_TSS.data.tab")
    threads: 16
    conda:
        "../envs/deeptools.yaml"
    log:
        "logs/deeptools_tss/tss.log"
    shell:
        """
        mkdir -p {DEEPTOOLS_DIR} logs/deeptools_tss
        computeMatrix reference-point \
            -p {threads} \
            --referencePoint TSS \
            -S {input.bigwigs} \
            -R {input.gtf} \
            -a 2000 -b 2000 \
            --skipZeros \
            -o {output.matrix} 2> {log}
        plotHeatmap \
            -m {output.matrix} \
            --dpi 300 \
            --zMin -3 --zMax 3 \
            --heatmapWidth 20 \
            -out {output.heatmap} \
            --plotFileFormat png \
            --sortUsing mean 2>> {log}
        plotProfile \
            -m {output.matrix} \
            --dpi 300 \
            -out {output.profile} \
            --plotFileFormat png \
            --outFileNameData {output.profiledata} 2>> {log}
        """


# 6b. NEW: downsampled TSS heatmap matrix (JSON) for the interactive report's
#     canvas heatmap; reuses the existing matrix.mat.gz (no computeMatrix rerun).
rule deeptools_tss_heatmap_downsample:
    input:
        matrix = os.path.join(DEEPTOOLS_DIR, "matrix.mat.gz")
    output:
        json = os.path.join(DEEPTOOLS_DIR, "tss_heatmap_downsampled.json")
    conda:
        "../envs/deeptools.yaml"
    log:
        "logs/deeptools_tss/downsample.log"
    shell:
        """
        python workflow/scripts/downsample_tss_matrix.py {input.matrix} \
            -o {output.json} --nrows 180 --ncols 80 > {log} 2>&1
        """


# 7. NEW: numeric TSS enrichment score per sample (from the profile-data table)
rule tss_enrichment_score:
    input:
        profile = os.path.join(DEEPTOOLS_DIR, "Profile_TSS.data.tab")
    output:
        tsv = os.path.join(QC_DIR, "tss_enrichment_scores.tsv"),
        mqc = os.path.join(QC_DIR, "tss_enrichment_mqc.txt")
    conda:
        "../envs/snakemake.yaml"
    log:
        "logs/tss_enrichment_score/tss.log"
    script:
        "../scripts/tss_score.py"


# ── IDR on within-condition replicate pairs (verbatim from ATAC) ─────────
# 10. IDR on within-condition replicate pairs (relaxed peaks)
rule idr:
    input:
        rep1 = lambda w: os.path.join(RELAXED_DIR, f"{w.rep1}_relaxed.{w.condition}"),
        rep2 = lambda w: os.path.join(RELAXED_DIR, f"{w.rep2}_relaxed.{w.condition}")
    output:
        peaks = os.path.join(IDR_DIR, "{group}--{rep1}--{rep2}--idr_peaks.{condition}.txt")
    conda:
        "../envs/idr.yaml"
    wildcard_constraints:
        condition = "narrowPeak|broadPeak"
    log:
        "logs/idr/{group}--{rep1}--{rep2}--idr_{condition}.log"
    shell:
        """
        mkdir -p {IDR_DIR} logs/idr
        idr --samples {input.rep1} {input.rep2} \
            --input-file-type {wildcards.condition} \
            --rank p.value \
            --output-file {output.peaks} \
            --plot \
            --log-output-file {log}
        """


# ── Library complexity (verbatim from ATAC) ──────────────────────────────
# 11. Library complexity (NRF / PBC1 / PBC2) on the pre-dedup filtered BAM
rule calculate_library_complexity:
    input:
        bam = os.path.join(FILTERED_DIR, "{sample}.sorted.filtered.bam")
    output:
        txt = os.path.join(COMPLEXITY_DIR, "{sample}_complexity.txt")
    conda:
        "../envs/bedtools.yaml"  # needs samtools, bedtools, bc
    threads: 8
    log:
        "logs/library_complexity/{sample}.log"
    params:
        temp_dir = os.path.join(TMP_DIR, "{sample}_complexity"),
        tmp = os.path.join(TMP_DIR, "{sample}_complexity.bed")
    shell:
        """
        mkdir -p {COMPLEXITY_DIR} logs/library_complexity {params.temp_dir}

        # Extract fragments (BEDPE) from the name-sorted BAM
        echo "Extracting fragments from BAM file..." >> {log}
        samtools sort -n -@ {threads} {input.bam} | \
        bedtools bamtobed -bedpe -i stdin > {params.tmp} 2>> {log}
        # position-only key (chrom, start, end) — NO read name: PCR duplicates share
        # coordinates but have different names, so the name must be dropped for
        # `uniq -c` to collapse them into one location with the right multiplicity.
        awk 'BEGIN {{OFS="\t"}} {{
            if ($1==$4) {{
                start = ($2 < $5) ? $2 : $5;
                end = ($3 > $6) ? $3 : $6;
                print $1, start, end;
            }}
        }}' {params.tmp} | \
        sort -k1,1 -k2,2n -k3,3n > {params.temp_dir}/fragments.bed

        fragment_count=$(wc -l < {params.temp_dir}/fragments.bed)
        echo "Total extracted fragments: $fragment_count" >> {log}

        # Fragment counts by genomic location (for PCR-duplicate stats)
        sort -k1,1 -k2,2n -k3,3n {params.temp_dir}/fragments.bed | \
        uniq -c > {params.temp_dir}/fragment_counts.txt

        unique=$(wc -l < {params.temp_dir}/fragment_counts.txt)
        one_read=$(awk '$1 == 1' {params.temp_dir}/fragment_counts.txt | wc -l)
        two_reads=$(awk '$1 == 2' {params.temp_dir}/fragment_counts.txt | wc -l)
        total_reads=$(samtools view -c {input.bam})

        echo "Unique locations: $unique" >> {log}
        echo "Locations with exactly one fragment: $one_read" >> {log}
        echo "Locations with exactly two fragments: $two_reads" >> {log}
        echo "Total mapped reads: $total_reads" >> {log}

        if [ "$one_read" -eq 0 ]; then
            echo "WARNING: No unique fragments found, setting to 1 to prevent division by zero" >> {log}
            one_read=1
        fi

        nrf=$(echo "scale=6; $unique / $fragment_count" | bc)
        pbc1=$(echo "scale=6; $one_read / $unique" | bc)
        if [ "$two_reads" -eq 0 ]; then
            two_reads=1
            echo "CRITICAL: two_reads is 0 before PBC2 calculation, forcing to 1" >> {log}
        fi
        pbc2=$(echo "scale=6; $one_read / $two_reads" | bc)

        echo "## Library Complexity Metrics for {wildcards.sample} ##" > {output.txt}
        echo -e "Total Reads\t$total_reads" >> {output.txt}
        echo -e "Total Fragments\t$fragment_count" >> {output.txt}
        echo -e "Distinct Fragment Locations (Nd)\t$unique" >> {output.txt}
        echo -e "Locations with 1 Fragment (N1)\t$one_read" >> {output.txt}
        echo -e "Locations with 2 Fragments (N2)\t$two_reads" >> {output.txt}
        echo -e "NRF (Nd/Total)\t$nrf" >> {output.txt}
        echo -e "PBC1 (N1/Nd)\t$pbc1" >> {output.txt}
        echo -e "PBC2 (N1/N2)\t$pbc2" >> {output.txt}

        echo -e "\n## Quality Assessment ##" >> {output.txt}
        if (( $(echo "$nrf > 0.9" | bc -l) )); then echo -e "NRF: $nrf - High complexity (>0.9)" >> {output.txt}
        elif (( $(echo "$nrf > 0.8" | bc -l) )); then echo -e "NRF: $nrf - Good complexity (0.8-0.9)" >> {output.txt}
        elif (( $(echo "$nrf > 0.7" | bc -l) )); then echo -e "NRF: $nrf - Moderate complexity (0.7-0.8)" >> {output.txt}
        else echo -e "NRF: $nrf - Low complexity (<0.7)" >> {output.txt}; fi

        if (( $(echo "$pbc1 > 0.9" | bc -l) )); then echo -e "PBC1: $pbc1 - Near ideal (>0.9)" >> {output.txt}
        elif (( $(echo "$pbc1 > 0.8" | bc -l) )); then echo -e "PBC1: $pbc1 - Good (0.8-0.9)" >> {output.txt}
        elif (( $(echo "$pbc1 > 0.7" | bc -l) )); then echo -e "PBC1: $pbc1 - Moderate (0.7-0.8)" >> {output.txt}
        else echo -e "PBC1: $pbc1 - Severe bottlenecking (<0.7)" >> {output.txt}; fi

        if (( $(echo "$pbc2 > 10" | bc -l) )); then echo -e "PBC2: $pbc2 - Near ideal (>10)" >> {output.txt}
        elif (( $(echo "$pbc2 > 3" | bc -l) )); then echo -e "PBC2: $pbc2 - Good (3-10)" >> {output.txt}
        elif (( $(echo "$pbc2 > 1" | bc -l) )); then echo -e "PBC2: $pbc2 - Moderate (1-3)" >> {output.txt}
        else echo -e "PBC2: $pbc2 - Severe bottlenecking (<1)" >> {output.txt}; fi

        rm -rf {params.temp_dir} {params.tmp}
        echo "Library complexity calculation completed for {wildcards.sample}" >> {log}
        """


# ── Reads in promoters vs enhancers (verbatim from ATAC) ─────────────────
# 13. NEW: reads in promoters vs enhancers (signal distribution)
rule reads_in_annotations:
    input:
        bams = expand(os.path.join(RMD_BAM_DIR, "{sample}.nobl.bam"), sample=SAMPLES),
        promoter = PROMOTER_BED,
        enhancer = ENHANCER_BED
    output:
        tsv = os.path.join(ANNOT_DIR, "reads_in_annotations.tsv"),
        mqc = os.path.join(ANNOT_DIR, "reads_in_annotations_mqc.txt")
    params:
        samples = SAMPLES,
        bamdir = RMD_BAM_DIR
    conda:
        "../envs/bedtools.yaml"
    log:
        "logs/reads_in_annotations/annot.log"
    shell:
        """
        mkdir -p {ANNOT_DIR} logs/reads_in_annotations
        printf "# id: reads_in_annotations\n# section_name: 'Reads in annotations'\n# description: 'Fraction of reads overlapping promoters / enhancers.'\n# plot_type: 'table'\nSample\tTotal\tIn promoter\tIn enhancer\tPromoter frac\tEnhancer frac\n" > {output.mqc}
        echo -e "sample\ttotal\tin_promoter\tin_enhancer\tpromoter_frac\tenhancer_frac" > {output.tsv}
        for s in {params.samples}; do
            bam={params.bamdir}/$s.nobl.bam
            total=$(samtools view -c $bam)
            prom=$(bedtools intersect -u -abam $bam -b {input.promoter} | samtools view -c)
            enh=$(bedtools intersect -u -abam $bam -b {input.enhancer} | samtools view -c)
            awk -v s=$s -v t=$total -v p=$prom -v e=$enh 'BEGIN{{
                pf=(t>0)?p/t:0; ef=(t>0)?e/t:0;
                printf "%s\\t%d\\t%d\\t%d\\t%.4f\\t%.4f\\n", s, t, p, e, pf, ef
            }}' | tee -a {output.tsv} >> {output.mqc}
        done 2> {log}
        """


# ── FastQC-only MultiQC (verbatim from ATAC) ─────────────────────────────
# 15. FastQC-only MultiQC report (everything else moves to the interactive report)
rule multiqc_fastqc:
    input:
        expand(os.path.join(FASTQC_DIR, "{sample}_R1_001_fastqc.html"), sample=SAMPLES),
        expand(os.path.join(FASTQC_DIR, "{sample}_R2_001_fastqc.html"), sample=SAMPLES)
    output:
        html = os.path.join(QC_DIR, "multiqc_fastqc.html")
    conda:
        "../envs/snakemake.yaml"
    log:
        "logs/multiqc_fastqc/multiqc.log"
    shell:
        """
        mkdir -p {QC_DIR} logs/multiqc_fastqc
        multiqc -f -m fastqc {FASTQC_DIR}/ \
            --outdir {QC_DIR} \
            --filename multiqc_fastqc.html > {log} 2>&1
        """


# ── QC relaxed peaks for IDR — narrow (all reps of narrow multi-rep conditions) ─
rule qc_relaxed_peaks_narrow:
    wildcard_constraints:
        sample = _alt(QC_IDR_NARROW_SAMPLES)
    input:
        bam = os.path.join(RMD_BAM_DIR, "{sample}.nobl.bam")
    output:
        peaks = os.path.join(RELAXED_DIR, "{sample}_relaxed.narrowPeak")
    params:
        outdir = RELAXED_DIR, name = "{sample}",
        genome = config["macs2_genome"], pvalue = config["idr_relaxed_pvalue"],
        top_n = config["idr_top_n_peaks"]
    conda:
        "../envs/macs2.yaml"
    log:
        "logs/qc_relaxed_peaks/{sample}.log"
    shell:
        r"""
        mkdir -p {params.outdir} logs/qc_relaxed_peaks
        macs2 callpeak -t {input.bam} -f BAMPE -g {params.genome} \
            --outdir {params.outdir} -n {params.name}_relaxedtmp \
            --nomodel -p {params.pvalue} > {log} 2>&1
        sort -k8,8gr {params.outdir}/{params.name}_relaxedtmp_peaks.narrowPeak \
            > {params.outdir}/{params.name}_relaxedtmp_sorted.narrowPeak
        head -n {params.top_n} {params.outdir}/{params.name}_relaxedtmp_sorted.narrowPeak > {output.peaks}
        rm -f {params.outdir}/{params.name}_relaxedtmp_peaks.narrowPeak \
              {params.outdir}/{params.name}_relaxedtmp_peaks.xls \
              {params.outdir}/{params.name}_relaxedtmp_summits.bed \
              {params.outdir}/{params.name}_relaxedtmp_sorted.narrowPeak
        """


# ── QC relaxed peaks for IDR — broad (all reps of broad multi-rep conditions) ─
rule qc_relaxed_peaks_broad:
    wildcard_constraints:
        sample = _alt(QC_IDR_BROAD_SAMPLES)
    input:
        bam = os.path.join(RMD_BAM_DIR, "{sample}.nobl.bam")
    output:
        peaks = os.path.join(RELAXED_DIR, "{sample}_relaxed.broadPeak")
    params:
        outdir = RELAXED_DIR, name = "{sample}",
        genome = config["macs2_genome"], pvalue = config["idr_relaxed_pvalue"],
        broad_cutoff = config["macs2_broad_cutoff"], top_n = config["idr_top_n_peaks"]
    conda:
        "../envs/macs2.yaml"
    log:
        "logs/qc_relaxed_peaks/{sample}.log"
    shell:
        r"""
        mkdir -p {params.outdir} logs/qc_relaxed_peaks
        macs2 callpeak -t {input.bam} -f BAMPE -g {params.genome} \
            --outdir {params.outdir} -n {params.name}_relaxedtmp \
            --nomodel --broad --broad-cutoff {params.broad_cutoff} -p {params.pvalue} > {log} 2>&1
        sort -k8,8gr {params.outdir}/{params.name}_relaxedtmp_peaks.broadPeak \
            > {params.outdir}/{params.name}_relaxedtmp_sorted.broadPeak
        head -n {params.top_n} {params.outdir}/{params.name}_relaxedtmp_sorted.broadPeak > {output.peaks}
        rm -f {params.outdir}/{params.name}_relaxedtmp_peaks.broadPeak \
              {params.outdir}/{params.name}_relaxedtmp_peaks.gappedPeak \
              {params.outdir}/{params.name}_relaxedtmp_peaks.xls \
              {params.outdir}/{params.name}_relaxedtmp_sorted.broadPeak
        """


# ── FRiP against the MACS2 peak set (per treatment sample) ────────────────
rule frip_macs2:
    input:
        bamfile = os.path.join(RMD_BAM_DIR, "{sample}.nobl.bam"),
        peakfile = lambda w: macs2_peak(w.sample)
    output:
        fripfile = os.path.join(FRIP_DIR, "{sample}.macs2.frip.txt")
    conda:
        "../envs/bedtools.yaml"
    log:
        "logs/FRiP/{sample}.macs2.log"
    shell:
        r"""
        mkdir -p {FRIP_DIR} logs/FRiP
        total=$(samtools view -c {input.bamfile})
        in_peaks=$(bedtools intersect -u -abam {input.bamfile} -b {input.peakfile} | samtools view -c) 2> {log}
        frip=$(echo "scale=4; $in_peaks / $total" | bc)
        echo -e "{wildcards.sample}\t$in_peaks\t$total\t$frip" > {output.fripfile}
        """


# ── FRiP against the SEACR peak set (per controlled treatment sample) ─────
rule frip_seacr:
    wildcard_constraints:
        sample = _alt(CONTROLLED_SAMPLES)
    input:
        bamfile = os.path.join(RMD_BAM_DIR, "{sample}.nobl.bam"),
        peakfile = lambda w: seacr_peak(w.sample)
    output:
        fripfile = os.path.join(FRIP_DIR, "{sample}.seacr.frip.txt")
    conda:
        "../envs/bedtools.yaml"
    log:
        "logs/FRiP/{sample}.seacr.log"
    shell:
        r"""
        mkdir -p {FRIP_DIR} logs/FRiP
        total=$(samtools view -c {input.bamfile})
        in_peaks=$(bedtools intersect -u -abam {input.bamfile} -b {input.peakfile} | samtools view -c) 2> {log}
        frip=$(echo "scale=4; $in_peaks / $total" | bc)
        echo -e "{wildcards.sample}\t$in_peaks\t$total\t$frip" > {output.fripfile}
        """


# ── Peak count / width summary + FRiP — MACS2 (per treatment sample) ──────
rule peak_summary_macs2:
    input:
        peaks = [macs2_peak(s) for s in TREATMENT_SAMPLES],
        frips = expand(os.path.join(FRIP_DIR, "{s}.macs2.frip.txt"), s=TREATMENT_SAMPLES)
    output:
        tsv = os.path.join(QC_DIR, "peak_summary_macs2.tsv"),
        mqc = os.path.join(QC_DIR, "peak_summary_macs2_mqc.txt")
    conda:
        "../envs/snakemake.yaml"
    log:
        "logs/peak_summary/macs2.log"
    shell:
        r"""
        mkdir -p {QC_DIR} logs/peak_summary
        printf "# id: peak_summary_macs2\n# section_name: 'MACS2 peak summary'\n# description: 'MACS2 peak count, width stats and FRiP per treatment sample.'\n# plot_type: 'table'\nSample\tPeaks\tMean width\tMin width\tMax width\tFRiP\n" > {output.mqc}
        echo -e "sample\tn_peaks\tmean_width\tmin_width\tmax_width\tFRiP" > {output.tsv}
        for pk in {input.peaks}; do
            s=$(basename "$pk"); s=${{s%_peaks.*}}
            frip=$(cut -f4 {FRIP_DIR}/$s.macs2.frip.txt)
            awk -v s=$s -v frip=$frip 'BEGIN{{mn=""}} {{
                w=$3-$2; n++; sum+=w; if(mn==""||w<mn)mn=w; if(w>mx)mx=w
            }} END{{
                printf "%s\t%d\t%.1f\t%d\t%d\t%s\n", s, n, (n>0?sum/n:0), (mn==""?0:mn), (mx==""?0:mx), frip
            }}' "$pk" | tee -a {output.tsv} >> {output.mqc}
        done 2> {log}
        """


# ── Peak count / width summary + FRiP — SEACR (per controlled sample) ─────
rule peak_summary_seacr:
    input:
        peaks = [seacr_peak(s) for s in CONTROLLED_SAMPLES],
        frips = expand(os.path.join(FRIP_DIR, "{s}.seacr.frip.txt"), s=CONTROLLED_SAMPLES)
    output:
        tsv = os.path.join(QC_DIR, "peak_summary_seacr.tsv"),
        mqc = os.path.join(QC_DIR, "peak_summary_seacr_mqc.txt")
    conda:
        "../envs/snakemake.yaml"
    log:
        "logs/peak_summary/seacr.log"
    shell:
        r"""
        mkdir -p {QC_DIR} logs/peak_summary
        printf "# id: peak_summary_seacr\n# section_name: 'SEACR peak summary'\n# description: 'SEACR peak count, width stats and FRiP per controlled sample.'\n# plot_type: 'table'\nSample\tPeaks\tMean width\tMin width\tMax width\tFRiP\n" > {output.mqc}
        echo -e "sample\tn_peaks\tmean_width\tmin_width\tmax_width\tFRiP" > {output.tsv}
        for pk in {input.peaks}; do
            s=$(basename "$pk"); s=${{s%.*.bed}}
            frip=$(cut -f4 {FRIP_DIR}/$s.seacr.frip.txt)
            awk -v s=$s -v frip=$frip 'BEGIN{{mn=""}} {{
                w=$3-$2; n++; sum+=w; if(mn==""||w<mn)mn=w; if(w>mx)mx=w
            }} END{{
                printf "%s\t%d\t%.1f\t%d\t%d\t%s\n", s, n, (n>0?sum/n:0), (mn==""?0:mn), (mx==""?0:mx), frip
            }}' "$pk" | tee -a {output.tsv} >> {output.mqc}
        done 2> {log}
        """


# ── Aggregate QC target (interactive report added in the report task) ─────
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
        [os.path.join(IDR_DIR, f"{g}--{r1}--{r2}--idr_peaks.{group_macs2_ext(g)}.txt")
         for g, r1, r2 in IDR_PAIRS],
        expand(os.path.join(COMPLEXITY_DIR, "{s}_complexity.txt"), s=SAMPLES),
        os.path.join(ANNOT_DIR, "reads_in_annotations_mqc.txt"),
        os.path.join(QC_DIR, "peak_summary_macs2_mqc.txt"),
        os.path.join(QC_DIR, "peak_summary_seacr_mqc.txt"),
        os.path.join(XCOR_DIR, "xcor_summary.tsv"),
        os.path.join(QC_DIR, "fingerprint_jsd.tsv"),
        os.path.join(QC_DIR, "multiqc_fastqc.html"),
        os.path.join(QC_DIR, "cutandrun_qc_report.html"),


# ── Interactive HTML QC report (spike-in removed, SEACR added) ────────────
rule qc_report:
    input:
        # numeric sources
        expand(os.path.join(ALIGN_DIR, "{s}.bowtie2.log"), s=SAMPLES),
        expand(os.path.join(FILTERED_DIR, "{s}.idxstats.txt"), s=SAMPLES),
        expand(os.path.join(DEDUP_DIR, "{s}.dedup.metrics.txt"), s=SAMPLES),
        expand(os.path.join(COMPLEXITY_DIR, "{s}_complexity.txt"), s=SAMPLES),
        os.path.join(QC_DIR, "peak_summary_macs2.tsv"),
        os.path.join(QC_DIR, "peak_summary_seacr.tsv"),
        os.path.join(QC_DIR, "tss_enrichment_scores.tsv"),
        os.path.join(QC_DIR, "blacklist_filtering_stats.txt"),
        os.path.join(ANNOT_DIR, "reads_in_annotations.tsv"),
        os.path.join(XCOR_DIR, "xcor_summary.tsv"),
        os.path.join(QC_DIR, "fingerprint_jsd.tsv"),
        os.path.join(CONSENSUS_DIR, "consensus_peaks.bed"),
        # embedded plots / data
        os.path.join(DEEPTOOLS_DIR, "fragmentSize.png"),
        os.path.join(DEEPTOOLS_DIR, "Heatmap_TSS.png"),
        os.path.join(DEEPTOOLS_DIR, "Profile_TSS.png"),
        os.path.join(DEEPTOOLS_DIR, "ATACseq_fingerprint.png"),
        os.path.join(DEEPTOOLS_DIR, "deeptools_heatmap.png"),
        os.path.join(DEEPTOOLS_DIR, "deeptools_PCA.png"),
        os.path.join(DEEPTOOLS_DIR, "deeptools_scatterplot.png"),
        expand(os.path.join(DEEPTOOLS_DIR, "{s}.gc_content.png"), s=SAMPLES),
        os.path.join(DEEPTOOLS_DIR, "fragment_lengths.txt"),
        os.path.join(DEEPTOOLS_DIR, "correlation_matrix.tab"),
        os.path.join(DEEPTOOLS_DIR, "tss_heatmap_downsampled.json"),
    output:
        html = os.path.join(QC_DIR, "cutandrun_qc_report.html")
    params:
        results = RESULT_DIR,
        samples = ",".join(SAMPLES),
    conda:
        "../envs/snakemake.yaml"
    log:
        "logs/qc_report/qc_report.log"
    shell:
        r"""
        mkdir -p {QC_DIR} logs/qc_report
        python workflow/scripts/build_qc_report.py \
            --results-dir {params.results} \
            --out {output.html} \
            --samples {params.samples} \
            --generated "$(date -u '+%Y-%m-%d %H:%M UTC')" > {log} 2>&1
        """


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
            awk -F'\t' -v s=$s 'BEGIN{{OFS="\t"}}
                {{split($3,a,","); print s, a[1], $9, $10, $11}}' \
                {params.xcordir}/$s.spp.out | tee -a {output.tsv} >> {output.mqc}
        done 2> {log}
        """


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
