# Primary CUT&RUN pipeline: FastQC/fastp → human-only Bowtie2 alignment →
# unique/proper-pair filter (+ mito-% QC) → Picard dedup → blacklist filter →
# RPGC + IgG-subtracted tracks → MACS2 & SEACR peaks → MACS2 fixed-width consensus
# + SEACR variable-width consensus → featureCounts matrices.
#
# Shared config, samples, directory constants and helpers live in common.smk.

rule cutandrun_all:
    input:
        # FastQC + fastp
        expand(f"{FASTQC_DIR}/{{s}}_R1_001_fastqc.html", s=SAMPLES),
        expand(f"{FASTQC_DIR}/{{s}}_R2_001_fastqc.html", s=SAMPLES),
        expand(f"{FASTP_DIR}/{{s}}.html", s=SAMPLES),
        # Filtered / dedup / blacklist-filtered BAMs
        expand(f"{FILTERED_DIR}/{{s}}.sorted.filtered.bam", s=SAMPLES),
        expand(f"{DEDUP_DIR}/{{s}}.dedup.bam", s=SAMPLES),
        expand(f"{BLACKLIST_FILTERED_DIR}/{{s}}.nobl.bam", s=SAMPLES),
        expand(f"{BLACKLIST_FILTERED_DIR}/{{s}}.nobl.bam.bai", s=SAMPLES),
        # Signal tracks
        expand(f"{BIGWIG_DIR}/{{s}}.bw", s=SAMPLES),
        expand(f"{LOG2_BIGWIG_DIR}/{{s}}.log2ratio.bw", s=CONTROLLED_SAMPLES),
        # MACS2 peaks (per-sample narrow/broad extension)
        [macs2_peak(s) for s in TREATMENT_SAMPLES],


# ── Genome chrom sizes (for SEACR bedgraph / genomecov) ──────────────────
rule genome_chrom_sizes:
    input:
        fasta = GENOME_FASTA
    output:
        sizes = CHROM_SIZES
    conda:
        "../envs/snakemake.yaml"
    log:
        "logs/genome_chrom_sizes/sizes.log"
    shell:
        r"""
        mkdir -p $(dirname {output.sizes}) logs/genome_chrom_sizes
        samtools faidx {input.fasta} 2> {log}
        cut -f1,2 {input.fasta}.fai > {output.sizes} 2>> {log}
        """


# ── Build the human-only Bowtie2 index (optionally chromosome-subset) ─────
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
    conda:
        "../envs/snakemake.yaml"
    log:
        "logs/build_genome_index/build.log"
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


# ── Bowtie2 alignment to the human index (CUT&RUN parameters) ─────────────
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
    conda:
        "../envs/snakemake.yaml"
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


# ── Filter: properly-paired + unique reads, drop orphans, mito-% QC, keep_chroms ─
rule samtools_sort_filter_index:
    input:
        f"{ALIGN_DIR}/{{sample}}.bam"
    output:
        bam = f"{FILTERED_DIR}/{{sample}}.sorted.filtered.bam",
        bai = f"{FILTERED_DIR}/{{sample}}.sorted.filtered.bam.bai",
        flagstat = f"{FILTERED_DIR}/{{sample}}_summary.txt",
        idxstats = f"{FILTERED_DIR}/{{sample}}.idxstats.txt"
    params:
        keep_chroms = config["keep_chroms"],
        prekeep     = f"{TMP_DIR}/{{sample}}.prekeep.bam"
    log:
        "logs/samtools/{sample}.log"
    threads: 20
    conda:
        "../envs/snakemake.yaml"
    shell:
        r"""
        mkdir -p {FILTERED_DIR} logs/samtools {TMP_DIR}
        # Raw flagstat of the aligned BAM (pre-filter QC)
        samtools flagstat {input} > {FILTERED_DIR}/{wildcards.sample}_raw_summary.txt 2>> {log}

        # Keep properly-paired, primary, mapped, UNIQUE reads (drop multi-mappers)
        samtools view -@ {threads} -hS -f 2 -F 2316 {input} | grep -v "XS:i:" \
            > {TMP_DIR}/temp_{wildcards.sample}.unsorted.sam 2>> {log}

        # Name-sort, then drop reads orphaned by filtering (keep complete pairs only)
        samtools sort -n -O sam -o {TMP_DIR}/temp_{wildcards.sample}.sorted.sam \
            {TMP_DIR}/temp_{wildcards.sample}.unsorted.sam 2>> {log}
        python workflow/scripts/process_sam.py {TMP_DIR}/temp_{wildcards.sample}.sorted.sam \
            {TMP_DIR}/temp_{wildcards.sample}.uc.unsorted.sam 2>> {log}

        # Coordinate-sort (incl chrM) to a BAM and index it
        samtools view -@ {threads} -bhS {TMP_DIR}/temp_{wildcards.sample}.uc.unsorted.sam | \
        samtools sort -@ {threads} -O bam -o {params.prekeep} 2>> {log}
        samtools index -@ {threads} {params.prekeep} 2>> {log}

        # Per-chromosome counts (mito-% QC) BEFORE dropping non-analysis chroms
        samtools idxstats {params.prekeep} > {output.idxstats} 2>> {log}

        # Keep only the analysis chromosomes (drop chrM/chrY/non-primary), index, flagstat
        samtools view -@ {threads} -b -o {output.bam} {params.prekeep} {params.keep_chroms} 2>> {log}
        samtools index -@ {threads} {output.bam} {output.bai} 2>> {log}
        samtools flagstat {output.bam} > {output.flagstat} 2>> {log}

        rm -f {TMP_DIR}/temp_{wildcards.sample}.unsorted.sam \
              {TMP_DIR}/temp_{wildcards.sample}.uc.unsorted.sam \
              {TMP_DIR}/temp_{wildcards.sample}.sorted.sam \
              {params.prekeep} {params.prekeep}.bai
        """


# ── Remove (or mark) duplicates with Picard ──────────────────────────────
rule remove_duplicates:
    input:
        filtered_bam = f"{FILTERED_DIR}/{{sample}}.sorted.filtered.bam"
    output:
        dedup_bam = f"{DEDUP_DIR}/{{sample}}.dedup.bam",
        metrics = f"{DEDUP_DIR}/{{sample}}.dedup.metrics.txt"
    params:
        remove_dups = "true" if config["remove_duplicates"] else "false"
    threads: 4
    conda:
        "../envs/snakemake.yaml"
    log:
        "logs/dedup/{sample}.log"
    shell:
        r"""
        mkdir -p {DEDUP_DIR}
        java -jar ref/picard.jar MarkDuplicates \
               INPUT={input.filtered_bam} \
               OUTPUT={output.dedup_bam} \
               METRICS_FILE={output.metrics} \
               REMOVE_DUPLICATES={params.remove_dups} \
               ASSUME_SORTED=true \
               VALIDATION_STRINGENCY=LENIENT \
               TMP_DIR=tmp 2> {log}
        samtools index {output.dedup_bam}
        """


# ── FastQC + fastp (verbatim from the ATAC workflow) ─────────────────────
# FastQC on raw reads
rule fastqc:
    input:
        r1 = "data/{sample}_R1_001.fastq.gz",
        r2 = "data/{sample}_R2_001.fastq.gz"
    output:
        r1_html = f"{FASTQC_DIR}/{{sample}}_R1_001_fastqc.html",
        r2_html = f"{FASTQC_DIR}/{{sample}}_R2_001_fastqc.html"
    params:
        outdir = FASTQC_DIR
    threads: 2
    conda:
        "../envs/snakemake.yaml"
    log:
        "logs/fastqc/{sample}.log"
    shell:
        """
        mkdir -p {params.outdir}
        fastqc -t {threads} -o {params.outdir} {input.r1} {input.r2} > {log} 2>&1
        """

# Fastp for read trimming and quality filtering
rule fastp:
    input:
        r1 = "data/{sample}_R1_001.fastq.gz",
        r2 = "data/{sample}_R2_001.fastq.gz"
    output:
        r1 = f"{FASTP_DIR}/{{sample}}_R1.trimmed.fastq.gz",
        r2 = f"{FASTP_DIR}/{{sample}}_R2.trimmed.fastq.gz",
        html = f"{FASTP_DIR}/{{sample}}.html",
        json = f"{FASTP_DIR}/{{sample}}.json"
    threads: 16
    conda:
        "../envs/snakemake.yaml"
    params:
        adapter_args = FASTP_ADAPTER_ARGS
    log:
        "logs/fastp/{sample}.log"
    shell:
        """
        mkdir -p {FASTP_DIR}
        fastp --in1 {input.r1} --in2 {input.r2} \
              --out1 {output.r1} --out2 {output.r2} \
              --thread {threads} \
              --html {output.html} \
              --json {output.json} \
              {params.adapter_args} \
              --trim_poly_g \
              --length_required 30 -p --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 > {log} 2>&1
        """


# ── Fragment-level ENCODE blacklist filtering (verbatim from ATAC) ───────
# Filter against blacklist regions
rule filter_blacklist:
    priority: 10
    input:
        bam = f"{DEDUP_DIR}/{{sample}}.dedup.bam",
        blacklist = config["blacklist"]
    output:
        filtered_bam = f"{BLACKLIST_FILTERED_DIR}/{{sample}}.nobl.bam",
        filtered_bai = f"{BLACKLIST_FILTERED_DIR}/{{sample}}.nobl.bam.bai",
        excluded_reads = f"{BLACKLIST_FILTERED_DIR}/{{sample}}.blacklisted.bam"
    params:
        temp_bedpe = f"{TMP_DIR}/{{sample}}.fragments.bedpe",
        temp_fragment_bed = f"{TMP_DIR}/{{sample}}.fragments.bed",
        temp_blacklist_fragments = f"{TMP_DIR}/{{sample}}.blacklisted.fragments.bed",
        temp_blacklist_ids = f"{TMP_DIR}/{{sample}}.blacklisted.ids.txt",
        temp_namesorted_bam = f"{TMP_DIR}/{{sample}}.namesorted.bam",
        temp_filtered_bam = f"{TMP_DIR}/{{sample}}.filtered.bam",
        temp_excluded_bam = f"{TMP_DIR}/{{sample}}.excluded.bam"
    threads: 8
    conda:
        "../envs/snakemake.yaml"
    log:
        "logs/blacklist_filter/{sample}.log"
    shell:
        """
        mkdir -p {BLACKLIST_FILTERED_DIR} {TMP_DIR}
        
        # Sort BAM by read name once and reuse for both BEDPE conversion and filtering
        samtools sort -n -@ {threads} -o {params.temp_namesorted_bam} {input.bam} 2> {log}

        # Convert name-sorted BAM to BEDPE format
        bedtools bamtobed -bedpe -i {params.temp_namesorted_bam} > {params.temp_bedpe} 2>> {log}

        # Convert BEDPE to fragment BED (one entry per fragment)
        # Extract the fragment coordinates (minimum start, maximum end)
        # and keep the read name for later filtering
        awk 'BEGIN {{OFS="\\t"}} {{if ($1==$4) print $1, ($2<$5?$2:$5), ($3>$6?$3:$6), $7, ".", ($9=="+"?"+":"-")}}' \
        {params.temp_bedpe} > {params.temp_fragment_bed} 2>> {log}

        # Find fragments that intersect with blacklisted regions
        bedtools intersect -a {params.temp_fragment_bed} -b {input.blacklist} -wa > {params.temp_blacklist_fragments} 2>> {log}

        # Extract read IDs from blacklisted fragments
        cut -f4 {params.temp_blacklist_fragments} | sort | uniq > {params.temp_blacklist_ids} 2>> {log}
        
        # Create properly paired BAMs - one with fragments that don't overlap blacklist
        
        # Filter out fragments overlapping blacklisted regions
        samtools view -@ {threads} -b -N ^{params.temp_blacklist_ids} \
            {params.temp_namesorted_bam} > {params.temp_filtered_bam} 2>> {log}
            
        # Extract fragments overlapping blacklisted regions
        samtools view -@ {threads} -b -N {params.temp_blacklist_ids} \
            {params.temp_namesorted_bam} > {params.temp_excluded_bam} 2>> {log}
            
        # Sort filtered BAM (non-blacklisted fragments) by coordinate for final output
        samtools sort -@ {threads} -o {output.filtered_bam} {params.temp_filtered_bam} 2>> {log}
        
        # Sort excluded BAM (blacklisted fragments) by coordinate for QC
        samtools sort -@ {threads} -o {output.excluded_reads} {params.temp_excluded_bam} 2>> {log}
        
        # Index the filtered BAM
        samtools index -@ {threads} {output.filtered_bam} {output.filtered_bai} 2>> {log}
        
        # Report stats (before cleanup so temp files are still available)
        echo "Blacklist filtering completed for {wildcards.sample}" >> {log}
        echo "$(wc -l < {params.temp_blacklist_fragments}) fragments overlap blacklisted regions" >> {log}
        echo "$(wc -l < {params.temp_blacklist_ids}) unique fragment IDs overlapping blacklisted regions" >> {log}
        echo "$(samtools view -c {output.excluded_reads}) total reads in excluded fragments" >> {log}
        echo "$(samtools view -c {output.filtered_bam}) total reads in filtered output" >> {log}

        # Clean up temporary files
        rm -f {params.temp_bedpe} {params.temp_fragment_bed} {params.temp_blacklist_fragments} \
            {params.temp_blacklist_ids} {params.temp_namesorted_bam} {params.temp_filtered_bam} \
            {params.temp_excluded_bam}
        """


# ── Depth-normalized (RPGC) bigWig per sample (verbatim from ATAC) ───────
# ── Module A: depth-normalized bigWig (before/after comparison) ─────────
rule create_bigwig:
    input:
        bam = f"{BLACKLIST_FILTERED_DIR}/{{sample}}.nobl.bam",
        bai = f"{BLACKLIST_FILTERED_DIR}/{{sample}}.nobl.bam.bai"
    output:
        bw = f"{BIGWIG_DIR}/{{sample}}.bw"
    params:
        egs       = config["effective_genome_size"],
        bin_size  = config["bin_size"],
        blacklist = config["blacklist"]
    threads: 8
    conda:
        "../envs/deeptools.yaml"
    log:
        "logs/bigwig/{sample}.log"
    shell:
        """
        mkdir -p {BIGWIG_DIR} logs/bigwig
        bamCoverage --bam {input.bam} \
            --normalizeUsing RPGC \
            --effectiveGenomeSize {params.egs} \
            --binSize {params.bin_size} \
            --numberOfProcessors {threads} \
            --extendReads \
            --blackListFileName {params.blacklist} \
            --outFileName {output.bw} > {log} 2>&1
        """


# ── IgG-subtracted log2(treatment/IgG) bigWig (deepTools bamCompare) ──────
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
    conda:
        "../envs/deeptools.yaml"
    log:
        "logs/log2ratio_bigwig/{sample}.log"
    shell:
        r"""
        mkdir -p {LOG2_BIGWIG_DIR} logs/log2ratio_bigwig
        bamCompare -b1 {input.treat} -b2 {input.ctrl} \
            --operation log2 --normalizeUsing CPM \
            --binSize {params.bin_size} --numberOfProcessors {threads} \
            --extendReads --blackListFileName {params.blacklist} \
            --outFileName {output.bw} > {log} 2>&1
        """


# ── MACS2 peak calling — narrow (per-sample IgG control if provided) ──────
rule call_peaks_macs2_narrow:
    wildcard_constraints:
        sample = _alt(NARROW_SAMPLES)
    input:
        treatment = f"{BLACKLIST_FILTERED_DIR}/{{sample}}.nobl.bam",
        control = lambda w: [igg_bam(w.sample)] if igg_bam(w.sample) else []
    output:
        peaks = f"{PEAKS_DIR}/{{sample}}_peaks.narrowPeak"
    params:
        outdir = PEAKS_DIR,
        name = "{sample}",
        genome = config["macs2_genome"],
        q = config["macs2_qvalue"],
        control_arg = lambda w: f"-c {igg_bam(w.sample)}" if igg_bam(w.sample) else ""
    conda:
        "../envs/macs2.yaml"
    log:
        "logs/macs2/{sample}.log"
    shell:
        r"""
        mkdir -p {params.outdir} logs/macs2
        macs2 callpeak -t {input.treatment} {params.control_arg} \
              -f BAMPE -g {params.genome} --outdir {params.outdir} \
              -n {params.name} --nomodel -q {params.q} > {log} 2>&1
        """


# ── MACS2 peak calling — broad (per-sample IgG control if provided) ───────
rule call_peaks_macs2_broad:
    wildcard_constraints:
        sample = _alt(BROAD_SAMPLES)
    input:
        treatment = f"{BLACKLIST_FILTERED_DIR}/{{sample}}.nobl.bam",
        control = lambda w: [igg_bam(w.sample)] if igg_bam(w.sample) else []
    output:
        peaks = f"{PEAKS_DIR}/{{sample}}_peaks.broadPeak"
    params:
        outdir = PEAKS_DIR,
        name = "{sample}",
        genome = config["macs2_genome"],
        q = config["macs2_qvalue"],
        broad_cutoff = config["macs2_broad_cutoff"],
        control_arg = lambda w: f"-c {igg_bam(w.sample)}" if igg_bam(w.sample) else ""
    conda:
        "../envs/macs2.yaml"
    log:
        "logs/macs2/{sample}.log"
    shell:
        r"""
        mkdir -p {params.outdir} logs/macs2
        macs2 callpeak -t {input.treatment} {params.control_arg} \
              -f BAMPE -g {params.genome} --outdir {params.outdir} \
              -n {params.name} --nomodel --broad --broad-cutoff {params.broad_cutoff} \
              -q {params.q} > {log} 2>&1
        """
