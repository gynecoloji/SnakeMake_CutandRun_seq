# Shared setup for cutandrun.smk and qc.smk: config validation, sample sheet
# (via workflow/scripts/samplesheet.py), directory constants, and helpers.
import os
import re
import sys
from snakemake.utils import validate

validate(config, "../schemas/config.schema.yaml")

sys.path.insert(0, os.path.join(workflow.basedir, "scripts"))
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

# IDR groups + samples split by peak mode (all replicates of a condition share one mode)
NARROW_IDR_GROUPS = [g for g in IDR_GROUPS if SS.peak_mode(GROUPS[g][0]) == "narrow"]
BROAD_IDR_GROUPS  = [g for g in IDR_GROUPS if SS.peak_mode(GROUPS[g][0]) == "broad"]
IDR_NARROW_SAMPLES = [s for s in IDR_SAMPLES if s in NARROW_SAMPLES]
IDR_BROAD_SAMPLES  = [s for s in IDR_SAMPLES if s in BROAD_SAMPLES]

# QC IDR runs on every within-condition replicate pair (any condition with >=2 reps)
QC_IDR_SAMPLES = sorted({s for (_g, a, b) in IDR_PAIRS for s in (a, b)})
QC_IDR_NARROW_SAMPLES = [s for s in QC_IDR_SAMPLES if s in NARROW_SAMPLES]
QC_IDR_BROAD_SAMPLES  = [s for s in QC_IDR_SAMPLES if s in BROAD_SAMPLES]

CONTROL_TYPE = config["control_type"]

def resolved_control(sample):
    return SS.resolved_control(sample, CONTROL_TYPE)

# Treatment samples that have an IgG/Input control (for log2 tracks + SEACR)
CONTROLLED_SAMPLES = [s for s in TREATMENT_SAMPLES if resolved_control(s)]
# SEACR is run for controlled treatment samples, split by stringency (static outputs)
SEACR_STRINGENT_SAMPLES = [s for s in CONTROLLED_SAMPLES
                           if SS.seacr_stringency(s, config) == "stringent"]
SEACR_RELAXED_SAMPLES   = [s for s in CONTROLLED_SAMPLES
                           if SS.seacr_stringency(s, config) == "relaxed"]

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
XCOR_DIR       = f"{RESULT_DIR}/xcor"
JSD_DIR        = f"{RESULT_DIR}/qc_fingerprint"

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

def control_bam(sample):
    """Blacklist-filtered BAM of a sample's resolved control ('' if none)."""
    c = resolved_control(sample)
    return f"{BLACKLIST_FILTERED_DIR}/{c}.nobl.bam" if c else ""

def macs2_peak(sample):
    """Per-sample MACS2 peak path with the correct narrow/broad extension."""
    return f"{PEAKS_DIR}/{sample}_peaks.{SS.macs2_ext(sample)}"

def group_macs2_ext(group):
    """narrowPeak/broadPeak extension for a condition (shared by its replicates)."""
    return SS.macs2_ext(GROUPS[group][0])

def seacr_peak(sample):
    """Per-sample SEACR output BED path."""
    return f"{SEACR_DIR}/{sample}.{SS.seacr_stringency(sample, config)}.bed"

def _group_relaxed_inputs(wildcards):
    ext = "broadPeak" if SS.peak_mode(GROUPS[wildcards.group][0]) == "broad" else "narrowPeak"
    return [f"{RELAXED_PEAKS_DIR}/{s}_relaxed.{ext}" for s in GROUPS[wildcards.group]]


# ── ENCODE IDR reproducibility (self-consistency + rescue ratios) ───────
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

# Downstream-analysis directories (workflow/rules/downstream.smk)
ANNOTATION_DIR = f"{RESULT_DIR}/annotation"   # ChIPseeker peak annotation + GO
MOTIF_DIR      = f"{RESULT_DIR}/motifs"        # HOMER motif enrichment
OVERLAP_DIR    = f"{RESULT_DIR}/peak_overlap"  # peak-set Jaccard / overlap

def peak_file(sample):
    """Alias used by the downstream rules/scripts: the MACS2 peak for an IP sample."""
    return macs2_peak(sample)

def all_peak_files():
    return [macs2_peak(s) for s in TREATMENT_SAMPLES]
