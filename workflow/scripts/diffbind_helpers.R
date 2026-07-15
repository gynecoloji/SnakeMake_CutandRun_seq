# Helpers for generic (no spike-in) DESeq2 differential binding (cutandrun_Dx.ipynb).

#' Read a featureCounts table -> coords + integer count matrix (cols = sample names).
read_featurecounts_matrix <- function(path) {
  df <- read.delim(path, comment.char = "#", check.names = FALSE, stringsAsFactors = FALSE)
  meta <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
  count_cols <- setdiff(colnames(df), meta)
  samples <- sub("\\.nobl\\.bam$", "", basename(count_cols))
  counts <- as.matrix(df[, count_cols, drop = FALSE])
  storage.mode(counts) <- "integer"          # featureCounts --countReadPairs -> integers
  colnames(counts) <- samples
  rownames(counts) <- df$Geneid
  list(coords = df[, c("Geneid", "Chr", "Start", "End")], counts = counts, samples = samples)
}

#' Read the CUT&RUN sample sheet -> treatment-sample metadata (peak_mode non-empty).
read_treatment_metadata <- function(samples_csv) {
  m <- read.csv(samples_csv, stringsAsFactors = FALSE, colClasses = "character")
  m$peak_mode <- trimws(ifelse(is.na(m$peak_mode), "", m$peak_mode))
  m <- m[m$peak_mode != "", c("sample_id", "condition", "replicate")]
  rownames(m) <- NULL
  m
}

#' TRUE where a peak overlaps any promoter interval from a BED file.
classify_promoter <- function(coords, promoter_bed_path) {
  pb <- read.delim(promoter_bed_path, header = FALSE, stringsAsFactors = FALSE)
  prom  <- GenomicRanges::GRanges(pb[[1]], IRanges::IRanges(pb[[2]] + 1L, pb[[3]]))  # BED 0-based -> 1-based
  peaks <- GenomicRanges::GRanges(coords$Chr, IRanges::IRanges(coords$Start, coords$End))
  IRanges::overlapsAny(peaks, prom)
}

#' DESeq2 differential binding for one pairwise contrast (cond_b vs cond_a) using
#' DESeq2 median-of-ratios normalization (no spike-in). A paired design
#' `~replicate + condition` is used automatically when the two conditions carry the
#' same set of (unique) replicate labels; otherwise `~condition`.
#' Positive log2FoldChange = higher signal in cond_b.
run_deseq2_contrast <- function(counts, meta, cond_a, cond_b, coords) {
  sub <- meta[meta$condition %in% c(cond_a, cond_b), , drop = FALSE]
  cols <- sub$sample_id
  cts <- counts[, cols, drop = FALSE]
  condition <- factor(sub$condition, levels = c(cond_a, cond_b))
  coldata <- data.frame(condition = condition, row.names = cols)
  reps_a <- sub$replicate[sub$condition == cond_a]
  reps_b <- sub$replicate[sub$condition == cond_b]
  paired <- length(reps_a) > 1 && length(reps_a) == length(reps_b) &&
            setequal(reps_a, reps_b) && !any(duplicated(reps_a)) && !any(duplicated(reps_b))
  if (paired) {
    coldata$replicate <- factor(sub$replicate)
    design <- ~replicate + condition
  } else {
    design <- ~condition
  }
  dds <- DESeq2::DESeqDataSetFromMatrix(cts, coldata, design = design)
  dds <- DESeq2::DESeq(dds)                                   # median-of-ratios normalization
  res <- DESeq2::results(dds, contrast = c("condition", cond_b, cond_a))
  coefname <- grep("^condition_", DESeq2::resultsNames(dds), value = TRUE)
  lfc <- res$log2FoldChange
  if (length(coefname) == 1) {
    lfc <- tryCatch(DESeq2::lfcShrink(dds, coef = coefname, type = "apeglm")$log2FoldChange,
                    error = function(e) res$log2FoldChange)
  }
  data.frame(coords,
             baseMean = res$baseMean, log2FoldChange = lfc,
             pvalue = res$pvalue, padj = res$padj,
             row.names = NULL, check.names = FALSE)
}

#' plotPCA data.frame (VST, blind) over a set of treatment samples.
vst_pca_data <- function(counts, meta) {
  condition <- factor(meta$condition)
  coldata <- data.frame(condition = condition, row.names = meta$sample_id)
  dds <- DESeq2::DESeqDataSetFromMatrix(counts[, meta$sample_id, drop = FALSE], coldata, design = ~condition)
  dds <- DESeq2::estimateSizeFactors(dds)
  vsd <- DESeq2::vst(dds, blind = TRUE)
  DESeq2::plotPCA(vsd, intgroup = "condition", returnData = TRUE)
}
