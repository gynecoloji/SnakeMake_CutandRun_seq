#!/usr/bin/env python3
"""Generate cutandrun_Dx.ipynb — a generic (no spike-in) DESeq2 differential-binding
notebook for the CUT&RUN pipeline.

Reads condition/replicate from config/samples.csv and a consensus count matrix,
runs DESeq2 (median-of-ratios) for every pairwise treatment-condition contrast,
split into promoter vs distal peaks, with MA/volcano/PCA and ChIPseeker annotation.

Usage:
    python workflow/scripts/build_diffbind_notebook.py [macs2|seacr|both]
(default: macs2). Runs on the `cutandrun_Dx` conda env (R / Bioconductor, ir kernel).
"""
import sys
import nbformat as nbf

MATRIX_PATHS = {
    "macs2": ("macs2", "results/consensus/consensus_counts.txt"),
    "seacr": ("seacr", "results/consensus_seacr/consensus_counts.txt"),
}

sel = sys.argv[1] if len(sys.argv) > 1 else "macs2"
if sel not in ("macs2", "seacr", "both"):
    sys.exit("argument must be one of: macs2, seacr, both")
matrices = [MATRIX_PATHS["macs2"], MATRIX_PATHS["seacr"]] if sel == "both" else [MATRIX_PATHS[sel]]
matrices_r = "list(\n" + ",\n".join(
    f"  list(label={lab!r}, path={path!r})".replace("'", '"') for lab, path in matrices
) + "\n)"

nb = nbf.v4.new_notebook()
md = nbf.v4.new_markdown_cell
code = nbf.v4.new_code_cell
cells = []

cells.append(md(
    "# CUT&RUN differential binding (DESeq2, no spike-in)\n\n"
    "DESeq2 (median-of-ratios normalization) on the workflow's consensus count matrix, run for "
    "**every pairwise treatment-condition contrast**, split into **promoter** and **distal** peaks.\n\n"
    "- Sample metadata (condition, replicate): `config/samples.csv` (treatment rows only)\n"
    "- Count matrix: MACS2 (`results/consensus/consensus_counts.txt`) and/or SEACR "
    "(`results/consensus_seacr/consensus_counts.txt`) — selected when the notebook was generated\n"
    "- Promoter set: `ref/promoter_chr1-22X.bed`; distal = non-overlapping\n"
    "- **Paired design** `~replicate + condition` is used automatically when both conditions share the "
    "same replicate labels; otherwise `~condition`\n"
    "- Contrast `B_vs_A`: **positive log2FC = higher signal in B**\n"
    "- Sample **PCA** on the VST count matrix; MA/volcano per contrast; ChIPseeker nearest-gene annotation\n"
    "- Runs on the `cutandrun_Dx` env (R kernel)."))

cells.append(code(
    "suppressMessages({library(DESeq2); library(GenomicRanges); library(ggplot2)})\n"
    "source('workflow/scripts/diffbind_helpers.R')\n"
    "outdir <- 'results/diff_region'; dir.create(outdir, recursive=TRUE, showWarnings=FALSE)\n"
    "promoter_bed <- 'ref/promoter_chr1-22X.bed'\n"
    "THRESH_PADJ <- 0.05; THRESH_LFC <- 1\n"
    "MATRICES <- " + matrices_r + "\n"
    "annotate_sig <- function(sig, out_tsv){\n"
    "  if(nrow(sig)==0){cat('  (no DB regions to annotate)\\n'); return(invisible())}\n"
    "  suppressMessages({library(ChIPseeker); library(TxDb.Hsapiens.UCSC.hg38.knownGene); library(org.Hs.eg.db)})\n"
    "  gr <- GRanges(sig$Chr, IRanges(sig$Start, sig$End))\n"
    "  an <- as.data.frame(annotatePeak(gr, TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene,\n"
    "                                   annoDb='org.Hs.eg.db', verbose=FALSE))\n"
    "  out <- cbind(sig, annotation=an$annotation, SYMBOL=an$SYMBOL, distanceToTSS=an$distanceToTSS)\n"
    "  write.table(out, out_tsv, sep='\\t', quote=FALSE, row.names=FALSE)\n"
    "  cat('  annotated', nrow(out), 'regions ->', out_tsv, '\\n')}"))

cells.append(code(
    "# Treatment-sample metadata + the pairwise contrasts that will be tested\n"
    "meta <- read_treatment_metadata('config/samples.csv')\n"
    "print(meta)\n"
    "conds <- unique(meta$condition)\n"
    "cat('conditions:', paste(conds, collapse=', '), '\\n')\n"
    "if(length(conds) >= 2) print(t(combn(conds, 2))) else cat('only one condition; nothing to contrast\\n')"))

cells.append(code(
    "# Analyze one count matrix: PCA + every pairwise contrast (promoter & distal)\n"
    "analyze_matrix <- function(mlabel, mpath){\n"
    "  cat('\\n==== matrix:', mlabel, '(', mpath, ') ====\\n')\n"
    "  fc <- read_featurecounts_matrix(mpath); counts <- fc$counts; coords <- fc$coords\n"
    "  meta_m <- meta[meta$sample_id %in% colnames(counts), , drop=FALSE]\n"
    "  counts <- counts[, meta_m$sample_id, drop=FALSE]\n"
    "  od <- file.path(outdir, mlabel); dir.create(od, recursive=TRUE, showWarnings=FALSE)\n"
    "  is_prom <- classify_promoter(coords, promoter_bed)\n"
    "  cat('peaks:', nrow(counts), ' promoter:', sum(is_prom), ' distal:', sum(!is_prom), '\\n')\n"
    "  if(ncol(counts) >= 2){\n"
    "    pca <- vst_pca_data(counts, meta_m); pv <- round(100*attr(pca,'percentVar'))\n"
    "    p <- ggplot(pca, aes(PC1, PC2, color=condition, label=name)) + geom_point(size=3) +\n"
    "      geom_text(vjust=-0.8, size=3, show.legend=FALSE) +\n"
    "      xlab(paste0('PC1: ',pv[1],'%')) + ylab(paste0('PC2: ',pv[2],'%')) +\n"
    "      ggtitle(paste('Sample PCA -', mlabel)) + theme_bw()\n"
    "    ggsave(file.path(od,'PCA_samples.png'), p, width=6, height=5, dpi=120); print(p)}\n"
    "  conds <- unique(meta_m$condition)\n"
    "  if(length(conds) < 2){cat('only one condition in this matrix; no contrasts\\n'); return(invisible())}\n"
    "  for(pr in combn(conds, 2, simplify=FALSE)){\n"
    "    a <- pr[1]; b <- pr[2]; tag <- paste0(b,'_vs_',a); cat('\\n-- contrast:', tag, '--\\n')\n"
    "    for(part in list(list('promoter', is_prom), list('distal', !is_prom))){\n"
    "      pname <- part[[1]]; sel <- part[[2]]\n"
    "      res <- run_deseq2_contrast(counts[sel,,drop=FALSE], meta_m, a, b, coords[sel,])\n"
    "      write.table(res, file.path(od, paste0(pname,'_',tag,'_deseq2.tsv')), sep='\\t', quote=FALSE, row.names=FALSE)\n"
    "      sig <- subset(res, !is.na(padj) & padj<THRESH_PADJ & abs(log2FoldChange)>THRESH_LFC)\n"
    "      write.table(sig, file.path(od, paste0(pname,'_',tag,'_sig.tsv')), sep='\\t', quote=FALSE, row.names=FALSE)\n"
    "      cat(sprintf('  %-9s tested=%d DB=%d up(%s)=%d down=%d\\n', pname, sum(!is.na(res$padj)),\n"
    "          nrow(sig), b, sum(sig$log2FoldChange>0), sum(sig$log2FoldChange<0)))\n"
    "      ma <- ggplot(res, aes(baseMean, log2FoldChange)) +\n"
    "        geom_point(aes(color=!is.na(padj)&padj<THRESH_PADJ), size=.4) + scale_x_log10() +\n"
    "        scale_color_manual(values=c('grey70','red'), guide='none') + labs(title=paste('MA',pname,tag)) + theme_bw()\n"
    "      vol <- ggplot(res, aes(log2FoldChange, -log10(padj))) +\n"
    "        geom_point(aes(color=!is.na(padj)&padj<THRESH_PADJ&abs(log2FoldChange)>THRESH_LFC), size=.4) +\n"
    "        scale_color_manual(values=c('grey70','red'), guide='none') +\n"
    "        geom_vline(xintercept=c(-THRESH_LFC,THRESH_LFC), lty=2) + labs(title=paste('Volcano',pname,tag)) + theme_bw()\n"
    "      ggsave(file.path(od, paste0('MA_',pname,'_',tag,'.png')), ma, width=5, height=4, dpi=120)\n"
    "      ggsave(file.path(od, paste0('volcano_',pname,'_',tag,'.png')), vol, width=5, height=4, dpi=120)\n"
    "      annotate_sig(sig, file.path(od, paste0(pname,'_',tag,'_sig_annotated.tsv')))\n"
    "    }}}\n"
    "for(m in MATRICES) analyze_matrix(m$label, m$path)\n"
    "cat('\\nDone. Results under', outdir, '/<matrix>/\\n')"))

nb.cells = cells
nb.metadata.kernelspec = {"name": "ir", "display_name": "R", "language": "R"}
nb.metadata.language_info = {"name": "R"}
nbf.write(nb, "cutandrun_Dx.ipynb")
print(f"wrote cutandrun_Dx.ipynb with {len(cells)} cells (ir kernel); matrices={sel}")
