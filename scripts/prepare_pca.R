suppressPackageStartupMessages(library(DESeq2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop('Usage: Rscript scripts/prepare_pca.R Input|Monosome|Polysome [output]')
}
fraction <- args[1]
prefixes <- c(Input = 'input', Monosome = 'mono', Polysome = 'poly')
if (!fraction %in% names(prefixes)) stop('Unknown RNA fraction.')
out <- if (length(args) > 1) args[2] else file.path('manuscript/work/pca', fraction)
dir.create(out, recursive = TRUE, showWarnings = FALSE)

meta <- read.delim('data/metadata.tsv', colClasses = 'character', check.names = FALSE)
path <- file.path('data', fraction, 'salmon.merged.gene_counts_length_scaled.tsv')
counts <- read.delim(path, check.names = FALSE)
selected <- paste0(prefixes[fraction], meta$sampleId)
stopifnot(all(selected %in% names(counts)), nrow(meta) == 6)
mat <- round(as.matrix(counts[, selected, drop = FALSE]))
rownames(mat) <- counts$gene_id
colnames(mat) <- meta$sampleId
mat <- mat[rowSums(mat) > 0, , drop = FALSE]
storage.mode(mat) <- 'integer'
rownames(meta) <- meta$sampleId

dds <- DESeqDataSetFromMatrix(countData = mat, colData = meta, design = ~1)
vsd <- varianceStabilizingTransformation(dds, blind = TRUE, fitType = 'parametric')
vst <- assay(vsd)
stopifnot(all(is.finite(vst)), ncol(vst) == 6)
write.csv(
    data.frame(geneId = rownames(vst), vst, check.names = FALSE),
    file.path(out, paste0(fraction, '_blind_VST_counts.csv')),
    row.names = FALSE, quote = FALSE
)
capture.output(sessionInfo(), file = file.path(out, 'R_session_info.txt'))
cat('Prepared', fraction, 'blind VST:', nrow(vst), 'genes,', ncol(vst), 'libraries\n')
