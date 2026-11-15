suppressPackageStartupMessages({
    library(BiocParallel)
    library(SummarizedExperiment)
    library(satuRn)
    library(stageR)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
    stop('Expected: counts samples analysis outputDirectory workers')
}

count_path <- args[[1]]
sample_path <- args[[2]]
analysis <- args[[3]]
output_dir <- args[[4]]
workers <- as.integer(args[[5]])

if (!analysis %in% c('input', 'monosome', 'polysome', 'redistribution')) {
    stop('Unknown analysis.')
}
if (is.na(workers) || workers < 1) {
    stop('workers must be a positive integer.')
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
counts <- read.delim(count_path, check.names = FALSE)
samples <- read.delim(sample_path, check.names = FALSE)

required_counts <- c('geneId', 'transcriptId')
required_samples <- c('sampleId', 'participantId', 'condition', 'fraction')

if (!all(required_counts %in% colnames(counts))) {
    stop('Count annotations are incomplete.')
}
if (!all(required_samples %in% colnames(samples))) {
    stop('Sample annotations are incomplete.')
}
if (anyDuplicated(counts$transcriptId) || anyDuplicated(samples$sampleId)) {
    stop('Transcript and sample identifiers must be unique.')
}
if (!identical(colnames(counts)[-(1:2)], samples$sampleId)) {
    stop('Count columns and samples are not identically ordered.')
}

counts <- counts[, c('geneId', 'transcriptId', samples$sampleId)]
count_matrix <- as.matrix(counts[, samples$sampleId])
storage.mode(count_matrix) <- 'double'
rownames(count_matrix) <- counts$transcriptId

if (any(!is.finite(count_matrix)) || any(count_matrix < 0)) {
    stop('Counts must be finite and nonnegative.')
}
if (any(table(counts$geneId) < 2)) {
    stop('Every tested gene must contain at least two transcripts.')
}

samples$condition <- factor(
    samples$condition,
    levels = c('Unaffected', 'ASD')
)
if (anyNA(samples$condition)) {
    stop('Unexpected condition label.')
}

if (analysis == 'redistribution') {
    samples$participantId <- factor(samples$participantId)
    samples$fraction <- factor(
        samples$fraction,
        levels = c('monosome', 'polysome')
    )
    samples$pmInteraction <- as.integer(
        samples$condition == 'ASD'
        & samples$fraction == 'polysome'
    )
    formula_full <- ~ participantId + fraction + pmInteraction
    formula_null <- ~ participantId + fraction
    coefficient <- 'pmInteraction'
} else {
    formula_full <- ~ condition
    formula_null <- ~ 1
    coefficient <- 'conditionASD'
}

design_full <- model.matrix(formula_full, data = samples)
design_null <- model.matrix(formula_null, data = samples)

if (qr(design_full)$rank != ncol(design_full)) {
    stop('The full design matrix is not full rank.')
}
if (qr(design_null)$rank != ncol(design_null)) {
    stop('The null design matrix is not full rank.')
}
if (ncol(design_full) - ncol(design_null) != 1) {
    stop('The requested test must have one model degree of freedom.')
}
if (nrow(design_full) - ncol(design_full) < 1) {
    stop('The full model has no residual degrees of freedom.')
}
if (!coefficient %in% colnames(design_full)) {
    stop('The requested coefficient is absent from the full design.')
}

contrast <- matrix(
    0,
    nrow = ncol(design_full),
    ncol = 1,
    dimnames = list(colnames(design_full), analysis)
)
contrast[coefficient, 1] <- 1

row_data <- S4Vectors::DataFrame(
    isoform_id = counts$transcriptId,
    gene_id = counts$geneId,
    row.names = counts$transcriptId
)
col_data <- S4Vectors::DataFrame(
    samples,
    row.names = samples$sampleId
)
object <- SummarizedExperiment(
    assays = list(counts = count_matrix),
    rowData = row_data,
    colData = col_data
)

parallel <- SnowParam(
    workers = workers,
    type = 'SOCK',
    RNGseed = 1077,
    progressbar = FALSE
)
parallel <- bpstart(parallel)
warning_messages <- character()
diagnostic_path <- file.path(
    output_dir,
    paste0(analysis, 'EmpiricalDiagnostics.pdf')
)

gene_totals <- rowsum(
    count_matrix,
    group = counts$geneId,
    reorder = FALSE
)
other_count <- gene_totals[
    match(counts$geneId, rownames(gene_totals)),
    ,
    drop = FALSE
] - count_matrix
fit_audit <- do.call(rbind, bplapply(
    seq_len(nrow(count_matrix)),
    function(index, count_matrix, other_count, design_full) {
        response <- cbind(
            count_matrix[index, ],
            other_count[index, ]
        ) + 1
        model <- suppressWarnings(try(
            glm(
                response ~ -1 + design_full,
                family = 'quasibinomial'
            ),
            silent = TRUE
        ))
        if (inherits(model, 'try-error')) {
            return(c(
                converged = FALSE,
                boundary = NA,
                finiteCoefficients = FALSE
            ))
        }
        c(
            converged = isTRUE(model$converged),
            boundary = (
                isTRUE(model$boundary)
                || any(fitted(model) < 10 * .Machine$double.eps)
                || any(fitted(model) > 1 - 10 * .Machine$double.eps)
            ),
            finiteCoefficients = all(is.finite(coef(model)))
        )
    },
    count_matrix = count_matrix,
    other_count = other_count,
    design_full = design_full,
    BPPARAM = parallel
))
fit_usable <- (
    fit_audit[, 'converged']
    & !fit_audit[, 'boundary']
    & fit_audit[, 'finiteCoefficients']
)
fit_usable[is.na(fit_usable)] <- FALSE

set.seed(1077)
pdf(diagnostic_path, width = 7, height = 5)
object <- tryCatch(
    withCallingHandlers(
        {
            fitted <- fitDTU(
                object,
                formula = formula_full,
                parallel = TRUE,
                BPPARAM = parallel,
                verbose = FALSE
            )
            fitted_models <- rowData(fitted)[['fitDTUModels']]
            fitted_models[!fit_usable] <- lapply(
                fitted_models[!fit_usable],
                function(model) StatModel(
                    type = 'fitError',
                    params = list(
                        coefficients = setNames(
                            rep(NA_real_, ncol(design_full)),
                            colnames(design_full)
                        ),
                        df.residual = NA_real_,
                        dispersion = NA_real_,
                        vcovUnsc = matrix(
                            NA_real_,
                            nrow = ncol(design_full),
                            ncol = ncol(design_full)
                        )
                    ),
                    varPosterior = NA_real_,
                    dfPosterior = NA_real_
                )
            )
            rowData(fitted)[['fitDTUModels']] <- fitted_models
            testDTU(
                fitted,
                contrasts = contrast,
                diagplot1 = TRUE,
                diagplot2 = TRUE,
                sort = FALSE
            )
        },
        warning = function(warning) {
            warning_messages <<- c(
                warning_messages,
                conditionMessage(warning)
            )
            invokeRestart('muffleWarning')
        }
    ),
    finally = {
        dev.off()
        bpstop(parallel)
    }
)

models <- rowData(object)[['fitDTUModels']]
model_type <- vapply(models, slot, character(1), name = 'type')
result <- as.data.frame(
    rowData(object)[[paste0('fitDTUResult_', analysis)]]
)

expected_result <- c(
    'estimates',
    'se',
    'df',
    't',
    'pval',
    'regular_FDR',
    'empirical_pval',
    'empirical_FDR'
)
if (!all(expected_result %in% colnames(result))) {
    stop('Unexpected satuRn result schema.')
}
if (!identical(rownames(result), counts$transcriptId)) {
    stop('satuRn results are not aligned to the input transcripts.')
}

primary_p <- result$empirical_pval
gene_factor <- factor(counts$geneId, levels = unique(counts$geneId))
gene_split <- split(seq_along(gene_factor), gene_factor)
finite_split <- lapply(
    gene_split,
    function(index) index[is.finite(primary_p[index])]
)
minimum_p <- vapply(
    finite_split,
    function(index) {
        if (length(index)) min(primary_p[index]) else NA_real_
    },
    numeric(1)
)
testable_gene <- is.finite(minimum_p)
screen_input <- rep(1, length(minimum_p))
names(screen_input) <- names(minimum_p)
screen_input[testable_gene] <- minimum_p[testable_gene]
theta <- unique(sort(minimum_p[testable_gene]))
features_per_gene <- lengths(finite_split)[testable_gene]
q_by_theta <- vapply(
    theta,
    function(threshold) {
        expected_false <- sum(
            1 - (1 - threshold)^features_per_gene
        )
        discoveries <- sum(minimum_p[testable_gene] <= threshold)
        expected_false / discoveries
    },
    numeric(1)
)
q_screen <- rep(1, length(minimum_p))
names(q_screen) <- names(gene_split)
q_screen[testable_gene] <- pmin(
    1,
    q_by_theta[match(minimum_p[testable_gene], theta)]
)

tx2gene <- data.frame(
    transcript = counts$transcriptId,
    gene = counts$geneId
)
p_confirmation <- matrix(
    primary_p,
    ncol = 1,
    dimnames = list(counts$transcriptId, 'transcript')
)
stage <- stageRTx(
    pScreen = q_screen,
    pConfirmation = p_confirmation,
    pScreenAdjusted = TRUE,
    tx2gene = tx2gene
)
stage <- stageWiseAdjustment(
    stage,
    method = 'dtu',
    alpha = 0.05,
    allowNA = TRUE
)
stage_adjusted <- getAdjustedPValues(
    stage,
    order = FALSE,
    onlySignificantGenes = FALSE
)
stage_adjusted <- data.frame(
    transcriptId = stage_adjusted[, 'txID'],
    stageGeneAdjustedPValue = stage_adjusted[, 'gene'],
    stageTranscriptAdjustedPValue = stage_adjusted[, 'transcript'],
    row.names = NULL,
    check.names = FALSE
)

if (anyDuplicated(stage_adjusted$transcriptId)) {
    stop('stageR returned duplicated transcript identifiers.')
}
if (!all(stage_adjusted$transcriptId %in% counts$transcriptId)) {
    stop('stageR returned unknown transcript identifiers.')
}

transcript_result <- data.frame(
    geneId = counts$geneId,
    transcriptId = counts$transcriptId,
    estimateLogOdds = result$estimates,
    standardError = result$se,
    degreesOfFreedom = result$df,
    statistic = result$t,
    pValue = result$empirical_pval,
    adjustedPValue = result$empirical_FDR,
    regularPValue = result$pval,
    regularAdjustedPValue = result$regular_FDR,
    fitStatus = unname(model_type),
    modelConverged = unname(fit_audit[, 'converged']),
    boundaryFit = unname(fit_audit[, 'boundary']),
    finiteCoefficients = unname(fit_audit[, 'finiteCoefficients']),
    fitUsable = unname(fit_usable),
    stringsAsFactors = FALSE
)
transcript_result <- merge(
    transcript_result,
    stage_adjusted,
    by = 'transcriptId',
    all.x = TRUE,
    sort = FALSE
)
transcript_result <- transcript_result[
    match(counts$transcriptId, transcript_result$transcriptId),
]

if (anyNA(transcript_result$geneId)) {
    stop('stageR results could not be aligned to satuRn transcripts.')
}
expected_stage_gene <- q_screen[
    match(transcript_result$geneId, names(q_screen))
]
if (!isTRUE(all.equal(
    transcript_result$stageGeneAdjustedPValue,
    unname(expected_stage_gene),
    tolerance = 1e-12,
    check.attributes = FALSE
))) {
    stop('stageR gene adjustments disagree with the screening q-values.')
}

fit_error <- model_type != 'glm'
gene_result <- data.frame(
    geneId = names(gene_split),
    minimumTranscriptPValue = minimum_p,
    screenInputPValue = screen_input,
    stageAdjustedPValue = q_screen,
    transcriptCount = lengths(gene_split),
    finiteTranscriptPValues = vapply(
        gene_split,
        function(index) sum(is.finite(primary_p[index])),
        integer(1)
    ),
    fitErrorTranscripts = vapply(
        gene_split,
        function(index) sum(fit_error[index]),
        integer(1)
    ),
    stringsAsFactors = FALSE
)

design_full_output <- data.frame(
    sampleId = samples$sampleId,
    design_full,
    row.names = NULL,
    check.names = FALSE
)
design_null_output <- data.frame(
    sampleId = samples$sampleId,
    design_null,
    row.names = NULL,
    check.names = FALSE
)
summary_output <- data.frame(
    analysis = analysis,
    testedGenes = nrow(gene_result),
    testedTranscripts = nrow(transcript_result),
    modelRank = qr(design_full)$rank,
    residualDegreesOfFreedom = nrow(design_full) - qr(design_full)$rank,
    testDegreesOfFreedom = ncol(design_full) - ncol(design_null),
    fitErrorTranscripts = sum(fit_error),
    nonConvergedTranscripts = sum(!fit_audit[, 'converged']),
    boundaryFitTranscripts = sum(fit_audit[, 'boundary'], na.rm = TRUE),
    nonfiniteCoefficientTranscripts = sum(
        !fit_audit[, 'finiteCoefficients']
    ),
    undefinedGeneScreens = sum(is.na(minimum_p)),
    undefinedTranscriptPValues = sum(is.na(transcript_result$pValue)),
    stageSignificantGenes = sum(
        gene_result$stageAdjustedPValue < 0.05,
        na.rm = TRUE
    ),
    stageSignificantTranscripts = sum(
        transcript_result$stageTranscriptAdjustedPValue < 0.05,
        na.rm = TRUE
    ),
    warningCount = length(warning_messages)
)
if (length(warning_messages)) {
    warning_output <- as.data.frame(table(warning_messages))
    colnames(warning_output) <- c('warningMessage', 'count')
} else {
    warning_output <- data.frame(
        warningMessage = character(),
        count = integer()
    )
}
versions <- data.frame(
    software = c(
        'R',
        'Bioconductor',
        'satuRn',
        'stageR',
        'BiocParallel'
    ),
    version = c(
        paste(R.version$major, R.version$minor, sep = '.'),
        as.character(BiocManager::version()),
        as.character(packageVersion('satuRn')),
        as.character(packageVersion('stageR')),
        as.character(packageVersion('BiocParallel'))
    )
)

write.table(
    gene_result,
    file.path(output_dir, paste0(analysis, 'Gene.tsv')),
    sep = '\t',
    quote = FALSE,
    row.names = FALSE,
    na = 'NA'
)
write.table(
    transcript_result,
    file.path(output_dir, paste0(analysis, 'Transcript.tsv')),
    sep = '\t',
    quote = FALSE,
    row.names = FALSE,
    na = 'NA'
)
write.table(
    counts[, c('geneId', 'transcriptId')],
    file.path(output_dir, paste0(analysis, 'Retained.tsv')),
    sep = '\t',
    quote = FALSE,
    row.names = FALSE
)
write.table(
    design_full_output,
    file.path(output_dir, paste0(analysis, 'DesignFull.tsv')),
    sep = '\t',
    quote = FALSE,
    row.names = FALSE
)
write.table(
    design_null_output,
    file.path(output_dir, paste0(analysis, 'DesignNull.tsv')),
    sep = '\t',
    quote = FALSE,
    row.names = FALSE
)
write.table(
    summary_output,
    file.path(output_dir, paste0(analysis, 'Summary.tsv')),
    sep = '\t',
    quote = FALSE,
    row.names = FALSE
)
write.table(
    warning_output,
    file.path(output_dir, paste0(analysis, 'Warnings.tsv')),
    sep = '\t',
    quote = FALSE,
    row.names = FALSE
)
write.table(
    versions,
    file.path(output_dir, paste0(analysis, 'Versions.tsv')),
    sep = '\t',
    quote = FALSE,
    row.names = FALSE
)
capture.output(
    sessionInfo(),
    file = file.path(output_dir, paste0(analysis, 'SessionInfo.txt'))
)
