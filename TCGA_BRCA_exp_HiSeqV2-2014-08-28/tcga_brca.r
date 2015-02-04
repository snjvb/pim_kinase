getExprs <- function(file) {
    e <- read.delim(file, as.is = T, fill = T, check.names = F, row.names = 1)

    e
}

getDesign <- function(file) {
    d <- read.delim(file, as.is = T, fill = T, check.names = F, row.names = 1)

    getDesiredSamples(d)
}

getDesiredSamples <- function(data) {
    .data <- data

    library(plyr)
    .data <- mutate(
        .data, 
        hr = ifelse(ER_Status_nature2012 == 'Positive' | 
                    PR_Status_nature2012 == 'Positive', 
                    'HR+', 'HR-'), 
        status = ifelse(
            HER2_Final_Status_nature2012 == 'Positive', 'HER2+',
            ifelse(
                ER_Status_nature2012 == 'Positive' | 
                PR_Status_nature2012 == 'Positive', 
                'HR+HER2-', 
                ifelse(ER_Status_nature2012 == 'Negative' & 
                       PR_Status_nature2012 == 'Negative' & 
                       HER2_Final_Status_nature2012 == 'Negative', 
                       'TN', NA))))

    subset(.data, 
           grepl('[Pp]rimary', .data$sample_type) & !is.na(.data$status))
}

RNASeqDataset <- function(exprsFile, designFile) {
    exprs <- getExprs(exprsFile)
    design <- getDesign(designFile)

    out <- list(exprs = exprs, design = design)
    class(out) <- "RNASeqDataset"

    out
}

syncData <- function(data) {
    .data <- data
    keep <- intersect(colnames(.data$exprs), rownames(.data$design))
    .data$exprs <- .data$exprs[, keep]
    .data$design <- .data$design[keep, ]

    .data
}

TCGA.BRCA <- function() {
    exprs.file <- "~/brca_datasets/TCGA_BRCA_exp_HiSeqV2-2014-08-28/genomicMatrix"
    design.file <- "~/brca_datasets/TCGA_BRCA_exp_HiSeqV2-2014-08-28/clinical_data"

    data <- RNASeqDataset(exprs.file, design.file)

    syncData(data)
}
