getExprs <- function(file) {
    e <- read.csv(file, as.is = T, fill = T, check.names = F, row.names = 1)

    e
}

getAnnot <- function(file) {
    a <- read.csv(file, as.is = T, fill = T, check.names = F, row.names = 1)

    a
}

getDesign <- function(file) {
    d <- read.csv(file, as.is = T, fill = T, check.names = F, row.names = 1)

    library('plyr')
    mutate(
        d, 
        hr = ifelse(er_status_ihc == 'P' | pr_status_ihc == 'P', 'HR+', 'HR-'), 
        status = ifelse(
            her2_status == 'P', 'HER2+',
            ifelse(
                er_status_ihc == 'P' | pr_status_ihc == 'P', 'HR+HER2-', 
                ifelse(er_status_ihc == 'N' & 
                       pr_status_ihc == 'N' & 
                       her2_status == 'N', 'TN', NA))), 
        time = drfs_t, 
        event = drfs_e)
}

MicroarrayDataset <- function(exprsFile, annotFile, designFile) {
    exprs <- getExprs(exprsFile)
    annot <- getAnnot(annotFile)
    design <- getDesign(designFile)

    out <- list(exprs = exprs, annot = annot, design = design)
    class(out) <- "MicroarrayDataset"

    out
}

collapseProbes <- function(arrayData, method = 'MaxMean') {
    array.data <- arrayData
    exprs <- array.data$exprs
    annot <- array.data$annot

    temp <- merge(exprs, annot['ENTREZ_GENE_ID'], by = 0)
    rownames(temp) <- row.names(temp)
    valid <- temp$ENTREZ_GENE_ID != '' & !is.na(temp$ENTREZ_GENE_ID)
    temp <- temp[valid, -1]

    library(WGCNA)
    rowGroup <- temp$ENTREZ_GENE_ID
    rowID <- rownames(temp)
    temp <- temp[, -ncol(temp)]
    temp <- collapseRows(temp, rowGroup = rowGroup, rowID = rowID, 
                         method = method)

    array.data$exprs <- as.data.frame(temp$datETcollapsed)
    array.data
}

syncArrayData <- function(arrayData) {
    array.data <- arrayData
    keep <- intersect(colnames(array.data$exprs), rownames(array.data$design))
    array.data$exprs <- array.data$exprs[, keep]
    array.data$design <- array.data$design[keep, ]

    array.data
}

GSE25066 <- function() {
    exprs.file <- "~/brca_datasets/gse25066_091114/raw_data/exprs.csv"
    annot.file <- "~/brca_datasets/gse25066_091114/raw_data/annot.csv"
    design.file <- "~/brca_datasets/gse25066_091114/raw_data/design.csv"

    data <- MicroarrayDataset(exprs.file, annot.file, design.file)
    data <- collapseProbes(data, method = 'MaxMean')

    data
}
