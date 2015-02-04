getExprs <- function(file) {
    e <- read.delim(file, as.is = T, fill = T, check.names = F, row.names = 1)

    e
}

getDesign <- function(file) {
    d <- read.delim(file, as.is = T, fill = T, check.names = F, row.names = 1)

    library('plyr')
    REGEX_PATTERN <- 'ER([a-z]+)_ERBB2([a-z]+)'
    er <- gsub(REGEX_PATTERN, '\\1', d$subGroup)
    her2 <- gsub(REGEX_PATTERN, '\\2', d$subGroup)

    mutate(d, 
           er_status = ifelse(er == 'pos', 'ER+', 'ER-'), 
           status = ifelse(her2 == 'pos', 'HER2+', 
                ifelse(er == 'pos' , 'ER+HER2-', 
                    ifelse(er == 'neg' & her2 == 'neg', 'ER-HER2-', NA))), 
           time = t_dmfs, 
           event = e_dmfs)
}

MicroarrayDataset <- function(exprsFile, designFile) {
    exprs <- getExprs(exprsFile)
    design <- getDesign(designFile)

    out <- list(exprs = exprs, design = design)
    class(out) <- "MicroarrayDataset"

    out
}

syncArrayData <- function(arrayData) {
    array.data <- arrayData
    keep <- intersect(colnames(array.data$exprs), rownames(array.data$design))
    array.data$exprs <- array.data$exprs[, keep]
    array.data$design <- array.data$design[keep, ]

    array.data
}

YauGeneExp <- function() {
    exprs.file <- "~/brca_datasets/YauGeneExp-2011-11-11/raw_data/genomicMatrix"
    design.file <- "~/brca_datasets/YauGeneExp-2011-11-11/raw_data/clinical_data"

    MicroarrayDataset(exprs.file, design.file)
}
