getExprs <- function(file) {
    e <- read.csv(file, as.is = T, fill = T, check.names = F, row.names = 1)

    e
}

getAnnot <- function(file) {
    a <- read.delim(file, as.is = T, fill = T, check.names = F, row.names = 1)

    a
}

getDesign <- function(file) {
    d <- read.csv(file, as.is = T, fill = T, check.names = F, row.names = 1)

    ## Convert days to years
    d$rfs.t <- d$rfs.t / 365

    library('plyr')
    mutate(
        d,
        hr = ifelse(er == 1 | pr == 1, 'HR+', 'HR-'),
        status = ifelse(
            her2 == 1, 'HER2+',
            ifelse(
                er == 1 | pr == 1, 'HR+HER2-',
                ifelse(er == 0 &
                       pr == 0 &
                       her2 == 0, 'TN', NA))),
        time = rfs.t,
        event = rfs.e)
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

    temp <- merge(exprs, annot['EntrezGeneID'], by = 0)
    rownames(temp) <- row.names(temp)
    valid <- temp$EntrezGeneID != '' & !is.na(temp$EntrezGeneID)
    temp <- temp[valid, -1]

    library(WGCNA)
    rowGroup <- temp$EntrezGeneID
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

ISPY1 <- function() {
    exprs.file <- "~/Workspace/pim_kinase/ispy1_082814/raw_data/exprs.csv"
    annot.file <- "~/Workspace/pim_kinase/ispy1_082814/raw_data/annot.txt"
    design.file <- "~/Workspace/pim_kinase/ispy1_082814/raw_data/design.csv"

    data <- MicroarrayDataset(exprs.file, annot.file, design.file)
    data <- collapseProbes(data, method = 'MaxMean')

    data
}
