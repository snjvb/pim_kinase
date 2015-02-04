library(FactoMineR)
library(GEOquery)
library(sva)

getExprs <- function(eset, probes) {
    ## replace missing values with the lowest possible value for the gene
    e <- t(apply(exprs(eset), 1, function(r) {
        r[is.na(r)] <- min(r, na.rm = T)
        r
    }))

    ## collapse data from duplicate probes
    features <- as.character(featureData(eset)@data$SPOT_ID)
    e <- aggregate(e, by = list(features), FUN = mean)

    ## return expression values for the annotated probeset
    subset(e, e[, 1] %in% probes)
}

getMetaData <- function(adf) {
    characteristics <- apply(adf[, 19:38], 1, paste, collapse = "||")

    getField <- function(field) {
        pattern <- sprintf(".*\\|\\|%s: (.+?)\\|\\|.*", field)

        a <- characteristics
        a[!grepl(pattern, a, perl = T)] <- NA

        gsub(pattern, "\\1", a, perl = T)
    }

    ispy.id <- getField("i-spy id")
    histologic.grade <- getField("histologic grade \\(1=grade i \\(low\\);  2= grade ii \\(intermediate\\); 3= grade iii \\(high\\); 4=indeterminate\\)")
    histology <- getField("histology \\(1=necrosis; 2=ductal carcinoma; 3=lobular; 4=mixed ductal/lobular carcinoma; 5=other 6=no invasive tumor present;\\)")
    clinical.tumor.size <- getField("clinical tumor size \\(pre-chemotherapy\\), cm \\(-1 = n/a\\)")
    clinical.tumor.stage <- getField("clinical t stage \\(pre-chemotherapy\\) \\(1 = <= 2cm; 2 = >2-5 cm; 3 = > 5cm\\)")
    er <- getField("er \\(0=negative; 1=positive;\\)")
    pr <- getField("pgr  \\(0=negative; 1=positive;\\)")
    her2 <- getField("her2 \\(0=negative; 1=positive;\\)")
    receptor <- apply(cbind(er, pr, her2), 1, function(r) {
        if (!any(is.na(r)) & all(r == "0"))
            "TN"
        else if (any(r == "1", na.rm = T))
            "NTN"
        else
            NA
    })
    subtype <- getField("intrinsic subtype by pam50")
    pcr <- getField("pathological complete response \\(pcr\\)")
    rfs.t <- getField("relapse-free survival time – time from chemo start date until earliest: local or distant progression or death \\(time unit is days\\)")
    rfs.e <- getField("relapse-free survival indicator \\(1=event; local or distant progression or death, 0=censor at last follow-up\\)")
    os.t <- getField("overall survival - time from study chemo start date to death or last follow-up; \\(time unit is days\\)")
    os.status <- getField("survival status \\(7=alive; 8=dead, 9=lost\\)")
    platform <- adf$platform_id

    data.frame(
        ispy.id, 
        histologic.grade, 
        histology, 
        clinical.tumor.size, 
        clinical.tumor.stage, 
        er, 
        pr, 
        her2, 
        receptor, 
        subtype, 
        pcr, 
        rfs.t, 
        rfs.e, 
        os.t, 
        os.status, 
        platform, 
        row.names = rownames(adf), 
        stringsAsFactors = F)
}

PCAClustering <- function(df, groups) {
    pca.exprs <- cbind(groups, as.data.frame(t(df)))
    pca.res <- PCA(pca.exprs, quali.sup = 1, ncp = 5, graph = F)

    concat <- cbind.data.frame(pca.exprs[, 1], pca.res$ind$coord)
    ellipse.coord <- coord.ellipse(concat, level.conf = 0.99, bary = T)

    plot.PCA(pca.res, axes=c(1, 2), choix = "ind", habillage = 1, 
        ellipse = ellipse.coord, label = "quali")
}

## extract the raw datasets and probe annotations
gpl1708 <- getGEO(filename = "GSE22226-GPL1708_series_matrix.txt")
gpl4133 <- getGEO(filename = "GSE22226-GPL4133_series_matrix.txt")
annot <- read.delim("annot.txt", as.is = T, fill = T, 
    check.names = F)

## prepare each dataset individually
## GSE22226-GPL1708
gpl1708.features <- featureData(gpl1708)@data
gpl1708.exprs <- getExprs(gpl1708, annot$ProbeID)
gpl1708.meta.data <- getMetaData(pData(gpl1708))

## GSE22226-GPL4133
gpl4133.features <- featureData(gpl4133)@data
gpl4133.exprs <- getExprs(gpl4133, annot$ProbeID)
gpl4133.meta.data <- getMetaData(pData(gpl4133))

## merge datasets and cleanup samples without annotations
exprs <- merge(gpl1708.exprs, gpl4133.exprs, by = 1)
rownames(exprs) <- exprs[, 1]
exprs <- data.matrix(exprs[, -1])
meta.data <- rbind(gpl1708.meta.data, gpl4133.meta.data)

good.samples <- rownames(meta.data)[!is.na(meta.data$ispy.id)]
exprs <- exprs[, colnames(exprs) %in% good.samples]
meta.data <- meta.data[rownames(meta.data) %in% good.samples, ]

## pre-ComBat PCA analysis
png("pre_combat_batch_clustering.png", 1024, 1024, "px")
PCAClustering(exprs, meta.data$platform)
dev.off()

## correcting for batch effects with ComBat
merged.exprs <- ComBat(exprs, batch = meta.data$platform, mod = NULL, 
    numCovs = NULL, par.prior = T, prior.plots = F)

## post-ComBat PCA analysis
png("post_combat_batch_clustering.png", 1024, 1024, "px")
PCAClustering(merged.exprs, meta.data$platform)
dev.off()

## write out exprs and design matrix
write.csv(merged.exprs, file = "exprs.csv")
write.csv(meta.data, file = "design.csv")
