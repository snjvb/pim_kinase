library(GEOquery)

getMetaData <- function(adf) {
    characteristics <- apply(adf[, 10:33], 1, paste, collapse = "||")

    getField <- function(field) {
        pattern <- sprintf(".*\\|\\|%s: (.+?)\\|\\|.*", field)

        a <- characteristics
        a[!grepl(pattern, a, perl = T)] <- NA

        gsub(pattern, "\\1", a, perl = T)
    }

    sample_id <- getField("sample id")
    source <- getField("source")
    age <- getField("age_years")
    er_status_ihc <- getField("er_status_ihc")
    pr_status_ihc <- getField("pr_status_ihc")
    her2_status <- getField("her2_status")
    clinical_t_stage <- getField("clinical_t_stage")
    clinical_nodal_status <- getField("clinical_nodal_status")
    clinical_ajcc_stage <- getField("clinical_ajcc_stage")
    grade <- getField("grade")
    pathologic_response_pcr_rd <- getField("pathologic_response_pcr_rd")
    pathologic_response_rcb_class <- getField("pathologic_response_rcb_class")
    drfs_e <- getField("drfs_1_event_0_censored")
    drfs_t <- getField("drfs_even_time_years")
    type_taxane <- getField("type_taxane")
    esr1_status <- getField("esr1_status")
    erbb2_status <- getField("erbb2_status")
    set_class <- getField("set_class")
    chemosensitivity_prediction <- getField("chemosensitivity_prediction")
    ggi_class <- getField("ggi_class")
    pam50_class <- getField("pam50_class")
    dlda30_prediction <- getField("dlda30_prediction")
    rcb_0_i_prediction <- getField("rcb_0_i_prediction")
    tissue <- getField("tissue")

    data.frame(
        sample_id, 
        source, 
        age, 
        er_status_ihc, 
        pr_status_ihc, 
        her2_status, 
        clinical_t_stage, 
        clinical_nodal_status, 
        clinical_ajcc_stage, 
        grade, 
        pathologic_response_pcr_rd, 
        pathologic_response_rcb_class, 
        drfs_e, 
        drfs_t, 
        type_taxane, 
        esr1_status, 
        erbb2_status, 
        set_class, 
        chemosensitivity_prediction, 
        ggi_class, 
        pam50_class, 
        dlda30_prediction, 
        rcb_0_i_prediction, 
        tissue, 
        row.names = rownames(adf), 
        stringsAsFactors = F)
}

## extract the raw datasets and probe annotations
gse25066 <- getGEO(filename = "GSE25066_series_matrix.txt")

## extract information from raw dataset
features <- featureData(gse25066)@data
exprs <- exprs(gse25066)
meta.data <- getMetaData(pData(gse25066))

## write out exprs and design matrix
write.csv(features, file = "annot.csv")
write.csv(exprs, file = "exprs.csv")
write.csv(meta.data, file = "design.csv")
