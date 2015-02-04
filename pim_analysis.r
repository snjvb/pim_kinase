setwd('~/brca_datasets/pim_analysis')

library(bear)
library(ggplot2)
library(grid)
library(limma)
library(survival)

TCGA.BRCA <- function() {
    source('~/brca_datasets/TCGA_BRCA_exp_HiSeqV2-2014-08-28/tcga_brca.r')

    TCGA.BRCA()
}

ISPY1_DATASET <- function() {
    source('~/brca_datasets/ispy1_082814/ispy1.r')

    ISPY1()
}

GSE25066_DATASET <- function() {
    source('~/brca_datasets/gse25066_091114/gse25066.r')

    ## remove ISPY1 samples from dataset
    data <- GSE25066()
    data$design <- subset(data$design, data$design$source != 'ISPY')
    
    syncArrayData(data)
}

YAU_DATASET <- function() {
    source('~/brca_datasets/YauGeneExp-2011-11-11/yau_gene_exp.r')

    ## remove samples without survival data from dataset
    data <- YauGeneExp()
    data$design <- subset(data$design, 
                          !is.na(data$design$e_dmfs) & 
                          !is.na(data$design$t_dmfs))

    syncArrayData(data)
}

discretizeSamplesByQuantile <- function(exprs, time, event, .quantile) {
    rank <- round(.quantile * length(exprs) / 100)
    groups <- discretizeSamplesByRank(exprs, rank)
    reg <- coxph(Surv(time, event) ~ groups)

    list(P = summary(reg)$sctest['pvalue'], 
         quantile = .quantile, 
         groups = groups)
}

discretizeSamplesByRank <- function(exprs, rank) {
    .exprs <- exprs[order(exprs)]
    cutoff <- .exprs[rank]

    discretizeSamples(exprs, cutoff)
}

discretizeSamples <- function(exprs, cutoff) {
    ifelse(exprs < cutoff, 1, 2)
}

getOptimalQuantile <- function(gene.exprs, time, event, min.quantile = 0.1) {
    start <- ceiling(min.quantile * length(gene.exprs))
    end <- floor((1 - min.quantile) * length(gene.exprs))

    p.dist <- lapply(start:end, function(rank) {
        groups <- discretizeSamplesByRank(gene.exprs, rank)
        reg <- coxph(Surv(time, event) ~ groups)

        list(P = summary(reg)$sctest['pvalue'], 
             quantile = rank / length(gene.exprs) * 100, 
             groups = groups)
    })

    optimal <- which.min(sapply(p.dist, function(x) x$P))
    p.dist[[optimal]]
}

drawKMPlot <- function(filename, title.prefix, gene.exprs, time, event, 
                       .quantile = NULL, ...) {
    temp <- if (is.null(.quantile)) {
        getOptimalQuantile(gene.exprs, time, event)
    } else {
        discretizeSamplesByQuantile(gene.exprs, time, event, .quantile)
    }
    fit <- survfit(Surv(time, event) ~ temp$groups)

    bitmap(filename, res = 600, type = 'jpeg', height = 9, width = 9, 
           units = 'in')
    col <- c('black', 'red')
    title <- sprintf('%s (n = %i)\nLog-rank p = %.3e', 
                     title.prefix, length(temp$groups), temp$P)
    par(mar = c(5, 6, 4, 2))
    plot(fit, col = col, main = title, ylim = c(0.3, 1), cex.lab = 2, 
         cex.axis = 1.7, cex.main = 1.5, las = 1, lwd = 3, font.lab = 2, ...)
    legend('topright', 
           legend = c(
                sprintf('< %.2f percentile (%i)', temp$quantile, 
                        sum(temp$groups == 1)), 
                sprintf('>= %.2f percentile (%i)', temp$quantile, 
                        sum(temp$groups == 2))), 
           col = col, 
           lty = c(1, 1), 
           lwd = c(5, 5), 
           cex = 1.7)
    dev.off()
}

generateDataForPlot <- function(gene.exprs, status) {
    .gene.exprs <- gene.exprs - median(gene.exprs)
    df <- data.frame(status = status, expr = .gene.exprs)
    df <- subset(df, !is.na(df$status))

    ## Perform pair-wise unpaired t-tests with pooled SD
    sig <- pairwise.t.test(df$expr, df$status, p.adj = 'none')$p.value

    .df <- summarySE(df, measurevar = 'expr', groupvars = 'status')
    .df$status <- as.factor(sprintf('%s\n(n=%i)', .df$status, .df$N))

    list(plot.data = .df, 
         significance = sig)
}

drawBarPlot <- function(filename, df) {
    plot <- ggplot(data = df, aes(x = status, y = expr)) + 
        geom_bar(aes(fill = status), position = position_dodge(), 
                 stat = 'identity', color = 'black', size = 1.3, 
                 width = 0.6) +
        geom_errorbar(aes(ymin = ifelse(expr < 0, expr - se, expr), 
                          ymax = ifelse(expr < 0, expr, expr + se)), 
                      width = 0.2, position = position_dodge(.7), size = 1) +
        geom_hline(y = 0, size = 1.5) + 
        scale_fill_manual(values = c('black', 'grey', 'white')) + 
        ylab('Median-centered log expression') + 
        theme_bw() + 
        theme(
            legend.position = 'none', 
            axis.line = element_line(size = 1.3), 
            axis.ticks = element_line(size = 1.5), 
            axis.ticks.length = unit(3, 'mm'), 
            axis.text.x = element_text(size = 20, face = 'bold', 
                                       angle = 45, hjust = 1, vjust = 1), 
            axis.text.y = element_text(size = 20), 
            axis.title.x = element_blank(), 
            axis.title.y = element_blank(), 
            panel.border = element_blank(), 
            panel.grid = element_blank())

    ggsave(filename = filename, plot = plot, height = 6.75, width = 4.5, 
           dpi = 600)
}

univariateAnalysis <- function(gene.exprs, group, time, event) {
    valid <- !is.na(group)
    .gene.exprs <- scale(gene.exprs[valid])
    .group <- group[valid]
    .time <- time[valid]
    .event <- event[valid]

    subsetAnalysis <- function(subset) {
        ..exprs <- .gene.exprs[subset]
        ..time <- .time[subset]
        ..event <- .event[subset]

        reg <- coxph(Surv(..time, ..event) ~ ..exprs)
        summ <- summary(reg)
        conf.int <- summ$conf.int
        p <- summ$logtest['pvalue']

        data.frame(N = length(..exprs), 
                   HazardRatio = sprintf('%.3f (%.3f - %.3f)',
                                         conf.int[, 1], 
                                         conf.int[, 3], 
                                         conf.int[, 4]), 
                   p = p, 
                   stringsAsFactors = F)
    }

    subsets <- cbind(
        rep(T, length(.gene.exprs)), 
        sapply(unique(.group), function(g) .group == g))
    
    res <- apply(subsets, 2, subsetAnalysis)
    comp.groups <- names(res)
    comp.groups[1] <- 'Overall'

    res <- rbind.fill(res)
    rownames(res) <- comp.groups

    res
}

##########
## Main ##
##########

## TCGA
tcga <- TCGA.BRCA()
tcga.df <- generateDataForPlot(as.numeric(tcga$exprs['PIM1', ]), 
                               factor(tcga$design$status, 
                                      levels = c('HER2+', 'HR+HER2-', 'TN')))
drawBarPlot('tcga_PIM1_exprs_barplot.jpg', tcga.df$plot.data)
write.csv(tcga.df$significance, file = 'tcga_PIM1_exprs_pairwise_t.csv')

## ISPY1
ispy1 <- ISPY1_DATASET()
sample.subset <- !is.na(ispy1$design$hr) & ispy1$design$hr == 'HR-'
drawKMPlot('ispy1_PIM1_kaplan_meier.jpg', 
           'I-SPY', 
           as.numeric(ispy1$exprs['5292', sample.subset]), 
           ispy1$design[sample.subset, 'rfs.t'], 
           ispy1$design[sample.subset, 'rfs.e'], 
           xlab = 'Time (years)')
ispy.df <- generateDataForPlot(as.numeric(ispy1$exprs['5292', ]), 
                               factor(ispy1$design$status, 
                                      levels = c('HER2+', 'HR+HER2-', 'TN')))
drawBarPlot('ispy1_PIM1_exprs_barplot.jpg', ispy.df$plot.data)
write.csv(ispy.df$significance, file = 'ispy1_PIM1_exprs_pairwise_t.csv')
ispy1.res <- univariateAnalysis(as.numeric(ispy1$exprs['5292', ]), 
                                ispy1$design$hr, 
                                ispy1$design$rfs.t, 
                                ispy1$design$rfs.e)
write.csv(ispy1.res, file = 'ispy1_univariate_analysis.csv')

## GSE25066
gse25066 <- GSE25066_DATASET()
sample.subset <- !is.na(gse25066$design$hr) & gse25066$design$hr == 'HR-'
drawKMPlot('gse25066_PIM1_kaplan_meier.jpg', 
           'Pooled Neoadjuvant Chemotherapy Treated', 
           as.numeric(gse25066$exprs['5292', sample.subset]), 
           gse25066$design[sample.subset, 'drfs_t'], 
           gse25066$design[sample.subset, 'drfs_e'], 
           xlab = 'Time (years)')
gse25066.df <- generateDataForPlot(as.numeric(gse25066$exprs['5292', ]), 
                                   factor(gse25066$design$status, 
                                          levels = c('HER2+', 
                                                     'HR+HER2-', 
                                                     'TN')))
drawBarPlot('gse25066_PIM1_exprs_barplot.jpg', gse25066.df$plot.data)
write.csv(gse25066.df$significance, file = 'gse25066_PIM1_exprs_pairwise_t.csv')
gse25066.res <- univariateAnalysis(as.numeric(gse25066$exprs['5292', ]), 
                                   gse25066$design$hr, 
                                   gse25066$design$drfs_t, 
                                   gse25066$design$drfs_e)
write.csv(gse25066.res, file = 'gse25066_univariate_analysis.csv')

## Yau dataset
yau <- YAU_DATASET()
sample.subset <- !is.na(yau$design$er_status) & yau$design$er_status == 'ER-'
drawKMPlot('yau_PIM1_kaplan_meier.jpg', 
           'Pooled Node-Negative Adjuvant Chemotherapy Naive', 
           as.numeric(yau$exprs['PIM1', sample.subset]), 
           time = yau$design[sample.subset, 't_dmfs'], 
           event = yau$design[sample.subset, 'e_dmfs'], 
           xlab = 'Time (years)')
yau.df <- generateDataForPlot(as.numeric(yau$exprs['PIM1', ]), 
                              factor(yau$design$status, 
                                     levels = c('HER2+', 
                                                'ER+HER2-', 
                                                'ER-HER2-')))
drawBarPlot('yau_PIM1_exprs_barplot.jpg', yau.df$plot.data)
write.csv(yau.df$significance, file = 'yau_PIM1_exprs_pairwise_t.csv')
yau.res <- univariateAnalysis(as.numeric(yau$exprs['PIM1', ]), 
                              yau$design$er_status, 
                              yau$design$t_dmfs, 
                              yau$design$e_dmfs)
write.csv(yau.res, file = 'yau_univariate_analysis.csv')
