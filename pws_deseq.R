# References
# Analyzing RNA-seq data wth DESeq2:
# 	http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#count-matrix-input
# Beginner's guide to DESeq2:
# 	https://bioc.ism.ac.jp/packages/2.14/bioc/vignettes/DESeq2/inst/doc/beginner.pdf

library(tidyverse)
#BiocManager::install("DESeq2")
library(DESeq2)

fc.1.2 <- read.table(file = "featureCountsoutput/1_2_featureCounts.output.txt",
                     sep = "\t",
                     header = TRUE)[c(1, 7)]
fc.3.4 <- read.table(file = "featureCountsoutput/3_4_featureCounts.output.txt",
                     sep = "\t",
                     header = TRUE)[c(1, 7)]
fc.5.6 <- read.table(file = "featureCountsoutput/5_6_featureCounts.output.txt",
                     sep = "\t",
                     header = TRUE)[c(1, 7)]
fc.7.8 <- read.table(file = "featureCountsoutput/7_8_featureCounts.output.txt",
                     sep = "\t",
                     header = TRUE)[c(1, 7)]
fc.9.10 <- read.table(file = "featureCountsoutput/9_10_featureCounts.output.txt",
                     sep = "\t",
                     header = TRUE)[c(1, 7)]
fc.11.12 <- read.table(file = "featureCountsoutput/11_12_featureCounts.output.txt",
                      sep = "\t",
                      header = TRUE)[c(1, 7)]
fc.13 <- read.table(file = "featureCountsoutput/13_featureCounts.output.txt",
                     sep = "\t",
                     header = TRUE)[c(1, 7)]
fc.14 <- read.table(file = "featureCountsoutput/14_featureCounts.output.txt",
                     sep = "\t",
                     header = TRUE)[c(1, 7)]

geneIdList <- c(fc.1.2$Geneid, fc.3.4$Geneid, fc.5.6$Geneid, fc.7.8$Geneid, 
                fc.9.10$Geneid, fc.11.12$Geneid, fc.13$Geneid, fc.14$Geneid) %>% unique()

# test
# tail(fc.14$Geneid)
# head(geneIdList)
# head(geneMatch.14)
# head(fc.14)

# gene matches
geneMatch.1.2 <- match(geneIdList, fc.1.2$Geneid)
geneMatch.3.4 <- match(geneIdList, fc.3.4$Geneid)
geneMatch.5.6 <- match(geneIdList, fc.5.6$Geneid)
geneMatch.7.8 <- match(geneIdList, fc.7.8$Geneid)
geneMatch.9.10 <- match(geneIdList, fc.9.10$Geneid)
geneMatch.11.12 <- match(geneIdList, fc.11.12$Geneid)
geneMatch.13 <- match(geneIdList, fc.13$Geneid)
geneMatch.14 <- match(geneIdList, fc.14$Geneid)

# gene counts
geneCount.1.2 <- fc.1.2[geneMatch.1.2, 2]
geneCount.3.4 <- fc.3.4[geneMatch.3.4, 2]
geneCount.5.6 <- fc.5.6[geneMatch.5.6, 2]
geneCount.7.8 <- fc.7.8[geneMatch.7.8, 2]
geneCount.9.10 <- fc.9.10[geneMatch.9.10, 2]
geneCount.11.12 <- fc.11.12[geneMatch.11.12, 2]
geneCount.13 <- fc.13[geneMatch.13, 2]
geneCount.14 <- fc.14[geneMatch.14, 2] 

countDf <- data.frame(
  gene.id = geneIdList,
  sample.1.2 = geneCount.1.2,
  sample.3.4 = geneCount.3.4,
  sample.5.6 = geneCount.5.6,
  sample.7.8 = geneCount.7.8,
  sample.9.10 = geneCount.9.10,
  sample.11.12 = geneCount.11.12,
  control.13 = geneCount.13,
  control.14 = geneCount.14
)

countDf[is.na(countDf)] <- 0
rownames(countDf) <- countDf$gene.id

#head(countDf)

countMat <- countDf[, -1] %>% as.matrix()

columnData <- data.frame(
  condition = c("sample", "sample", "sample", "sample", "sample", "sample",
                "control", "control"),
  type = c("paired-end", "paired-end", "paired-end", "paired-end", 
           "paired-end", "paired-end", "paired-end", "paired-end")
)
rownames(columnData) <- c(
  "sample.1.2",
  "sample.3.4",
  "sample.5.6",
  "sample.7.8",
  "sample.9.10",
  "sample.11.12",
  "control.13",
  "control.14"
)

ddsFromMatrix <- DESeqDataSetFromMatrix(countData = countMat,
                              colData = columnData,
                              design = ~ condition)


dds <- DESeq(ddsFromMatrix)
ddsResults <- results(dds)

# Error for vsd: Error in vst(dds, blind = FALSE) : 
# less than 'nsub' rows with mean normalized count > 5, 
# it is recommended to use varianceStabilizingTransformation directly
vsd <- vst(dds, blind=FALSE) 
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
head(assay(rld), 3)

# Section: Effects of transformations on the variance
ntd <- normTransform(dds)
# BiocManager::install("vsn")
library(vsn)
# meanSdPlot(assay(ntd))
# Warning message:
# Computation failed in `stat_binhex()`:
# Package `hexbin` required for `stat_binhex`.
# Please install and try again.


# Heatmap of the count matrix
# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#count-matrix-input

# install.packages("pheatmap")
library(pheatmap)
select <- order(rowMeans(counts(dds, normalized = TRUE)), decreasing = TRUE)

df <- as.data.frame(colData(dds)[,c("condition", "type")])

# heatmap - cluster rows and cols
pdf("heatmap_ntd_clustered_50.pdf")
pheatmap(assay(ntd)[select[1:50],], cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=FALSE, annotation_col=df, main = "Heatmap of count matrix of 50 genes")
dev.off()

# heatmap - cluster rows only
pdf("heatmap_ntd_clustered_rows_50.pdf")
pheatmap(assay(ntd)[select[1:50],], cluster_rows=TRUE, cluster_cols=FALSE, show_rownames=FALSE, annotation_col=df, main = "Heatmap of count matrix of 50 genes")
dev.off()

# heatmap - no clustering
pdf("heatmap_ntd_not_clustered_50.pdf")
pheatmap(assay(ntd)[select[1:50],], cluster_rows=FALSE, cluster_cols=FALSE, show_rownames=FALSE, annotation_col=df, main = "Heatmap of count matrix of 50 genes")
dev.off()



