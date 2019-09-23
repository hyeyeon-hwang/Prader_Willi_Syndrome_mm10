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

fc.output <- fc.1.2

# test
tail(fc.14$Geneid)
head(geneIdList)
head(geneMatch.14)
head(fc.14)

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
  s.1.2 = geneCount.1.2,
  s.3.4 = geneCount.3.4,
  s.5.6 = geneCount.5.6,
  s.7.8 = geneCount.7.8,
  s.9.10 = geneCount.9.10,
  s.11.12 = geneCount.11.12,
  s.13 = geneCount.13,
  s.14 = geneCount.14
)

countDf[is.na(countDf)] <- 0

rownames(countDf) <- countDf$gene.id
countDf

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


dds <- DESeqDataSetFromMatrix(countData = countMat,
                              colData = columnData,
                              design = ~ condition)



