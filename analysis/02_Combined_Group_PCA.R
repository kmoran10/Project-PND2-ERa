
# PCA


library(tidyverse)
library(DESeq2)


rawcounts <- read.csv("rawdata/pnd2era_counts.csv")
rawcounts <- rawcounts[-1,]
rownames(rawcounts) <- NULL
cts <- rawcounts %>% column_to_rownames(.,var = "X")

id <- read.csv("rawdata/id.csv")
coldata1 <- id %>% column_to_rownames(.,var = "id")

coldata1$group <- paste0(coldata1$Sex, "_",coldata1$Treatment)

all(rownames(coldata1) == colnames(cts)) #check they match
coldata <- coldata1

### 


dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design= ~ group)


dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients


#removing 0s
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

###


vsd <- vst(dds, blind=FALSE)

library("vsn")

meanSdPlot(assay(vsd))

library("RColorBrewer")


plotPCA(vsd, intgroup=c("group"))
plotPCA(vsd, intgroup=c("Treatment"))
plotPCA(vsd, intgroup=c("Sex"))




pcaData <- plotPCA(vsd, intgroup=c("group","Sex","Treatment"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=group.1, label = rownames(pcaData))) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + geom_text() + labs(title="Combined Treatment")




