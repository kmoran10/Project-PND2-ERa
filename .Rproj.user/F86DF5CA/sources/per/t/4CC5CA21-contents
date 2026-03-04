

# female GO analysis

library(limma)
library(edgeR)
library(Mus.musculus)
organism = 'org.Mm.eg.db'
library(organism, character.only = TRUE)
library(biomaRt)
library(AnnotationDbi)
library(pheatmap)
library(annotables)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(tidyverse)
mgenes <- grcm38 # mouse genes


source("functions/gettop10GO.R")

my_logFC_threshold = 0.2

Fem <- readRDS("results/Female_limma_results.RDS")

gettop10GO(Fem, my_showCategory) %>% 
  mutate(comparison = "Oil - TP") -> top10_GOterms_Female


# indiv gene lookup

Fem %>% filter(symbol == "Esr1")






