

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


write.csv(top10_GOterms_Female,"results/top100_GOterms_Female.csv", row.names = F)


# indiv gene lookup

Fem %>% filter(symbol == "Esr1")




#### searching GO terms for common denominator genes 

treatgo <- read.csv("results/top100_GOterms_Female.csv")

# First, process the treatgo data as before
gene_counts <- treatgo %>%
  separate_rows(geneID, sep = "/") %>%
  count(direction, geneID, name = "gene_count") %>%
  arrange(direction, desc(gene_count)) %>%
  group_by(direction) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# Join with mgenes to get the gene descriptions
result_with_desc <- gene_counts %>%
  # Join with mgenes using geneID = symbol
  left_join(
    mgenes %>% select(symbol, gene_description = description),
    by = c("geneID" = "symbol")
  ) %>%
  # Now join with original data to get the treatgo descriptions
  left_join(
    treatgo %>%
      separate_rows(geneID, sep = "/") %>%
      select(direction, geneID, treatgo_description = Description),
    by = c("direction", "geneID")
  ) %>%
  # Group and combine the treatgo descriptions
  group_by(direction, geneID, gene_count, rank, gene_description) %>%
  summarise(
    treatgo_descriptions = paste(unique(treatgo_description), collapse = "; "),
    .groups = "drop"
  ) %>%
  # Keep only top genes
  filter(rank <= 10) %>%
  arrange(direction, desc(gene_count))

result_with_desc



GO_top_genes <- result_with_desc %>%
  mutate(
    treat_data = map(geneID, ~ Fem %>% filter(symbol == .x))
  )%>%
  unnest(treat_data, keep_empty = TRUE)

## NONE of the top genes in the GO terms are significantly up or down regulated themselves.

GO_top_genes_table_Fem <- GO_top_genes %>% 
  select(1,2,3,4,5,6,8,9,10)



write.csv(GO_top_genes_table_Fem,"results/GO_top_genes_table_Fem.csv", row.names = F)

