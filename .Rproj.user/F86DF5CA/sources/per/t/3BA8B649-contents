
# volcano plot Females

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
library(organism, character.only = TRUE)
library(DOSE)
library(EnhancedVolcano)
library(tidyverse)
grcm38 # mouse genes


my_logFC_threshold = 0.2




y1a <- readRDS("results/Female_limma_results.RDS")


dc <- y1a %>% mutate(contrast = "Oil vs. TP") %>% mutate(log10 = -log10(P.Value))


dc$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
dc$diffexpressed[dc$logFC < 0.2 & dc$P.Value < 0.05] <- "UP"

# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
dc$diffexpressed[dc$logFC > -0.2 & dc$P.Value < 0.05] <- "DOWN"
## REMEMBER, NEGATIVE LOGFC = HIGHER IN EXPT


dcx <- dc %>% filter(.,logFC >= .4)%>% filter(P.Value < 0.05)
dcxx <- dc %>% filter(.,logFC <= -.45) %>% filter(P.Value < 0.05)

dcxxx<- dc %>%
  filter(logFC >= 0.2 | logFC <= -0.2) %>%
  filter(P.Value < 0.05) %>% 
  arrange(-abs(logFC), P.Value) %>% 
  slice_head(n = 10)

dc$log10 <- ifelse(dc$log10 == Inf, 4,dc$log10)


vp_Fem <- ggplot(data = dc, 
                         aes(x = logFC, 
                             y = log10, 
                             colour=diffexpressed)) +
  geom_point(alpha=0.25, size=3.5) +
  scale_color_manual(values=c("#B00B69", "grey","#EF7627"))+
  xlim(c(-1.5, 1.5)) +
  ylim(0,4)+
  geom_vline(xintercept=c(-.2,.2),lty=4,col="black",lwd=0.8) +
  #  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 1.301,lty=4,col="black",lwd=0.8) +
  #  geom_text_repel(data = dcx, aes(x = logFC, y = -log10(P.Value),label = symbol), color = "black",vjust =1, hjust =.45)+
  #  geom_text_repel(data = dcxx, aes(x = logFC, y = -log10(P.Value),label = symbol), color = "black", hjust = 1, vjust = -.55)+
  geom_text_repel(data = dcxxx, aes(x = logFC, y = -log10(P.Value),label = symbol), color = "black", hjust = 0, vjust = .5, max.overlaps = 7)+
  labs(title = "Differential Gene Expression in MBH - Females",
       x="log2 Fold Change",
       y=bquote(~-Log[10]~italic(eFDR)))+
  theme_bw() +
  annotate(geom="text", x=1, y=2.75, label=paste0(" ", "\n", "Lower in TP"),
           color="black", size = 5)+
  annotate(geom="text", x=-1, y=2.75, label=paste0(" ", "\n", "Higher in TP"),
           color="black", size = 5)+
  scale_x_continuous(limits = c(-1.5,1.5),breaks = c(-1.5,-1,-.5,0,.5,1,1.5))+
  theme(axis.text.x = element_text(vjust = 1,size = 20),
        # axis.ticks = element_blank(),
        axis.text.y = element_text(hjust = 0.5,size = 20),
        axis.text = element_text(color="#3C3C3C",size = 20),
        axis.title = element_text(size = 20),   
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 20),
        text = element_text(size = 20)
        
  )

vp_Fem #875 x 550


top50sigdiff_genes_Female <- y1a %>% 
  filter(logFC >= 0.2 | logFC <= -0.2) %>%
  filter(P.Value < 0.05) %>% 
  arrange(-abs(logFC), P.Value) %>% 
  select(symbol, logFC, P.Value, chr, description) %>% 
  head(50)


write.csv(top50sigdiff_genes_Female,"results/top50sigdiff_genes_Female.csv", row.names = F)


# Log2FC of 0.5 as cutoff - 19 total
y1a %>% 
  filter(logFC >= 0.5 | logFC <= -0.5) %>%
  filter(P.Value < 0.05) %>% 
  arrange(-abs(logFC), P.Value) %>% 
  select(symbol, logFC, P.Value, chr, description)

y1a %>% filter(chr == "X") %>% arrange(-abs(logFC), P.Value)

y1a %>% filter(symbol == "Sts")

cts %>% filter(X == "Sts")
