
# initial cleaning and checks

library(limma)
library(tidyverse)
library(scales)


library(tidyverse)
library(DESeq2)
library(pheatmap)
library(annotables)
library(ggpubr)
library(lsr)

#source("functions/geom_boxjitter.R")


# load data 

rawcounts <- read.csv("rawdata/pnd2era_counts.csv")
rawcounts <- rawcounts[-1,]
rownames(rawcounts) <- NULL
bcx1 <- rawcounts %>% column_to_rownames(.,var = "X")

id <- read.csv("rawdata/id.csv")
idxx <- id %>% column_to_rownames(.,var = "id")

all(rownames(idxx) == colnames(bcx1)) #check they match

#### count checks

colSums(bcx1[,]) %>% 
  as_tibble_row() %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'id') %>% 
  dplyr::rename(genecounts = V1) -> genecounts

genecounts <- genecounts %>%
  full_join(id) 


#########

genecounts %>% 
  ggplot(aes(genecounts)) +
  geom_histogram(bins = 54,alpha =0.5,color = 'grey') +
  theme_classic() 

counts <- genecounts %>%
  ggplot(aes(id,genecounts))+
  geom_bar(stat = 'identity')+
  coord_flip()

##

# genecounts %>%
#   ggplot(aes(Treatment,genecounts, color = Treatment))+
#   geom_boxjitter(outlier.color = NA, jitter.shape = 21, 
#                  alpha = 1,
#                  width = 0.5,
#                  jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
#                  position = position_dodge(0.6)) +
#   theme_classic()+
#   theme(legend.position = "none")


genecounts %>%
  ggplot(aes(Treatment,genecounts, color = Treatment))+
  geom_boxplot() +
  geom_point()+
  theme_classic()




##

t.test(genecounts~Treatment, genecounts)

genecounts %>%
  group_by(Treatment) %>%
  summarise(sd_var1 = sd(genecounts, na.rm=TRUE))
cohensD(genecounts~Treatment, genecounts)








