setwd("/Volumes/Samsung_T5/IECOH_SGP/Correlation_new")
library(readxl)
clinical_data_BL <- read_excel("clinical_data_BL.xlsx")
#View(clinical_data_BL)    
clinical_data_FU <- read_excel("clinical_data_FU.xlsx")
#View(clinical_data_FU)   

#to match subject ID with SGP sample ID
mapping_v1 <- read_excel("mapping_v1.xlsx")
View(mapping_v1)                                                                                                                            
library(readxl)
mapping_v2 <- read_excel("mapping_v2.xlsx")

clinical_data_BL_merged <- merge(clinical_data_BL, mapping_v1, by = "SubjectID")
clinical_data_FU_merged <- merge(clinical_data_FU, mapping_v2, by = "SubjectID")

write.csv(clinical_data_BL_merged, file = "clinical_data_BL_merged.csv")
write.csv(clinical_data_FU_merged, file = "clinical_data_FU_merged.csv")


################################################################################
#import clinical data
clinical_data <- read_excel("clinical_data.xlsx")
#import taxa abundance data
level6 <-read_excel("level6.xlsx")
#transpose dataset and maintain the first column as heading
level6_trans <-t(level6)
level6_trans <- as.data.frame(t(level6[,-1]))
colnames(level6_trans) <- level6$SampleID

#calculate relative abundance
level6_trans_relative = apply(level6_trans,2,function(x){(x/sum(x))})

#transpose back
level6_relative_t <-as.data.frame(t(level6_trans_relative))
#row name to column
library(dplyr)
library(tidyverse)
level6_relative_t  <- tibble::rownames_to_column(level6_relative_t , "SampleID")
#import MSD data
MSD_data <- read_excel("SGP_MSD_mapping_updated.xlsx")

#merge datasets
master_data_1 <-merge(level6_relative_t , MSD_data, by ="SampleID")
master_data <-merge(master_data_1, clinical_data, by = "SampleID")
#remove columns not used
master_data <- master_data[ -c(145:147, 164:166) ]
master_data <- master_data[ -c(155) ]#remove subjectID
#column to row names
master_data<- master_data %>% remove_rownames %>% column_to_rownames(var="SampleID")

write.csv(master_data, file="master_data.csv")
#########################################################################################
#Computing correlation matrix
master_data_cor<-cor(master_data, use = "complete.obs")#default is pearson method!!!
master_data_cor <-as.data.frame(master_data_cor)
write.csv(master_data_cor, file = "master_data_cor.csv")
#load function
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
#test
library("Hmisc")
res2 <- rcorr(as.matrix(master_data), type = "pearson")
res2
sig <- flattenCorrMatrix(res2$r, res2$P)
write.csv(sig, file="cor_significance_pearson.csv")
#########################################################################################
#clean the data
cor_data <- master_data_cor[ -c(1:143) ]
#select top20 genus only
#import names of the top20 genus
top20 <- read_excel("top20_genus.xlsx")
#column to row names
top20<- top20 %>% remove_rownames %>% column_to_rownames(var="SampleID")
#merge cor table with top20 genus list
top20_cor_data<-merge(top20, cor_data, by = "row.names")
top20_cor_data <-top20_cor_data[-c(1) ]
top20_cor_data<- top20_cor_data %>% remove_rownames %>% column_to_rownames(var="Genus")

write.csv(top20_cor_data, file = "top20_cor_data.csv")
#import clean data
top20_cor_data <- read_excel("top20_cor_data_plus.xlsx")
top20_cor_data<- top20_cor_data %>% remove_rownames %>% column_to_rownames(var="...1")

top20_cor_neg_data <- read_excel("top20_cor_data_minus.xlsx")
top20_cor_neg_data<- top20_cor_neg_data %>% remove_rownames %>% column_to_rownames(var="...1")
#circlize plot
library(circlize)
grid_col <-  c(Actinomyces = '#E6E6FA',
               Rothia = '#D8BFD8',
               Corynebacterium = '#DDA0DD',
               Porphyromonas = '#EE82EE',
               Alloprevotella = '#DA70D6',
               Prevotella = '#BA55D3',
               Capnocytophaga = '#9370DB',
               Streptococcus = '#8A2BE2',
               Selenomonas = '#BA55D3',
               Dialister = '#FFA07A',
               Veillonella = '#FA8072',
               Fusobacterium = '#DC143C',
               Leptotrichia = '#FF1493',
               Lautropia = '#4B0082',
               Neisseria = '#FFC0CB',
               Campylobacter = '#FFB6C1',
               Haemophilus = '#DB7093',
               'Saccharibacteria_(TM7)_[G-1]' = '#C71585',
               'Saccharibacteria_(TM7)_[G-5]' = '#FF69B4',
               Treponema = '#FF0000',
               'IFN-γ'= '#333300',
               'IL-10'='#666600',
               'IL-12p70'='#999900',
               'IL-13'='#CCCC00',
               'IL-1β'='#FFFF00',
               'IL-2'='#FFFF33',
               'IL-4'='#FFFF66',
               'IL-6'='#FFFF99',
               'IL-8'='#FFA500',
               'TNF-α'='#FF8C00',
               PD.avg.site = '#00008B',
               saliva.rate = '#2F4F4F',
               BoP.rate = '#1E90FF',
               cm_level = '#5F9EA0',
               cotinine = '#778899'
               )

chordDiagram(t(top20_cor_data), grid.col = grid_col)
circos.clear()


chordDiagram(t(top20_cor_neg_data), grid.col = grid_col)
circos.clear()




