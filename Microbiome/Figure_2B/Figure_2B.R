#Author: Ziyan Lin
#Description: this is the script for plotting Figure 2B
######################
#library package
library(phyloseq)
library(biomformat)
library(ggplot2)
library(vegan)
library(ggsignif)
library(ggpubr)
library(plotly)
library("RColorBrewer")
library(dplyr)
library(picante)
library(rstatix)
library(tidyr)
library(gridExtra)

options(scipen = 999)

############################################
#Figure 2B
############################################

V1 <- read.csv("Alp-Groups-V1-data.csv",row.names=1)
V2 <- read.csv("Alp-Groups-V2-data.csv",row.names=1)
meta <- read.delim("metadata.tsv")

merge.df1 <- merge(V1,meta,by="SampleID")
merge.df2 <- merge(V2,meta,by="SampleID")
df1 <- subset(merge.df1,select=c("SampleID","Visit.x","Group.x","PatientID","variable","value"))
colnames(df1) <- c("SampleID","Visit","Group","PatientID","variable","value")
df1 <- df1[order(df1$PatientID,df1$variable),]

df2 <- subset(merge.df2,select=c("SampleID","Visit.x","Group.x","PatientID","variable","value"))
colnames(df2) <- c("SampleID","Visit","Group","PatientID","variable","value")
df2 <- df2[order(df2$PatientID,df2$variable),]

df <- subset(df1,select=c("PatientID","Group","variable"))
df$SampleID_V1 <- df1$SampleID
df$SampleID_V2 <- df2$SampleID
df$value_V1 <- df1$value
df$value_V2 <- df2$value

df$value <- df$value_V2 - df$value_V1

myColor <- c("firebrick3","goldenrod2","dodgerblue2")

set.seed(123)
#genere box plot
p<-NULL
p <- ggplot(df,aes(x=Group, y=value))+
facet_wrap(~variable,scales = "free")+
geom_boxplot(width=0.8,outlier.shape = NA) +
geom_jitter(aes(color=Group),width=.15) +
theme_classic() +
xlab("")+
ylab("Alpha Diversity Measure")+
theme(axis.text.x=element_text(angle=0,size=12,color="black"),
      axis.text.y=element_text(size=10,color="black"))+
theme(strip.background = element_rect(fill=NULL,linetype="blank"), #facet background
      strip.text.x = element_text(size=10,face = "bold.italic")) #facet font size
p <- p + scale_color_manual(values=myColor,
                            name="Group")

#prepare for p-value (wilcox.test)
lev <- levels(as.factor(df$Group)) # get the variables
L.pairs <- combn(seq_along(lev), 2, simplify = FALSE, FUN = function(i) lev[i])# make a pairwise list that we want to compare.

#add p value(significant level)
#p <- p + stat_compare_means(
#  comparisons = L.pairs,
#  hide.ns = T,
#  label = "p.signif",#"p.format"
#  method = "wilcox.test",
#  paired = FALSE
#)
            

pdf(paste0("Alp-Group-V2-V1.pdf"),height=6,width=6)
print(p)
dev.off()

#output alpha p-value
pval_tbl <- compare_means(formula=value ~ Group,data=df,method="wilcox.test",group.by="variable")
write.csv(as.data.frame(pval_tbl),paste0("Alp-Group-pval-V2-V1.csv"))
