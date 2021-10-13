#Author: Ziyan Lin
#Description: this is the script for plotting Figure 2D
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

options(scipen = 999)

load("my.physeq.Robj")

############################################
#Figure 2D
############################################
physeq_norm  = physeq #rarefied data
myColor <- c("firebrick3","goldenrod2","dodgerblue2")

plot_title = "Group"
d_method = "bray"
px =  0.45
py =  0.45

#calculate distance
iDist <- distance(physeq_norm, method=d_method)
   
#get coordinate
iPCoA <- ordinate(physeq_norm, "PCoA", distance=iDist)
   
#adonis(p-value)
set.seed(123)
stats<- adonis(formula = iDist ~ sample_data(physeq_norm)$Group)
p_val<-stats$aov.tab$"Pr(>F)"[1]
   
#generate PCoA plot
p <- NULL
p <- plot_ordination(physeq_norm, iPCoA,color="Group") #shape=

pc1 <- gsub("Axis.1  ","PCo1",p$labels$x)
pc2 <- gsub("Axis.2  ","PCo2",p$labels$y)

df <- p$data
p <- ggplot(df,aes(x=Axis.1,y=Axis.2,color=Group,fill=Group)) +
         geom_point(size=3) +
         ggtitle(paste0("PCoA (", d_method,") --- ",plot_title)) +
         xlab(paste0(pc1)) +
         ylab(paste0(pc2)) +
         stat_ellipse(geom = "polygon",type = "norm",level = 0.95,size=1.2,alpha=0.1) + theme_bw()+ #stat_ellipse(type = "t")
         annotate("text", x = px, y = py, label = paste0("p = ",p_val),size=5,fontface = 'italic') +
         theme_classic()+
         theme(axis.text.x=element_text(angle=0,size=13,color="black"),
               axis.text.y=element_text(size=13,color="black"),
               axis.title.x=element_text(size=15),
               axis.title.y=element_text(size=15))

p <- p + scale_color_manual(values=myColor,
                            name  ="Group") +
         scale_fill_manual(values=alpha(myColor),
                            name  ="Group")

pdf(paste0("PCoA-",d_method,"-Groups.pdf"),height=5,width=6)
print(p)
dev.off()

#output PCoA p-value
write.csv(as.matrix(stats$aov.tab),paste0("PCoA-",d_method,"-Groups-pval.csv"))

###############################
#Pairwise
group_list <- c("ES","NS") #replace with each pairwise: c("CS","ES"), c("CS","NS"),c("ES","NS")
plot_title = "ES_vs_NS"
d_method = "bray"

myColor <- c("firebrick3","goldenrod2","dodgerblue2")
myColor <- myColor[-1]

px =  0.45
py =  0.45

sub_physeq <- subset_samples(physeq_norm,Group==group_list[1]|Group==group_list[2]) #select interested group

#calculate distance
iDist <- distance(sub_physeq, method=d_method)
   
#get coordinate
iPCoA <- ordinate(sub_physeq, "PCoA", distance=iDist)
   
#adonis(p-value)
set.seed(123)
stats<- adonis(formula = iDist ~ sample_data(sub_physeq)$Group)
p_val<-stats$aov.tab$"Pr(>F)"[1]
   
#generate PCoA plot
p <- NULL
p <- plot_ordination(sub_physeq, iPCoA,color="Group") #shape=

pc1 <- gsub("Axis.1  ","PCo1",p$labels$x)
pc2 <- gsub("Axis.2  ","PCo2",p$labels$y)

df <- p$data
p <- ggplot(df,aes(x=Axis.1,y=Axis.2,color=Group,fill=Group))+
      geom_point(size=3) +
      xlab(paste0(pc1)) +
      ylab(paste0(pc2)) +
         ggtitle(paste0("PCoA (", d_method,") --- ",plot_title)) +
         stat_ellipse(geom = "polygon",type = "norm",level = 0.95,size=1.2,alpha=0.1)  +
         theme_bw()+ #stat_ellipse(type = "t")
         annotate("text",x = px , y = py, label = paste0("p = ",p_val),size=5,fontface = 'italic') +
         theme_classic()+
         theme(axis.text.x=element_text(angle=0,size=13,color="black"),
               axis.text.y=element_text(size=13,color="black"),
               axis.title.x=element_text(size=15),
               axis.title.y=element_text(size=15))

p <- p + scale_color_manual(values=myColor,
                            name  ="Group") +
         scale_fill_manual(values=alpha(myColor),
                            name  ="Group")

pdf(paste0("PCoA-",d_method,"-",plot_title,".pdf"),height=5,width=6)
print(p)
dev.off()

#output PCoA p-value
write.csv(as.matrix(stats$aov.tab),paste0("PCoA-",d_method,"-",plot_title,"-pval.csv"))




