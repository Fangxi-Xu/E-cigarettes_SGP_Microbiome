#Author: Ziyan Lin
#Description: this is the script for plotting Fig. S2 B-C
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
#Fig. S2 B-C
############################################
physeq_norm  = physeq #rarefied data
myColor <- c("firebrick3","goldenrod2","dodgerblue2")
group_list <- c("V1","V2")

pval_pos <- list(c(0.45,0.45),c(0.45,0.45))

for (i in 1:length(group_list)){
    group <- group_list[i]
    plot_title <- group
    
    sub_physeq <- subset_samples(physeq_norm,Visit==group) #select interested group
    
    d_method = "bray"
    
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
    p <- ggplot(df,aes(x=Axis.1,y=Axis.2,color=Group,fill=Group)) +
         geom_point(size=3) +
         xlab(paste0(pc1)) +
         ylab(paste0(pc2)) +
             ggtitle(paste0("PCoA (", d_method,") --- ",plot_title)) +
             stat_ellipse(geom = "polygon",type = "norm",level = 0.95,size=1.2,alpha=0.1) + theme_bw()+ #stat_ellipse(type = "t")
             annotate("text", x = pval_pos[[i]][1] , y = pval_pos[[i]][2], label = paste0("p = ",p_val),size=5,fontface = 'italic') +
             theme_classic()+
             theme(axis.text.x=element_text(angle=0,size=13,color="black"),
                   axis.text.y=element_text(size=13,color="black"),
                   axis.title.x=element_text(size=15),
                   axis.title.y=element_text(size=15))

    p <- p + scale_color_manual(values=myColor,
                                name  ="Group") +
             scale_fill_manual(values=alpha(myColor),
                                name  ="Group")

    pdf(paste0("PCoA-",d_method,"-Groups-",group,".pdf"),height=5,width=6)
    print(p)
    dev.off()
    
    #output PCoA p-value
    write.csv(as.matrix(stats$aov.tab),paste0("PCoA-",d_method,"-Groups-",group,"-pval.csv"))
    
}

############################################
#Pairwise
myVisit <- "V2" #V1/V2
my_physeq <- subset_samples(physeq_norm,Visit==myVisit) #select interested group

group_list <- c("CS","NS") #replace with each pairwise: c("CS","ES"), c("CS","NS"),c("ES","NS")
plot_title = "CS_vs_NS"
d_method = "bray" #"wunifrac"/"unifrac"/"bray"
myColor <- c("firebrick3","goldenrod2","dodgerblue2")
myColor <- myColor[-2]
px =  0.45
py =  0.45

sub_physeq <- subset_samples(my_physeq,Group==group_list[1]|Group==group_list[2]) #select interested group

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

pdf(paste0("PCoA-",d_method,"-",plot_title,"-",myVisit,".pdf"),height=5,width=6)
print(p)
dev.off()

#output PCoA p-value
write.csv(as.matrix(stats$aov.tab),paste0("PCoA-",d_method,"-",plot_title,"-pval-",myVisit,".csv"))


