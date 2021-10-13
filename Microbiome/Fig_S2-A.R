#Author: Ziyan Lin
#Description: this is the script for plotting Fig. S2 A
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
#Fig. S2 A
############################################
physeq_norm  = physeq #rarefied data
myColor <- c("firebrick3","dodgerblue2","goldenrod2","tomato","mediumpurple","forestgreen")

plot_title = "All"
d_method = "bray"
px =  0.45
py =  0.45

#calculate distance
iDist <- distance(physeq_norm, method=d_method)
   
#get coordinate
iPCoA <- ordinate(physeq_norm, "PCoA", distance=iDist)

#generate PCoA plot
p <- NULL
p <- plot_ordination(physeq_norm, iPCoA,color="Group",shape="Visit") #shape=
df <- p$data
df$Condition <- paste0(df$Visit,":",df$Group)

#adonis(p-value)
set.seed(123)
stats<- adonis(formula = iDist ~ df$Condition)
p_val<-stats$aov.tab$"Pr(>F)"[1]

pc1 <- gsub("Axis.1  ","PCo1",p$labels$x)
pc2 <- gsub("Axis.2  ","PCo2",p$labels$y)

p <- NULL
p <- ggplot(df,aes(x=Axis.1,y=Axis.2,color=Condition,fill=Condition)) + geom_point(size=3) +
         ggtitle(paste0("PCoA (", d_method,") --- ",plot_title)) +
         xlab(paste0(pc1)) +
         ylab(paste0(pc2)) +
         stat_ellipse(geom = "polygon",type = "norm",level = 0.95,size=1.2,alpha=0.1)  +#stat_ellipse(type = "t")
         annotate("text", x = px, y = py, label = paste0("p = ",p_val),size=5,fontface = 'italic') +
         theme_classic()+
         theme(axis.text.x=element_text(angle=0,size=13,color="black"),
               axis.text.y=element_text(size=13,color="black"),
               axis.title.x=element_text(size=15),
               axis.title.y=element_text(size=15))

p <- p + scale_color_manual(values=myColor,
                            name  ="Condition") +
         scale_fill_manual(values=alpha(myColor),
                            name  ="Condition")

pdf(paste0("PCoA-",d_method,"-All.pdf"),height=4,width=6)
print(p)
dev.off()

#output PCoA p-value
write.csv(as.matrix(stats$aov.tab),paste0("PCoA-",d_method,"-All-pval.csv"))

###########################################################
#Baseline vs Follow-up visits within cohort
myColor <- c("tomato","mediumpurple")
group_list <- c("CS","ES","NS")

pval_pos <- list(c(0.45,0.45),c(0.45,0.45),c(0.45,0.45)) #bray

for (i in 1:length(group_list)){
    group <- group_list[i]
    plot_title <- group
    
    sub_physeq <- subset_samples(physeq_norm,Group==group) #select interested group
    
    d_method = "bray"
    
    #calculate distance
    iDist <- distance(sub_physeq, method=d_method)
    
    #get coordinate
    iPCoA <- ordinate(sub_physeq, "PCoA", distance=iDist)
    
    #adonis(p-value)
    set.seed(123)
    stats<- adonis(formula = iDist ~ sample_data(sub_physeq)$Visit)
    p_val<-stats$aov.tab$"Pr(>F)"[1]
    
    #generate PCoA plot
    p <- NULL
    p <- plot_ordination(sub_physeq, iPCoA,color="Visit") #shape=
    
    pc1 <- gsub("Axis.1  ","PCo1",p$labels$x)
    pc2 <- gsub("Axis.2  ","PCo2",p$labels$y)
    
    df <- p$data
    p <- ggplot(df,aes(x=Axis.1,y=Axis.2,color=Visit,fill=Visit)) +
         geom_point(size=3) +
         xlab(paste0(pc1)) +
         ylab(paste0(pc2)) +
             ggtitle(paste0("PCoA (", d_method,") --- ",plot_title)) +
             stat_ellipse(geom = "polygon",type = "norm",level = 0.95,size=1.2,alpha=0.1) + theme_bw()+ #stat_ellipse(type = "t")
             annotate("text",  x = pval_pos[[i]][1] , y = pval_pos[[i]][2], label = paste0("p = ",p_val),size=5,fontface = 'italic') +
             theme_classic()+
             theme(axis.text.x=element_text(angle=0,size=13,color="black"),
                   axis.text.y=element_text(size=13,color="black"),
                   axis.title.x=element_text(size=15),
                   axis.title.y=element_text(size=15))

    p <- p + scale_color_manual(values=myColor,
                                name  ="Visit") +
             scale_fill_manual(values=alpha(myColor),
                                name  ="Visit")

    pdf(paste0("PCoA-",d_method,"-Visits-",group,".pdf"),height=5,width=6)
    print(p)
    dev.off()
    
    #output PCoA p-value
    write.csv(as.matrix(stats$aov.tab),paste0("PCoA-",d_method,"-Visits-",group,"-pval.csv"))
    
}
