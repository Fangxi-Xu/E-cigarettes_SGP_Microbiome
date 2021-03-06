#Author: Ziyan Lin
#Description: this is the script for plotting Fig.S1 A-B
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

load("my.physeq.Robj")

############################################
#Fig.S1 A-B
############################################
myColor <- c("firebrick3","goldenrod2","dodgerblue2")
group_list <- c("V1","V2")

for (i in 1:length(group_list)){
    group <- group_list[i]
    plot_title <- group
    sub_physeq <- subset_samples(physeq,Visit==group) #select interested group
    phy_alp  <- prune_taxa(taxa_sums(sub_physeq) > 0, sub_physeq) #prune OTUs that are not present in any of the samples
    p <- NULL
    p <- plot_richness(phy_alp,measures=c("Shannon","Observed","Chao1"))
    df1 <- subset(p$data,select=c("SampleID","Visit","Group","variable","value","se"))

    #Calculate PD
    phy_alp_otu <- as.data.frame(phy_alp@otu_table)
    phy_alp_tree <- phy_alp@phy_tree# it is a rooted tree
    df.pd.temp <- pd(t(phy_alp_otu), phy_alp_tree,include.root=T) # t(ou_table) transposes the table for use in picante
    
    #format PD result
    df.pd <- subset(phy_alp@sam_data,select=c("SampleID","Visit","Group"))
    df.pd$variable <- "PD"
    df.pd$value <- df.pd.temp$PD
    df.pd$se <- "NA"

    #combine all method
    alp_data <- rbind(df1,df.pd)
    alp_data$variable <- factor(alp_data$variable,levels=c("Observed","Shannon","Chao1","PD"))#customize facet order

    set.seed(123)
    #genere box plot
    p<-NULL
    p <- ggplot(alp_data,aes(x=Group, y=value))+
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
    lev <- levels(as.factor(alp_data$Group)) # get the variables
    L.pairs <- combn(seq_along(lev), 2, simplify = FALSE, FUN = function(i) lev[i])# make a pairwise list that we want to compare.

    #add p value(significant level)
    #p <- p + stat_compare_means(
    #  comparisons = L.pairs,
    #  hide.ns = T,
    #  label = "p.signif",#"p.format"
    #  method = "wilcox.test",
    #  paired = FALSE
    #  )

    pdf(paste0("Alp-Groups-",plot_title,".pdf"),height=6,width=6)
    print(p)
    dev.off()

    #output alpha diversity info
    write.csv(as.data.frame(alp_data),paste0("Alp-Groups-",plot_title,"-data.csv"))
    
    #output alpha p-value
    pval_tbl <- compare_means(formula=value ~ Group,data=alp_data,method="wilcox.test",group.by="variable")
    write.csv(as.data.frame(pval_tbl),paste0("Alp-Groups-",plot_title,"-pval.csv"))

    }
