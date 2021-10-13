#Author: Ziyan Lin
#Description: this is the script for plotting Figure 2A
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
#Figure 2A
############################################
phy_alp  <- prune_taxa(taxa_sums(physeq) > 0, physeq) #prune OTUs that are not present in any of the samples
p <- NULL
p <- plot_richness(phy_alp,measures=c("Shannon","Observed","Chao1"))#color="Group"
df1 <- subset(p$data,select=c("SampleID","Visit","Group","PatientID","variable","value","se"))

#Calculate PD
phy_alp_otu <- as.data.frame(phy_alp@otu_table)
phy_alp_tree <- phy_alp@phy_tree# it is a rooted tree
df.pd.temp <- pd(t(phy_alp_otu), phy_alp_tree,include.root=T) # t(ou_table) transposes the table for use in picante

#format PD result
df.pd <- subset(phy_alp@sam_data,select=c("SampleID","Visit","Group","PatientID"))
df.pd$variable <- "PD"
df.pd$value <- df.pd.temp$PD
df.pd$se <- "NA"

#combine all method
alp_data <- rbind(df1,df.pd)
alp_data$variable <- factor(alp_data$variable,levels=c("Observed","Shannon","Chao1","PD"))#customize facet order

alp_data <- alp_data[order(alp_data$variable,alp_data$PatientID),]

unique(alp_data$variable) #
#Observed Shannon Chao1 PD

AlphaDodge <- function(variable){
    myColor <- c("black","gray")
    myColor2 <- c("firebrick3","dodgerblue2","goldenrod2")
    set.seed(123)
    variable <- variable
    subset_df <- alp_data[alp_data$variable==variable,]
    #genere box plot
    p<-NULL
    p <- ggplot(subset_df,aes(x=Group, y=value,color=Visit,fill=Group))+
        geom_boxplot(position=position_dodge(width=0.9),outlier.shape = NA) +
        geom_point(position = position_jitterdodge(dodge.width=0.9)) +
        theme_classic() +
        ggtitle(variable) +
        xlab("")+
        ylab("Alpha Diversity Measure")+
        theme(plot.title = element_text(face = "bold.italic",hjust = 0.5,size=12),
            axis.text.x=element_text(angle=0,size=12,color="black"),
            axis.text.y=element_text(size=10,color="black"))

        p <- p + scale_color_manual(values=myColor,
                            name="Visit") +
                scale_fill_manual(values=myColor2,
                                    name="Group")

        stat.test <- subset_df %>%
                group_by(Group) %>%
                wilcox_test(value ~ Visit,paired = TRUE) %>%
                add_significance("p")

        # Add p-values onto the box plots
        stat.test <- stat.test %>%
                add_xy_position(x = "Group", dodge = 0.8)

        p <- p + stat_pvalue_manual(stat.test,label = "p.signif",hide.ns = TRUE)#tip.length = 0

       return(p)
       }

A <- AlphaDodge("Observed")
B <- AlphaDodge("Shannon")
C <- AlphaDodge("Chao1")
D <- AlphaDodge("PD")

pdf(paste0("Figure2A-Alp-All.pdf"),height=6,width=6)
grid.arrange(A,B,C,D)
dev.off()
