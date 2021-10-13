#Author: Ziyan Lin
#Description: this is the script for plotting Figure 3A
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
#Figure 3A
############################################
#normalize per sample (persentage)
physeq_norm  = transform_sample_counts(physeq, function(x) x / sum(x)) #calculate each proportion for each sample

phy_taxo <- tax_glom(physeq_norm,"Class")

set.seed(123)
myColor <- c(palette(rainbow(19)),"violetred1","springgreen4","purple3","yellow1","firebrick1","turquoise3","darkorange2","gray48")

df <- psmelt(phy_taxo)
df$Abundance <- df$Abundance*100 #change to percentage(%)
#calculate standard error

taxo_filter <- df %>%
           group_by(Class,Visit) %>%  # the grouping variable
           summarise(mean_Abun = mean(Abundance),
           n_Abun = n()) %>%
           filter(mean_Abun < 0.1) %>% #obtain phlyums <1%
           summarise(num = sum(n_Abun), .groups = 'drop') %>%
           filter(num >= nrow(sample_data(phy_taxo))) #keep phylum <1% in all conditions

df$Class[which(df$Class %in% taxo_filter$Class)] <- "Others" #assign <1% Taxo to others

df.stats <- df %>%
            group_by(Class) %>% summarise(mean_Abun = mean(Abundance))

taxo_order <- df.stats$Class[order(df.stats$mean_Abun,decreasing=T)]
taxo_order #move others to bottom if not
taxo_name <- c(gsub("c__","",taxo_order[-length(taxo_order)]),"Others") #formatted name for Taxo
names(taxo_name) <- taxo_order

df_summary <- df %>%
              group_by(Group,Class,Visit)  %>% # the grouping variable
              summarise(mean_Abun = mean(Abundance))

df_summary$Class <- factor(df_summary$Class,levels=taxo_order,labels=taxo_name) # set facet order

p <- NULL
p <- ggplot(df_summary, aes(fill=Class, x=Visit, y=mean_Abun)) +
              facet_wrap(~Group,scales = "free") +
geom_bar(aes(color=Class,fill=Class),stat="identity",color="black",width=0.65) +
              ylab("Relative Abundance(%)") +
              xlab("") +
              theme_classic()+
              theme(axis.text.x=element_text(angle=0,size=12,color="black"),
                    axis.text.y=element_text(color="black"))+
              theme(axis.title.y=element_text(size=20,color="black"))+
              theme(strip.background = element_rect(fill=NULL,linetype="blank"), #facet background
                    strip.text.x = element_text(size=15,face="bold")) #facet font size

myColor.plot <- c(myColor[1:length(taxo_order)-1],myColor[length(myColor)]) #select colors

p <- p + scale_fill_manual(values=myColor.plot,
                           breaks=taxo_name)

pdf(paste0("Barplot-Class-Visits-Stack.pdf"),height=6,width=12)
print(p)
dev.off()
