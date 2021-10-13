#Author: Ziyan Lin
#Description: this is the script for plotting Figure 3B
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

phy_taxo <- tax_glom(physeq_norm,"Genus")

set.seed(123)
myColor <- c(palette(rainbow(60)),"violetred1","springgreen4","purple3","yellow1","firebrick1","turquoise3","darkorange2","gray48")


df <- psmelt(phy_taxo)
df$Abundance <- df$Abundance*100 #change to percentage(%)

taxo_filter <- df %>%
           group_by(Genus,Group) %>%  # the grouping variable
           summarise(mean_Abun = mean(Abundance),
           n_Abun = n()) %>%
           filter(mean_Abun < 0.1) %>% #obtain phlyums <1%
           summarise(num = sum(n_Abun), .groups = 'drop') %>%
           filter(num >= nrow(sample_data(phy_taxo))) #keep phylum <1% in all conditions

df$Genus[which(df$Genus %in% taxo_filter$Genus)] <- "Others" #assign <1% Taxo to others

df.stats <- df %>%
            group_by(Genus) %>% summarise(mean_Abun = mean(Abundance))

taxo_order <- df.stats$Genus[order(df.stats$mean_Abun,decreasing=T)]
taxo_order #move others to bottom if not
taxo_name <- c(gsub("g__","",taxo_order[-length(taxo_order)]),"Others") #formatted name for Taxo
names(taxo_name) <- taxo_order

#calculate standard error
df_summary <- df %>%
              group_by(Group,Genus) %>%  # the grouping variable
              summarise(mean_Abun = mean(Abundance),# calculates the mean of each group
              sd_Abun = sd(Abundance), # calculates the standard deviation of each group
              n_Abun = n(),  # calculates the sample size per group
              SE_Abun = sd(Abundance)/sqrt(n()))# calculates the standard error of each group

df_summary$Genus <- factor(df_summary$Genus,levels=taxo_order,labels=taxo_name) # set facet order

p<-NULL
p <- ggplot(df_summary, aes(fill=Genus, x=Group, y=mean_Abun)) +
              facet_wrap(~Genus,scales = "free") +
geom_bar(aes(color=Genus,fill=Genus),stat="identity",color="black",width=0.55) +
              ylab("Relative Abundance(%)") +
              xlab("") +
              geom_errorbar(aes(ymin=mean_Abun,ymax=mean_Abun+SE_Abun),width=.2) +#one side
              theme_classic()+
              theme(axis.text.x=element_text(size=12,color="black"),
                    axis.text.y=element_text(color="black"))+
              theme(axis.title.y=element_text(size=18,color="black"))+
              theme(strip.background = element_rect(fill=NULL,linetype="blank"), #facet background
                    strip.text.x = element_text(size=8,face = "bold.italic")) #facet font size

myColor.plot <- c(myColor[1:length(taxo_order)-1],myColor[length(myColor)]) #select colors

p <- p + scale_fill_manual(values=myColor.plot,
                           breaks=taxo_name)

pdf("Barplot-Genus-Groups.pdf",height=18,width=22)
print(p)
dev.off()

#output data in csv format
write.csv(as.data.frame(df_summary),"Barplot-Genus-Groups-data.csv")

#calculate p-value using Kruskal-Wallis
pval_tbl <- compare_means(formula=Abundance ~ Group,data=df,method="kruskal.test",group.by="Genus")

#output p-value
write.csv(as.data.frame(pval_tbl),"Barplot-Genus-Groups-pval.csv")

#calculate p-value using wilcox
pval_tbl_pairwise <- compare_means(formula=Abundance ~ Group,data=df,method="wilcox.test",group.by="Genus")

#output p-value(pairwise)
write.csv(as.data.frame(pval_tbl_pairwise),"Barplot-Genus-Groups-pval-pairwise.csv")
