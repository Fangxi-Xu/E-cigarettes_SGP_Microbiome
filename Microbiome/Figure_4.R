#Author: Ziyan Lin
#Description: this is the script for plotting Figure 4
######################
#library package
library(phyloseq)
library(pheatmap)
library(dplyr)

load("my.physeq.Robj")

############################################
#Figure 4
############################################
#(top 20 genera taxa)
physeq_norm  = transform_sample_counts(physeq, function(x) x / sum(x)) #calculate each proportion(Relative abundance) for each sample
phy_genus <- tax_glom(physeq_norm,"Genus") #add up the otu proportion based on selected taxa (group otu)
top20_taxa <- prune_taxa(names(sort(taxa_sums(phy_genus),TRUE)[1:20]), phy_genus)

df <- psmelt(top20_taxa)
df_summary <- df %>%
              group_by(Family,Genus,Group) %>%
              summarise(mean_Abun = mean(Abundance))
              
df_summary$myName <- paste0(df_summary$Genus," (",df_summary$Family,")")

######################
#formatting the plot data
tbl <- matrix(ncol=length(unique(df_summary$Group)),nrow=length(unique(df_summary$myName)))
colnames(tbl) <- unique(df_summary$Group)
rownames(tbl) <- unique(df_summary$myName)

df_summary <- df_summary[order(df_summary$Group,df_summary$myName),]
tbl <- tbl[order(rownames(tbl)),]

tbl[rownames(tbl),colnames(tbl)] <- df_summary$mean_Abun #assign the Abundance value
#make sure the assign is correct
head(tbl)
head(df_summary)

######################
#generate heatmap
cell_colors = colorRampPalette(c("#043177", "#244B88", "#FAFAFA", "#C62E2E", "#BF0F0F"))(50)

#file name, plot title
plot_title <- "Heatmap"
file_name <- paste0("Heatmap")

#set parameters
show_rowname = TRUE #TRUE/FALSE
fontsize_row = 10 #row font
fontsize_col = 12 #colname font
cluster_rows = TRUE #TRUE/FALSE
cluster_cols = FALSE #TRUE/FALSE
scale = "row" #"none","row","column"

#plot heatmap in pdf
pheatmap(as.matrix(tbl) , color = cell_colors, border_color = NA, scale = scale, cluster_rows =cluster_rows, cluster_cols = cluster_cols, main = plot_title, fontsize_row = fontsize_row, fontsize_col = fontsize_col, show_rownames = show_rowname,filename = paste0(file_name,".pdf"), width =7, height = 6)
