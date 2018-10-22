#################################################################################
## Exploratory graphs for the RAPID analysis
#
# Created by Rebecca Maher
# Created on 10/9/18
# Edited on 
#################################################################################

## clear workspace------------------------
rm(list=ls())

# load libraries
library('ggplot2')
library('phyloseq')

# Plot 1
# Description: boxplot of community beta diversity distances by time.
sps_palette <- c("#56B4E9", "#D55E00", "#CC79A7")
og_names <- c("1", "2", "3", "4", "5")
new_names <- c("Jan 16", "March 16", "May 16", "July 16", "Jan 17")
p1 <- ggplot(data = distances, aes(x=time, y=distbc, color=species)) +
  geom_boxplot(aes(color=species), alpha=0.5) + 
  scale_color_manual(values = sps_palette) +
  scale_x_discrete(breaks = og_names, labels = new_names) +
  xlab("") + ylab("Bray Curtis distance") +
  theme(panel.background = element_rect(fill="white"),
        legend.background = element_rect(fill="white"),
        legend.key = element_rect(fill="white"),
        panel.grid.major = element_line(colour="white"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black"),
        legend.title=element_blank(),
        axis.text = element_text(size =12),
        axis.text.x = element_text(angle = 35, hjust = 1),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 10),
#        legend.position="none"
  )

p1 <- ggplot(data = distances, aes(x=time, y=distbc)) +
  geom_boxplot(alpha=0.5) + 
  scale_color_manual(values = sps_palette) +
  scale_x_discrete(breaks = og_names, labels = new_names) +
  xlab("") + ylab("Bray Curtis distance") +
  theme(panel.background = element_rect(fill="white"),
        legend.background = element_rect(fill="white"),
        legend.key = element_rect(fill="white"),
        panel.grid.major = element_line(colour="white"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black"),
        legend.title=element_blank(),
        axis.text = element_text(size =12),
        axis.text.x = element_text(angle = 35, hjust = 1),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.position="none"
  )
p1


# PCoA plots 

plot_ordination(qd, ord_un, type = "split", color = "Phylum", shape = "time")

plot_heatmap(qd, "MDS", "unifrac", "time", "Family", weighted = TRUE)

# Family table made from summarize_taxa.py in Qiime
fam_table <- read.csv(file = "/Users/Becca/Box Sync/Rapid/family-table/map_nona_L5.csv")
color_palette <- colorRampPalette(c("#e0ecf4", "#9ebcda", "#8856a7"))

ggplot(data = fam_table, aes(x=time, y=species)) +
  geom_tile(aes(fill = f__Endozoicimonaceae), color = "white") +
  scale_fill_gradient(low = "white", high = "mediumpurple3") +
  ggtitle("Endozoicimonaceae") + 
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(), legend.position = "none",
        plot.title = element_text(hjust = 0.5, vjust = 0))

ggplot(data = fam_table, aes(x=time, y=species)) +
  geom_tile(aes(fill = f__.Amoebophilaceae.), color = "white") +
  scale_fill_gradient(low = "white", high = "mediumpurple3") +
  ggtitle("[Amoebophilaceae]") + 
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, vjust = 0))

ggplot(data = fam_table, aes(x=time, y=species)) +
  geom_tile(aes(fill = f__Vibrionaceae), color = "white") +
  scale_fill_gradient(low = "white", high = "mediumpurple3") +
  ggtitle("Vibrionaceae") + 
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, vjust = 0))
