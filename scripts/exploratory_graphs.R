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
