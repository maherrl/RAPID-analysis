#######################################################
## Figures for RAPID manuscript
#######################################################

library("ggplot2")
library("phyloseq")
library("microbiome")
library("ggthemes")
library(ggnetwork)
library(cowplot)

# Figure 1 - Taxa plot
#######################################################
rm(list=ls())
# load the rarefied data table with all species
load("~/Box Sync/RAPID/RAPID-analysis/data/qd_881.RData")
qd
# transform to relative abundance
qd_rel <- transform(qd, "compositional")
# melt the data at the Genus level
qd_rel_genus_melt <- qd_rel %>%
  tax_glom(taxrank = "Genus") %>%
  psmelt()
head(qd_rel_genus_melt)

# Get genera with mean realtive abundance >0.01 across all samples 
genus_sum <- qd_rel_genus_melt %>% group_by(Genus) %>% dplyr::summarise(Aver = mean(Abundance))
genus_sub <- genus_sum[which(genus_sum$Aver > 0.01),]
names <- genus_sub$Genus
# Replace genera with <0.01 abundance with "NA"
qd_rel_genus_melt$genus <- qd_rel_genus_melt$Genus

qd_rel_genus_melt$genus[qd_rel_genus_melt$genus != "g__Endozoicomonas" & 
                          qd_rel_genus_melt$genus != "g__Vibrio" &
                          qd_rel_genus_melt$genus != "g__Acinetobacter" &
                          qd_rel_genus_melt$genus != "g__Candidatus Amoebophilus" &
                          qd_rel_genus_melt$genus != "g__Corynebacterium 1" &
                          qd_rel_genus_melt$genus != "g__Halobacteriovorax" &
                          qd_rel_genus_melt$genus != "g__Halomonas" &
                          qd_rel_genus_melt$genus != "g__Paenibacillus" &
                          qd_rel_genus_melt$genus != "g__Pseudomonas" &
                          qd_rel_genus_melt$genus != "g__Spiroplasma" &
                          qd_rel_genus_melt$genus != "g__Staphylococcus" &
                          qd_rel_genus_melt$genus != "g__Streptococcus" &
                          qd_rel_genus_melt$genus != "g__Tenacibaculum"] <- NA



bar_species = ggplot(qd_rel_genus_melt, aes(x = treatment, y=Abundance)) + 
  geom_bar(stat="identity", position="fill", aes(fill = genus)) +
#  scale_fill_colorblind() +
  scale_fill_manual(values=c("#56B4E9","#CBD588","#5F7FC7", "orange","#DA5724","#CD9BCD",
                             "gray80", "#AD6F3B", "#673770","#D14285", "#652926","#8569D5",
                             "#5E738F","#D1A33D", "#8A7C64","lightsalmon","aquamarine4",
                             "lightblue4", "lightpink", "ivory4","royalblue4", "darkorchid",
                             "palevioletred1", "#56B4E9","#CBD588","yellow2","#5F7FC7", "orange","#DA5724",
                             "#CD9BCD", "gray80",
                             "#AD6F3B", "#673770","#D14285", "#652926","#8569D5", "#5E738F",
                             "#56B4E9","#CBD588","#5F7FC7", "orange","#DA5724","#CD9BCD")) +
  facet_grid(species~time) +
  theme_bw() +
  ylab("Relative Abundance")
bar_species


# Figure 2
######################################################
rm(list=ls())
alphadiv <- read.csv(file = "./data/alphadiv.csv")

breakss <- c("T0","T1","T2","T3","T6")
labelss <- c("Jan 16","March 16", "May 16", "July 16", "Jan 17")

A <- ggplot(alphadiv, aes(x=time, y=richness)) +
  geom_boxplot(outlier.shape = NA, color = "gray35") +
  geom_point(aes(color = species), position = position_jitter(width = .2, height = 0)) +
  scale_x_discrete(breaks=breakss, labels=labelss) +
  ylab("Chao1 Index") + xlab("Time") +
  scale_colour_colorblind() + theme_bw()
A
B <- ggplot(alphadiv, aes(x=time, y=evenness)) +
  geom_boxplot(outlier.shape = NA, color = "gray35") +
  geom_point(aes(color = species), position = position_jitter(width = .2, height = 0)) +
  scale_x_discrete(breaks=breakss, labels=labelss) +
  ylab("Simpson's Diversity") + xlab("Time") +
  scale_colour_colorblind() + theme_bw()
B
C <- ggplot(alphadiv, aes(x=time, y=faithPD)) +
  geom_boxplot(outlier.shape = NA, color = "gray35") +
  geom_point(aes(color = species), position = position_jitter(width = .2, height = 0)) +
  scale_x_discrete(breaks=breakss, labels=labelss) +
  ylab("Faith's Phylogenetic Diversity") + xlab("Time") + 
  scale_colour_colorblind() + theme_bw()
C

# Figure 3
###############################################################
rm(list=ls())
load("~/Box Sync/RAPID/RAPID-analysis/data/qd_881.RData")

# transform the phyloseq object with log
otus_log <- as(otu_table(qd), "matrix")
otus_log <- decostand(otus_log, method = "log")
OTU_log <- otu_table(otus_log, taxa_are_rows = TRUE)
otu_table(qd) <- OTU_log

qd_wu <- phyloseq::distance(qd, method = "wunifrac")
nmds_wu <- metaMDSiter(qd_wu, k=2, trymax = 10000, maxit = 10000, autotransform=FALSE)

a <- plot_ordination(qd, nmds_wu, color = "time") + geom_point(size = 2) + theme_classic() +
  theme(legend.position = "left") + scale_colour_colorblind()

# Extract axis coordinates
mds1 <- as.data.frame(nmds_wu[["points"]])
mds1 <- cbind(mds1, sample_data(qd))

p <- position_dodge(0.8)
b <- ggplot(mds1, aes(time, MDS1)) + 
  geom_boxplot(aes(color = time, fill = time), outlier.colour = NULL, position = p) + 
  stat_summary(geom = "crossbar", width = 0.7, fatten=0, color="white", position = p, 
               fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))}) + 
  theme_classic() + coord_flip() +
  theme(legend.position = "left") + 
  scale_colour_colorblind() +
  scale_fill_colorblind() +
  scale_x_discrete(limits = rev(levels(mds1$time)))
b  
plot_grid(a,b, nrow = 2, rel_heights = c(3,1))


# Figure 4
##############################################################
rm(list = ls())

library("ggplot2")
library("phyloseq")
library("microbiome")
library("ggthemes")
library(ggnetwork)
library(cowplot)
library("dplyr")

# load the rarefied data table with all species
load("~/Box Sync/RAPID/RAPID-analysis/data/qd_881.RData")
qd
qd <- subset_samples(qd, species == "ACR")
qd
# transform to relative abundance
qd_rel <- transform(qd, "compositional")
# melt the data at the Genus level
qd_rel_genus_melt <- qd_rel %>%
  tax_glom(taxrank = "Genus") %>%
  psmelt()
head(qd_rel_genus_melt)

# Get genera with mean realtive abundance >0.01 across all samples 
genus_sum <- qd_rel_genus_melt %>% group_by(Genus) %>% dplyr::summarise(Aver = mean(Abundance))
genus_sub <- genus_sum[which(genus_sum$Aver > 0.005),]
names <- genus_sub$Genus
# Replace genera with <0.01 abundance with "NA"
qd_rel_genus_melt$genus <- qd_rel_genus_melt$Genus

qd_rel_genus_melt$genus[qd_rel_genus_melt$genus != "g__Endozoicomonas" & 
                          qd_rel_genus_melt$genus != "g__Vibrio" &
                          qd_rel_genus_melt$genus != "g__Acinetobacter" &
                          qd_rel_genus_melt$genus != "g__Alteromonas" &
                          qd_rel_genus_melt$genus != "g__Corynebacterium 1" &
                          qd_rel_genus_melt$genus != "g__Halobacteriovorax" &
                          qd_rel_genus_melt$genus != "g__Halomonas" &
                          qd_rel_genus_melt$genus != "g__Paenibacillus" &
                          qd_rel_genus_melt$genus != "g__Pseudomonas" &
                          qd_rel_genus_melt$genus != "g__Staphylococcus" &
                          qd_rel_genus_melt$genus != "g__Streptococcus" &
                          qd_rel_genus_melt$genus != "g__Psychrobacter" &
                          qd_rel_genus_melt$genus != "g__Reyranella"] <- NA
head(qd_rel_genus_melt)

breakss <- c("T0","T1","T2","T3","T6")
labelss <- c("Jan 16","March 16", "May 16", "July 16","Jan 17")


bar_species = ggplot(qd_rel_genus_melt, aes(x = treatment, y=Abundance)) + 
  geom_bar(stat="identity", position="fill", aes(fill = genus)) +
  #  scale_fill_colorblind() +
  scale_fill_manual(values=c("#56B4E9","#CBD588","#5F7FC7", "orange","#DA5724","#CD9BCD",
                             "gray80", "#AD6F3B", "#673770","#D14285", "#652926","#8569D5",
                             "#5E738F","#D1A33D", "#8A7C64","lightsalmon","aquamarine4",
                             "lightblue4", "lightpink", "ivory4","royalblue4", "darkorchid",
                             "palevioletred1", "#56B4E9","#CBD588","yellow2","#5F7FC7", "orange","#DA5724",
                             "#CD9BCD", "gray80",
                             "#AD6F3B", "#673770","#D14285", "#652926","#8569D5", "#5E738F",
                             "#56B4E9","#CBD588","#5F7FC7", "orange","#DA5724","#CD9BCD")) +
  facet_grid(~qd_rel_genus_melt$time) +
  theme_bw() +
  ylab("Relative Abundance")
bar_species

# Figure 5
######################################################
rm(list=ls())
alphadiv <- read.csv(file = "./data/alphadiv.csv")
alphadiv <- alphadiv[which(alphadiv$species =="ACR"),]

breakss <- c("T0","T1","T2","T3","T6")
labelss <- c("Jan 16","March 16", "May 16", "July 16", "Jan 17")

A <- ggplot(alphadiv, aes(x=time, y=richness)) +
  geom_boxplot(outlier.shape = NA, color = "gray35") +
  geom_point(aes(color = treatment), position = position_jitter(width = .2, height = 0)) +
  scale_x_discrete(breaks=breakss, labels=labelss) +
  ylab("Chao1 Index") + xlab("Time") +
  scale_colour_colorblind() + theme_bw()
A
B <- ggplot(alphadiv, aes(x=time, y=evenness)) +
  geom_boxplot(outlier.shape = NA, color = "gray35") +
  geom_point(aes(color = treatment), position = position_jitter(width = .2, height = 0)) +
  scale_x_discrete(breaks=breakss, labels=labelss) +
  ylab("Simpson's Diversity") + xlab("Time") +
  scale_colour_colorblind() + theme_bw()
B
C <- ggplot(alphadiv, aes(x=time, y=faithPD)) +
  geom_boxplot(outlier.shape = NA, color = "gray35") +
  geom_point(aes(color = treatment), position = position_jitter(width = .2, height = 0)) +
  scale_x_discrete(breaks=breakss, labels=labelss) +
  ylab("Faith's Phylogenetic Diversity") + xlab("Time") + 
  scale_colour_colorblind() + theme_bw()
C



# Figure 6
###############################################################
rm(list=ls())

library("ggplot2")
library("phyloseq")
library("microbiome")
library("ggthemes")
library(ggnetwork)
library(cowplot)
library("vegan")

load("~/Box Sync/RAPID/RAPID-analysis/data/qd_881.RData")
qd <- subset_samples(qd, species == "ACR")

# transform the phyloseq object with log
otus_log <- as(otu_table(qd), "matrix")
otus_log <- decostand(otus_log, method = "log")
OTU_log <- otu_table(otus_log, taxa_are_rows = TRUE)
otu_table(qd) <- OTU_log

qd_wu <- phyloseq::distance(qd, method = "wunifrac")
nmds_wu <- metaMDSiter(qd_wu, k=2, trymax = 10000, maxit = 10000, autotransform=FALSE)

a <- plot_ordination(qd, nmds_wu, color = "time") + geom_point(size = 2) + theme_classic() +
  theme(legend.position = "left") + scale_colour_colorblind()

# Extract axis coordinates
mds1 <- as.data.frame(nmds_wu[["points"]])
mds1 <- cbind(mds1, sample_data(qd))

# extract distance to centroid
sampledf <- data.frame(sample_data(qd))
disp <- betadisper(qd_wu, sampledf$time, bias.adjust = TRUE)
dispd <- as.data.frame(disp$distances)
dispd <- cbind(dispd, sample_data(qd))
colnames(dispd)[1] <- "distance"


p <- position_dodge(0.8)
b <- ggplot(mds1, aes(time, MDS1)) + 
  geom_boxplot(aes(color = time, fill = time), outlier.colour = NULL, position = p) + 
  stat_summary(geom = "crossbar", width = 0.7, fatten=0, color="white", position = p, 
               fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))}) + 
  theme_classic() + coord_flip() +
  theme(legend.position = "left") + 
  scale_colour_colorblind() +
  scale_fill_colorblind() +
  scale_x_discrete(limits = rev(levels(mds1$time)))
b  

p <- position_dodge(0.8)
c <- ggplot(dispd, aes(time, distance)) + 
  geom_boxplot(aes(color = time, fill = time), outlier.colour = NULL, position = p) + 
  stat_summary(geom = "crossbar", width = 0.7, fatten=0, color="white", position = p, 
               fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))}) + 
  theme_classic() + 
  #  stat_summary(aes(color = time), geom = 'text', label = sigwudisp, fun.y = max, vjust = -1, size = 3) +
  theme(legend.position = "none") + labs(y = "Distance to Centroid")+ 
  scale_colour_colorblind() +
  scale_fill_colorblind()
c
plot_grid(a,c,b, nrow = 2, rel_heights = c(3,1), rel_widths = c(3,1))


