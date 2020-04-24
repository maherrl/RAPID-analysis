#######################################################
## Figures for RAPID manuscript
#######################################################

library("ggplot2")
library("phyloseq")
library("microbiome")
library("ggthemes")
library(ggnetwork)
library(cowplot)
library("vegan")

# Figure 1 - Taxa plot
#######################################################
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

supp.labs <- c("Jan 16","March 16", "May 16", "July 16", "Jan 17")
names(supp.labs) <- c("T0","T1","T2","T3","T6")

bar_species = ggplot(qd_rel_genus_melt, aes(x = treatment, y=Abundance)) + 
  geom_bar(stat="identity", position="fill", aes(fill = genus)) +
  scale_fill_manual(values=c("#56B4E9","#CBD588","#5F7FC7", "orange","#DA5724","#CD9BCD",
                             "gray80", "#AD6F3B", "#673770","#D14285", "#652926","#8569D5",
                             "#5E738F","#D1A33D", "#8A7C64","lightsalmon","aquamarine4",
                             "lightblue4", "lightpink", "ivory4","royalblue4", "darkorchid",
                             "palevioletred1", "#56B4E9","#CBD588","yellow2","#5F7FC7", "orange","#DA5724",
                             "#CD9BCD", "gray80",
                             "#AD6F3B", "#673770","#D14285", "#652926","#8569D5", "#5E738F",
                             "#56B4E9","#CBD588","#5F7FC7", "orange","#DA5724","#CD9BCD")) +
  facet_grid(species~time, labeller = labeller(time = supp.labs)) +
  theme_bw() +
  ylab("Relative Abundance")
bar_species


# Figure 2 Alpha diversity
######################################################
rm(list=ls())
alphadiv <- read.csv(file = "./data/alphadiv.csv")

breakss <- c("T0","T1","T2","T3","T6")
labelss <- c("Jan 16","March 16", "May 16", "July 16", "Jan 17")

A <- ggplot(alphadiv, aes(x=time, y=richness, color = time)) +
  geom_boxplot(outlier.shape = NA, color = "gray35") +
  facet_grid(~species) + 
  geom_point(aes(color = time), position = position_jitter(width = .2, height = 0)) +
  scale_x_discrete(breaks=breakss, labels=labelss) +
  ylab("Chao1 Index") +
  scale_colour_colorblind() + 
  stat_summary(geom = 'text', 
               label = c("ac","ab","b","ab","c", "a","a","b","a","b","ab","b","ab","a"), 
               fun.y = max, vjust = -1, size = 4, col = "black") +
  scale_y_continuous(expand = expand_scale(mult = c(0, .1))) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 35, hjust = 1))

B <- ggplot(alphadiv, aes(x=time, y=evenness, color = time)) +
  geom_boxplot(outlier.shape = NA, color = "gray35") +
  facet_grid(~species) + 
  geom_point(aes(color = time), position = position_jitter(width = .2, height = 0)) +
  scale_x_discrete(breaks=breakss, labels=labelss) +
  ylab("Simpson's Diversity Index") +
  scale_y_continuous(expand = expand_scale(mult = c(0, .1))) +
  scale_colour_colorblind() + 
  stat_summary(geom = 'text', 
               label = c("ab","b","b","b","a","abc", "b","c","b","d","ab","ab","a","b"), 
               fun.y = max, vjust = -1, size = 4, col = "black") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 35, hjust = 1))

C <- ggplot(alphadiv, aes(x=time, y=faithPD, color = time)) +
  geom_boxplot(outlier.shape = NA, color = "gray35") +
  facet_grid(~species) + 
  geom_point(aes(color = time), position = position_jitter(width = .2, height = 0)) +
  scale_x_discrete(breaks=breakss, labels=labelss) +
  ylab("Faith's Phylogenetic Diversity") +
  scale_y_continuous(expand = expand_scale(mult = c(0, .1))) +
  scale_colour_colorblind() +
  stat_summary(geom = 'text', 
               label = c("a","a","a","a","b","abc","b", "c","b","d","ab","a","ab","b"), 
               fun.y = max, vjust = -1, size = 4, col = "black") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 35, hjust = 1))

plot_grid(A,B,C, nrow = 3, labels = c("A","B","C"))



# Figure 3 Ordinations
###############################################################
rm(list=ls())
load("~/Box Sync/RAPID/RAPID-analysis/data/qd_881.RData")

# transform the phyloseq object with log
otus_log <- as(otu_table(qd), "matrix")
otus_log <- decostand(otus_log, method = "log")
OTU_log <- otu_table(otus_log, taxa_are_rows = TRUE)
otu_table(qd) <- OTU_log

# all corals by species
qd_wu <- phyloseq::distance(qd, method = "wunifrac")
nmds_wu <- metaMDSiter(qd_wu, k=2, trymax = 10000, maxit = 10000, autotransform=FALSE)

# A
a <- plot_ordination(qd, nmds_wu, color = "species") + 
  geom_point(size = 2) + 
  theme_classic() +
  theme(legend.position = "left", legend.title = element_blank()) + 
  scale_colour_manual(values = c("#0072B2", "#D55E00", "#CC79A7")) +
  ggtitle("All corals") +
  annotate(geom = "text", label = "stress = 0.17", x = 1.5, y = -3.2)

# B
# acr by time
qd_acr <- subset_samples(qd, species == "ACR")
qd_acr_wu <- phyloseq::distance(qd_acr, method = "wunifrac")
nmds_acr_wu <- metaMDSiter(qd_acr_wu, k=2, trymax = 10000, maxit = 10000, autotransform=FALSE)

b <- plot_ordination(qd_acr, nmds_acr_wu, color = "time") + 
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442"), 
                    name=element_blank(),
                    breaks=c("T0", "T1", "T2", "T3", "T6"),
                    labels=c("Jan 16", "March 16", "May 16", "July 16", "Jan 17")) +
  theme_classic() +
  theme(legend.position = "left") +
  ggtitle("Acropora") +
  annotate(geom = "text", label = "stress = 0.16", x = 1.1, y = -1.9)

# C
qd_poc <- subset_samples(qd, species == "POC")
qd_poc_wu <- phyloseq::distance(qd_poc, method = "wunifrac")
nmds_poc_wu <- metaMDSiter(qd_poc_wu, k=2, trymax = 10000, maxit = 10000, autotransform=FALSE)

c <- plot_ordination(qd, nmds_poc_wu, color = "time") + 
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442"), 
                     name=element_blank(),
                     breaks=c("T0", "T1", "T2", "T3", "T6"),
                     labels=c("Jan 16", "March 16", "May 16", "July 16", "Jan 17")) +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("Pocillopora") +
  annotate(geom = "text", label = "stress = 0.14", x = 1.2, y = -1.9)

# D
qd_por <- subset_samples(qd, species == "POR")
qd_por_wu <- phyloseq::distance(qd_por, method = "wunifrac")
nmds_por_wu <- metaMDSiter(qd_por_wu, k=2, trymax = 10000, maxit = 10000, autotransform=FALSE)

d <- plot_ordination(qd, nmds_por_wu, color = "time") + 
  geom_point(size = 2) + 
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), 
                     name=element_blank(),
                     breaks=c("T1", "T2", "T3", "T6"),
                     labels=c("March 16", "May 16", "July 16", "Jan 17")) +
  theme_classic() +
  theme(legend.position = "none")+
  ggtitle("Porites") +
  annotate(geom = "text", label = "stress = 0.17", x = .8, y = -1.5)


plot_grid(a,c,b,d, nrow = 2, labels = c("A","C","B","D"), rel_widths = c(1,.8))



# Figure 4 Dispersion
###############################################################
rm(list=ls())
load("~/Box Sync/RAPID/RAPID-analysis/data/qd_881.RData")

# transform the phyloseq object with log
otus_log <- as(otu_table(qd), "matrix")
otus_log <- decostand(otus_log, method = "log")
OTU_log <- otu_table(otus_log, taxa_are_rows = TRUE)
otu_table(qd) <- OTU_log

# all corals by species
qd_wu <- phyloseq::distance(qd, method = "wunifrac")

# A
# extract distance to centroid
sampledf <- data.frame(sample_data(qd))
disp <- betadisper(qd_wu, sampledf$species, bias.adjust = TRUE)
disp <- as.data.frame(disp$distances)
disp <- cbind(disp, sample_data(qd))
colnames(disp)[1] <- "distance"

p <- position_dodge(0.8)
a <- ggplot(disp, aes(species, distance)) + 
  geom_boxplot(outlier.shape = NA, color = "gray35") +
  geom_point(aes(color = species), position = position_jitter(width = .2, height = 0)) +
  theme_classic() + 
  stat_summary(geom = 'text', 
               label = c("a","b","b"), 
               fun.y = max, vjust = -1, size = 4, col = "black") +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 35, hjust = 1)) +
#  scale_y_continuous(expand = expand_scale(mult = c(0, .1))) +
  labs(y = "Distance to Centroid") + 
  scale_colour_manual(values = c("#0072B2", "#D55E00", "#CC79A7")) +
  scale_fill_colorblind() +
  ggtitle("All corals") +
  ylim(0,0.4)

# B

qd_acr <- subset_samples(qd, species == "ACR")
qd_acr_wu <- phyloseq::distance(qd_acr, method = "wunifrac")

# extract distance to centroid
sampledf <- data.frame(sample_data(qd_acr))
dispa <- betadisper(qd_acr_wu, sampledf$time, bias.adjust = TRUE)
dispa <- as.data.frame(dispa$distances)
dispa <- cbind(dispa, sample_data(qd_acr))
colnames(dispa)[1] <- "distance"

p <- position_dodge(0.8)
b <- ggplot(dispa, aes(time, distance)) + 
  geom_boxplot(outlier.shape = NA, color = "gray35") +
  geom_point(aes(color = time), position = position_jitter(width = .2, height = 0)) +
  theme_classic() + 
  stat_summary(geom = 'text', 
               label = c("a","ab","b","ab","c"), 
               fun.y = max, vjust = -1, size = 4, col = "black") +
  theme(legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 35, hjust = 1),
        axis.title.y = element_blank()) + 
  labs(y = "Distance to Centroid") + 
#  scale_y_continuous(expand = expand_scale(mult = c(0, .1))) +
  scale_colour_colorblind() +
  scale_fill_colorblind() +
  ggtitle("Acropora") +
  ylim(0,0.4) +
  scale_x_discrete( breaks=c("T0", "T1", "T2", "T3", "T6"),
                    labels=c("Jan 16", "March 16", "May 16", "July 16", "Jan 17"))

# C

qd_poc <- subset_samples(qd, species == "POC")
qd_poc_wu <- phyloseq::distance(qd_poc, method = "wunifrac")

# extract distance to centroid
sampledf <- data.frame(sample_data(qd_poc))
dispc <- betadisper(qd_poc_wu, sampledf$time, bias.adjust = TRUE)
dispc <- as.data.frame(dispc$distances)
dispc <- cbind(dispc, sample_data(qd_poc))
colnames(dispc)[1] <- "distance"

p <- position_dodge(0.8)
c <- ggplot(dispc, aes(time, distance)) + 
  geom_boxplot(outlier.shape = NA, color = "gray35") +
  geom_point(aes(color = time), position = position_jitter(width = .2, height = 0)) +
  theme_classic() + 
 # stat_summary(aes(color = time), geom = 'text', label = c("a","ab","b","ab","c"), fun.y = max, vjust = -1, size = 3) +
  theme(legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 35, hjust = 1),
        axis.title.y = element_blank()) + 
  labs(y = "Distance to Centroid") + 
  #scale_y_continuous(expand = expand_scale(mult = c(0, .1))) +
  scale_colour_colorblind() +
  scale_fill_colorblind() +
  ggtitle("Pocillopora") +
  ylim(0,0.4) +
  scale_x_discrete( breaks=c("T0", "T1", "T2", "T3", "T6"),
                    labels=c("Jan 16", "March 16", "May 16", "July 16", "Jan 17"))

# D

qd_por <- subset_samples(qd, species == "POR")
qd_por_wu <- phyloseq::distance(qd_por, method = "wunifrac")

# extract distance to centroid
sampledf <- data.frame(sample_data(qd_por))
dispr <- betadisper(qd_por_wu, sampledf$time, bias.adjust = TRUE)
dispr <- as.data.frame(dispr$distances)
dispr <- cbind(dispr, sample_data(qd_por))
colnames(dispr)[1] <- "distance"

p <- position_dodge(0.8)
d <- ggplot(dispr, aes(time, distance)) + 
  geom_boxplot(outlier.shape = NA, color = "gray35") +
  geom_point(aes(color = time), position = position_jitter(width = .2, height = 0)) +
  theme_classic() + 
  #stat_summary(aes(color = time), geom = 'text', label = c("a","ab","b","ab","c"), fun.y = max, vjust = -1, size = 3) +
  theme(legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 35, hjust = 1),
        axis.title.y = element_blank()) + 
  labs(y = "Distance to Centroid") + 
  #scale_y_continuous(expand = expand_scale(mult = c(0, .1))) +
  scale_colour_manual(values = c("#E69F00", "#56B4E9", "#009E73","#F0E442")) +
  scale_fill_colorblind() +
  ggtitle("Porites") +
  ylim(0,0.4) +
  scale_x_discrete( breaks=c("T1", "T2", "T3", "T6"),
                    labels=c("March 16", "May 16", "July 16", "Jan 17"))

plot_grid(a,b,c,d, nrow = 1, labels = c("A","B","C","D"), rel_widths = c(4/5, 1,1, 4/5),
          axis = c("b"))
