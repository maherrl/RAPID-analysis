###################################################
## Code for exploratory figures for early RAPID ##
## analysis
## Created by Rebecca Maher
## Created on Jan 24, 2019
## Edited on
##################################################

# load librarys
library(cowplot)
library(ggplot2)
library(ggnetwork)
library(phyloseq)
library(ape)
library(reshape)


#------------------------------------------------------------------
## Beta diversity plots

# Making beta diversity figure. Uses commands from beta_div analysis.R to get distances
# NMDS ordination
nmds_wu <- metaMDSiter(qd_wu, k=2, trymax = 1000, maxit = 1000, autotransform=FALSE)
nmds_un <- metaMDSiter(qd_un, k=2, trymax = 1000, maxit = 1000, autotransform=FALSE)

# plot
a <- plot_ordination(qd, nmds_wu, color = "time") + geom_point(size = 4) + theme_classic() +
  theme(legend.position = "left")
  scale_color_brewer(palette = "Set2") 

# Log transformation
otus_log <- as(otu_table(qd), "matrix")
otus_log <- decostand(otus_log, method = "log")
OTU_log <- otu_table(otus_log, taxa_are_rows = TRUE)
otu_table(qd) <- OTU_log

# Max transformation
otus_max <- as(otu_table(qd), "matrix")
otus_max <- decostand(otus_max, method = "max", MARGIN = 1)
OTU_max <- otu_table(otus_max, taxa_are_rows = TRUE)
otu_table(qd) <- OTU_max

# Extract axis coordinates
mds1 <- as.data.frame(nmds_wu[["points"]])
mds1 <- cbind(mds1, sample_data(qd))

# extract distance to centroid
disp <- betadisper(qd_wu, sampledf$time, bias.adjust = TRUE)
dispd <- as.data.frame(disp$distances)
dispd <- cbind(dispd, sample_data(qd))
colnames(dispd)[1] <- "distance"

# Plot boxplot for underneath ordination
sigwuman <- c("A","A","B","C","D") # ACR-only, max transformed, WU
sigwudisp <- c("A","B","B","AB","B") # ACR-only, max transformed, WU

p <- position_dodge(0.8)
b <- ggplot(mds1, aes(time, MDS1)) + 
  geom_boxplot(aes(color = time, fill = time), outlier.colour = NULL, position = p) + 
  stat_summary(geom = "crossbar", width = 0.7, fatten=0, color="white", position = p, 
               fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))}) + 
  theme_classic() + coord_flip() +
  stat_summary(aes(color = time), geom = 'text', label = sigwuman, fun.y = max, hjust = -1, size = 3) +
  theme(legend.position = "left")
b  

p <- position_dodge(0.8)
c <- ggplot(dispd, aes(time, distance)) + 
  geom_boxplot(aes(color = time, fill = time), outlier.colour = NULL, position = p) + 
  stat_summary(geom = "crossbar", width = 0.7, fatten=0, color="white", position = p, 
               fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))}) + 
  theme_classic() + stat_summary(aes(color = time), geom = 'text', label = sigwudisp, fun.y = max, vjust = -1, size = 3) +
  theme(legend.position = "none") + labs(y = "Distance to Centroid")
c

plot_grid(a,c,b, nrow = 2, rel_heights = c(3,1), rel_widths = c(3,1))

#--------------------------------------------------------------------------------------------
# Making ANCOM plot
# open unrarefied qd object
qd <- prune_samples(sample_sums(qd) > 1000, qd)
qd <- subset_samples(qd, species =="ACR")
qd = transform_sample_counts(qd, function(x) x/sum(x)) # turn into relative abundances
qd <- tax_glom(qd, taxrank = rank_names(qd)[5], NArm = TRUE)

### Reextract your relative abundance data as community structure and taxonomic annotation
comm = as.data.frame(as(object = phyloseq::otu_table(qd), Class = "matrix"))
tax = as.data.frame(as(object = phyloseq::tax_table(qd), Class = "matrix"))
meta = as.data.frame(as(object = phyloseq::sample_data(qd), Class = "matrix"))
data = cbind(tax,comm)
 
# subset data
ancom <- c("Corynebacteriaceae","Paenibacillaceae","Staphylococcaceae","Rhodospirillaceae","Sphingomonadaceae",
           "Alteromonadaceae","Endozoicimonaceae","Moraxellaceae","Pseudomonadaceae","Vibrionaceae","Xanthomonadaceae")
data_melt = melt(data[,c(5,8:76)], id = c("Family"))  # Melt data by OTU and all samples
data_melt$Family <- gsub("f__", "", data_melt$Family) # get rid of the f__ before each family name
data_melt2 <- data_melt[data_melt$Family %in% ancom,] # subset for only rows with specified Family name

meta_t <- meta[,c(1,9:10)]
colnames(meta_t)[1] <- "variable"
temp <- merge(data_melt2, meta_t)


# plot
og_names <- c("T0", "T1", "T2", "T3", "T4")
new_names <- c("Jan 16", "March 16", "May 16", "July 16", "Jan 17")
temp$time <- factor(temp$time, labels = c("Jan 16", "March 16", "May 16", "July 16", "Jan 17"))
colorblind_pallette = c("#999999", "#E69F00", "#56B4E9", "white", "#009E73", "#660066", "#FFFF00", "#0072B2", "#CC3300", "#CC79A7","#000000")
p <- ggplot(temp,aes(temp$variable,temp$Family)) + geom_point(aes(size = temp$value, fill = temp$Family), shape = 21) + 
        facet_grid(cols = vars(temp$time), rows = vars(temp$nutrient)) +
  theme_facet() +
  scale_fill_manual(values = colorblind_pallette) +
  scale_x_discrete(breaks = og_names, labels = new_names) +
  guides(fill = guide_legend(override.aes = list(size=4))) + 
  labs(fill = "Family", size = "Relative abundance")
p


#--------------------------------------------------------------
# Alpha plots of Chao1 index by time and treatment
# Colorfull plots
og_names <- c("T0", "T1", "T2", "T3", "T4")
new_names <- c("Jan 16", "March 16", "May 16", "July 16", "Jan 17")


p <- position_dodge(0.8)
sigsim <- c("B","AB","AB","AB","A")
a <- ggplot(alphadiv, aes(nutrient, richness)) + 
  geom_boxplot(aes(color = nutrient, fill = nutrient), outlier.colour = NULL, position = p) + 
  stat_summary(geom = "crossbar", width = 0.7, fatten=0, color="white", position = p, 
               fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))}) + 
  theme_classic() + 
  #  stat_summary(aes(color = time_count), geom = 'text', label = sigsim, fun.y = max, vjust = -1, size = 3) +
  theme(legend.position = "none") + labs(y = "Chao1 Index", x = "Nutrient") +
  #  scale_x_discrete(breaks = og_names, labels = new_names) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  
# geom_point(aes(y = mmm*30.5), color = "red") +
#scale_y_continuous(sec.axis = sec_axis(~./30.5, name = "Mean Monthly Max Temp"))


sigchao <- c("AB","B","C","BC","A")

ggplot(alphadiv, aes(time, richness)) + 
  geom_boxplot(aes(color = time, fill = time), outlier.colour = NULL, position = p) + 
  stat_summary(geom = "crossbar", width = 0.7, fatten=0, color="white", position = p, 
               fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))}) + 
  theme_classic() + 
  stat_summary(aes(color = time), geom = 'text', label = sigchao, fun.y = max, vjust = -1, size = 3) +
  theme(legend.position = "none") + labs(y = "Chao1") +  
  scale_x_discrete(breaks = og_names, labels = new_names) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
  ylab("Chao1") + xlab("Month")
