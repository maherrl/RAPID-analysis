
#################################################################################
# This script is the alpha diversity analysis for the RAPID Project:
# Here, I calculate statistics for community richness, evenness,
# and phylogenetic diversity.
#
# Created by Rebecca Maher
# Created on 10/2/18
# Edited on 10/4/18
#################################################################################

## clear workspace------------------------
rm(list=ls())

# load libraries
library('phyloseq'); packageVersion('phyloseq')
library('vegan'); packageVersion('vegan')
library('picante')
library('cowplot')

# set working directory-------------------
#setwd("~/Box Sync/RAPID-analysis/")

## functions (if any)----------------------
estimate_pd <- function(phylo){
  
  # Error if input is not of class phylo
  if(class(phylo) != "phyloseq"){
    stop("Input file is not of class 'phyloseq'.")
  }
  
  # Error if no class phy_tree
  if(!(.hasSlot(phylo, "phy_tree"))){
    stop("Could not find tree slot in phylo object.")
  }
  
  # Transpose if needed
  # Adapted from phyloseq/vegan import
  OTU <- phyloseq::otu_table(phylo)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  
  # Get matrix version of OTU table
  otutable <- as(OTU, "matrix")
  
  # Get phylogenetic tree from pyloseq object
  tree <- phyloseq::phy_tree(phylo)
  
  # Print status message
  message("Calculating Faiths PD-index...")
  
  # If object is greater than 10mb, then print status message
  if(object.size(otutable) > 10000000){
    message("This is a large object, it may take awhile...")
  }
  
  # Calculate Faith's PD-index
  pdtable <- picante::pd(otutable, tree, include.root = F)
  
  # Return data frame of results
  return(pdtable)
}

## Data Analysis----------------------------

# load data
load(file = "~/Box Sync/RAPID/RAPID-analysis/data/qd_rare.RData")

# Initialize matrices to store richness and evenness estimates
nsamp = nsamples(qd)

richness <- matrix(nrow = nsamp)
row.names(richness) <- sample_names(qd)

evenness <- matrix(nrow =nsamp)
row.names(evenness) <- sample_names(qd)

faithPD <- matrix(nrow = nsamp)
row.names(faithPD) <- sample_names(qd)
##


##
# Options for measures = ("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")

# Calculate richness
rich <- as.numeric(as.matrix(subset(estimate_richness(qd, measures = "Chao1"), select = c(1))))
richness[ ,] <- rich
colnames(richness) [1] <- "richness"

# Calculate evenness
even <- as.numeric(as.matrix(estimate_richness(qd, measures = "Simpson")))
evenness[ ,] <- even
colnames(evenness) [1] <- "evenness"

# Calculate Faith's PD
faith <- as.numeric(as.matrix(subset(estimate_pd(qd), select = c(1))))  # estimate_pd is a function assigned at the end of the script
faithPD[ ,] <- faith
colnames(faithPD) [1] <- "faithPD"

# Included the subset in "rich" and "faith" because the Chao1 and Faith's PD measurement outputs two measures per sample (Chao1 and se.chao1)
# and we only want Chao1, so we select for the first column
##


##
# Combine our estimates for richness, evenness, and faith's PD into one dataframe
alpha <- cbind(richness, evenness[,1][match(rownames(richness), rownames(evenness))],
               faithPD[,1][match(rownames(richness), rownames(faithPD))])
colnames(alpha) [2] <- "evenness"
colnames(alpha) [3] <- "faithPD"

# Add the sample metadata into this dataframe
# Definitely not the cleanest way to do this...
colnames(sample_data(qd))[1] <- "SampleID"
s <- data.frame(sample_data(qd))
rownames(s) <- NULL
alpha <- cbind(SampleID = rownames(alpha), alpha)
rownames(alpha) <- NULL
alphadiv <- cbind(alpha, s)
alphadiv <- alphadiv[,-5]
##



# Exploratory plots
theme_set(theme_bw())
# Rename time factors
og_names <- c("1", "2", "3", "4", "5")
new_names <- c("Jan 16", "March 16", "May 16", "July 16", "Jan 17")
# Make factors into numerics
alphadiv$richness <- as.numeric(alphadiv$richness)
alphadiv$evenness <- as.numeric(alphadiv$evenness)
alphadiv$faithPD <- as.numeric(alphadiv$faithPD)

a <- ggplot(alphadiv, aes(x = time, y = richness)) +
  geom_boxplot() + labs(y = "Chao1 Index") + 
  scale_x_discrete(breaks = og_names, labels = new_names) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
b <- ggplot(alphadiv, aes(x = time, y = evenness)) +
  geom_boxplot() + labs(y = "Simpson's Index") + 
  scale_x_discrete(breaks = og_names, labels = new_names) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
c <-ggplot(alphadiv, aes(x = time, y = faithPD)) +
  geom_boxplot() + labs(y = "Faith's PD") +  
  scale_x_discrete(breaks = og_names, labels = new_names) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

d <-ggplot(alphadiv, aes(x = species, y = richness)) +
  geom_boxplot() + labs(y = "Chao1 Index")
e <- ggplot(alphadiv, aes(x = species, y = evenness)) +
  geom_boxplot() + labs(y = "Simpson's Index")
f <- ggplot(alphadiv, aes(x = species, y = faithPD)) +
  geom_boxplot() + labs(y = "Faith's PD")

g <- ggplot(alphadiv, aes(x = treatment, y = richness)) +
  geom_boxplot() + labs(y = "Chao1 Index")
h <- ggplot(alphadiv, aes(x = treatment, y = evenness)) +
  geom_boxplot() + labs(y = "Simpson's Index")
i <- ggplot(alphadiv, aes(x = treatment, y = faithPD)) +
  geom_boxplot() + labs(y = "Faith's PD")

plot_grid(a,b,c,d,e,f,g,h,i, ncol = 3)

# Exploratory stats
pairwise.t.test(alphadiv$richness, alphadiv$time, p.adjust.method = "fdr")
pairwise.t.test(alphadiv$evenness, alphadiv$time, p.adjust.method = "fdr")
pairwise.t.test(alphadiv$faithPD, alphadiv$time, p.adjust.method = "fdr")

pairwise.t.test(alphadiv$richness, alphadiv$species, p.adjust.method = "fdr")
pairwise.t.test(alphadiv$evenness, alphadiv$species, p.adjust.method = "fdr")
pairwise.t.test(alphadiv$faithPD, alphadiv$species, p.adjust.method = "fdr")

pairwise.t.test(alphadiv$richness, alphadiv$treatment, p.adjust.method = "fdr")
pairwise.t.test(alphadiv$evenness, alphadiv$treatment, p.adjust.method = "fdr")
pairwise.t.test(alphadiv$faithPD, alphadiv$treatment, p.adjust.method = "fdr")
