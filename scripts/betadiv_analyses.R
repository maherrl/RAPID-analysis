
#################################################################################
# This script is the beta diversity analysis for the RAPID Project:
# Here, I calculate statistics for PERMANOVA and PERMDISP tests.
#
# Created by Rebecca Maher
# Created on 10/2/18
# Edited on 
#################################################################################

## clear workspace------------------------
rm(list=ls())

# load libraries
library('phyloseq'); packageVersion('phyloseq')
library('vegan'); packageVersion('vegan')

# set working directory-------------------
#setwd("~/Box Sync/RAPID-analysis/")

## functions (if any)----------------------

## Data Analysis----------------------------

# First import qd object, this is the phyloseq object from the import_qiime_to_phyloseq.R
#import(qd, file = "~/Box Sync/RAPID/RAPID-analysis/data/qd.RData")

# Calculate all four distance matrices: weighted unifrac, unweighted unifrac, bray curtis, binary jaccard
qd_wu <- phyloseq::distance(qd, method = "wunifrac")
qd_un <- phyloseq::distance(qd, method = "unifrac")
qd_bc <- phyloseq::distance(qd, method = "bray")
qd_bj <- distance(qd, method = "jaccard", binary =TRUE)

# PERMANOVA's with Adonis -  Permutational Multivariate Analasis of Variance
# Make a data frame from the sample_data
sampledf <- data.frame(sample_data(qd))

# Adonis individual tests
adonis(qd_bc ~ time, data = sampledf)
adonis(qd_bj ~ time, data = sampledf)
adonis(qd_wu ~ time, data = sampledf)
adonis(qd_un ~ time, data = sampledf)

adonis(qd_bc ~ treatment, data = sampledf)
adonis(qd_bj ~ treatment, data = sampledf)
adonis(qd_wu ~ treatment, data = sampledf)
adonis(qd_un ~ treatment, data = sampledf)

adonis(qd_bc ~ species, data = sampledf)
adonis(qd_bj ~ species, data = sampledf)
adonis(qd_wu ~ species, data = sampledf)
adonis(qd_un ~ species, data = sampledf)

# Adonis test for between group diversity, with full formula
adonis(qd_bc ~ time*species*treatment, data = sampledf)
adonis(qd_bj ~ time*species*treatment, data = sampledf) 
adonis(qd_wu ~ time*species*treatment, data = sampledf) 
adonis(qd_un ~ time*species*treatment, data = sampledf) 


# PERMDISP with betadisper - Multivariate Homogeneity of group dispersions
anova(betadisper(qd_bc, sampledf$species, bias.adjust = TRUE))
anova(betadisper(qd_bj, sampledf$species, bias.adjust = TRUE))
anova(betadisper(qd_wu, sampledf$species, bias.adjust = TRUE))
anova(betadisper(qd_un, sampledf$species, bias.adjust = TRUE))

anova(betadisper(qd_bc, sampledf$treatment, bias.adjust = TRUE))
anova(betadisper(qd_bj, sampledf$treatment, bias.adjust = TRUE))
anova(betadisper(qd_wu, sampledf$treatment, bias.adjust = TRUE))
anova(betadisper(qd_un, sampledf$treatment, bias.adjust = TRUE))

anova(betadisper(qd_bc, sampledf$time, bias.adjust = TRUE))
anova(betadisper(qd_bj, sampledf$time, bias.adjust = TRUE))
anova(betadisper(qd_wu, sampledf$time, bias.adjust = TRUE))
anova(betadisper(qd_un, sampledf$time, bias.adjust = TRUE))

