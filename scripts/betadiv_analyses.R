
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

## functions----------------------
pairwise.adonis.dm <- function(x,factors,stratum=NULL,p.adjust.m="bonferroni",perm=999){
  
  library(vegan)
  
  if(class(x) != "dist") stop("x must be a dissimilarity matrix (dist object)")
  
  co = as.matrix(combn(unique(factors),2))
  
  pairs = c()
  
  F.Model =c()
  
  R2 = c()
  
  p.value = c()
  
  
  
  for(elem in 1:ncol(co)){
    
    sub_inds <- factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))
    
    resp <- as.matrix(x)[sub_inds,sub_inds]
    
    ad = adonis(as.dist(resp) ~
                  
                  factors[sub_inds], strata=stratum[sub_inds], permutations=perm);
    
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    
    R2 = c(R2,ad$aov.tab[1,5]);
    
    p.value = c(p.value,ad$aov.tab[1,6])
    
  }
  
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  
  return(pairw.res)
  }

## Data Analysis----------------------------

# First import qd object, this is the phyloseq object from the import_qiime_to_phyloseq.R
import(qd, file = "~/Box Sync/RAPID/RAPID-analysis/data/qd_rare.RData")

# Calculate all four distance matrices: weighted unifrac, unweighted unifrac, bray curtis, binary jaccard
qd_wu <- phyloseq::distance(qd, method = "wunifrac")
qd_un <- phyloseq::distance(qd, method = "unifrac")
qd_bc <- phyloseq::distance(qd, method = "bray")
qd_bj <- distance(qd, method = "jaccard", binary =TRUE)

ord_wu <- ordinate(qd, "PCoA", distance = qd_wu)
ord_un <- ordinate(qd, "PCoA", distance = qd_un)
ord_bc <- ordinate(qd, "PCoA", distance = qd_bc)
ord_bj <- ordinate(qd, "PCoA", distance = qd_bj)

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

## Add distances into a mapping file for plotting
wu_means <- data.frame(rowMeans(as.matrix(qd_wu)))
un_means <- data.frame(rowMeans(as.matrix(qd_un)))
bc_means <- data.frame(rowMeans(as.matrix(qd_bc)))
bj_means <- data.frame(rowMeans(as.matrix(qd_bj)))
# Make a dataframe
# Extracts the SampleIDs from the dataframe
SampleID <- row.names(wu_means) 
# makes a new dataframe with the specified columnes
raw_distances <- data.frame(SampleID, wu_means, un_means, bc_means,bj_means) 
# Makes metadata into a df to work with
s <- data.frame(sample_data(qd)) 
# Change first column title to SampleID to match distances dataframe
colnames(s)[1] <- "SampleID" 
# merges metadata df and distances df
distances <- merge(raw_distances, s, by = "SampleID") 
colnames(distances)[2:5] <- c("distwu", "distun", "distbc","distbj")

# Pairwise adonis
pairwise.adonis.dm(qd_bc, sample_data(qd)$treatment, p.adjust.m = "fdr")
pairwise.adonis.dm(qd_bc, sample_data(qd)$time, p.adjust.m = "fdr")
pairwise.adonis.dm(qd_bc, sample_data(qd)$species, p.adjust.m = "fdr")
# Pairwise betadisper with fdr correction
p.adjust(permutest(betadisper(qd_un, sampledf$treatment, 
                              bias.adjust = TRUE), 
                   pairwise=TRUE)$pairwise$permuted, method = 'fdr')
p.adjust(permutest(betadisper(qd_un, sampledf$species, 
                              bias.adjust = TRUE), 
                   pairwise=TRUE)$pairwise$permuted, method = 'fdr')
p.adjust(permutest(betadisper(qd_un, sampledf$time, 
                              bias.adjust = TRUE), 
                   pairwise=TRUE)$pairwise$permuted, method = 'fdr')


