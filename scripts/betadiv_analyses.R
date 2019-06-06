
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
import(file = "~/Box Sync/RAPID/RAPID-analysis/data/qd_rare.RData")

# subset by species
qd <- subset_samples(qd, species =="ACR")
#qd <- subset_samples(qd, species =="POC")
#qd <- subset_samples(qd, time == "T2")

# Log-transform OTU table
otus_log <- as(otu_table(qd), "matrix")
otus_log <- decostand(otus_log, method = "log")
OTU_log <- otu_table(otus_log, taxa_are_rows = TRUE)
otu_table(qd) <- OTU_log

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

adonis(qd_bc ~ time*treatment, data = sampledf)
adonis(qd_bj ~ time*treatment, data = sampledf) 
adonis(qd_wu ~ time*treatment, data = sampledf) 
adonis(qd_un ~ time*treatment, data = sampledf) 

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

## Add between group distances into a mapping file for plotting
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
bet_distances <- merge(raw_distances, s, by = "SampleID") 
colnames(bet_distances)[2:5] <- c("betdistwu", "betdistun", "betdistbc","betdistbj")

### Within group distances

# Here I define groups as trt.time by making a new column in the mapping file with treatment 
# and time merged into one column so that groups are for ex: May 2016.N or Jan 2017.C
# Run betadisper on the time_trt combo
betabc <- betadisper(qd_bc, sampledf$trt.time, bias.adjust = TRUE)
betabj <- betadisper(qd_bj, sampledf$trt.time, bias.adjust = TRUE)
betaun <- betadisper(qd_un, sampledf$trt.time, bias.adjust = TRUE)
betawu <- betadisper(qd_wu, sampledf$trt.time, bias.adjust = TRUE)

# Extract within group distances and export as .csv
# Gives the distances of each group to its centroid by interaction
withdistbc <- betabc$distances
withdistbj <- betabj$distances
withdistun <- betaun$distances
withdistwu <- betawu$distances

# IMPORTANT STEP: Merge the distance to centroid for each distance measure by SampleID
# Put them into the map file for use in the graphs

# Makes a dataframe of the distance to centroid of one distance measure
# Just necessary to get the SampleID labels
df_bc <- data.frame(betabc$distances) 
# Extracts the SampleIDs from the dataframe
SampleID <- row.names(df_bc) 
# makes a new dataframe with the specified columnes
with_distances <- data.frame(SampleID, withdistbc, withdistbj, withdistun, withdistwu) 
# merges metadata df and distances df
withdistmax <- merge(bet_distances, with_distances, by = "SampleID") 
write.csv(withdistmax, file = "/Users/Becca/Box Sync/RAPID/RAPID-analysis/data/beta_group_distances_rarefied.csv")


# Pairwise adonis
pairwise.adonis.dm(qd_bc, sample_data(qd)$treatment, p.adjust.m = "fdr")
pairwise.adonis.dm(qd_wu, sample_data(qd)$time, p.adjust.m = "fdr")
pairwise.adonis.dm(qd_bc, sample_data(qd)$species, p.adjust.m = "fdr")
# Pairwise betadisper with fdr correction
p.adjust(permutest(betadisper(qd_bc, sampledf$treatment, 
                              bias.adjust = TRUE), 
                   pairwise=TRUE)$pairwise$permuted, method = 'fdr')
p.adjust(permutest(betadisper(qd_un, sampledf$species, 
                              bias.adjust = TRUE), 
                   pairwise=TRUE)$pairwise$permuted, method = 'fdr')
p.adjust(permutest(betadisper(qd_wu, sampledf$time, 
                              bias.adjust = TRUE), 
                   pairwise=TRUE)$pairwise$permuted, method = 'none')


load("/Users/Becca/Box Sync/RAPID/RAPID-analysis/data/qdrar1778.RData") # RData object was created and saved using the import_qiime_to_phyloseq.R script and loaded here as "qd"
# Loaded full rarefied phyloseq object again
qd <- subset_samples(qd, species =="ACR")
qd <- subset_samples(qd, time == "T1" | time == "T2")
summary(sample_data(qd))


#Calculate the dissimilarity matrics

qd_wu <- phyloseq::distance(qd, method = "wunifrac")
qd_un <- phyloseq::distance(qd, method = "unifrac")
qd_bc <- phyloseq::distance(qd, method = "bray")
qd_bj <- distance(qd, method = "jaccard", binary =TRUE)

#Make a sample dataframe.

sampledf <- data.frame(sample_data(qd))

#Adonis tests.

adonis(qd_bc ~ treatment, data = sampledf)
adonis(qd_bj ~ treatment, data = sampledf)
adonis(qd_wu ~ treatment, data = sampledf)
adonis(qd_un ~ treatment, data = sampledf)

adonis(qd_bc ~ time, data = sampledf)
adonis(qd_bj ~ time, data = sampledf)
adonis(qd_wu ~ time, data = sampledf)
adonis(qd_un ~ time, data = sampledf)

adonis(qd_bc ~ nutrient, data = sampledf)
adonis(qd_bj ~ nutrient, data = sampledf)
adonis(qd_wu ~ nutrient, data = sampledf)
adonis(qd_un ~ nutrient, data = sampledf)

adonis(qd_bc ~ time*nutrient, data = sampledf)
adonis(qd_bj ~ time*nutrient, data = sampledf)
adonis(qd_wu ~ time*nutrient, data = sampledf)
adonis(qd_un ~ time*nutrient, data = sampledf)

adonis(qd_bc ~ time*treatment, data = sampledf)
adonis(qd_bj ~ time*treatment, data = sampledf) 
adonis(qd_wu ~ time*treatment, data = sampledf) 
adonis(qd_un ~ time*treatment, data = sampledf) 

#Betadisper tests.

anova(betadisper(qd_bc, sampledf$treatment, bias.adjust = TRUE))
anova(betadisper(qd_bj, sampledf$treatment, bias.adjust = TRUE))
anova(betadisper(qd_wu, sampledf$treatment, bias.adjust = TRUE))
anova(betadisper(qd_un, sampledf$treatment, bias.adjust = TRUE))

anova(betadisper(qd_bc, sampledf$time, bias.adjust = TRUE))
anova(betadisper(qd_bj, sampledf$time, bias.adjust = TRUE))
anova(betadisper(qd_wu, sampledf$time, bias.adjust = TRUE))
anova(betadisper(qd_un, sampledf$time, bias.adjust = TRUE))

anova(betadisper(qd_bc, sampledf$nutrient, bias.adjust = TRUE))
anova(betadisper(qd_bj, sampledf$nutrient, bias.adjust = TRUE))
anova(betadisper(qd_wu, sampledf$nutrient, bias.adjust = TRUE))
anova(betadisper(qd_un, sampledf$nutrient, bias.adjust = TRUE))

betabc <- betadisper(qd_bc, sampledf$treatment, bias.adjust = TRUE)
distbc <- betabc$distances

anova(betadisper(qd_bc, sampledf$trt.time, bias.adjust = TRUE))
anova(betadisper(qd_bj, sampledf$trt.time, bias.adjust = TRUE))
anova(betadisper(qd_wu, sampledf$trt.time, bias.adjust = TRUE))
anova(betadisper(qd_un, sampledf$trt.time, bias.adjust = TRUE))

# FDR corrected gave no significant results so trying with no correction
pairwise.dips.un.padj <- p.adjust(permutest(betadisper(qd_un, sampledf$trt.time, 
                              bias.adjust = TRUE), 
                   pairwise=TRUE)$pairwise$permuted, method = 'none')

sig.pairwise.dips.un.padj <- which(pairwise.dips.un.padj < 0.05)
# Try without pval correction
pairwise.disp.un <- permutest(betadisper(qd_un, sampledf$trt.time, 
                     bias.adjust = TRUE), 
          pairwise=TRUE)
