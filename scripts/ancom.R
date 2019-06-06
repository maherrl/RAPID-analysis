########################################################
# This script is for performing an ANCOM (Analysis of
# Composition of Microbiomes (ANCOM) to detect 
# differentially abundant taxa in microbial surveys.
# Created by Rebecca Maher
# Created on 10/26/18
# Edited on 01/17/19
########################################################

# clear workspace-----------------------------
rm(list=ls())

# load libraries
library(exactRankTests)
library(nlme)
library(stats)
library(ggplot2)

# load data
load("~/Box Sync/RAPID/RAPID-analysis/data/qdunrar1000.RData")

# Subset to only acropora samples
qd_acr <- subset_samples(qd, species == "ACR")
# Agglomerate samples to family ignoring, but still including taxa not identified to family
qd_acr_fam <- tax_glom(qd_acr, taxrank = "Family", bad_empty = c(""))

OTUdf <- as.data.frame(t(otu_table(qd_acr_fam)))
OTUdf <- cbind(Sample.ID = rownames(OTUdf), OTUdf)
rownames(OTUdf) <- NULL

# OLD CODE
# OTU data or taxa data: This should be a data frame with each
# sample in rows and OTUs (or taxa) in columns. The first 
# column should be the sample identifier with column name
# "Sample.ID"

#OTUdat <- read.csv(file = "/Users/Becca/Box Sync/RAPID/family-table/otu-table-no-mitochondria-no-chloroplast-min2-names-wtax_L5.csv") # by family

# Metadata: Dataframe with the first columns being the sample 
# identifier with column name "Sample.ID"

Vardat <- read.csv(file = "/Users/Becca/Box Sync/RAPID/RAPID-analysis/data/map.csv")
colnames(Vardat)[1] <- "Sample.ID"

# Subset map to only includ ACR samples from March and May with no NAs
Vardat <- Vardat[which(Vardat$species == "ACR" & Vardat$Sample.ID != "March2016.acr35"),]


# ANCOM test
# First must run the function ANCOM.main from the ANCOM_updated_code.R files

longitudinal_comp_test = ANCOM.main(OTUdat = OTUdf, 
                                    Vardat = Vardat, 
                                    adjusted = F,
                                    repeated = T,
                                    main.var = "treatment",
                                    adj.formula = NULL,
                                    repeat.var = "time",
                                    longitudinal = T,
                                    random.formula = "~1|indiv",
                                    multcorr = 2,
                                    sig=0.05,
                                    prev.cut = 0.90)

longitudinal_comp_test$W.taxa

