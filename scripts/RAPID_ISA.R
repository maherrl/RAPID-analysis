#####################################################
## Indicator Species Analysis for RAPID ##
## September 12, 2019
## By Rebecca Maher
####################################################

rm(list=ls())

library("phyloseq")
library("data.table")
library("plyr")
library("dplyr")
library("ggplot2")
library("reshape2")
library("indicspecies")
library("ggnetwork")
library("ape")

load("~/Box Sync/RAPID/RAPID-analysis/data/qd_881.RData")
qd
qd <- subset_samples(qd, species =="ACR")

qd <- tax_glom(qd, taxrank = "Genus", bad_empty = c(NA, "", " ", "\t"))
qd
tab <- as.data.frame(t(otu_table(qd)))


# Making the named vector for the cluster option
v = sample_data(qd)$time
names(v) = rownames(sample_data(qd))
levels(v)
# Run indicator species analysis
vals <- multipatt(tab, cluster = v, func = "IndVal.g", control = how(nperm =999), duleg = TRUE)
summary(vals)

