########################################################
# This script is for importing microbiome data from Qiime2
# to phyloseq for the RAPID data
# Created by Rebecca Maher
# Created on 9/20/18
# Edited on 9/27/18
########################################################

# clear workspace-----------------------------
rm(list=ls())

# load libraries
library('phyloseq')
packageVersion('phyloseq')
library(ape)
library("biomformat");packageVersion("biomformat")

setwd("~/Box Sync/RAPID-analysis/data")

##
# Importing data and transformations ##
# Import qiime mapping file, biom otu table, and tree
mapfile = "~/Box Sync/RAPID/RAPID-analysis/data/map.txt"
map = import_qiime_sample_data(mapfile)
class(map)

# Import tree
tree = read_tree("~/Box Sync/RAPID/RAPID-analysis/data/tree.nwk")

# In qiime2, I already removed mitochondria and chloroplasts, removed otus with < 10 occurrences and
# present in less than 1 sample
biomfile = "/Users/Becca/Box Sync/RAPID/RAPID-analysis/data/otu-table-no-mitochondria-no-chloroplast-min2-names-wtax.biom"
biom = import_biom(biomfile, parseFunction = parse_taxonomy_default)

qd = merge_phyloseq(map,tree,biom) 

# Changing the rank names of the phyloseq object
colnames(tax_table(qd)) = c(k="Kingdom", p="Phylum", c="Class", o="Order",f="Family", g="Genus", s="Species")
rank_names(qd)

# Normalization technique using Naive Proportions from the "Waste Not, Want Not"
# paper
# Define the naive (simple proportion) normalization function.
# Normalize total sequences represented
normf = function(x, tot=max(sample_sums(qd))){tot*x/sum(x)}
qd = transform_sample_counts(qd, normf)
# Scale by dividing each variable by its standard deviation
qd = transform_sample_counts(qd,function(x) x/sd(x))
# Center by subtracting the median
# qd = transform_sample_counts(qd, function(x) x - median(x)) # didnt do anything?

# Export it to Qiime2
otu <- as(otu_table(qd), "matrix")
otu_biom <- make_biom(data=otu)
write_biom(otu_biom, "/Users/Becca/Box Sync/RAPID/RAPID-analysis/data/otu-table-prop.biom")

# Make a numeric factor into a categorical
# In this case, i had temperature as 26 or 29, R considers this numerical but I want to consider it categorical
sample_data(qd)$temp=factor(get_variable(qd,"temp"))
##

##
# Save the Formal class phyloseq to an external file to load in other scripts (qd = qiimedata)
save(qd, file = "/filepath/qd.RData")
# The phyloseq object qd is now ready to use in diversity analyses