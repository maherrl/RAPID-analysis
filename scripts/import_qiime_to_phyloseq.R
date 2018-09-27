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

qd1 = merge_phyloseq(map,tree,biom) 

# Changing the rank names of the phyloseq object
colnames(tax_table(qd)) = c(k="Kingdom", p="Phylum", c="Class", o="Order",f="Family", g="Genus", s="Species")
rank_names(qd)

# Normalization technique using Naive Proportions from the "Waste Not, Want Not"
# paper

# I could not get this code to work. I was able to remove singletons, but not
# remove OTUs appearing in at least 3 samples, error in prune_taxa (wh1, qd)
# But not necessary because I did this in qiime.

# Remove any samples or OTUs that don't have any reads.
# Also remove singletons (only 1 count as sample or OTU)
# remove_singletons_empties = function(qd){
#   require("phyloseq")
#     qd = prune_taxa(taxa_sums(qd) > 2.5, qd)
#     qd = prune_samples(sample_sums(qd) > 2.5, qd)
#   # Remove OTUs not appearing in at least 3 samples
#   if( taxa_are_rows(qd)){
#     y=otu_table(qd)
#     wh1 = apply(apply(y, 1, function(x){x>=1}), MARGIN = 1, sum) >= 3
#   } else {
#     y = as(t(otu_table(qd)), "matrix")
#     wh1 = apply(apply(y, 1, function(x){x>=1}), MARGIN = 1, sum) >= 3
#   }
#   qd = prune_taxa(wh1, qd)
#     return(qd)
# }

# remove_singletons_empties(qd)

# Define the naive (simple proportion) normalization function.

proportion = function(qd){
  # Normalize total sequences represented
  normf = function(x, tot=max(sample_sums(qd))){tot*x/sum(x)}
  qd = transform_sample_counts(qd, normf)
  # Scale by dividing each variable by its standard deviation
  qd = transform_sample_counts(qd,function(x) x/sd(x))
  # Center by subtracting the median
  qd = transform_sample_counts(qd, function(x) (x-median(x)))
  return(qd)
}

proportion(qd)

qd1000 <- prune_samples(sample_sums(qd)>=1000, qd)
qd1000_acr <- subset_samples(qd1000, species == "ACR")




# Make a numeric factor into a categorical
# In this case, i had temperature as 26 or 29, R considers this numerical but I want to consider it categorical
sample_data(qd)$temp=factor(get_variable(qd,"temp"))
##

##
# Save the Formal class phyloseq to an external file to load in other scripts (qd = qiimedata)
save(qd, file = "/filepath/qd.RData")
# The phyloseq object qd is now ready to use in diversity analyses