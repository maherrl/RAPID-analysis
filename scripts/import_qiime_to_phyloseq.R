########################################################
# This script is for importing microbiome data from Qiime2
# to phyloseq for the RAPID data
# Created by Rebecca Maher
# Created on 9/20/18
# Edited on 10/3/18
########################################################

# clear workspace-----------------------------
rm(list=ls())

# load libraries
library('phyloseq')
packageVersion('phyloseq')
library(ape)
library("biomformat");packageVersion("biomformat")


##
# Importing data and transformations ##
# Import qiime mapping file, biom otu table, and tree
mapfile = "~/Box Sync/RAPID/RAPID-analysis/data/map.txt"
map = import_qiime_sample_data(mapfile)
class(map)

# Import tree
tree = read_tree("~/Box Sync/RAPID/RAPID-analysis/data/filt_tree.nwk")

# In qiime2, I already removed mitochondria and chloroplasts, removed otus present in only a single sample
# I also annotated the taxonomy with the Silva database
biomfile = "~/Box Sync/RAPID/RAPID-analysis/data/otu-table-no-mitochondria-no-chloroplast-min2-names-wtax-silva.biom"
biom = import_biom(biomfile, parseFunction = parse_taxonomy_default)

qd = merge_phyloseq(map,tree,biom) 
print(qd)
##phyloseq-class experiment-level object
##otu_table()   OTU Table:         [ 898 taxa and 280 samples ]
##sample_data() Sample Data:       [ 280 samples by 12 sample variables ]
##tax_table()   Taxonomy Table:    [ 898 taxa by 7 taxonomic ranks ]
##phy_tree()    Phylogenetic Tree: [ 898 tips and 897 internal nodes ]


# Changing the rank names of the phyloseq object
colnames(tax_table(qd)) = c(k="Kingdom", p="Phylum", c="Class", o="Order",f="Family", g="Genus", s="Species")
rank_names(qd)

# Make a numeric factor into a categorical
sample_data(qd)$time=factor(get_variable(qd,"time"))

# prune NA samples
samples_wna <- c("July2016.228.2", "July2016.228.1","July2016.251.2", "July2016.251.1", 
                 "July2016.220.2" ,"July2016.220.1" ,"March2016.acr35", "positive",
                 "July2016.294.1", "July2016.294.2")
samples <- sample_names(qd)
samples_nona <- setdiff(samples, samples_wna)
qd <- prune_samples(samples_nona, qd)
qd

# rarefy
set.seed(400)
qd <- rarefy_even_depth(qd, sample.size = 881, verbose = FALSE, replace = TRUE)
qd

# save rarefied phyloseq object
save(qd, file = "~/Box Sync/RAPID/RAPID-analysis/data/qd_881.RData")


# For differential abundance analyses (instead of rarefying)
# Remove samples with less than 881 reads
qd <- prune_samples(sample_sums(qd) > 880, qd)
qd
save(qd, file = "~/Box Sync/RAPID/RAPID-analysis/data/qd_DA.RData")

