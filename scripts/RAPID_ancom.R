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
library(dplyr)

# function
ANCOM.main = function(OTUdat,Vardat,
                      adjusted,repeated,
                      main.var,adj.formula,
                      repeat.var,longitudinal,
                      random.formula,
                      multcorr,sig,
                      prev.cut){
  
  p.zeroes=apply(OTUdat[,-1],2,function(x){
    s=length(which(x==0))/length(x)
  })
  
  zeroes.dist=data.frame(colnames(OTUdat)[-1],p.zeroes,row.names=NULL)
  colnames(zeroes.dist)=c("Taxon","Proportion_zero")
  
  zero.plot = ggplot(zeroes.dist, aes(x=Proportion_zero)) + 
    geom_histogram(binwidth=0.1,colour="black",fill="white") + 
    xlab("Proportion of zeroes") + ylab("Number of taxa") +
    theme_bw()
  
  #print(zero.plot)
  
  OTUdat.thinned=OTUdat
  OTUdat.thinned=OTUdat.thinned[,c(1,1+which(p.zeroes<prev.cut))]
  
  otu.names=colnames(OTUdat.thinned)[-1]
  
  W.detected   <- ancom.W(OTUdat.thinned,Vardat,
                          adjusted,repeated,
                          main.var,adj.formula,
                          repeat.var,longitudinal,random.formula,
                          multcorr,sig)
  
  W_stat       <- W.detected
  
  
  ### Bubble plot
  
  W_frame = data.frame(otu.names,W_stat,row.names=NULL)
  W_frame = W_frame[order(-W_frame$W_stat),]
  
  W_frame$detected_0.9=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.8=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.7=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.6=rep(FALSE,dim(W_frame)[1])
  
  W_frame$detected_0.9[which(W_frame$W_stat>0.9*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.8[which(W_frame$W_stat>0.8*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.7[which(W_frame$W_stat>0.7*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.6[which(W_frame$W_stat>0.6*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  
  final_results=list(W_frame,zero.plot)
  names(final_results)=c("W.taxa","PLot.zeroes")
  return(final_results)
}

# Ancom function
ancom.W = function(otu_data,var_data,
                   adjusted,repeated,
                   main.var,adj.formula,
                   repeat.var,long,rand.formula,
                   multcorr,sig){
  
  n_otu=dim(otu_data)[2]-1
  
  otu_ids=colnames(otu_data)[-1]
  
  if(repeated==F){
    data_comp=data.frame(merge(otu_data,var_data,by="Sample.ID",all.y=T),row.names=NULL)
    #data_comp=data.frame(merge(otu_data,var_data[,c("Sample.ID",main.var)],by="Sample.ID",all.y=T),row.names=NULL)
  }else if(repeated==T){
    data_comp=data.frame(merge(otu_data,var_data,by="Sample.ID"),row.names=NULL)
    # data_comp=data.frame(merge(otu_data,var_data[,c("Sample.ID",main.var,repeat.var)],by="Sample.ID"),row.names=NULL)
  }
  
  base.formula = paste0("lr ~ ",main.var)
  if(repeated==T){
    repeat.formula = paste0(base.formula," | ", repeat.var)
  }
  if(adjusted==T){
    adjusted.formula = paste0(base.formula," + ", adj.formula)
  }
  
  if( adjusted == F & repeated == F ){
    fformula  <- formula(base.formula)
  } else if( adjusted == F & repeated == T & long == T ){
    fformula  <- formula(base.formula)   
  }else if( adjusted == F & repeated == T & long == F ){
    fformula  <- formula(repeat.formula)   
  }else if( adjusted == T & repeated == F  ){
    fformula  <- formula(adjusted.formula)   
  }else if( adjusted == T & repeated == T  ){
    fformula  <- formula(adjusted.formula)   
  }else{
    stop("Problem with data. Dataset should contain OTU abundances, groups, 
         and optionally an ID for repeated measures.")
  }
  
  
  
  if( repeated==FALSE & adjusted == FALSE){
    if( length(unique(data_comp[,which(colnames(data_comp)==main.var)]))==2 ){
      tfun <- exactRankTests::wilcox.exact
    } else{
      tfun <- stats::kruskal.test
    }
  }else if( repeated==FALSE & adjusted == TRUE){
    tfun <- stats::aov
  }else if( repeated== TRUE & adjusted == FALSE & long == FALSE){
    tfun <- stats::friedman.test
  }else if( repeated== TRUE & adjusted == FALSE & long == TRUE){
    tfun <- nlme::lme
  }else if( repeated== TRUE & adjusted == TRUE){
    tfun <- nlme::lme
  }
  
  logratio.mat <- matrix(NA, nrow=n_otu, ncol=n_otu)
  for(ii in 1:(n_otu-1)){
    for(jj in (ii+1):n_otu){
      data.pair <- data_comp[,which(colnames(data_comp)%in%otu_ids[c(ii,jj)])]
      lr <- log((1+as.numeric(data.pair[,1]))/(1+as.numeric(data.pair[,2])))
      
      lr_dat <- data.frame( lr=lr, data_comp,row.names=NULL )
      
      if(adjusted==FALSE&repeated==FALSE){  ## Wilcox, Kruskal Wallis
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
      }else if(adjusted==FALSE&repeated==TRUE&long==FALSE){ ## Friedman's 
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
      }else if(adjusted==TRUE&repeated==FALSE){ ## ANOVA
        model=tfun(formula=fformula, data = lr_dat,na.action=na.omit)   
        picker=which(gsub(" ","",row.names(summary(model)[[1]]))==main.var)  
        logratio.mat[ii,jj] <- summary(model)[[1]][["Pr(>F)"]][picker]
      }else if(repeated==TRUE&long==TRUE){ ## GEE
        model=tfun(fixed=fformula,data = lr_dat,
                   random = formula(rand.formula),
                   correlation=corAR1(),
                   na.action=na.omit)   
        picker=which(gsub(" ","",row.names(anova(model)))==main.var)
        logratio.mat[ii,jj] <- anova(model)[["p-value"]][picker]
      }
      
    }
  } 
  
  ind <- lower.tri(logratio.mat)
  logratio.mat[ind] <- t(logratio.mat)[ind]
  
  
  logratio.mat[which(is.finite(logratio.mat)==FALSE)] <- 1
  
  mc.pval <- t(apply(logratio.mat,1,function(x){
    s <- p.adjust(x, method = "BH")
    return(s)
  }))
  
  a <- logratio.mat[upper.tri(logratio.mat,diag=FALSE)==TRUE]
  
  b <- matrix(0,ncol=n_otu,nrow=n_otu)
  b[upper.tri(b)==T] <- p.adjust(a, method = "BH")
  diag(b)  <- NA
  ind.1    <- lower.tri(b)
  b[ind.1] <- t(b)[ind.1]
  
  #########################################
  ### Code to extract surrogate p-value
  surr.pval <- apply(mc.pval,1,function(x){
    s0=quantile(x[which(as.numeric(as.character(x))<sig)],0.95)
    # s0=max(x[which(as.numeric(as.character(x))<alpha)])
    return(s0)
  })
  #########################################
  ### Conservative
  if(multcorr==1){
    W <- apply(b,1,function(x){
      subp <- length(which(x<sig))
    })
    ### Moderate
  } else if(multcorr==2){
    W <- apply(mc.pval,1,function(x){
      subp <- length(which(x<sig))
    })
    ### No correction
  } else if(multcorr==3){
    W <- apply(logratio.mat,1,function(x){
      subp <- length(which(x<sig))
    })
  }
  
  return(W)
  }


# load data for differential abundance analysis
load("~/Box Sync/RAPID/RAPID-analysis/data/qd_DA.RData")

qd


# Agglomerate samples to family ignoring, but still including taxa not identified to family
qd <- tax_glom(qd, taxrank = "Genus", bad_empty = c(NA, "", " ", "\t"))

qd_acr <- subset_samples(qd, species == "ACR")
qd_poc <- subset_samples(qd, species == "POC")
qd_por <- subset_samples(qd, species == "POR")

OTUdf <- as.data.frame(t(otu_table(qd)))
OTUdf <- cbind(Sample.ID = rownames(OTUdf), OTUdf)
rownames(OTUdf) <- NULL

OTUdf_acr <- as.data.frame(t(otu_table(qd_acr)))
OTUdf_acr <- cbind(Sample.ID = rownames(OTUdf_acr), OTUdf_acr)
rownames(OTUdf_acr) <- NULL

OTUdf_poc <- as.data.frame(t(otu_table(qd_poc)))
OTUdf_poc <- cbind(Sample.ID = rownames(OTUdf_poc), OTUdf_poc)
rownames(OTUdf_poc) <- NULL

OTUdf_por <- as.data.frame(t(otu_table(qd_por)))
OTUdf_por <- cbind(Sample.ID = rownames(OTUdf_por), OTUdf_por)
rownames(OTUdf_por) <- NULL


# OLD CODE
# OTU data or taxa data: This should be a data frame with each
# sample in rows and OTUs (or taxa) in columns. The first 
# column should be the sample identifier with column name
# "Sample.ID"

#OTUdat <- read.csv(file = "/Users/Becca/Box Sync/RAPID/family-table/otu-table-no-mitochondria-no-chloroplast-min2-names-wtax_L5.csv") # by family

# Metadata: Dataframe with the first columns being the sample 
# identifier with column name "Sample.ID"

Vardat <- read.csv(file = "./data/map.csv")
colnames(Vardat)[1] <- "Sample.ID"

# Subset map to only includ ACR samples
Vardat_acr <- Vardat[which(Vardat$species == "ACR"),]
Vardat_poc <- Vardat[which(Vardat$species == "POC"),]
Vardat_por <- Vardat[which(Vardat$species == "POR"),]

# ANCOM test
# First must run the function ANCOM.main from the ANCOM_updated_code.R files

# ACR
comp_test_time = ANCOM.main(OTUdat = OTUdf_acr, 
                       Vardat = Vardat_acr, 
                       adjusted = F,
                       repeated = F,
                       main.var = "time",
                       adj.formula = NULL,
                       repeat.var = NULL,
                       longitudinal = F,
                       random.formula = "~1|indiv",
                       multcorr = 2,
                       sig=0.05,
                       prev.cut = 0.90)

resA <- comp_test_time$W.taxa
head(resA)
resA <- resA[which(resA$detected_0.9 == "TRUE"),]
tax <- as.data.frame(qd_acr@tax_table@.Data)
tax$otu.names <- rownames(tax)
resA_tax <- merge(resA, tax, by = "otu.names")
write.csv(resA_tax, file = "~/Box Sync/RAPID/RAPID-analysis/results/ancom_acr.csv")

# plotting
otu.names <- as.vector(resA_tax$otu.names)

qdrr  = transform_sample_counts(qd_acr, function(x) x / sum(x) )
qdrr
sample_sums(qdrr)
qdrr <- prune_taxa(otu.names, qdrr)
qdrr

taxa.sum <- as.data.frame(sort(taxa_sums(qdrr), decreasing = F))
otu.names.sorted <- as.vector(rownames(taxa.sum))

plot_heatmap(qdrr, sample.label = "time", taxa.order = otu.names.sorted,
             taxa.label = "Genus", sample.order = "time",
             low="grey", high="darkblue", na.value = "white")


# POC
comp_test_time = ANCOM.main(OTUdat = OTUdf_poc, 
                            Vardat = Vardat_poc, 
                            adjusted = F,
                            repeated = F,
                            main.var = "time",
                            adj.formula = NULL,
                            repeat.var = NULL,
                            longitudinal = F,
                            random.formula = "~1|indiv",
                            multcorr = 2,
                            sig=0.05,
                            prev.cut = 0.90)

resP <- comp_test_time$W.taxa
resP <- resP[which(resP$detected_0.9 == "TRUE"),]
tax <- as.data.frame(qd_poc@tax_table@.Data)
tax$otu.names <- rownames(tax)
resP_tax <- merge(resP, tax, by = "otu.names")
write.csv(resP_tax, file = "~/Box Sync/RAPID/RAPID-analysis/results/ancom_poc.csv")

# plotting
otu.names <- as.vector(resP_tax$otu.names)

qdrr  = transform_sample_counts(qd_poc, function(x) x / sum(x) )
qdrr
sample_sums(qdrr)
qdrr <- prune_taxa(otu.names, qdrr)
qdrr

taxa.sum <- as.data.frame(sort(taxa_sums(qdrr), decreasing = F))
otu.names.sorted <- as.vector(rownames(taxa.sum))

plot_heatmap(qdrr, sample.label = "time", taxa.order = otu.names.sorted,
             taxa.label = "Genus", sample.order = "time",
             low="grey", high="darkblue", na.value = "white")

# POR
comp_test_time = ANCOM.main(OTUdat = OTUdf_por, 
                            Vardat = Vardat_por, 
                            adjusted = F,
                            repeated = F,
                            main.var = "time",
                            adj.formula = NULL,
                            repeat.var = NULL,
                            longitudinal = F,
                            random.formula = "~1|indiv",
                            multcorr = 2,
                            sig=0.05,
                            prev.cut = 0.90)

resR <- comp_test_time$W.taxa
resR <- resR[which(resR$detected_0.9 == "TRUE"),]
tax <- as.data.frame(qd_por@tax_table@.Data)
tax$otu.names <- rownames(tax)
resR_tax <- merge(resR, tax, by = "otu.names")
write.csv(resR_tax, file = "~/Box Sync/RAPID/RAPID-analysis/results/ancom_por.csv")

# plotting
otu.names <- as.vector(resR_tax$otu.names)

qdrr  = transform_sample_counts(qd_por, function(x) x / sum(x) )
qdrr
sample_sums(qdrr)
qdrr <- prune_taxa(otu.names, qdrr)
qdrr

taxa.sum <- as.data.frame(sort(taxa_sums(qdrr), decreasing = F))
otu.names.sorted <- as.vector(rownames(taxa.sum))

plot_heatmap(qdrr, sample.label = "time", taxa.order = otu.names.sorted,
             taxa.label = "Genus", sample.order = "time",
             low="grey", high="darkblue", na.value = "white")

# ALL
comp_test_time = ANCOM.main(OTUdat = OTUdf, 
                            Vardat = Vardat, 
                            adjusted = F,
                            repeated = F,
                            main.var = "time",
                            adj.formula = NULL,
                            repeat.var = NULL,
                            longitudinal = F,
                            random.formula = "~1|indiv",
                            multcorr = 2,
                            sig=0.05,
                            prev.cut = 0.90)

res <- comp_test_time$W.taxa
res <- res[which(res$detected_0.9 == "TRUE"),]
tax <- as.data.frame(qd@tax_table@.Data)
tax$otu.names <- rownames(tax)
res_tax <- merge(res, tax, by = "otu.names")
write.csv(res_tax, file = "~/Box Sync/RAPID/RAPID-analysis/results/ancom_all.csv")

# plotting
otu.names <- as.vector(res_tax$otu.names)

qdrr  = transform_sample_counts(qd, function(x) x / sum(x) )
qdrr
sample_sums(qdrr)
qdrr <- prune_taxa(otu.names, qdrr)
qdrr

taxa.sum <- as.data.frame(sort(taxa_sums(qdrr), decreasing = F))
otu.names.sorted <- as.vector(rownames(taxa.sum))

plot_heatmap(qdrr, sample.label = "time", taxa.order = otu.names.sorted,
             taxa.label = "Genus", sample.order = "time",
             low="grey", high="darkblue", na.value = "white")
