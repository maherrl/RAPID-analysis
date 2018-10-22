otutable <- data.frame(t(otu_table(qd)))
sampledf <- as.data.frame(sample_data(qd))
#merge(otu, sampledf, by.x = rownames(otu), by.y = 'sample.id')
# Why can I never get this to work!!!!
# I just did it in excell
write.csv(otutable, file = "/Users/Becca/Documents/otu.csv")
write.csv(sampledf, file = "/Users/Becca/Documents/sampledf.csv")

rapid <- read.csv(file = "/Users/Becca/Documents/otu.csv")
colnames(rapid)[1] <- "sample"

# Cutoff for a fraction of total counts the OTU has to represent (otu.cut).
# OTUs are removed if they are observed in less than specified fraction 
# of all samples (zero.cut)
goods=purgeOutliers(rapid, count.columns = 6:903, otu.cut = 0.0001, zero.cut = 0)
# Try this again with a zero.cut = 0.2

# what is the proportion of samples with data for these OTUs?
withData=apply(goods[,6:length(goods[1,])],2,function(x){sum(x>0)/length(x)})
hist(withData)

# what percentage of total counts each OTU represents?
props=apply(goods[,6:length(goods[1,])],2,function(x){sum(x)/sum(goods[,6:length(goods[1,])])})
barplot(sort(props,decreasing=T),xaxt="n",log="y")

# stacking the data; adjust otu.columns and condition.columns values for your data
ggs=otuStack(goods, count.columns = c(6:length(goods[1,])), condition.columns = c(1:5))

# fitting the model.
mm=mcmc.otu(
       fixed = "time+species+time:species",
       random = "indiv",
       data = ggs,
       nitt=55000,thin=50,burnin=5000
       )

# Selecting the OTUs that were modeled reliably
#(OTUs that are too rare for confident parameter estimates are discarded)
acpass=otuByAutocorr(mm,gss)

# calculating differences and p-values between all pairs of factor combinations
smm0=OTUsummary(mm,gss,otus = acpass, summ.plot = FALSE)
