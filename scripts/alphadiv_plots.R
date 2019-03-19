###################################################
## Alpha Diversity plots with time
## RAPID data
## Created by Becca Maher
## March 19, 2019
###################################################

alphadiv <- read.csv(file = "/Users/Becca/Box Sync/RAPID/RAPID-analysis/data/ACR_alphadiv.csv")
levels(alphadiv$time)

#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summariezed
# groupnames : vector of column names to be used as
# grouping variables
sderr <- function(x) {sd(x)/sqrt(length(x))}

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sderr(x[[col]]))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

df4 <- data_summary(alphadiv, varname="faithPD", 
                    groupnames=c("treatment", "time"))

breakss <- c("T0","T1","T2","T3","T4","T5","T6")
labelss <- c("Jan 16","March 16", "May 16", "July 16", "Sept 16","Nov 16","Jan 17")

p <- ggplot(df2, aes(x=time, y=richness, group=treatment, color = treatment)) + 
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin=richness-sd, ymax=richness+sd), width=.2,
                position = position_dodge(0.05)) +
  scale_x_discrete(breaks=breakss, labels=labelss) +
  ylab("Chao1 Index") + ggtitle("ACR Chao1 Index by Time")+ xlab("Time")
p
