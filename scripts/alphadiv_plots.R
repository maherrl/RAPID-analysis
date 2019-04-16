###################################################
## Alpha and Beta Diversity plots with time
## RAPID data
## Created by Becca Maher
## March 19, 2019
###################################################

library("ggplot2")
library("ggthemes")

alphadiv <- read.csv(file = "/Users/Becca/Box Sync/RAPID/RAPID-analysis/data/all_alphadiv.csv")
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

df_fpd <- data_summary(alphadiv, varname="faithPD", groupnames=c("treatment", "time"))
df_simp <- data_summary(alphadiv, varname="evenness", groupnames=c("treatment", "time"))
df_chao <- data_summary(alphadiv, varname="richness", groupnames=c("treatment", "time"))

breakss <- c("T0","T1","T2","T3","T4","T5","T6")
labelss <- c("Jan 16","March 16", "May 16", "July 16", "Sept 16","Nov 16","Jan 17")

p <- ggplot(df_chao, aes(x=time, y=richness, group=treatment, color = treatment)) + 
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin=richness-sd, ymax=richness+sd), width=.2,
                position = position_dodge(0.05)) +
  scale_x_discrete(breaks=breakss, labels=labelss) +
  ylab("Chao1 Index") + ggtitle("All Species Chao1 Index by Time") + xlab("Time")
p + scale_colour_colorblind() + theme_bw()

p <- ggplot(df_simp, aes(x=time, y=evenness, group=treatment, color = treatment)) + 
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin=evenness-sd, ymax=evenness+sd), width=.2,
                position = position_dodge(0.05)) +
  scale_x_discrete(breaks=breakss, labels=labelss) +
  ylab("Simpson's Index") + ggtitle("All Species Simpson's Index by Time") + xlab("Time")
p + scale_colour_colorblind() + theme_bw()

p <- ggplot(df_fpd, aes(x=time, y=faithPD, group=treatment, color = treatment)) + 
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin=faithPD-sd, ymax=faithPD+sd), width=.2,
                position = position_dodge(0.05)) +
  scale_x_discrete(breaks=breakss, labels=labelss) +
  ylab("Faith's Phylogenetic Distance") + ggtitle("All Species Faith's PD by Time") + xlab("Time")
p + scale_colour_colorblind() + theme_bw()

## Beta plots
betadiv <- read.csv(file = "/Users/Becca/Box Sync/RAPID/RAPID-analysis/data/beta_group_distances_rarefied.csv")
df_bet_wu <- data_summary(betadiv, varname="betdistwu", groupnames=c("treatment", "time"))
df_bet_bc <- data_summary(betadiv, varname="betdistbc", groupnames=c("treatment", "time"))
df_bet_bj <- data_summary(betadiv, varname="betdistbj", groupnames=c("treatment", "time"))
df_bet_un <- data_summary(betadiv, varname="betdistun", groupnames=c("treatment", "time"))
df_with_wu <- data_summary(betadiv, varname="withdistwu", groupnames=c("treatment", "time"))
df_with_bc <- data_summary(betadiv, varname="withdistbc", groupnames=c("treatment", "time"))
df_with_bj <- data_summary(betadiv, varname="withdistbj", groupnames=c("treatment", "time"))
df_with_un <- data_summary(betadiv, varname="withdistun", groupnames=c("treatment", "time"))


p <- ggplot(df_with_un, aes(x=time, y=withdistun, group=treatment, color = treatment)) + 
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin=withdistun-sd, ymax=withdistun+sd), width=.2,
                position = position_dodge(0.05)) +
  scale_x_discrete(breaks=breakss, labels=labelss) +
  ylab("Between Group Distance") + ggtitle("All Species - UN") + xlab("Time")
p + scale_colour_colorblind() + theme_bw()
