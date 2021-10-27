setwd("/cloud-data/snf-mgln-dds/AIDA/Bioinformatics/i0439277/RTcure/Benchmark/Pilot1_Sanger_TT/plot_time_execution/")
library(ggplot2)
library(tidyr)
#install.packages("plyr")
library(plyr)
library(ggplot2)

time_exec <- read.table("Time_Execution_Methods_TT_plus_Linda.txt",header=T)
#names(time_exec)[8] <- "cov1mln.250k"
long_df <- time_exec
#long_df <- time_exec %>%  gather(coverage, time, cov50k:cov1mln.250k)
#long_df <- time_exec %>%  gather(coverage, time, Method)

head(long_df)
long_df$Time <- as.numeric(long_df$Time)
long_df$Coverage <- factor(long_df$Coverage, levels=c("50k", "100k", "250k", "500k", "750k", "1mln", "1m250k"))
long_df$Read_Length <- factor(long_df$Read_Length, levels=c("25bp", "50bp", "75bp", "100bp"))
#now lets plot using the standatd deviation

#1) first create a function for this

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

#2) compute the standard deviation
df2 <- data_summary(long_df, varname="Time", 
                    groupnames=c("Coverage", "Method","Read_Length"))
#quality check
head(df2)

#3) plot the standard deviation and the median of each method 
p <- ggplot(df2,aes(Coverage,Time,colour=Method,group=Read_Length))
pd<-position_dodge(.9)

p + 
  geom_errorbar(aes(ymin=Time-sd,ymax=Time+sd),width=.1,position=pd,colour="black") +
  geom_point(position=pd,size=4) + geom_line(position=pd) + 
  theme_bw() +  
  facet_grid(Method~Read_Length)+theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) +
  scale_y_continuous(trans = log2_trans(),breaks = trans_breaks("log2", function(x) 2^x),labels = trans_format("log2", math_format(2^.x)))+
  scale_y_log10()+
  annotation_logticks()+
  ggtitle("Time of execution \n different coverage and reads length internal dataset")+ theme(plot.title = element_text(hjust = 0.5))+theme(text=element_text(size=15))+ xlab("Coverage") + ylab("Time (Seconds)")









########################################################
########################################################
rm(list = ls())
library(data.table)
library(reshape2)

time_exec <- read.table("Time_Execution_Methods_TT_plus_Linda.txt",header=T)
#names(time_exec)[8] <- "cov1mln.250k"
long_df <- time_exec
#long_df <- time_exec %>%  gather(coverage, time, cov50k:cov1mln.250k)
#long_df <- time_exec %>%  gather(coverage, time, Method)

head(long_df)
long_df$Time <- as.numeric(long_df$Time)
long_df$Coverage <- factor(long_df$Coverage, levels=c("50k", "100k", "250k", "500k", "750k", "1mln", "1m250k"))
long_df$Read_Length <- factor(long_df$Read_Length, levels=c("25bp", "50bp", "75bp", "100bp"))

long_df
#long_df <- subset(long_df, Coverage == "1m250k")
long_df <- subset(long_df, Read_Length == "25bp")
p <- ggplot(data=long_df,
       aes(x=Coverage, y=Time, colour=Method)) +
        geom_boxplot()

p + scale_y_continuous(trans = log2_trans(),
                       breaks = trans_breaks("log2", function(x) 2^x),
                       labels = trans_format("log2", math_format(2^.x)))+
  ggtitle("Time of execution \n different coverage 25bp internal dataset")+ theme(plot.title = element_text(hjust = 0.5))+theme(text=element_text(size=15))+ xlab("Coverage") + ylab("Time (Seconds)")



