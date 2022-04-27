#Reproduction figure 6

library(ggplot2)
library(tidyr)
library(plyr)
library(ggplot2)
library(scales)

#input file is within this folder
time_exec <- read.table("Time_Execution_Methods_TT_plus_Linda.txt",header=T)
head(time_exec)
#names(time_exec)[8] <- "cov1mln.250k"
long_df <- time_exec
#long_df <- time_exec %>%  gather(coverage, time, cov50k:cov1mln.250k)
#long_df <- time_exec %>%  gather(coverage, time, Method)

head(long_df)
long_df$Time <- as.numeric(long_df$Time)
long_df$Coverage <- factor(long_df$Coverage, levels=c("50K", "100K", "250K", "500K", "750K", "1M", "1M250K"))
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
long_df <- as.data.frame(long_df)
df2 <- data_summary(long_df, varname="Time", 
                    groupnames=c("Coverage", "Method","Read_Length"))
#quality check
head(df2)

#3) plot the standard deviation and the median of each method 
p <- ggplot(df2,aes(Coverage,Time,colour=Method,group=Read_Length))
pd<-position_dodge(.9)

p +   geom_errorbar(aes(ymin=Time-sd,ymax=Time+sd),width=.1,position=pd,colour="black") +
  geom_point(position=pd,size=4) + geom_line(position=pd) + 
  theme_bw() +  
  facet_grid(Method~Read_Length)+theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) +
  scale_y_continuous(trans = log2_trans(),breaks = trans_breaks("log2", function(x) 2^x),labels = trans_format("log2", math_format(2^.x)))+
  scale_y_log10()+
  annotation_logticks()+
  ggtitle("")+ theme(plot.title = element_text(hjust = 0.5))+theme(text=element_text(size=15))+ xlab("Coverage") + ylab("Time (Seconds)")+
  theme(text=element_text(size=15))+theme(axis.text.x = element_text(face="bold", color="#000001")) + theme(axis.text.y = element_text(face="bold", color="#000001"))+theme(text = element_text(size = 18, face="bold"))+theme(axis.title.y = element_text(size = 18))+theme(legend.text = element_text(size = 18,face="plain"))+theme(legend.title = element_text(size = 18))  


#p+axis.title=element_text(size=18,face="bold")+theme(plot.title = element_text(hjust = 0.5))+theme(axis.text=element_text(size=14),axis.title=element_text(size=14,face="bold"))+ggtitle("")+ xlab("") + ylab("Time (Seconds)") +





