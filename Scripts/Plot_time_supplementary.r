setwd("/cloud-data/snf-mgln-dds/AIDA/Bioinformatics/i0439277/RTcure/Benchmark/Pilot1_Sanger_TT/plot_time_execution/")
library(ggplot2)
library(tidyr)
library(plyr)
library(ggplot2)

time_exec <- read.table("Time_Execution_Methods_TT_plus_Linda.txt",header=T)
#names(time_exec)[8] <- "cov1mln.250k"
long_df <- time_exec
#long_df <- time_exec %>%  gather(coverage, time, cov50k:cov1mln.250k)
#long_df <- time_exec %>%  gather(coverage, time, Method)

head(long_df)
long_df$Time <- as.numeric(long_df$Time)
long_df$Coverage <- factor(long_df$Coverage, levels=c("50K", "100K", "250K", "500K", "750K", "1M", "1M250K"))
long_df$Read_Length <- factor(long_df$Read_Length, levels=c("25bp", "50bp", "75bp", "100bp"))
head(long_df)

long_df_25bp <- subset(long_df, Read_Length  == "100bp")
head(long_df_25bp)

p<-ggplot(long_df_25bp, aes(x=Coverage, y=Time, color=Method)) +
  geom_boxplot()

p1 <- p+scale_y_continuous(trans = log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x)))

p1+ theme(axis.text=element_text(size=18),
          axis.title=element_text(size=18,face="bold"))+  theme(plot.title = element_text(hjust = 0.5))+theme(axis.text=element_text(size=14),axis.title=element_text(size=14,face="bold"))+ggtitle("")+ xlab("") + ylab("Time (Seconds)") +theme(text=element_text(size=15))+theme(axis.text.x = element_text(face="bold", color="#000001")) + theme(axis.text.y = element_text(face="bold", color="#000001"))+theme(text = element_text(size = 18, face="bold"))+theme(axis.title.y = element_text(size = 18))+theme(legend.text = element_text(size = 18,face="plain"))+theme(legend.title = element_text(size = 18))


