rm(list = ls())
library(ggplot2)
library(tidyr)
library(dplyr)
library("RColorBrewer")

sensitivity <- read.table("/root/cloud-data/snf-mgln-dds/AIDA/Bioinformatics/i0439277/RTcure/Benchmark/plot_results/BALDR_BASI_THIS_work_Datasets.tsv",header=T,sep="\t")
sensitivity <- sensitivity %>% mutate(Coverage = str_replace(Coverage, "1mln.250k", "1m250k"))

colnames(sensitivity)[2] <- "newname2"
sensitivity_Canzar <- subset(sensitivity, Dataset == "Canzar_2016")
sensitivity_Canzar <- subset(sensitivity_Canzar, Read_Length > 25)

sensitivity_Upd <- subset(sensitivity, Dataset == "Upadhyay_2018")
sensitivity_Upd <- subset(sensitivity_Upd, Read_Length == 100)

sensitivity_This_Work <- subset(sensitivity, Dataset == "This_Work")
sensitivity_This_Work <- subset(sensitivity_This_Work, Read_Length == 100)



sensitivity <- rbind(sensitivity_Canzar,sensitivity_Upd,sensitivity_This_Work)

sensitivity$Coverage <- factor(sensitivity$Coverage, levels=c("50k", "100k", "250k", "500k", "750k", "1mln", "1m250k"))
sensitivity$Dataset <- factor(sensitivity$Dataset, levels=c("Upadhyay_2018", "Canzar_2016","This_Work"))
sensitivity$Method <- factor(sensitivity$Method, levels=c("BASIC", "MIXCR", "BALDR", "VDJpuzzle","BRACER","TRUST4"))
sensitivity$Chain <- factor(sensitivity$Chain, levels=c("Heavy", "Light"))

sensitivity$Percentage_Assembled <- (sensitivity$Percentage_Assembled)*100
sensitivity$Percentage_Productive <- (sensitivity$Percentage_Productive)*100
sensitivity$Percentage_Assembled <- lapply(sensitivity$Percentage_Assembled, as.integer)
sensitivity$Percentage_Productive <- lapply(sensitivity$Percentage_Productive,as.integer)

sensitivity$Percentage_Assembled <- as.numeric(sensitivity$Percentage_Assembled)
sensitivity$Percentage_Productive <- as.numeric(sensitivity$Percentage_Productive)
sensitivity$Sensitivity <- (sensitivity$Sensitivity)*100



sensitivity_this_work <- subset(sensitivity, Dataset == "Leiden")
sensitivity_Canzar <- subset(sensitivity, Dataset == "Canzar_2016")
sensitivity_Upadhyay <- subset(sensitivity, Dataset == "Upadhyay_2018")


## Plot Fig 2 A
p1 <- ggplot(data = sensitivity_Canzar,
             aes(x = Method, y = Percentage_Assembled,size=Percentage_Productive,shape=Sensitivity)) +
  geom_point(aes(size = Percentage_Productive, fill = Sensitivity), shape = 21) +
  scale_fill_viridis_c()+  scale_size_continuous(range = c(1, 4))


# parsing ok but facet labels are on two different rows
p2 <- p1 +
  facet_grid(Coverage ~ Chain)+ 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))+ylim(c(10,100))+ggtitle(label="Performance tools on heavy and light chains on Canzar Dataset")+theme(plot.title = element_text(hjust = 0.5))+theme(text=element_text(size=15))

p2 +theme(axis.text=element_text(size=14),
          axis.title=element_text(size=14,face="bold"))+  theme(plot.title = element_text(hjust = 0.5))+theme(axis.text=element_text(size=14),axis.title=element_text(size=14,face="bold"))+ggtitle("Canzar Dataset")+ xlab("Methods") + ylab("% Assembled") +theme(text=element_text(size=15))


## Plot Fig 2 B
p3 <- ggplot(data = sensitivity_Upadhyay,
             aes(x = Method, y = Percentage_Assembled,size=Percentage_Productive,shape=Sensitivity)) +
  geom_point(aes(size = Percentage_Productive, fill = Sensitivity), shape = 21) +
  scale_fill_viridis_c()+  scale_size_continuous(range = c(1, 4))

p4 <- p3 +
  facet_grid(Coverage ~ Chain)+ 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))+ylim(c(10,100))+ggtitle(label="Performance tools on heavy and light chains on Upadhyay Dataset")+theme(plot.title = element_text(hjust = 0.5))+theme(text=element_text(size=15))

p4 +theme(axis.text=element_text(size=14),
          axis.title=element_text(size=14,face="bold"))+  theme(plot.title = element_text(hjust = 0.5))+theme(axis.text=element_text(size=14),axis.title=element_text(size=14,face="bold"))+ggtitle("Upadhyay Dataset")+ xlab("Methods") + ylab("% Assembled") +theme(text=element_text(size=15))


## Plot Fig 2 C
p5 <- ggplot(data = sensitivity_this_work,
             aes(x = Method, y = Percentage_Assembled,size=Percentage_Productive,shape=Sensitivity)) +
  geom_point(aes(size = Percentage_Productive, fill = Sensitivity), shape = 21) +
  scale_fill_viridis_c()+  scale_size_continuous(range = c(1, 4))

p6 <- p5 +
  facet_grid(Coverage ~ Chain)+ 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))+ylim(c(10,100))+ggtitle(label="Performance tools on heavy and light chains on our Dataset")+theme(plot.title = element_text(hjust = 0.5))+theme(text=element_text(size=15))

p6 +theme(axis.text=element_text(size=14),
          axis.title=element_text(size=14,face="bold"))+  theme(plot.title = element_text(hjust = 0.5))+theme(axis.text=element_text(size=14),axis.title=element_text(size=14,face="bold"))+ggtitle("Leiden Dataset")+ xlab("Methods") + ylab("% Assembled") +theme(text=element_text(size=15))


##########


