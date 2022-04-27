#Reproduction figure 2

rm(list = ls())
library(ggplot2)
library(tidyr)
library(dplyr)
library(RColorBrewer)


#input file is within this folder
sensitivity <- read.table("BALDR_BASIC_LEIDEN.txt",header=T,sep="\t")
head(sensitivity)

colnames(sensitivity)[2] <- "newname2"
sensitivity_Canzar <- subset(sensitivity, Dataset == "Canzar")
sensitivity_Canzar <- subset(sensitivity_Canzar, Read_Length > 25)

sensitivity_Upd <- subset(sensitivity, Dataset == "Upadhyay")
sensitivity_Upd <- subset(sensitivity_Upd, Read_Length == 100)

sensitivity_Leiden <- subset(sensitivity, Dataset == "Leiden")
sensitivity_Leiden <- subset(sensitivity_Leiden, Read_Length == 100)



sensitivity <- rbind(sensitivity_Canzar,sensitivity_Upd,sensitivity_Leiden)

sensitivity$Coverage <- factor(sensitivity$Coverage, levels=c("1M250K", "1M","750K","500K", "250K","100K","50K"))
sensitivity$Dataset <- factor(sensitivity$Dataset, levels=c("Upadhyay", "Canzar","Leiden"))
sensitivity$Method <- factor(sensitivity$Method, levels=c("BASIC", "MIXCR", "BALDR", "VDJpuzzle","BRACER","TRUST4"))
sensitivity$Chain <- factor(sensitivity$Chain, levels=c("Heavy", "Light"))

sensitivity$Percentage_Assembled <- (sensitivity$Percentage_Assembled)*100
sensitivity$Percentage_Productive <- (sensitivity$Percentage_Productive)*100
sensitivity$Percentage_Assembled <- lapply(sensitivity$Percentage_Assembled, as.integer)
sensitivity$Percentage_Productive <- lapply(sensitivity$Percentage_Productive,as.integer)

sensitivity$Percentage_Assembled <- as.numeric(sensitivity$Percentage_Assembled)
sensitivity$Percentage_Productive <- as.numeric(sensitivity$Percentage_Productive)
sensitivity$Sensitivity <- (sensitivity$Sensitivity)*100
head(sensitivity)

names(sensitivity)[3] <- "Assembled_2"
names(sensitivity)[4] <- "Assembled"
names(sensitivity)[5] <- "Productive_2"
names(sensitivity)[6] <- "Productive"
head(sensitivity)

sensitivity_Leiden <- subset(sensitivity, Dataset == "Leiden")
sensitivity_Canzar <- subset(sensitivity, Dataset == "Canzar")
sensitivity_Upadhyay <- subset(sensitivity, Dataset == "Upadhyay")


## Plot Fig 2 A
p1 <- ggplot(data = sensitivity_Canzar,
             aes(x = Method, y = Assembled,size=Productive,shape=Sensitivity)) +
  geom_point(aes(size = Productive, fill = Sensitivity), shape = 21) +
  scale_fill_viridis_c()+  scale_size_continuous(range = c(1, 4))
p1 + theme(axis.text.x = element_text(face="bold", color="#000001", size=14, angle=45))
p1

p2 <- p1 +
  facet_grid(Coverage ~ Chain)+ 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))+ylim(c(25,100))+ggtitle(label="Performance tools on heavy and light chains on Canzar Dataset")+theme(plot.title = element_text(hjust = 0.5))+theme(text=element_text(size=15))



p2 +theme(axis.text=element_text(size=18),
          axis.title=element_text(size=18,face="bold"))+  theme(plot.title = element_text(hjust = 0.5))+theme(axis.text=element_text(size=14),axis.title=element_text(size=14,face="bold"))+ggtitle("Canzar Dataset")+ xlab("") + ylab(" Assembled") +theme(text=element_text(size=15))+theme(axis.text.x = element_text(face="bold", color="#000001")) + theme(axis.text.y = element_text(face="bold", color="#000001"))+theme(text = element_text(size = 18, face="bold"))+theme(axis.title.y = element_text(size = 18))+theme(legend.text = element_text(size = 18,face="plain"))+theme(legend.title = element_text(size = 18))



## Plot Fig 2 B

p3 <- ggplot(data = sensitivity_Upadhyay,
             aes(x = Method, y = Assembled,size=Productive,shape=Sensitivity)) +
  geom_point(aes(size = Productive, fill = Sensitivity), shape = 21) +
  scale_fill_viridis_c()+  scale_size_continuous(range = c(1, 4))

p4 <- p3 +
  facet_grid(Coverage ~ Chain)+ 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))+ylim(c(25,100))+ggtitle(label="Performance tools on heavy and light chains on Upadhyay Dataset")+theme(plot.title = element_text(hjust = 0.5))+theme(text=element_text(size=15))

p4 +theme(axis.text=element_text(size=14),
          axis.title=element_text(size=18,face="bold"))+  theme(plot.title = element_text(hjust = 0.5))+theme(axis.text=element_text(size=14),axis.title=element_text(size=14,face="bold"))+ggtitle("Upadhyay Dataset")+ xlab("") + ylab(" Assembled") +theme(text=element_text(size=15))+theme(axis.text.x = element_text(face="bold", color="#000001")) + theme(axis.text.y = element_text(face="bold", color="#000001"))+theme(text = element_text(size = 18, face="bold"))+theme(axis.title.y = element_text(size = 18))+theme(legend.text = element_text(size = 18,face="plain"))+theme(legend.title = element_text(size = 18)) 


## Plot Fig 2 C
p5 <- ggplot(data = sensitivity_Leiden,
             aes(x = Method, y = Assembled,size=Productive,shape=Sensitivity)) +
  geom_point(aes(size = Productive, fill = Sensitivity), shape = 21) +
  scale_fill_viridis_c()+  scale_size_continuous(range = c(1, 4))

p6 <- p5 +
  facet_grid(Coverage ~ Chain)+ 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))+ylim(c(25,100))+ggtitle(label="Performance tools on heavy and light chains on our Dataset")+theme(plot.title = element_text(hjust = 0.5))+theme(text=element_text(size=15))

p6 +theme(axis.text=element_text(size=14),
          axis.title=element_text(size=18,face="bold"))+  theme(plot.title = element_text(hjust = 0.5))+theme(axis.text=element_text(size=14),axis.title=element_text(size=14,face="bold"))+ggtitle("Leiden Dataset")+ xlab("") + ylab(" Assembled") +theme(text=element_text(size=15))+theme(axis.text.x = element_text(face="bold", color="#000001")) + theme(axis.text.y = element_text(face="bold", color="#000001"))+theme(text = element_text(size = 18, face="bold"))+theme(axis.title.y = element_text(size = 18))+theme(legend.text = element_text(size = 18,face="plain"))+theme(legend.title = element_text(size = 18)) 


