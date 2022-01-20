rm(list = ls())
library(ggplot2)
library(tidyr)
library(dplyr)
library("RColorBrewer")

#input file is within this folder
accuracy <- read.table("Master_Table_shm.tsv",header=T,sep="\t")
head(accuracy)
colnames(accuracy)[1] <- "SHMs"
accuracy$SHMs <- factor(accuracy$SHMs, levels=c("SHM15", "SHM30", "SHM45", "SHM60"))
accuracy$Chain <- factor(accuracy$Chain, levels=c("Heavy", "Light Kappa","Light Lambda"))
accuracy$Method <- factor(accuracy$Method, levels=c("BASIC", "MIXCR", "BALDR", "VDJpuzzle","BRACER","TRUST4"))
accuracy$Percentage_Assembled <- (accuracy$Assembled/accuracy$Total_chains)*100
accuracy$Percentage_Productive <- (accuracy$Productive/accuracy$Total_chains)*100
accuracy$Accuracy <- as.numeric(accuracy$Accuracy)*100
accuracy$Percentage_Assembled <- as.numeric(accuracy$Percentage_Assembled)
accuracy$Percentage_Productive <- as.numeric(accuracy$Percentage_Productive)

names(accuracy)[3] <- "Assembled_2"
names(accuracy)[4] <- "Assembled"
names(accuracy)[5] <- "Productive_2"
names(accuracy)[6] <- "Productive"
p1 <- ggplot(data = accuracy,
             aes(x = Method, y = Assembled,size=Productive,shape=Accuracy)) +
  geom_point(aes(size = Productive, fill = Accuracy), shape = 21) +
  scale_fill_viridis_c() +
  scale_size_continuous(range = c(3, 3))

p1
# parsing ok but facet labels are on two different rows
p2 <- p1 +
  facet_grid(SHMs ~ Chain)+ 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))+ggtitle(label="Simulated Chains with different levels of SHMs")+theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(size = Assembled, fill = Accuracy), shape = 21) +
  scale_fill_viridis_c()+  scale_size_continuous(range = c(1, 4))


p2 +theme(axis.text=element_text(size=14),
          axis.title=element_text(size=14,face="bold"))+  theme(plot.title = element_text(hjust = 0.5))+theme(axis.text=element_text(size=14),
                                                                                                              axis.title=element_text(size=14,face="bold"))+ggtitle("Simulated Chains with different SHMs")+theme(text=element_text(size=20)) 
p2 + theme(legend.title = element_text(size = 8))
p2 + theme(legend.key.size = unit(0.2, "cm"))

p2 +theme(axis.text=element_text(size=14),
    axis.title=element_text(size=18,face="bold"))+  theme(plot.title = element_text(hjust = 0.5))+theme(axis.text=element_text(size=14),axis.title=element_text(size=14,face="bold"))+ggtitle("SHMs Dataset")+ xlab("") + ylab(" Assembled") +theme(text=element_text(size=15))+theme(axis.text.x = element_text(face="bold", color="#000001")) + theme(axis.text.y = element_text(face="bold", color="#000001"))+theme(text = element_text(size = 18, face="bold"))+theme(axis.title.y = element_text(size = 18))+theme(legend.text = element_text(size = 18,face="plain"))+theme(legend.title = element_text(size = 18)) 




