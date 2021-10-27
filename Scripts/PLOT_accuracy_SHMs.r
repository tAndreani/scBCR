rm(list = ls())
library(ggplot2)
library(tidyr)
library(dplyr)
library("RColorBrewer")

accuracy <- read.table("/root/cloud-data/snf-mgln-dds/AIDA/Bioinformatics/i0439277/RTcure/Benchmark/simulation_mutations/statistical_analysis_corrected/plot/Master_Table_shm.tsv",header=T,sep="\t")
head(accuracy)
colnames(accuracy)[1] <- "SHMs"
accuracy$SHMs <- factor(accuracy$SHMs, levels=c("shm15", "shm30", "shm45", "shm60"))
accuracy$Chain <- factor(accuracy$Chain, levels=c("Heavy", "Light_Kappa","Light_Lambda"))
accuracy$Method <- factor(accuracy$Method, levels=c("BASIC", "MIXCR", "BALDR", "vdjpuzzle","BRACER","TRUST4"))
accuracy$Percentage_Assembled <- (accuracy$Assembled/accuracy$Total_chains)*100
accuracy$Percentage_Productive <- (accuracy$Productive/accuracy$Total_chains)*100
accuracy$Accuracy <- as.numeric(accuracy$Accuracy)*100
accuracy$Percentage_Assembled <- as.numeric(accuracy$Percentage_Assembled)
accuracy$Percentage_Productive <- as.numeric(accuracy$Percentage_Productive)

p1 <- ggplot(data = accuracy,
             aes(x = Method, y = Percentage_Assembled,size=Percentage_Productive,shape=Accuracy)) +
  geom_point(aes(size = Percentage_Productive, fill = Accuracy), shape = 21) +
  scale_fill_viridis_c() +
  scale_size_continuous(range = c(3, 3))

p1
# parsing ok but facet labels are on two different rows
p2 <- p1 +
  facet_grid(SHMs ~ Chain)+ 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))+ggtitle(label="Simulated Chains with different levels of SHMs")+theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(size = Percentage_Assembled, fill = Accuracy), shape = 21) +
  scale_fill_viridis_c()+  scale_size_continuous(range = c(1, 4))


p2 +theme(axis.text=element_text(size=14),
         axis.title=element_text(size=14,face="bold"))+  theme(plot.title = element_text(hjust = 0.5))+theme(axis.text=element_text(size=14),
                                                                                                             axis.title=element_text(size=14,face="bold"))+ggtitle("Simulated Chains with different SHMs")+theme(text=element_text(size=20)) 
p2 + theme(legend.title = element_text(size = 8))
p2 + theme(legend.key.size = unit(0.2, "cm"))

p2 +theme(axis.text=element_text(size=14),
          axis.title=element_text(size=14,face="bold"))+  theme(plot.title = element_text(hjust = 0.5))+theme(axis.text=element_text(size=14),axis.title=element_text(size=14,face="bold"))+ggtitle("Simulated Chains with different SHMs")+ ylab("% Assembled") +theme(text=element_text(size=15))




