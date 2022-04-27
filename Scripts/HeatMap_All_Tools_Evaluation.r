##Reproduction Figure 7

getwd()
library(viridis)


#input file is within this folder
a <- read.table("Final_Ensembled_Average.tsv",header = T, sep="\t")
a
a$Dataset <- factor(a$Dataset, levels=c("Leiden","Canzar", "Upadhyay", "SHMs","Average"))
a$Method <- factor(a$Method, levels=c("BASIC","BALDR","BRACER","TRUST4","MiXCR","VDJpuzzle"))
a$Method <- factor(a$Method, levels=c("VDJpuzzle","MiXCR","TRUST4","BRACER","BALDR","BASIC"))

head(a)
ggplot(a, aes(Dataset,Method, fill=Score )) + geom_tile() +  scale_fill_viridis(discrete=FALSE) + geom_text(aes(label=round(Score,2)), size = 7) +  ggtitle(label="Methods performance")+theme(plot.title = element_text(hjust = 0.5))+theme(text=element_text(size=25))+  theme(axis.text = element_text(size = 20)) +
  labs(x = NULL, y = NULL)+ theme(plot.title = element_text(size=35))+ theme(text=element_text(size=15))+theme(axis.text.x = element_text(face="bold", color="#000001")) + theme(axis.text.y = element_text(face="bold", color="#000001"))+theme(text = element_text(size = 18, face="bold"))+theme(axis.title.y = element_text(size = 18))+theme(legend.text = element_text(size = 18,face="plain"))+theme(legend.title = element_text(size = 18))  


#modality 2
ggplot(a, aes(Dataset,Method, fill=Score )) + geom_tile() +  scale_fill_viridis(discrete=FALSE) + geom_text(aes(label=round(Score,2)), size = 5) +  ggtitle(label=" ")+theme(plot.title = element_text(hjust = 0.5))+theme(text=element_text(size=10))+  theme(axis.text = element_text(size = 14)) +
  labs(x = NULL, y = NULL)+ theme(plot.title = element_text(size=10))+ theme(text=element_text(size=10))+theme(axis.text.x = element_text(face="bold", color="#000001")) + theme(axis.text.y = element_text(face="bold", color="#000001"))+theme(text = element_text(size = 14, face="bold"))+theme(axis.title.y = element_text(size = 14))+theme(legend.text = element_text(size = 14,face="plain"))+theme(legend.title = element_text(size = 14))  


