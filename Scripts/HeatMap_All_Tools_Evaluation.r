getwd()

#install.packages("viridis")
library(viridis)

a <- read.table("Final_Ensembled_Average.tsv",header = T)
a$Dataset <- factor(a$Dataset, levels=c("Canzar", "Upadhyay", "This_Work", "Somatic_Hypermutations","Average"))
a$Method <- factor(a$Method, levels=c("VDJpuzzle", "MiXCR", "TRUST4","BRACER","BALDR","BASIC"))

ggplot(a, aes(Dataset,Method, fill=Score )) + geom_tile() +  scale_fill_viridis(discrete=FALSE) + geom_text(aes(label = round(Score, 2))) +  ggtitle("BCR tool performance")+ xlab(" ") + ylab("Method") +theme(text=element_text(size=15))+theme(plot.title = element_text(hjust = 0.5))+theme(text=element_text(size=15))+ theme(axis.text.x = element_text(angle = 15, hjust = 1, vjust = 0.5))+theme(axis.text=element_text(size=14),
                                                                                                                                                                                                                                                                                                                                                axis.title=element_text(size=14,face="bold"))

