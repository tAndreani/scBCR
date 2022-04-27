rm(list=ls())
library(reshape2)
library(viridis)
library(ggplot2)

#example for HC - J genes (you can change for any gene type)
#install.packages("reshape2")
matrix <- read.table("Canzar/HC_J.txt")
matrix <- as.matrix(matrix)
melted_cormat <- melt(matrix)
melted_cormat$Var1 <- factor(melted_cormat$Var1, levels=c("VDJpuzzle","MIXCR","TRUST4","BRACER","BALDR","BASIC"))
melted_cormat$Var2 <- factor(melted_cormat$Var2, levels=c("VDJpuzzle","MIXCR","TRUST4","BRACER","BALDR","BASIC"))



# Get lower triangle of the correlation matrix
get_lower_tri<-function(matrix){
  matrix[upper.tri(matrix)] <- NA
  return(matrix)
}


lower_tri <- get_lower_tri(matrix)
lower_tri

# Melt the correlation matrix
library(reshape2)
melted_cormat <- melt(lower_tri, na.rm = TRUE)
head(melted_cormat)
melted_cormat$Var1 <- factor(melted_cormat$Var1, levels=c("VDJpuzzle","MIXCR","TRUST4","BRACER","BALDR","BASIC"))
melted_cormat$Var2 <- factor(melted_cormat$Var2, levels=c("VDJpuzzle","MIXCR","TRUST4","BRACER","BALDR","BASIC"))

# Heatmap
library(ggplot2)

# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_viridis(discrete=FALSE)+
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()

# Print the heatmap
print(ggheatmap)



ggheatmap + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank())+
   guides(fill = guide_colorbar(barwidth = 7, barheight = ,
                               title.position = "top", title.hjust = 0.5))+ theme(plot.title = element_text(size=10))+ theme(text=element_text(size=10))+theme(axis.text.x = element_text(face="bold", color="#000001")) + theme(axis.text.y = element_text(face="bold", color="#000001"))+theme(text = element_text(size = 14, face="bold"))+theme(axis.title.y = element_text(size = 14))+theme(legend.text = element_text(size = 14,face="plain"))+theme(legend.title = element_text(size = 14))+scale_y_discrete(position = "right")

