
setwd("~/Bioinformatics Work/Meth & RNA/Meta-analysis WGCNA")

library(dplyr)
library(ggpubr)

write.table(geneModuleMembership1, "MM_TCGA.txt", sep = "\t")
write.table(geneModuleMembership2, "MM_METABRIC.txt", sep = "\t")


MM_both <- read.table("MM_both.txt", sep = "\t", header = TRUE, row.names = 1)


#Blue
plot <- ggplot(MM_both, aes(x = Blue_METABRIC, y = Blue_TCGA)) 

plot + 
  theme_bw(base_size = 26) + 
  geom_point(shape = 16, colour = "blue3") + 
  geom_smooth(method=lm, colour = "grey") +
  stat_cor(method = "pearson", label.x = 0.2, label.y = 1.5, size = 8) +
  ggtitle("Mod-blue") +
  xlab("METABRIC") + ylab("TCGA")
       
cor.test(MM_both)

res <- cor.test(MM_both$Blue_METABRIC, MM_both$Blue_TCGA, 
                method = "pearson")
res

#yellow
plot <- ggplot(MM_both, aes(x = Yellow_METABRIC, y = Yellow_TCGA)) 

plot + 
  theme_bw(base_size = 26) + 
  geom_point(shape = 16, colour = "yellow3") + 
  geom_smooth(method=lm, colour = "grey") +
  stat_cor(method = "pearson", label.x = 0.2, label.y = 1.5, size = 8) +
  ggtitle("Mod-yellow") +
  xlab("METABRIC") + ylab("TCGA")

#Green
plot <- ggplot(MM_both, aes(x = Green_METABRIC, y = Green_TCGA)) 

plot + 
  theme_bw(base_size = 26) + 
  geom_point(shape = 16, colour = "forestgreen") + 
  geom_smooth(method=lm, colour = "grey") +
  stat_cor(method = "pearson", label.x = 0.2, label.y = 1.5, size = 8) +
  ggtitle("Mod-green") +
  xlab("METABRIC") + ylab("TCGA")


# Black
plot <- ggplot(MM_both, aes(x = Black_METABRIC, y = Black_TCGA)) 

plot + 
  theme_bw(base_size = 26) + 
  geom_point(shape = 16, colour = "black") + 
  geom_smooth(method=lm, colour = "grey") +
  stat_cor(method = "pearson", label.x = 0.2, label.y = 1.5, size = 8) +
  ggtitle("Mod-black") +
  xlab("METABRIC") + ylab("TCGA")

#Browm
plot <- ggplot(MM_both, aes(x = Brown_METABRIC, y = Brown_TCGA)) 

plot + 
  theme_bw(base_size = 26) + 
  geom_point(shape = 16, colour = "brown") + 
  geom_smooth(method=lm, colour = "grey") +
  stat_cor(method = "pearson", label.x = 0.2, label.y = 1.5, size = 8) +
  ggtitle("Mod-brown") +
  xlab("METABRIC") + ylab("TCGA")


#Red
plot <- ggplot(MM_both, aes(x = Red_METABRIC, y = Red_TCGA)) 

plot + 
  theme_bw(base_size = 26) + 
  geom_point(shape = 16, colour = "red3") + 
  geom_smooth(method=lm, colour = "grey") +
  stat_cor(method = "pearson", label.x = 0.2, label.y = 1.5, size = 8) +
  ggtitle("Mod-red") +
  xlab("METABRIC") + ylab("TCGA")

#tan
plot <- ggplot(MM_both, aes(x = Tan_METABRIC, y = Tan_TCGA)) 

plot + 
  theme_bw(base_size = 26) + 
  geom_point(shape = 16, colour = "tan3") + 
  geom_smooth(method=lm, colour = "grey") +
  stat_cor(method = "pearson", label.x = 0.2, label.y = 1.5, size = 8) +
  ggtitle("Mod-tan") +
  xlab("METABRIC") + ylab("TCGA")

#magenta

plot <- ggplot(MM_both, aes(x = Magenta_METABRIC, y = Magenta_TCGA)) 

plot + 
  theme_bw(base_size = 26) + 
  geom_point(shape = 16, colour = "magenta3") + 
  geom_smooth(method=lm, colour = "grey") +
  stat_cor(method = "pearson", label.x = 0.2, label.y = 1.5, size = 8) +
  ggtitle("Mod-magenta") +
  xlab("METABRIC") + ylab("TCGA")

