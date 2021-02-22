### 

### kwithin v mm

setwd("~/Bioinformatics Work/Meth & RNA/Meta-analysis WGCNA")

library(WGCNA)
library(plyr)
library(flashClust)
library(ggplot2)
library(gplots)
library(lattice)
library(extrafont)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(ggpubr)

loadfonts()
##############


load(file = "MetaAnalysis_trimmed_input.RData")
load(file = "Modules_DS0.RData")
colorsA1 = names(table(modules1))


#####

# load the kwithin vs MM table (Metabric)
kw_mm <- read.table


datKME=signedKME(datExpr1, datME, outputColumnName="MM.")

mm <- read.table("MM_KW_META_BLUE.txt", sep ="\t", header = TRUE, row.names = 1)

####
#Blue
mm <- read.table("MM_KW_META_Blue.txt", sep ="\t", header = TRUE, row.names = 1)
plot <- ggplot(mm, aes(x = (MM_Blue_METABRIC ^ 6), y = Kwithin_METABRIC)) 

plot + 
  theme_bw(base_size = 26) + 
  geom_point(shape = 16, colour = "blue3", size = 4) + 
  geom_smooth(method=lm, colour = "grey", se = TRUE) +
  stat_cor(method = "pearson", label.x = 0.3, label.y = 0.0, size = 8) +
  ggtitle("Blue Module") +
  xlab("Intramodular connectivity") + ylab("Module Membership ^ 6")



#yellow
mm <- read.table("MM_KW_META_Yellow.txt", sep ="\t", header = TRUE, row.names = 1)
plot <- ggplot(mm, aes(x = (MM_Yellow_METABRIC ^ 6), y = Kwithin_METABRIC)) 

plot + 
  theme_bw(base_size = 26) + 
  geom_point(shape = 16, colour = "yellow3", size = 4) + 
  geom_smooth(method=lm, colour = "grey") +
  stat_cor(method = "pearson", label.x = 0.3, label.y = 0.0, size = 8) +
  ggtitle("Yellow Module") +
  xlab("Intramodular Connectivity") + ylab("Module Membership ^ 6")

#Green
mm <- read.table("MM_KW_META_Green.txt", sep ="\t", header = TRUE, row.names = 1)
plot <- ggplot(mm, aes(x = (MM_Green_METABRIC ^ 6), y = Kwithin_METABRIC)) 


plot + 
  theme_bw(base_size = 26) + 
  geom_point(shape = 16, colour = "forestgreen", size = 4) + 
  geom_smooth(method=lm, colour = "grey") +
  stat_cor(method = "pearson", label.x = 0.2, label.y = 0.0, size = 8) +
  ggtitle("Green Module") +
  xlab("Intramodular Connectivity") + ylab("Module Membership ^ 6")


# Black
mm <- read.table("MM_KW_META_Black.txt", sep ="\t", header = TRUE, row.names = 1)
plot <- ggplot(mm, aes(x = (MM_Black_METABRIC ^ 6), y = Kwithin_METABRIC)) 


plot + 
  theme_bw(base_size = 26) + 
  geom_point(shape = 16, colour = "black", size = 4) + 
  geom_smooth(method=lm, colour = "grey") +
  stat_cor(method = "pearson", label.x = 0.3, label.y = 0.0, size = 8) +
  ggtitle("Black Module") +
  xlab("Intramodular Connectivity") + ylab("Module Membership ^ 6")

#Brown
mm <- read.table("MM_KW_META_Brown.txt", sep ="\t", header = TRUE, row.names = 1)
plot <- ggplot(mm, aes(x = (MM_Brown_METABRIC ^ 6), y = Kwithin_METABRIC)) 


plot + 
  theme_bw(base_size = 26) + 
  geom_point(shape = 16, colour = "brown", size = 4) + 
  geom_smooth(method=lm, colour = "grey") +
  stat_cor(method = "pearson", label.x = 0.25, label.y = 0.0, size = 8) +
  ggtitle("Brown Module") +
  xlab("Intramodular Connectivity") + ylab("Module Membership ^ 6")


#Red
mm <- read.table("MM_KW_META_Red.txt", sep ="\t", header = TRUE, row.names = 1)
plot <- ggplot(mm, aes(x = (MM_Red_METABRIC ^ 6), y = Kwithin_METABRIC)) 


plot + 
  theme_bw(base_size = 26) + 
  geom_point(shape = 16, colour = "red3", size = 4) + 
  geom_smooth(method=lm, colour = "grey") +
  stat_cor(method = "pearson", label.x = 0.25, label.y = 0.0, size = 8) +
  ggtitle("Red Module") +
  xlab("Intramodular Connectivity") + ylab("Module Membership ^ 6")

#tan
mm <- read.table("MM_KW_META_Tan.txt", sep ="\t", header = TRUE, row.names = 1)
plot <- ggplot(mm, aes(x = (MM_Tan_METABRIC ^ 6), y = Kwithin_METABRIC)) 

plot + 
  theme_bw(base_size = 26) + 
  geom_point(shape = 16, colour = "tan3", size = 4) + 
  geom_smooth(method=lm, colour = "grey") +
  stat_cor(method = "pearson", label.x = 0.3, label.y = -0.5, size = 8) +
  ggtitle("Tan Module") +
  xlab("Intramodular Connectivity") + ylab("Module Membership ^ 6")

#magenta
mm <- read.table("MM_KW_META_Magenta.txt", sep ="\t", header = TRUE, row.names = 1)
plot <- ggplot(mm, aes(x = (MM_Magenta_METABRIC ^ 6), y = Kwithin_METABRIC)) 


plot + 
  theme_bw(base_size = 26) + 
  geom_point(shape = 16, colour = "magenta3", size = 4) + 
  geom_smooth(method=lm, colour = "grey") +
  stat_cor(method = "pearson", label.x = 0.3, label.y = 0.0, size = 8) +
  ggtitle("Magenta Module") +
  xlab("Intramodular Connectivity") + ylab("Module Membership ^ 6")

