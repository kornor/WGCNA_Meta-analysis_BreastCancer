
############################
###########   Script to calculate centrality measures for between module communciation

setwd("~/Bioinformatics Work/Meth & RNA/Meta-analysis WGCNA")

library(WGCNA)
library(plyr)
library(flashClust)
library(ggplot2)
library(gplots)
library(lattice)
library(extrafont)
library(reshape2)
library(dplyr)
library(broom)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(igraph)
loadfonts()
##############


load(file = "MetaAnalysis_trimmed_input.RData")
load(file = "Modules_DS0.RData")
colorsA1 = names(table(modules1))

### add in the rownames for the ME_ data frames

rownames(ME_1A) <- colnames(datExpr1)
rownames(ME_2A) <- colnames(datExpr2)


## need to set the rownames for ME_1A
rownames(ME_1A) <- colnames(datExpr1)


clin1 <- read.table("Clinical_final_TCGA.txt", sep = "\t",
                   header = TRUE, row.names = 1)


#### read in the clinical file for metabric
clin2 <- read.table("Clinical_final_METABRIC.txt", sep = "\t",
                   header = TRUE, row.names = 1)




datTraits1 <- read.table("Traits_trimmed_TCGA.txt", sep = "\t",
                         header = TRUE, row.names = 1)

list <- intersect(rownames(ME_1A), rownames(datTraits1))

datTraits1 <- datTraits1[list,]
ME_1trim <- ME_1A[list,]



# Define numbers of genes and samples
nGenes = ncol(datExpr1);
nSamples = nrow(datExpr1);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr1, modules1)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

##
moduleTraitCor = cor(ME_1trim, datTraits1, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits1),
               yLabels = names(ME_1trim),
               ySymbols = names(ME_1trim),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


# Define variable weight containing the weight column of datTrait
CIN = as.data.frame(datTraits1$CIN.metagene);
names(CIN) = "CIN Metagene"

datExpr <- as.data.frame(t(datExpr1))
datExpr <- datExpr[list,]
# names (colors) of the modules
modNames = substring(names(ME_1trim), 3)
geneModuleMembership = as.data.frame(cor(datExpr, ME_1trim, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));


names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, , use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(CIN), sep="");
names(GSPvalue) = paste("p.GS.", names(CIN), sep="");

module = "green"
column = match(module, modNames);
moduleGenes = modules1==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for CIN Metagene",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module = "blue"
column = match(module, modNames);
moduleGenes = modules1==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for CIN Metagene",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
