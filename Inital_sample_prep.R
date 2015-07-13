setwd("~/Bioinformatics Work/Meth & RNA/WGCNA_meta")

### This analysis includes all the genes for mRNA analysis of samples
## METABRIC VALIDATION SET

### Load packages
library(WGCNA)
library(flashClust)
library(permute)
library(lattice)
library(vegan)
library(scatterplot3d)
library(MASS)
library(nlme)
library(mgcv)
library(cluster)
library(caret)
library(ape)
library(sparcl)
library(plyr)
library(ggplot2)
library(gplots)
library(RColorBrewer)

## #########################################  This can be skipped - move down to load prepped file
### Load exp file & traits & clinical

exp <- read.table("Exp_final_METABRIC.txt", sep = "\t", header = TRUE, row.names = 1)
exp <- as.data.frame(t(exp))

clin <- read.table("Clinical_METABRIC.txt", sep = "\t", header = TRUE, row.names = 1)

traits <- read.table("Traits_METABRIC.txt", sep = "\t", header = TRUE, row.names = 1)

## look for genes with missing values
##Exclude genes with results for less than 400 samples; exclude samples with 
#results for less than 20 000 genes. 
gsg = goodSamplesGenes(exp,verbose = 5);
gsg$allOK


## if return is true, all good to go
### Otherwise

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(exp)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(exp)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  exp = exp[gsg$goodSamples, gsg$goodGenes]
}


## Create a dendro of the samples to look for outliers

sampleTree = flashClust(dist(exp), method = "average");

#### Make factor labels for the subtypes
Pam50 <- as.factor(clin$Pam50_subtype)

# Relabel blockwise modules

PamColours = labels2colors(Pam50)
count(PamColours)

# Plot the sample tree

sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, 
     main = "Sample clustering to detect outliers", 
     sub="", xlab="", 
     cex.lab = 1.5,
     cex.axis = 1.5, 
     cex.main = 2, labels = FALSE)

ColorDendrogram(sampleTree, y = PamColours, 
                main = "Sample clustering to detect outliers", 
                sub="", xlab="", 
                labels = FALSE, branchlength = 35)

legend("bottomright",legend=levels(factor(Pam50)),
       fill = (c("turquoise","blue", "brown", "yellow", "green" )), cex = 1)
##########   Different way?

plotDendroAndColors(sampleTree, PamColours,
                    groupLabels = "Pam50 subtype",
                    main = "Sample dendrogram and subtype", dendroLabels = FALSE)
legend("right",legend=levels(factor(Pam50)),
       fill = (c("turquoise","blue", "brown", "yellow", "green" )), cex = 1)



## Might remove the 4 outliers on the right
## Plot a line to show the cut

abline(h = 210, col = "red");

# Determine clusters under the line (use h from above)

#clust = cutreeStatic(sampleTree, cutHeight = 300, minSize = 10)
clust = cutreeStatic(sampleTree, cutHeight = 210, minSize = 10)
table(clust)

# clust 1 contains the samples we want to keep.

keepSamples = (clust==1)
datExpr = exp[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


##Traits

Samples = rownames(datExpr);
traitRows = match(Samples, rownames(traits));
datTraits = methTraits[traitRows,];


datClin = clin[traitRows,]

collectGarbage();



###################

# Re-cluster samples to look for underlying patterns

sampleTree2 = flashClust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, 
#grey means missing entry
# need to consider how to do this for discrete variables?

traitColors = numbers2colors(datTraits, 
                             signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap", dendroLabels = FALSE)



#### Note that the side panel - the "basal" group, is hypomethylated across all (comparatively)


### Save out the data and proceed to soft thresholding script ("SingleBlock_moduleCreation.R")

save(datExpr, datTraits, datClin, file = "METABRIC_trimmed_input.RData")



