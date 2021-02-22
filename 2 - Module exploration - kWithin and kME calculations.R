
### This script involves calculating the Module Membership values (kME) 
### as well as the intramodular connectivity (kWithin) 
### And then plotting kME vs kWithin for examination



setwd("~/Bioinformatics Work/Meth & RNA/Meta-analysis WGCNA")

library(WGCNA)
library(plyr)
library(flashClust)
library(ggplot2)
library(gplots)
library(lattice)
library(extrafont)
loadfonts()
##############


load(file = "MetaAnalysis_trimmed_input.RData")
load(file = "Modules_DS0.RData")
colorsA1 = names(table(modules1))

### add in the rownames for the ME_ data frames

rownames(ME_1A) <- colnames(datExpr1)
rownames(ME_2A) <- colnames(datExpr2)


### kME - module membership values

geneModuleMembership1 = signedKME(t(datExpr1), ME_1A)
colnames(geneModuleMembership1)=paste("PC",colorsA1,".cor",sep=""); 

MMPvalue1=corPvalueStudent(as.matrix(geneModuleMembership1),dim(datExpr1)[[2]]); 
colnames(MMPvalue1)=paste("PC",colorsA1,".pval",sep="");

Gene       = rownames(datExpr1)
kMEtable1  = cbind(Gene,Gene,modules1)
for (i in 1:length(colorsA1))
  kMEtable1 = cbind(kMEtable1, geneModuleMembership1[,i], MMPvalue1[,i])
colnames(kMEtable1)=c("PSID","Gene","Module",
                      sort(c(colnames(geneModuleMembership1), colnames(MMPvalue1))))

#write.csv(kMEtable1,"kMEtable1.csv",row.names=FALSE)

###  Now repeat for METABRIC, using the module assignments from TCGA to determine kME values.

# First calculate MEs for A2, since we haven't done that yet
PCs2A = moduleEigengenes(t(datExpr2),  colors=modules1) 
ME_2A = PCs2A$eigengenes

geneModuleMembership2 = signedKME(t(datExpr2), ME_2A)
colnames(geneModuleMembership1)=paste("PC",colorsA1,".cor",sep=""); 

MMPvalue2=corPvalueStudent(as.matrix(geneModuleMembership2),
                           dim(datExpr2)[[2]]); 
colnames(MMPvalue2)=paste("PC",colorsA1,".pval",sep="");

kMEtable2  = cbind(Gene,Gene,modules1)
for (i in 1:length(colorsA1))
  kMEtable2 = cbind(kMEtable2, geneModuleMembership2[,i], MMPvalue2[,i])
colnames(kMEtable2)=colnames(kMEtable1)

#write.csv(kMEtable2,"kMEtable2.csv",row.names=FALSE)

#pdf("all_kMEtable2_vs_kMEtable1.pdf",height=8,width=8)
#for (c in 1:length(colorsA1)){
#  verboseScatterplot(geneModuleMembership2[,c],geneModuleMembership1[,c],main=colorsA1[c],
#                     xlab="kME in A2",ylab="kME in A1")
#}; dev.off()


###############  with some changes for colour etc etc


pdf("inMod_kMEtable2_vs_kMEtable1.pdf",height=8,width=8)
par(mar = c(5,5,5,5))
for (c in 1:length(colorsA1)){
  inMod = modules1== colorsA1[c]
  verboseScatterplot(geneModuleMembership2[inMod,c],geneModuleMembership1[inMod,c],
                     main=colorsA1[c],
                     xlab="Module Membership in METABRIC",
                     ylab="Module Membership in TCGA",
                     col = colorsA1[c],
                     abline = TRUE,
                     pch = 19)

}; dev.off()

mm_scatter <- ggplot(data(y = (geneModuleMembership2[inMod,c] x = geneModuleMembership1[inMod,c]),
                     aes(main=colorsA1[c]))

                     
                     
#### beautiful!!! 

### save out some stuff
#save(geneTree1, geneTree2, modules1, ME_1A, ME_2A, file = "Modules_DS0.RData")

####  NOW determine which genes are hubs in both networks

### rank the module memberships so that LOWEST RANK = HIGHEST MM
### => numbered 1 - 12588
### cbind the two ranked lists side by side
### then find the max value between the two 
## and rank the genes by the max value
## then sort the genes into their modules
## and you'll have a list of the highest MM genes (picked from both datasets)

topGenesKME = NULL
for (c in 1:length(colorsA1)){
  kMErank1    = rank(-geneModuleMembership1[,c])
  kMErank2    = rank(-geneModuleMembership2[,c])
  maxKMErank  = rank(apply(cbind(kMErank1,kMErank2+.00001),1,max))
  topGenesKME = cbind(topGenesKME,Gene[maxKMErank<=100])
}; colnames(topGenesKME) = colorsA1
topGenesKME

#write.table(topGenesKME, "Top100_all_modules_both_networks.txt", sep = "\t")

### This is the top 100 most interconnected genes by module (in both networks)



#################

### Calculate adjacency of genes, and their intramodular connectivity (kWithin)

ADJ1=abs(cor(t(datExpr2),use="p"))^9
Alldegrees1=intramodularConnectivity(ADJ1, modules1)
head(Alldegrees1)

intra <- Alldegrees1
intra$Module <- modules1

### write out the kWithin table

write.table(intra, "kWithin for all modules_TCGA.txt", sep = "\t")


modcol <- as.data.frame(modules1)
rownames(modcol) <- rownames(datExpr1)

###########  this is a script to get MM (signedKME) and intramodular connectivity
#### and plot against each other

sizeGrWindow(8,6)
par(mfrow=c(2,2))
# We choose 4 modules to plot: yellow, blue, brown, green.
# For simplicity we write the code out explicitly for each module.

#### set the genes to the colour module want
which.col = as.factor(modcol$modules1);
whichGenes <- subset(modcol, which.col == "yellow")
restrictGenes <- rownames(whichGenes)

### now set up the kME for this
## Didn't set a GS1 ??
FilterGenes= abs(GS1)> .2 & abs(kMEtable1$PCgreen.cor)>.8
table(FilterGenes)