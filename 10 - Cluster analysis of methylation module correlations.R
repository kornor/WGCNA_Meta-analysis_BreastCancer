##script to create correlation table of MEs (from meta-analysis) 
##vs promoter methylation

setwd("~/Bioinformatics Work/Meth & RNA/Meta-analysis WGCNA")

library(WGCNA)
library(plyr)
library(dplyr)
library(flashClust)
library(ggplot2)
library(gplots)
library(ggdendro)
library(pvclust)
library(lattice)
library(extrafont)
library(tidyverse)
library(reshape)
loadfonts()
##############


load(file = "MetaAnalysis_trimmed_input.RData")
load(file = "Modules_DS0.RData")
colorsA1 = names(table(modules1))

### add in the rownames for the ME_ data frames

rownames(ME_1A) <- colnames(datExpr1)
rownames(ME_2A) <- colnames(datExpr2)

## ME_1A is for TCGA
## can remove the METABRIC ones just for clearing the table
##load in the promoter methylation data 
### Read in the meth data (promoters only, all CpG regions)
meth <- read.table("Meth_promoters_final.txt", sep = "\t", header = TRUE, row.names = 1)

## Trim and sort the file to the exp

list <- intersect(rownames(meth), rownames(ME_1A))

meth <- meth[list,]
datExpr1 <- datExpr1[,list]

Me_3 <- ME_1A[list,]


######################################## Files are now ready for work

#################################################  Start with the total sample set
### Perform the Pearson correlation calculations

nGenes = ncol(Me_3);
nSamples = nrow(Me_3);
# Recalculate MEs with color labels
#MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
#MEs = orderMEs(MEs0)

moduleTraitCor = cor(Me_3, meth, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

#########################  Transpose for dendros etc
t_moduleTraitCor <- as.data.frame(t(moduleTraitCor))
t_moduleTraitPvalue <- as.data.frame(t(moduleTraitPvalue))

write.table(t_moduleTraitCor, "All_meth_trait_corr.txt", sep = "\t")
write.table(t_moduleTraitPvalue, "All_meth_trait_pValue.txt", sep = "\t")

#####################################   Creating heatmaps #####################################
sizeGrWindow(10,6)

par(mar = c(6, 8.5, 3, 3));

# Display the values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = NULL,
               yLabels = names(Me_3),
               ySymbols = names(Me_3),
               colorLabels = FALSE,
               colors = blueWhiteRed(10),
               #textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.85,
               zlim = c(-1,1),
               main = paste("High module expression vs trait genes"))




##cormat2 <- read.table("All_meth_trait_corr.txt", sep = "\t", header = TRUE, row.names = 5)
#### 
cormat2 <- moduleTraitCor


pdf("All_promoter_meth_traits.pdf")
par(oma= c(3,1,1,1))

### Info about heatmap.2
### must be done on a matrix not a dataframe
### dendro can be set to none or some
### key gives you a colour key
### srtCol/ srtRow gives an angle for the text
##sepCol, etc gives a cell outline
### RowSideColors - can set colours for marking the cols or rows



### colorpanel in gplots
###redgreen(n)
###redblue(n)
###bluered(n)
###colorpanel(n, low, mid, high)
## or can set the colour panel before
##eg palette <- colorRampPalette(c('#f0f3ff','#0033BB'))(256)


heatmap.2(as.matrix(cormat2),dendrogram="none", 
          notecol="black",col=bluered(200),scale="none",
          labRow = NULL, 
          key=TRUE, keysize=5.5,key.title = NA, key.xlab = "Correlation",
          density.info="none", trace="none",
          cexRow=5.2,cexCol=5.2, srtCol = 75,
          main = "Traits expression \n by module in TCGA dataset")

library(pheatmap)
library(RColorBrewer)

pdf("Meth_Module_correlation_all.pdf",height=8,width=10)
fontsize = 10


pheatmap(t(cormat2), color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(100), 
         kmeans_k = NA, breaks = NA, border_color = "grey60",
         cellwidth = NA, cellheight = NA, scale = "none", 
         cluster_rows = TRUE, cluster_cols = F, fontsize_row = 10,
         annotation_legend = TRUE, drop_levels = TRUE, 
         show_rownames = F,show_colnames = T, 
         treeheight_col = 0,
         main = "Correlation of methylation by module")
dev.off()

################################################################

##save as file that can be directly opened####

save(datExpr1, ME_1A, meth, moduleTraitCor,moduleTraitPvalue, t_moduleTraitCor, t_moduleTraitPvalue,
     file = "Meth_ME_correlations_all_TCGA.RData")

load(file = "Meth_ME_correlations_all_TCGA.RData")


#### Perform hierarchichal clustering on the dataset
## use: dist method euclidean
## clustering method complete

methDist <- dist(t(moduleTraitCor))
methClust <- hclust(methDist, "complete")

plot(methClust,labels = FALSE, main='Methylation of correlations dendro')
## save this dendro file 

#### want 10 clusters
groups_x <- cutree(methClust, k = 10)
count(groups_x)  ## good, there is 10
factor_x <- as.factor(groups_x)


mat_cormat3 <- as.data.frame(t_cormat3)
mat_cormat3$Cluster_allocation <- groups_x
table(mat_cormat3$Cluster_allocation)

write.table(mat_cormat3, "Meth_Mod_Cluster_Correlations.txt", sep = "\t")
######## This splits the data into an array of dataframes based on clusters
cluster_array <- array(split(mat_cormat3), groups_x)
dim(cluster_array)
### write out this?

for (i in 1:dim(cluster_array)){ 
  write.csv(cluster_array[[i]], paste("Clusters_whole_correlation_", i, ".csv", sep="")) 
} 

t_cormat3 <- t(cormat3)

### create a summary of the values from each cluster  in green?
summary_array<-tapply(t_cormat3$MEgreen,factor_x,summary)
complete_summary_array <- do.call(rbind, summary_array)
write.table(complete_summary_array, file = "Cluster_summary_GreenModule.txt", sep = "\t")

### now save this data all as an RData file

save(methClust, methDist, cluster_array, complete_summary_array, file = "Clustering_all_correlations.RData")

save(methClust, methDist, cluster_array, file = "Meth_Dendro_Clusters.RData")

############# validating clusters

library(NbClust)
nb <- NbClust(t_cormat3, distance = "euclidean", min.nc = 2,
              max.nc = 10, method = "complete", index ="all")
# Visualize the result
library(factoextra)
fviz_nbclust(nb) + theme_minimal()

# K-means clustering (k = number of clusters you want)
km.res <- eclust(t_cormat3, "kmeans", k = 10,
                 nstart = 25, graph = FALSE)
# k-means group number of each observation
km.res$cluster
# Visualize k-means clusters
fviz_cluster(km.res, geom = "point", frame.type = "norm")

# PAM clustering
pam.res <- eclust(t_cormat3, "pam", k = 3, graph = FALSE)
pam.res$cluster
# Visualize pam clusters
fviz_cluster(pam.res, geom = "point", frame.type = "norm")

test_clust <- km.res$cluster

frame <- as.matrix(test_clust)
colnames(frame) <- "Cluster"
frame <- as.data.frame(frame)

write.table(frame, "K_means clustering.txt", sep = "\t")

factor <- as.factor(frame$Cluster)

list_c <- read.table("Cluster3_kmeans.txt", sep = "\t", header = TRUE, row.names = 1)

comp <- intersect(rownames(clust4), rownames(list_c)) 

library(cluster)
sil <- silhouette(km.res$cluster, dist(t_cormat3))
head(sil[, 1:3], 10)

# Silhouette plot
plot(sil)

# Summary of silhouette analysis
si.sum <- summary(sil)
# Average silhouette width of each cluster
si.sum$clus.avg.widths
fviz_silhouette(sil)
fviz_silhouette(km.res, palette = "jco", ggtheme = theme_classic())

viz_sil <- fviz_silhouette(km.res, palette = "jco", ggtheme = theme_classic())
viz_sil$layers[[2]]$aes_params$colour <- NA

viz_sil + scale_y_continuous(name="Silhouette width", limits=c(-0.15, 0.80), 
                             breaks=c(0.00,0.20, 0.40, 0.60, 0.80),
                             labels=c("0.00" = "0.00", "0.2" = "0.25","0.4" = "0.50", "0.6" = "0.75", "0.8" = "1.0")) +
    labs(title="Silhouette values for methylation clusters", x = "Cluster number") +
    theme(axis.text.y = element_text(size = 18), plot.title = element_text(size=22, face = "bold"),
     axis.title.y = element_text(size=18, face="bold"), axis.title.x = element_text(size=18, face="bold"), legend.position = "none")
    #geom_text(aes(label = viz_sil$data$cluster))



#################################################################################################
#######       Clustering on the TNBC only!       #########


## Load the TNBC list
tnbc <- read.table("TNBC_jodi.txt", sep = "\t", header = TRUE, row.names = 1)

## Trim the methylation and MEs to this
keeps <- intersect(rownames(tnbc), rownames(meth))
#count = 79
keeps2 <- intersect(rownames(tnbc), rownames(ME_1A))
# count = 167
# discrepancy between numbers : outliers removed during WGNCA, no meth info available

### subset MEs and meth

meth_tnbc <- meth[keeps,]
ME_tnbc <- ME_1A[keeps,]

### Save this for final open as you need

save(meth_tnbc, ME_tnbc, file = "Methylation_modules_TNBC_only.RData")
load( file = "Methylation_modules_TNBC_only.RData")

### correlations
### Perform the Pearson correlation calculations

nGenes = ncol(ME_tnbc);
nSamples = nrow(ME_tnbc);
# Recalculate MEs with color labels
#MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
#MEs = orderMEs(MEs0)

moduleTraitCor = cor(ME_tnbc, meth_tnbc, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(10,6)
par(mar = c(6, 8.5, 3, 3));

# Display the values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = NULL,
               yLabels = names(Me_3),
               ySymbols = names(Me_3),
               colorLabels = FALSE,
               colors = blueWhiteRed(10),
               #textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.85,
               zlim = c(-1,1),
               main = paste("High module expression vs trait genes"))
######
t_moduleTraitCorTNBC <- as.data.frame(t(moduleTraitCor))
t_moduleTraitPvalue_TNBC <- as.data.frame(t(moduleTraitPvalue))

write.table(t_moduleTraitCorTNBC, "TNBC_meth_ME_corr.txt", sep = "\t")
write.table(t_moduleTraitPvalue_TNBC, "TNBC_meth_ME_pValue.txt", sep = "\t")

##cormat2 <- read.table("All_meth_trait_corr.txt", sep = "\t", header = TRUE, row.names = 5)
#### 
cormat2 <- moduleTraitCor


pdf("Methylation_module_correlation_TNBC.pdf")
par(oma= c(3,1,1,1))

### Info about heatmap.2
### must be done on a matrix not a dataframe
### dendro can be set to none or some
### key gives you a colour key
### srtCol/ srtRow gives an angle for the text
##sepCol, etc gives a cell outline
### RowSideColors - can set colours for marking the cols or rows



### colorpanel in gplots
###redgreen(n)
###redblue(n)
###bluered(n)
###colorpanel(n, low, mid, high)
## or can set the colour panel before
##eg palette <- colorRampPalette(c('#f0f3ff','#0033BB'))(256)


heatmap.2(as.matrix(cormat2),dendrogram="none", 
          notecol="black",col=bluered(200),scale="none",
          labRow = NULL, 
          key=TRUE, keysize=5.5,key.title = NA, key.xlab = "Correlation",
          density.info="none", trace="none",
          cexRow=5.2,cexCol=5.2, srtCol = 75,
          main = "Traits expression \n by module in TCGA dataset")

dev.off()
library(pheatmap)
library(RColorBrewer)

pdf("Methylation_module_correlation_TNBC.pdf",height=8,width=10)
fontsize = 10


pheatmap(t(cormat2), color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(100), 
         kmeans_k = NA, breaks = NA, border_color = "grey60",
         cellwidth = NA, cellheight = NA, scale = "none", 
         cluster_rows = TRUE, cluster_cols = F, fontsize_row = 10,
         annotation_legend = TRUE, drop_levels = TRUE, 
         show_rownames = F,show_colnames = T, 
         treeheight_col = 0,
         main = "Methylation of gene promoters by module")
dev.off()



#### Perform hierarchichal clustering on the dataset
## use: dist method euclidean
## clustering method complete

methDist <- dist(t_moduleTraitCorTNBC)
methClust <- hclust(methDist, "complete")

plot(methClust,labels = FALSE, main='Methylation module correlations dendro TNBC')
## save this dendro file 

#### want 14-16 clusters
groups_x <- cutree(methClust, 14)
count(groups_x)  ## good, there is 10
factor_x <- as.factor(groups_x)
count(factor_x)

######## This splits the data into an array of dataframes based on clusters
cluster_array <- array(split(t_moduleTraitCorTNBC, factor_x))
dim(cluster_array)
### write out this?

for (i in 1:dim(cluster_array)){ 
  write.csv(cluster_array[[i]], paste("Clusters_TNBC_correlation_", i, ".csv", sep="")) 
} 


### create a summary of the values from each cluster  in green?
summary_array<-tapply(t_moduleTraitCorTNBC$MEgreen,factor_x,summary)
complete_summary_array <- do.call(rbind, summary_array)
write.table(complete_summary_array, file = "Cluster_summary_GreenModule_TNBC.txt", sep = "\t")

### Sd?

sd_array <- tapply(t_moduleTraitCorTNBC$MEgreen, factor_x, sd)
sd_array
### now save this data all as an RData file

save(methClust, methDist, cluster_array, complete_summary_array, file = "Clustering_TNBC_correlations.RData")



############# validating clusters

library(NbClust)
nb <- NbClust(t_moduleTraitCorTNBC, distance = "euclidean", min.nc = 2,
              max.nc = 10, method = "complete", index ="all")
# Visualize the result
library(factoextra)
fviz_nbclust(nb) + theme_minimal()

# K-means clustering (k = number of clusters you want)
km.res <- eclust(t_moduleTraitCor, "kmeans", k = 2,
                 nstart = 25, graph = FALSE)
# k-means group number of each observation
km.res$cluster
# Visualize k-means clusters
fviz_cluster(km.res, geom = "point", frame.type = "norm")

# PAM clustering
pam.res <- eclust(t_moduleTraitCor, "pam", k = 5, graph = FALSE)
pam.res$cluster
# Visualize pam clusters
fviz_cluster(pam.res, geom = "point", frame.type = "norm")

test_clust <- km.res$cluster

frame <- as.matrix(test_clust)
colnames(frame) <- "Cluster"
frame <- as.data.frame(frame)

write.table(frame, "K_means clustering.txt", sep = "\t")

factor <- as.factor(frame$Cluster)

list_c <- read.table("Cluster7_kmeans.txt", sep = "\t", header = TRUE, row.names = 1)

comp <- intersect(rownames(clust4), rownames(list_c)) 

########### different heir

res.hc <- eclust(t_moduleTraitCor, "hclust", k = 3,
                 method = "complete", graph = FALSE) 
head(res.hc$cluster, 15)
fviz_dend(res.hc, rect = TRUE, show_labels = FALSE) 


#######silhouette analysis

library(cluster)
sil <- silhouette(km.res$cluster, dist(t_cormat3))
head(sil[, 1:3], 10)

# Silhouette plot
plot(sil, main ="Silhouette plot - K-means")

# Summary of silhouette analysis
si.sum <- summary(sil)
# Average silhouette width of each cluster
si.sum$clus.avg.widths


#######################################################################
#### try scaling the data first

medians = apply(t(meth_tnbc),2,median)
mads = apply(t(meth_tnbc),2,mad)
cors.use = scale(t(meth_tnbc),center=medians,scale=mads)


methDist <- dist(cors.use)
methClust <- hclust(methDist, "complete")

plot(methClust, labels = FALSE, main='Methylation of promoters dendro')

### cut at 50?


methDist1 <- dist(t(cors.use))
methClust1 <- hclust(methDist1, "complete")

plot(methClust1, main='Methylation of modules dendro')


#####
pdf("TNBC_Meth_scaled.pdf",height=8,width=10)
fontsize = 10


pheatmap(cors.use, color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(100), 
         kmeans_k = NA, breaks = NA, border_color = "grey60",
         cellwidth = NA, cellheight = NA, scale = "none", 
         cluster_rows = F, cluster_cols = T, fontsize_row = 10,
         annotation_legend = TRUE, drop_levels = TRUE, 
         show_rownames = F,show_colnames = F, 
         treeheight_col = 0,
         main = "Methylation of gene promoters by module")
dev.off()


############# validating clusters

library(NbClust)
nb <- NbClust(cors.use, distance = "euclidean", min.nc = 2,
              max.nc = 10, method = "complete", index ="all")
# Visualize the result
library(factoextra)
fviz_nbclust(nb) + theme_minimal()

# K-means clustering (k = number of clusters you want)
km.res <- eclust(cors.use, "kmeans", k = 7,
                 nstart = 25, graph = FALSE)
# k-means group number of each observation
km.res$cluster
# Visualize k-means clusters
fviz_cluster(km.res, geom = "point", frame.type = "norm")

# PAM clustering
pam.res <- eclust(t_moduleTraitCor, "pam", k = 5, graph = FALSE)
pam.res$cluster
# Visualize pam clusters
fviz_cluster(pam.res, geom = "point", frame.type = "norm")

test_clust <- km.res$cluster

frame <- as.matrix(test_clust)
colnames(frame) <- "Cluster"
frame <- as.data.frame(frame)

write.table(frame, "K_means clustering.txt", sep = "\t")



################ cutting the hclust on the scaled methylation at 50

##cut at 2.5 - two groups

groups_2 <- cutree(methClust, 2.5)
count(groups_2)

list <- rownames(cors.use)[groups_2 == 2]
frame <- as.data.frame(list)
write.table(frame, "List from methylation not cor cluster 2.txt", sep = "\t")

clust2 <- as.data.frame(cors.use[list,])
t_clust2 <- as.data.frame(t(clust2))

summary(t_clust2$SOX10)
sd(t_clust2$SOX10)

summary(t_clust2$PAX8)
sd(t_clust2$PAX8)


clust2$mean <- rowMeans(clust2)
clust2$SD <- 


list2 <- rownames(cors.use)[groups_2 == 1]
frame <- as.data.frame(list2)
write.table(frame, "List from methylation not cor cluster 1.txt", sep = "\t")

clust1 <- cors.use[list2,]
t_clust1 <- as.data.frame(t(clust1))

summary(t_clust1$RANGRF)
sd(t_clust1$RANGRF)

summary(t_clust2$PAX8)
sd(t_clust2$PAX8)


