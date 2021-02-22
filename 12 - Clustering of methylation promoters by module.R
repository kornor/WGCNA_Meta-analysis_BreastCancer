##script to cluster the beta-value methylation of promoters by module


setwd("~/Bioinformatics Work/Meth & RNA/Meta-analysis WGCNA")

library(WGCNA)
library(plyr)
library(dplyr)
library(flashClust)
library(ggplot2)
library(gplots)
library(pvclust)
library(lattice)
library(extrafont)
loadfonts()
##############


load(file = "MetaAnalysis_trimmed_input.RData")

### we are going to use the methylation beta values
load(file = "Meth_ME_correlations_all_TCGA.RData")
### Read in the meth data (promoters only, all CpG regions)


#### scaling the data first

medians <- apply(t(meth),2,median)
mads <- apply(t(meth), 2, mad)
meth.all.scaled <- as.data.frame(scale(t(meth), center = medians, scale = mads))

methDist <- dist(meth.all.scaled)
methClust <- hclust(methDist, "complete")
sizeGrWindow(10,6)
plot(methClust, labels = FALSE, main='Methylation of promoters dendro')


#### heatmap - this is a mess
par(mar = c(6, 8.5, 3, 3));

# Display the values within a heatmap plot
labeledHeatmap(Matrix = meth.all.scaled,
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



pdf("Meth_betavalues_cluster_all.pdf",height=8,width=10)
fontsize = 10


pheatmap(t(meth.all.scaled), color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(100), 
         kmeans_k = NA, breaks = NA, border_color = "grey60",
         cellwidth = NA, cellheight = NA, scale = "none", 
         cluster_rows = F, cluster_cols = F, fontsize_row = 10,
         annotation_legend = TRUE, drop_levels = TRUE, 
         show_rownames = F,show_colnames = F, 
         treeheight_col = 0,
         main = "Methylation of promoters by module")
dev.off()

################################################################



#### want ? clusters
groups_x <- cutree(methClust, 10.5)
count(groups_x)  ## good, there is 10
factor_x <- as.factor(groups_x)

######## This splits the data into an array of dataframes based on clusters
cluster_array <- array(split(t_moduleTraitCor, factor_x))
dim(cluster_array)
### write out this?

for (i in 1:dim(cluster_array)){ 
  write.csv(cluster_array[[i]], paste("Clusters_whole_correlation_", i, ".csv", sep="")) 
} 


### create a summary of the values from each cluster  in green?
summary_array<-tapply(t_moduleTraitCor$MEgreen,factor_x,summary)
complete_summary_array <- do.call(rbind, summary_array)
write.table(complete_summary_array, file = "Cluster_summary_GreenModule.txt", sep = "\t")

### now save this data all as an RData file

##save as file that can be directly opened####

save(meth.all.scaled, methDist, methClust, cluster_array, complete_summary_array,
     file = "Meth_promoter_clusters_all_TCGA.RData")




################################## for TNBC only


medians = apply(t(meth_tnbc),2,median)
mads = apply(t(meth_tnbc),2,mad)
cors.use = scale(t(meth_tnbc),center=medians,scale=mads)


methDist <- dist(cors.use)
methClust <- hclust(methDist, "complete")

plot(methClust, labels = FALSE, main='Methylation of promoters TNBC dendro')

### cut at 50?


methDist1 <- dist(t(cors.use))
methClust1 <- hclust(methDist1, "complete")

plot(methClust1, main='Methylation of TNBC samples dendro')


#####
pdf("TNBC_Meth_scaled.pdf",height=8,width=10)
fontsize = 10


pheatmap(cors.use, color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(100), 
         kmeans_k = NA, breaks = NA, border_color = "grey60",
         cellwidth = NA, cellheight = NA, scale = "none", 
         cluster_rows = T, cluster_cols = T, fontsize_row = 10,
         annotation_legend = TRUE, drop_levels = TRUE, 
         show_rownames = F,show_colnames = T, 
         treeheight_col = 0,
         main = "Methylation of gene promoters TNBC dendro")
dev.off()



#### want 14-16 clusters
groups_x <- cutree(methClust, 2)
count(groups_x)  ## good, there is 10
factor_x <- as.factor(groups_x)
levels(factor_x)
cors.df %>% 
  group_by(factor_x) %>%
  summarise(no_rows = length(factor_x))

######## This splits the data into an array of dataframes based on clusters
cors.df <- as.data.frame(cors.use)
cluster_array <- array(split(cors.df, factor_x))
dim(cluster_array)
### write out this?

for (i in 1:dim(cluster_array)){ 
  write.csv(cluster_array[[i]], paste("Clusters_TNBC_promoters_", i, ".csv", sep="")) 
} 

### let's get some SD values?? summary values??


### create a summary of the values from each cluster  in green?
summary_array<-tapply($MEgreen,factor_x,summary)
complete_summary_array <- do.call(rbind, summary_array)
write.table(complete_summary_array, file = "Cluster_summary_GreenModule_TNBC.txt", sep = "\t")

summary_info <- cors.df
summary_info$mean <- rowMeans(summary_info)
summary_info$SD <- apply(summary_info,1, sd)

##now trim that by cluster?
cluster_stat <- array(split(summary_info, factor_x))
for (i in 1:dim(cluster_stat)){ 
  write.csv(cluster_stat[[i]], paste("Clusters_TNBC_stats_", i, ".csv", sep="")) 
} 



##Validating this cluster set:

############# validating clusters

library(NbClust)
nb <- NbClust(t_moduleTraitCor, distance = "euclidean", min.nc = 2,
              max.nc = 15, method = "complete", index ="all")
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
pam.res <- eclust(t_moduleTraitCor, "pam", k = 2, graph = FALSE)
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


######silhouette analysis

library("cluster")
sil <- silhouette(km.res$cluster, dist(t_moduleTraitCor))
head(sil[, 1:3], 10)

# Silhouette plot
plot(sil, main ="Silhouette plot - K-means")

# Summary of silhouette analysis
si.sum <- summary(sil)
# Average silhouette width of each cluster
si.sum$clus.avg.widths

fviz_silhouette(km.res, palette = "jco", 
                ggtheme = theme_classic())

