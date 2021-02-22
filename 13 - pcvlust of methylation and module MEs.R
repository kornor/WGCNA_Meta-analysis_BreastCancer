

#### pvclust uses columns not rows, so don't need the transposed version of the df
### have to set the method.hclust as complete, and the method.dist as euclidean.
## might take a while :/

methall.pv <- pvclust(moduleTraitCor, method.hclust ='complete',
                      method.dist = 'euclidean', parallel = FALSE)

plot(methall.pv)


### matches this one:
methDist <- dist(t_moduleTraitCor)
methClust <- hclust(methDist, "complete")

plot(methClust, main='Methylation of promoters dendro')

### back to pvclust

pvrect(methall.pv, alpha = 0.95) ### highlights clusters with >95% likihood of being real



