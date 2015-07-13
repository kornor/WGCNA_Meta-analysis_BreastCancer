### This is code to set up the traits file for METABRIC set

setwd("~/Bioinformatics Work/Meth & RNA/WGCNA_meta")

## Load in the exp file

exp <- read.table("Exp_final_METABRIC.txt", sep = "\t", header = TRUE, row.names = 1)
exp <- as.data.frame(t(exp))

## Load in the CIN metagene list, 
##  Load in the ER metagene list

cin <- read.table("CIN_metagene.txt", sep = "\t", header = TRUE, row.names = 1)
er <- read.table("ER_metagene.txt", sep = "\t", header = TRUE)

list <- which(duplicated(er[,1]))
er <- er[-list,]
### Get subsets of exp with only those genes 

list <- intersect(colnames(exp), rownames(cin))
cin_exp <- subset(exp[,list])


list <- intersect(colnames(exp), er)
er_exp <- subset(exp[,list])

### Mean-centred average for those
# center with 'colMeans()'
center_colmeans <- function(x) {
  xcenter = colMeans(x)
  x - rep(xcenter, rep.int(nrow(x), ncol(x)))
}

# apply it
centre_cin <- center_colmeans(cin_exp)

cin_ave <- as.data.frame(rowMeans(centre_cin))

rownames(cin_ave) <- rownames(cin_exp) 
colnames(cin_ave)[1] <- "CIN metagene"
write.table(cin_ave, "CIN_centred_average.txt", sep = "\t")

### for ER metagene

# apply it
centre_er <- center_colmeans(er_exp)

er_ave <- as.data.frame(rowMeans(centre_er))

rownames(er_ave) <- rownames(er_exp) 
colnames(er_ave)[1] <- "ER metagene"
write.table(er_ave, "ER_centred_average.txt", sep = "\t")


## Add in the other trait elements

traits <- c("ESR1", "PGR","AR", "ERBB2", "FOXM1", "HOTAIR")

list <- intersect(traits, colnames(exp))

## Create single file

traits_exp <- subset(exp[,list])
traits_exp$CIN <- cin_ave$`CIN metagene`
traits_exp$ER <- er_ave$`ER metagene`


### Export traits file

write.table(traits_exp, "Traits_METABRIC.txt", sep ="\t")
