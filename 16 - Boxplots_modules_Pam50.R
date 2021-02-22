


setwd("~/Bioinformatics Work/Meth & RNA/Meta-analysis WGCNA")

library(WGCNA)
library(plyr)
library(dplyr)
library(flashClust)
library(ggplot2)
library(gplots)
library(pvclust)
library(lattice)

###
load(file = "MetaAnalysis_trimmed_input.RData")
load(file = "Modules_DS0.RData")
colorsA1 = names(table(modules1))


#### read in the clinical file for metabric
clin <- read.table("Clinical_final_METABRIC.txt", sep = "\t",
                   header = TRUE, row.names = 1)
#### Intersect this with the trimmed sample list for METABRIC
list <- intersect(colnames(datExpr2), rownames(clin))

clin <- clin[list,]
#### Now that you have trimmed the list, set the factors you would like
Pam50 <- as.factor(clin$Pam50_subtype)
ER <- as.factor(clin$ESR1_status)
grade <- as.factor(clin$grade)
site <- as.factor(clin$metastasis_site)
mets <- as.factor(clin$MET_Cat)
type <- as.factor(clin$histological_type)



## this does multiple boxplots by colour then by subtype
pdf("Module_Pam50.pdf",height=8,width=12)
for (c in 1:length(colorsA1)){
  inMod = modules1 == colorsA1[c]
  boxplot(ME_2A[,c] ~ Pam50,
          main=colorsA1[c],
          xlab="Pam 50",ylab="ME in METABRIC")
}; dev.off()


# try with ggplot2

write.table(ME_2A, "ME_Metabric.txt", sep = "\t")
ME_meta <- read.table("ME_Metabric.txt", sep = "\t", header = TRUE)

#rechape data into one table
library(reshape2)
dat.m <- melt(ME_meta,id.vars= "Sample", 
              measure.vars=c("Black","Blue","Brown","Green", "Magenta","Red","Tan", "Yellow"))
#add column with pam50
dat.m$pam50 <- as.factor(Pam50)

library(ggplot2)
plot <- ggplot(dat.m, aes(x=pam50, y=value, fill = pam50)) 

pdf("Module_Pam50x.pdf",height=12,width=8)
plot +
  geom_boxplot() +
  facet_wrap(~as.factor(variable), nrow = 4) +
  theme_bw() +
  theme(legend.position="none") +
  labs(title="Module Expression by Pam50 Subtype",
       x ="", y = "Module Eigengene") +
  theme(plot.title = element_text(hjust = 0.5))
  #scale_fill_brewer(palette="BuPu")
  
dev.off()


##########################
plot <- ggplot(dat.m, aes(x=as.factor(variable), y=value, fill = variable)) 
#wrap all these as a single theme you can call, make it easier 
themeMain <- theme(strip.text.x = element_text(size = 14, face = "bold"), #font info for facet headings
                   plot.title = element_text(hjust = 0.5, size = 22, face = "bold"), #font info for main title
                   axis.text.x = element_blank(), #no x label
                   axis.text=element_text(size=12), #axis mark font info, 
                   axis.ticks = element_blank(),#get rid of axis ticks
                   axis.title=element_text(size=14,face="bold")) # axis label font info



pdf("Pam50_Module_no legend.pdf",height=8,width=12)

plot +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~as.factor(pam50), nrow = 2) +
  theme_bw() +
  theme(legend.position="none") + #toggle the legend
  labs(fill = "Module") + #leged title
  labs(title="Module Expression by Pam50 Subtype",x ="", y = "Module Eigengene") + #set the axis labels
  scale_fill_manual(values = c("black","blue3","chocolate4","forestgreen", "magenta3","red3","tan", "yellow3")) +
  themeMain


dev.off()




 


