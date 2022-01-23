#######################################

### Mutational signatures vs MEs

setwd("~/Bioinformatics Work/Meth & RNA/Meta-analysis WGCNA")

library(WGCNA)
library(plyr)
library(flashClust)
library(ggplot2)
library(gplots)
library(Hmisc)
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

## read in the mutational signatures and MEs
dat <- read.table("ME_MutSig.txt", sep = "\t", header = TRUE)
dat <- read.table("ME_MutSig_all.txt", sep = "\t", header = TRUE)
### Function for flat correlation matrix using Hmisc
# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

res2<-rcorr(as.matrix(dat[,3:27]), type = "spearman")
Cors <- flattenCorrMatrix(res2$r, res2$P)

res <- cor.test(dat$Green.ME, dat$sig.3, methods = "spearman", exact = TRUE)
res$p.value

## write the correlation matriz
write.table(Cors, "ME_MutSig_correlations_all.txt", sep = "\t")


### reload with unneeded rows trimmed
sig <- read.table("ME_MutSig_correlations.txt", sep = "\t", header = TRUE)
#sig <- read.table("ME_MutSig_correlations_Pearson.txt", sep = "\t", header = TRUE)
# log 2 the p values
sig$log_p <- -log10(sig$p)

# trim to the green ME cors 
sig_g <- sig %>% filter(row == "Green.ME")

sig_g <- read.table("sig_g.txt", sep = "\t", header = TRUE)

# plots
p <- ggplot(sig_g, aes(x=log_p, y=cor)) +
  geom_rect(aes(xmin = -Inf, xmax = 1.30, ymin = -Inf, ymax = Inf), ## ad a rectangle of gry 
            fill = "gray88", alpha = 0.03) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0.0), ## ad a rectangle of gry 
            fill = "gray88", alpha = 0.03) +
  geom_point(data=subset(sig_g, sig_g$log_p <= 1.3),
             color = 'darkgrey', size=5) +
  geom_point(data=subset(sig_g, sig_g$log_p >= 1.3),
             size= 9, color = 'forestgreen', alpha = 0.4) +
  #geom_point(size= 8, color = 'forestgreen', alpha = 0.6) + ## scatterplot, colour fill is green, but transparent
  theme_classic(base_size = 22) + # just the bars, no inner markings / colours on graph 
  labs(x = "-Log10 p value", y = "Correlation")

p

sig_g$text <- if_else(sig_g$log_p <= 1.3 | is.na(sig_g$log_p), 6, 9)

## Ggrepel for the labels
#require("ggrepel")
set.seed(50)
p + geom_text_repel(aes(label = sig_g$column),# add the labetls, repelled from the dots
                    size = sig_g$text, colour = "black")


###### blue
# trim to the green ME cors 
sig_b <- sig %>% filter(row == "Blue.ME")
sig_b <- read.table("sig_b.txt", sep = "\t", header = TRUE)


# plots
p <- ggplot(sig_b, aes(x=log_p, y=cor)) +
  geom_rect(aes(xmin = -Inf, xmax = 1.30, ymin = -Inf, ymax = Inf), ## ad a rectangle of gry 
            fill = "gray88", alpha = 0.03) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0.0), ## ad a rectangle of gry 
            fill = "gray88", alpha = 0.03) +
  geom_point(data=subset(sig_b, sig_b$log_p <= 1.3),
             color = 'darkgrey', size=5) +
  geom_point(data=subset(sig_b, sig_b$log_p >= 1.3),
             size= 8, color = 'blue3', alpha = 0.4) +
  #geom_point(size= 8, color = 'blue3', alpha = 0.6) + ## scatterplot, colour fill is green, but transparent
  theme_classic(base_size = 22) + # just the bars, no inner markings / colours on graph 
  labs(x = "-Log10 p value", y = "Correlation")

p


## create a tezxt size column for plotting
sig_b$text <- if_else(sig_b$log_p <= 1.3 | is.na(sig_b$log_p), 6, 8)

## Ggrepel for the labels
#require("ggrepel")
set.seed(42)
p + geom_text_repel(aes(label = sig_b$column),# add the labetls, repelled from the dots
                     size = sig_b$text, colour = "black")
                    
                    

###########################################
####### creating 

dat1 <- read.table("ME_MutSig_Quart.txt", sep = "\t", header = TRUE, row.names = 1)
quart <- dat1 %>% mutate(quantile = ntile(Total, 4 ))

quart <- dat1 %>% select_if(is.numeric) %>% mutate(quantile = ntile(., 4 ))
#quart <- dat1  %>% mutate(is.numeric, ntile(., 4 ))
quart <- dat1 %>% mutate(across(where(is.numeric), .quantile = ntile(.col, 4 )))


quart$quantile <- as.factor(quart$quantile)
levels(quart)

quart$Group <- ifelse(quart$quantile == 3 | quart$quantile == 4, "High", ifelse(quart$quantile == 1, "Low", NA))

fin <- quart %>% filter(!is.na(Group))

write.table(fin, "Quartiles_MutSig.txt", sep = "\t")

grp <- as.factor(fin$Group)

es <- fin %>% 
  select_if(is.numeric) %>%
  map_df(~ broom::tidy(t.test(. ~ grp)), .id = 'var')

es$log <- -log10(es$p.value)




# The arguments to gather():
# - data: Data object
# - key: Name of new key column (made from names of data columns)
# - value: Name of new value column
# - ...: Names of source columns that contain values
# - factor_key: Treat the new key column as a factor (instead of character vector)
data_long <- gather(dat1, signature, measurement, sig.1:HRDetect.p, factor_key=TRUE)
data_long$signature <- as.factor(data_long$signature)
##### fuck it lets make the data long and just do it by groups that way 


## now lets quartile by group
fuk <- data_long %>%
  group_by(signature) %>%
  mutate(quartile = ntile(measurement, 4))


## add the quartile grouping
fuk$Group <- ifelse(fuk$quartile == 3 | fuk$quartile == 4, "High", ifelse(fuk$quartile == 1, "Low", NA))

grp <- as.factor(fuk$Group)


wat <- fuk  %>% 
  group_by(signature) %>%
  map_df(~ broom::tidy(t.test(fuk$Green.ME ~ grp)), .id = 'var')


write.table(fuk, "Fukit.txt", sep = "\t")

fukit <- read.table("Fukit.txt", sep = "\t", header = TRUE)

### means by group

dt_mean <- fukit %>% group_by(signature, Group) %>% summarise_each(funs(mean, sd))
write.table(dt_mean, "Means tables.txt", sep = "\t")


### t test by signature by quartile 
dt_result_1s <- fukit %>% group_by(signature) %>% do(tidy(t.test(Green.ME~Group, data=.,alternative = "greater")))
dt_result$log_p <- -log10(dt_result$p.value)
dt_result2_1s <- fukit %>% group_by(signature) %>% do(tidy(t.test(Blue.ME~Group, data=.)))
dt_result2$log_p <- -log10(dt_result2$p.value)


write.table(dt_result, "Quartiles_Indiv_Green_Pvalues.txt", sep = "\t" )
write.table(dt_result2, "Quartiles_Indiv_Blue_Pvalues.txt", sep = "\t" )


###############


########## Cross ref sig graph


xref <- read.table("Xref_Sig.txt", sep = "\t", header = TRUE)

p <- xref %>% filter(Module == "Green") %>% ggplot(aes(x= Group, y= Mean)) + 
  geom_bar(aes(fill = Group), stat="identity", position="dodge", width=.5) +
  facet_wrap(~Signature) +
  geom_signif(comparisons = list(c("High", "Low")), map_signif_level=TRUE)

p


p <- xref %>% filter(Module == "Blue") %>% ggplot(aes(x= Group, y= Mean)) + 
  geom_bar(aes(fill = Group), stat="identity", position="dodge", width=.5) +
  facet_wrap(~Signature) +
  geom_signif(comparisons = list(c("High", "Low")), map_signif_level=TRUE)

p



###############
fukit <- read.table("Fukit2.txt", sep = "\t", header = TRUE)


themeMain <- theme(strip.text.x = element_blank(), #font info for facet headings
                   plot.title = element_text(hjust = 0.5, size = 22, face = "bold"), #font info for main title
                   axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5), #angle of x lab
                   axis.text=element_text(size=16), #axis mark font info, 
                   axis.ticks = element_blank(),#get rid of axis ticks
                   axis.title=element_text(size=18,face="bold"),
                   plot.margin = margin(1, 1.5, 1, 1, "cm")) # increase plot margins to fit the 

library(lemon)
p <- fukit %>% 
  filter(signature %in% c("substitutions", "insertion.deletions", "sig.3", "sig.13", "RS2")) %>% #just the ones we want
  filter(!is.na(Group)) %>% # drop the NAs
  ggplot(aes(x= Group, y= ME)) + ## plotting Blue ME by group
  geom_boxplot(aes(fill = Group)) + ## it's a boxplot
  scale_fill_manual(values = c( "#CC79A7","#999999")) +# here;s the colour fil
  #facet_wrap(~signature, ncol = 3) + # wrapping by mutsig
  lemon::facet_rep_grid(Module ~ signature) + 
  coord_capped_cart(bottom='both', left='both') +
  theme_bw() +
  themeMain +
  xlab("") +
  ylab("Module Eigengene") +
  geom_signif(comparisons = list(c("High", "Low")), 
              test="t.test", 
              test.args=list(alternative = "greater", var.equal = FALSE, paired=FALSE), 
              map_signif_level=TRUE, textsize=8) +
  coord_cartesian(ylim = c(-0.05,0.18))
p




p <- fukit %>% 
  #filter(signature %in% c("substitutions", "insertion.deletions", "sig.3", "sig.13")) %>% #just the ones we want
  filter(!is.na(Group)) %>% # drop the NAs
  ggplot(aes(x= Group, y= Green.ME)) + ## plotting Blue ME by group
  geom_boxplot(aes(colour = Group)) + ## it's a boxplot
  ylim(-0.08, 0.18) +
  facet_wrap(~signature) + # wrapping by mutsig
  theme_classic() +
  themeMain +
  geom_signif(comparisons = list(c("High", "Low")), 
              test="t.test", 
              test.args=list(alternative = "greater", var.equal = FALSE, paired=FALSE), 
              map_signif_level=TRUE)
p











######################################

#####################      MEs vs DMM proteins




#### First lets set up the DMM protein list
DMM <- read.table("DMM_short.txt", sep = "\t", header = TRUE)
DMM <- read.table("Methyl Gene List.txt", sep = "\t", header = TRUE)
list3 <- intersect(rownames(datExpr1),DMM$Gene)
## and link it to Modules
dmm1 <- modcall %>% filter(Gene %in% list3 )
write.table(dmm1, "Meth_genes_modules.txt", sep = "\t")

## we'll also do the editors even though not invluded in the thesis
modcall <- read.table("All_genes_ModuleCall.txt", sep = "\t", header = TRUE)
call1 <- modcall %>% filter(Gene %in% list3)
write.table(call1, "DMM_Protein_modulecall_all.txt", sep = "\t")


##### We will create expression dataframe of these DMM proteins and add in the MEs
dmm_exp <- datExpr1[list3,]
dmm_exp <- as.data.frame(t(dmm_exp))
dmm_exp$Green.ME <- ME_1A$MEgreen
dmm_exp$Blue.ME <- ME_1A$MEblue
dmm_exp$Yellow.ME <- ME_1A$MEyellow

##### do the correlations on all samples
dmm2 <-rcorr(as.matrix(dmm_exp), type = "spearman")
Cors <- flattenCorrMatrix(dmm2$r, dmm2$P)
write.table(Cors, "DMM_ME_Cor.txt", sep = "\t")

### Now we'll just do it on the TNBC samples (n = 67)
clin1$TNBC <-  ifelse(clin1$ER_Status_nature2012=="Negative" & clin1$PR_Status_nature2012=="Negative" & clin1$HER2_Final_Status_nature2012=="Negative", "Y", "N")

clin1 <- clin1 %>% filter (TNBC == "Y")
list <- intersect(rownames(clin1), rownames(dmm_exp))
dmm_tnbc <- dmm_exp[list,]


dmm2 <-rcorr(as.matrix(dmm_tnbc), type = "spearman")
Cors <- flattenCorrMatrix(dmm2$r, dmm2$P)
write.table(Cors, "DMM_ME_Cor_TNBC.txt", sep = "\t")

### plot this
me_dmm <- read.table("ME_DMM.txt", sep = "\t", header = TRUE)


sig_g <- me_dmm %>% filter(Module == "Green.ME")
sig_g$log_p <- -log10(sig_g$P.value)

# plots
p <- ggplot(sig_g, aes(x=log_p, y=Cor)) +
  geom_rect(aes(xmin = -Inf, xmax = 1.30, ymin = -Inf, ymax = Inf), ## ad a rectangle of gry 
            fill = "gray88", alpha = 0.03) +
  #geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0.0), ## ad a rectangle of gry 
  #          fill = "gray88", alpha = 0.03) +
  geom_hline(yintercept = 0.0,color="red", linetype="dashed", size=1) +
  geom_point(data=subset(sig_g, sig_g$log_p <= 1.3),
             color = 'darkgrey', size=5) +
  geom_point(data=subset(sig_g, sig_g$log_p >= 1.3),
             size= 9, color = 'forestgreen', alpha = 0.4) +
  #geom_point(size= 8, color = 'forestgreen', alpha = 0.6) + ## scatterplot, colour fill is green, but transparent
  theme_classic(base_size = 22) + # just the bars, no inner markings / colours on graph 
  labs(x = "-Log10 p value", y = "Correlation")


sig_g$text <- if_else(sig_g$log_p <= 1.3 | is.na(sig_g$log_p), 0, 8)

## Ggrepel for the labels
require("ggrepel")
set.seed(42)
p + geom_text_repel(aes(label = sig_g$Gene),# add the labetls, repelled from the dots
                    size = sig_g$text, colour = "black")


#######

sig_b <- me_dmm %>% filter(Module == "Blue.ME")
sig_b$log_p <- -log10(sig_b$P.value)


# plots
p <- ggplot(sig_b, aes(x=log_p, y=Cor)) +
  geom_rect(aes(xmin = -Inf, xmax = 1.30, ymin = -Inf, ymax = Inf), ## ad a rectangle of gry 
            fill = "gray88", alpha = 0.03) +
  #geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0.0), ## ad a rectangle of gry 
   #         fill = "gray88", alpha = 0.03) +
  geom_hline(yintercept = 0.0,color="red", linetype="dashed", size=1) +
  geom_point(data=subset(sig_b, sig_b$log_p <= 1.3),
             color = 'darkgrey', size=5) +
  geom_point(data=subset(sig_b, sig_b$log_p >= 1.3),
             size= 8, color = 'blue3', alpha = 0.4) +
  #geom_point(size= 8, color = 'blue3', alpha = 0.6) + ## scatterplot, colour fill is green, but transparent
  theme_classic(base_size = 22) + # just the bars, no inner markings / colours on graph 
  labs(x = "-Log10 p value", y = "Correlation")

p


## create a tezxt size column for plotting
sig_b$text <- if_else(sig_b$log_p <= 1.3 | is.na(sig_b$log_p), "", 8)

## Ggrepel for the labels
#require("ggrepel")
set.seed(42)
p + geom_text_repel(aes(label = sig_b$Gene),# add the labetls, repelled from the dots
                    size = sig_b$text, colour = "black")





##############

###############   Can we try for FC vs p value?


### create a table with samples 

dmm_tnbc$Sample <- rownames(dmm_tnbc)
molten.data <- melt(dmm_tnbc, id = c("Sample","Green.ME","Blue.ME", "Yellow.ME"))


# The arguments to gather():
# - data: Data object
# - key: Name of new key column (made from names of data columns)
# - value: Name of new value column
# - ...: Names of source columns that contain values
# - factor_key: Treat the new key column as a factor (instead of character vector)
data_long <- gather(me_dmm, signature, measurement, sig.1:HRDetect.p, factor_key=TRUE)
data_long$signature <- as.factor(data_long$signature)
##### fuck it lets make the data long and just do it by groups that way 


## now lets quartile by group
fuk <- molten.data %>%
  group_by(variable) %>%
  mutate(quartile = ntile(value, 4))


## add the quartile grouping
fuk$Group <- ifelse(fuk$quartile == 4, "High", ifelse(fuk$quartile == 1, "Low", NA))

grp <- as.factor(fuk$Group)
gene <- as.factor(fuk$variable)


fukit <- fuk[,-c(1,7)]

fuk  %>% 
group_by(variable) %>%
  t.test(fuk$Green.ME ~ grp)

wat <- fuk  %>% 
  group_by(variable) %>%
  map_df(~ broom::tidy(t.test(fuk$Green.ME ~ grp)), .id = '')

fukit <- fukit %>% filter(!is.na(Group)) 

dt_mean <- fukit %>% group_by(variable, Group) %>% summarise_each(funs(mean, sd))
write.table(dt_mean, "Means tables_ME_DMM.txt", sep = "\t")
dt_mean$FC <- dt_mean %>% group_by(variable) %>% mutate(select(Group == "High", Green.ME) - select(Group == "Low", Green.ME))


df$freq[df$x==TRUE]
dt_mean %>% group_by(variable) %>% dt_mean$Green.ME_mean[dt_mean$Group == "High"]
dt_mean$Green.ME_mean[dt_mean$Group == "High"] - dt_mean$Green.ME_mean[dt_mean$Group == "Low"]


fc <- log2(dt_mean$value_mean[dt_mean$Group == "High"] - dt_mean$value_mean[dt_mean$Group == "Low"])
fc <- dt_mean$Blue.ME_mean[dt_mean$Group == "High"] - dt_mean$Blue.ME_mean[dt_mean$Group == "Low"]
log_fc <- log2(fc)

#### let's do it as a volcano plot
### load in the cool volcano plot library
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)

EnhancedVolcano(df,
                lab = rownames(df),
                x = 'log2FoldChange',
                y = 'pvalue')
