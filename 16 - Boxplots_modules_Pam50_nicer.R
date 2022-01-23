


setwd("~/Bioinformatics Work/Meth & RNA/Meta-analysis WGCNA")

library(WGCNA)
library(plyr)
library(dplyr)
library(flashClust)
library(ggplot2)
library(gplots)
library(pvclust)
library(lattice)
library(ggsignif)
library(ggpubr)
library(grid)
library(gridExtra)

###
load(file = "MetaAnalysis_trimmed_input.RData")
load(file = "Modules_DS0.RData")
colorsA1 = names(table(modules1))


#### read in the clinical file for metabric
clin <- read.table("Clinical_final_METABRIC.txt", sep = "\t",
                  header = TRUE, row.names = 1)

clin <- read.table("Clinical_final_TCGA.txt", sep = "\t",
                   header = TRUE, row.names = 1)

norm <- read.table("Normals_list.txt", sep = "\t", header = TRUE, row.names = 1)

#### Intersect this with the trimmed sample list for METABRIC
list <- intersect(colnames(datExpr2), rownames(clin))

#### Intersect this with the trimmed sample list for TCGA
list <- intersect(colnames(datExpr1), rownames(clin))

clin <- clin[list,]
#### Now that you have trimmed the list, set the factors you would like
Pam50 <- as.factor(clin$PAM50Call_RNAseq)
ER <- as.factor(clin$ER_Status_nature2012)
grade <- as.factor(clin$Tumor_nature2012)
HER <- as.factor(clin$HER2_Final_Status_nature2012)
mets <- as.factor(clin$Metastasis_Coded_nature2012)
#site <- as.factor(clin$metastasis_site)
type <- as.factor(clin$histological_type)
Node <- as.factor(clin$Node_Coded_nature2012)


## for metabric
Pam50 <- as.factor(clin$Pam50_subtype)

## this does multiple boxplots by colour then by subtype
pdf("Module_Pam50.pdf",height=8,width=12)
for (c in 1:length(colorsA1)){
  inMod = modules1 == colorsA1[c]
  boxplot(ME_2A[,c] ~ Pam50,
          main=colorsA1[c],
          xlab="Pam 50",ylab="ME in TCGA Dataset")
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
  labs(title="Module Expression by Pam50 Subtype - METABRIC",x ="", y = "Module Eigengene") + #set the axis labels
  scale_fill_manual(values = c("black","blue3","chocolate4","forestgreen", "magenta3","red3","tan", "yellow3")) +
  themeMain


dev.off()

plot +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~as.factor(pam50), nrow = 2) +
  theme_bw() +
  theme(legend.position="none") + #toggle the legend
  labs(fill = "Module") + #leged title
  labs(title="Module Expression by Pam50 Subtype",x ="", y = "Module Eigengene") + #set the axis labels
  scale_fill_manual(values = c("black","blue3","chocolate4","forestgreen", "magenta3","red3","tan", "yellow3")) +
  geom_signif(comparisons = list(c("", "")), map_signif_level=TRUE) +
  themeMain


dev.off()

############# TCGA

ME_tcga <- ME_1A
rownames(ME_tcga) <- colnames(datExpr1)
ME_tcga <- ME_tcga[list,]

write.table(ME_tcga, "ME_tcga.txt", sep = "\t")
ME_tcga <-read.table("ME_tcga.txt", sep = "\t", header = TRUE)
#ME_tcga <- ME_tcga[list,]


#rechape data into one table
library(reshape2)
dat.m <- melt(ME_tcga,id.vars= "Sample", 
              measure.vars=c("Black","Blue","Brown","Green", "Magenta","Red","Tan", "Yellow"))
#add column with pam50
levels(Pam50) <- c("Basal", "Her2 Enriched", "Luminal A", "Luminal B", "Normal-like" )
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
  geom_signif(stat = , map_signif_level=TRUE) +
  theme(plot.title = element_text(hjust = 0.5))
#scale_fill_brewer(palette="BuPu")

dev.off()



##########################
plot <- ggplot(dat.t, aes(x=as.factor(Module), y=ME, fill = Module)) 
#wrap all these as a single theme you can call, make it easier 
themeMain <- theme(strip.text.x = element_text(size = 14, face = "bold"), #font info for facet headings
                   plot.title = element_text(hjust = 0.5, size = 22, face = "bold"), #font info for main title
                   axis.text.x = element_blank(), #no x label
                   axis.text=element_text(size=12), #axis mark font info, 
                   axis.ticks = element_blank(),#get rid of axis ticks
                   axis.title=element_text(size=14,face="bold")) # axis label font info



pdf("TCGA_Pam50_Module_legend.pdf",height=8,width=12)

plot +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~as.factor(Pam50), nrow = 2) +
  theme_bw() +
  theme(legend.position="left") + #toggle the legend
  labs(fill = "Module") + #legend title
  labs(title="Module Expression by Pam50 Subtype - TCGA",x ="", y = "Module Eigengene") + #set the axis labels
  scale_fill_manual(values = c("black","blue3","chocolate4","forestgreen", "magenta3","red3","tan", "yellow3")) +
  themeMain


dev.off()

################### Module expression by Pam50, faceted by subtype


themeMain2 <- theme(strip.text.x = element_text(size = 14, face = "bold"),#font info for facet headings
                    strip.background = element_rect(fill = "white", color = "grey80", size = 1), # facet header colour
                   plot.title = element_text(hjust = 0.5, size = 22, face = "bold"), #font info for main title
                   axis.text.x = element_text(size=12, angle = 45, hjust = 1), #x label for this version
                   axis.text=element_text(size=12), #axis mark font info, 
                   axis.ticks = element_blank(),#get rid of axis ticks
                   axis.title=element_text(size=14,face="bold")) # axis label font info


plot <- ggplot(dat.m, aes(x=as.factor(pam50), y=value, fill = pam50)) 

plot +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(dat.m$variable, nrow = 2) +
  theme_minimal() +
  theme(legend.position="none") + #toggle the legend
  labs(fill = "Module") + #leged title
  labs(title="Module Expression by Pam50 Subtype",x ="", y = "Module Eigengene") + #set the axis labels +
  themeMain2


compare_means(value ~ pam50,  data = dat.m, ref.group = ".all.",method = "t.test") 
compare_means(value ~ pam50,  data = dat.m, group.by = "variable", method = "wilcox.test", ) 

table <- compare_means(value ~ pam50, group.by = "variable", ref.group = "variable", data = dat.m, method = "wilcox.test", ) 

# Visualize
plot +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(dat.m$variable, nrow = 2) +
  theme_bw() +
  theme(legend.position="none") + #toggle the legend
  labs(fill = "Module") + #leged title
  labs(title="Module Expression by Pam50 Subtype",x ="", y = "Module Eigengene") + #set the axis labels
  themeMain2 +
  stat_compare_means(method = "anova", label.y = .15)+      # Add global p-value, this is comparing the ME to the MEAN ME for all groups
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = ".all.", hide.ns = TRUE)                  # Pairwise comparison against all, hiding the non-sig



###########################

plot +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(dat.m$variable, nrow = 2) +
  theme_minimal() +
  theme(legend.position="none") + #toggle the legend
  labs(fill = "Module") + #leged title
  labs(title="Module Expression by Pam50 Subtype",x ="", y = "Module Eigengene") + #set the axis labels
  themeMain2 +
  #stat_compare_means(comparisons=levels(dat.m$pam50), method="wilcox.test", label="p.signif", color="red")
  geom_signif(test="wilcox.test", comparisons = combn(levels(dat.m$pam50),2, simplify = F)[-4], step_increase = 0.2, map_signif_level=c("***"=0.001,"**"=0.01, "*"=0.05, " "=2))
  
  
library(tidyverse)
  
stat_pvalue <- dat.m %>% rstatix::wilcox_test(value ~ pam50) %>% filter(p < 0.05) %>% 
    rstatix::add_significance("p") %>% 
    rstatix::add_y_position() %>% 
    mutate(y.position = seq(min(y.position), max(y.position),length.out = n()))
           
ggplot(dat.m, aes(x=pam50, y=value)) + geom_boxplot() +
ggpubr::stat_pvalue_manual(stat_pvalue, label = "p.signif") +
theme_bw(base_size = 16)


plot +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(dat.m$variable, nrow = 2) +
  theme_minimal() +
  theme(legend.position="none") + #toggle the legend
  labs(fill = "Module") + #leged title
  labs(title="Module Expression by Pam50 Subtype",x ="", y = "Module Eigengene") + #set the axis labels
  themeMain2 +
  ggpubr::stat_pvalue_manual(stat_pvalue, label = "p.signif")



###########################################
## Function for  Pam50 by module
# Function to generate and lay out the plots
signif_plot = function(signif.cutoff=0.05, height.factor=0.23) {
  
  # Get full range of y-values
  y_rng = range(dat.m$value)
  
  # Generate a list of three plots, one for each Species (these are the facets)
  plot_list = lapply(split(dat.m, dat.m$variable), function(d) {
    
    # Get pairs of x-values for current facet
    pairs = combn(levels(d$pam50), 2, simplify=FALSE)
    
    # Run wilcox test on every pair
    w.tst =  pairs %>% 
      map_df(function(lv) { 
        p.value = wilcox.test(d$value[d$pam50==lv[1]], d$value[d$pam50==lv[2]])$p.value
        data.frame(levs=paste(lv, collapse=" "), p.value)
      })
    
    # Record number of significant p.values. We'll use this later to adjust the top of the
    # y-range of the plots
    num_signif = sum(w.tst$p.value <= signif.cutoff)
    
    # Plot significance levels only for combinations with p <= signif.cutoff
    p = ggplot(d, aes(x=pam50, y=value, fill = pam50)) +
      geom_boxplot() +
      facet_grid(~variable, scales="free", space="free_x") +
      geom_signif(test="wilcox.test", comparisons = pairs[which(w.tst$p.value <= signif.cutoff)], 
                map_signif_level = c("****"=0.0001,"***"=0.001,"**" = 0.01 ,"*"=0.05, " "=2),            
                vjust=0,
                textsize=3,
                size=0.5,
                step_increase = .5) +
      theme_bw() +
      theme(axis.title=element_blank(),
            axis.text.x = element_text(angle=45, hjust=1), legend.position = "none",
            strip.background = element_rect(fill = "white", color = "grey80", size = 1))
    
    # Return the plot and the number of significant p-values
    return(list(num_signif, p))
  })
  
  # Get the highest number of significant p-values across all three "facets"
 max_signif = max(sapply(plot_list, function(x) x[[1]]))
  
  # Lay out the three plots as facets (one for each Species), but adjust so that y-range is same
  # for each facet. Top of y-range is adjusted using max_signif.
  grid.arrange(grobs=lapply(plot_list, function(x) x[[2]] + 
  scale_y_continuous(limits=c(y_rng[1], y_rng[2] + 0.5))), 
  ncol=3, left="Module Eigengene")
}

signif_plot(0.05) 
signif_plot(0.0001)


for (var in unique(mydata$Variety)) {
  dev.new()
  print( ggplot(mydata[mydata$Variety==var,], aes(Var1, Var2)) + geom_point() )
}



#################################################


dat.m <- read.table("METABRIC_ME_subtype_long.txt", sep = "\t", header = TRUE)
dat.t <- read.table("TCGA_ME_subtype long.txt", sep = "\t", header = TRUE)

levels(dat.t$Pam50) <- c("Basal", "Her2 Enriched", "Luminal A", "Luminal B", "Normal-like" )


plot <- ggplot(dat.m, aes(x=as.factor(Pam50), y=ME))
#wrap all these as a single theme you can call, make it easier 
themeMain <- theme(strip.text.x = element_text(size = 18, face = "bold"), #font info for facet headings
                   plot.title = element_text(hjust = 0.5, size = 22, face = "bold"), #font info for main title
                   axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5), #angle of x lab
                   axis.text=element_text(size=16), #axis mark font info, 
                   axis.ticks = element_blank(),#get rid of axis ticks
                   axis.title=element_text(size=18,face="bold"),
                   plot.margin = margin(1, 1.5, 1, 1, "cm")) # increase plot margins to fit the 



pdf("META_Module_Pam50_no legend_test.pdf",height=8,width=12)

plot +
  geom_boxplot(aes(fill = Module, alpha = Pam50)) +
  facet_wrap(~as.factor(Module), nrow = 2) +
  theme_bw() +
  theme(legend.position="none") + #toggle the legend
  labs(fill = "Module") + #leged title
  labs(title="Module Expression by Pam50 Subtype - METABRIC",x ="", y = "Module Eigengene") + #set the axis labels
  themeMain +
  scale_alpha_discrete(
    range = c(0.3, 0.9),
    guide = guide_legend(override.aes = list(fill = "black"))) +
    scale_fill_manual(values = c("black","blue3","chocolate4","forestgreen", "magenta3","red3","tan", "yellow3")) +
  #stat_compare_means(method = "anova", label.y = -0.06) +        # Add global annova p-value
  #stat_compare_means(label = "p.signif", method = "t.test",
  #                   ref.group = ".all.", label.y = 0.10)                      # Pairwise comparison against all




dev.off()

###################### TCGA

plot <- ggplot(dat.t, aes(x=as.factor(Pam50), y=ME))
#wrap all these as a single theme you can call, make it easier 
themeMain <- theme(strip.text.x = element_text(size = 18, face = "bold"), #font info for facet headings
                   plot.title = element_text(hjust = 0.5, size = 22, face = "bold"), #font info for main title
                   axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5), #angle of x lab
                   axis.text=element_text(size=16), #axis mark font info, 
                   axis.ticks = element_blank(),#get rid of axis ticks
                   axis.title=element_text(size=18,face="bold"),
                   plot.margin = margin(1, 1.5, 1, 1, "cm")) # increase plot margins to fit the 



pdf("TCGA_Module_Pam50_no legend.pdf",height=8,width=12)

plot +
  geom_boxplot(aes(fill = Module, alpha = Pam50)) +
  facet_wrap(~as.factor(Module), nrow = 2) +
  theme_bw() +
  theme(legend.position="none") + #toggle the legend
  labs(fill = "Module") + #leged title
  labs(title="Module Expression by Pam50 Subtype - TCGA",x ="", y = "Module Eigengene") + #set the axis labels
  themeMain +
  scale_alpha_discrete(
    range = c(0.3, 0.9),
    guide = guide_legend(override.aes = list(fill = "black"))) +
  scale_fill_manual(values = c("black","blue3","chocolate4","forestgreen", "magenta3","red3","tan", "yellow3")) 



dev.off()


