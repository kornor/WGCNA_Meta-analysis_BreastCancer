####
setwd("~/Bioinformatics Work/Meth & RNA/Meta-analysis WGCNA")

library(survival)
library(ggplot2)
library(survminer)
library(survival)
library(tidyverse)
library(tidytidbits)
library(survivalAnalysis)

######### This is a script to see about surivival just using the 
### modules in both TCGA and metabric

### load in the files as needed

load(file = "MetaAnalysis_trimmed_input.RData")
load(file = "Modules_DS0.RData")
colorsA1 = names(table(modules1))

### add in the rownames for the ME_ data frames

rownames(ME_1A) <- colnames(datExpr1)
rownames(ME_2A) <- colnames(datExpr2)
#### do binning of surivival by modules
### just for now, just do within the basals and use a 1/2 (??1/3)


### you're going to need the clinical info for both sets

######## start with TCGA

clin1 <- read.table("Clinical_final_TCGA.txt", sep = "\t",
                    header = TRUE, row.names = 1)
#### Intersect this with the trimmed sample list for METABRIC
list <- intersect(colnames(datExpr1), rownames(clin1))

clin1 <- clin1[list,]
datExpr1 <- datExpr1[list,]
#### Now that you have trimmed the list, set the factors you would like
Pam50_1 <- as.factor(clin1$PAM50Call_RNAseq)


clin_bas_tcga <- clin1[Pam50_1 == "Basal",]
datExpr1_bas <- datExpr1[,Pam50_1 == "Basal"]

list <- intersect(rownames(clin_bas_tcga), rownames(ME_1A))
mod_bas_tcga <- ME_1A[list,]


###now bin

##############  Trying to make a for loop over a data frame

iterations = 76  ## number of samples (rows) in the data frame
variables = 13  ## number of variables/genes (cols) in the data frame

output <- matrix(ncol=variables, nrow=iterations)


for(i in 1:ncol(mod_bas_tcga))
{
  col = mod_bas_tcga[,i]
  brk = unique(quantile(col, probs=c(0, 0.25, 0.75, 1)))
  
  if ( length(brk) < 3)
  { 
    for(j in 1:nrow(mod_bas_tcga))
      output[,i] = 'NA'
  }
  else
  {
    output[,i] =  as.factor(cut(col,
                                breaks=brk,
                                labels = FALSE, 
                                include.lowest=TRUE))
  }
}


output <- data.frame(output)
rownames(output) <- rownames(mod_bas_tcga)
colnames(output) <- colnames(mod_bas_tcga)

bads <- names(which(sapply(output, function(x) any(x == "NA"))))

for(i in 1:length(bads))
{
  output[,bads[i]] <- NULL
}


#write.table(output, "Bins_.txt", sep = "\t")


##################  Then do survival testing (Anova / cox) and order by which is sig. 
output <- read.table("Bins_exp_basal.txt", sep = "\t", row.names =1, header = TRUE)
#check if clinical file is there.  Toggle values
### To do the tops and tails, will need to exclude all the samples with a "level2" 
### ie, the middle group - they are unhelpful. 
### or do the anova one again??

library(survival)

SurvTest <- output

attach(SurvTest)

aov.out <- matrix(ncol = 5, nrow = 13)

## How to use only the top and bottom factors, and exlclude the others
#call <- factor(SurvTest$SOX10, levels = c(1, 3))

for(i in 1:ncol(SurvTest)) {
  
  #level = factor(SurvTest[,i])
  level = factor(SurvTest[,i], levels = c(1,3))
  survival_time = clin_bas_tcga$OS_Time_nature2012
  survival_event = clin_bas_tcga$OS_event_nature2012
  survival = Surv(survival_time,survival_event == 1) ~ level
  cox = coxph(survival, data = SurvTest)
  coeffs <- as.matrix(coef(summary(cox)))  ## coeffs has more info
  #aov = anova(cox)  ### anovas, for eg more than 2 groups
  aov.out[i,] <- coeffs
  #aov.out[i,] = (aov$`Pr(>|Chi|)`)[2]  ### output anovas from above
  
}

#### PUT IN GENE NAMES AS ROWNAMES HERE
rownames(aov.out) <- colnames(SurvTest)
colnames(aov.out) <- c("log Hazard Ratio", "Hazard Ratio", "SE Hazard Ratio", "Z score", "p.value" )

write.table(aov.out, "Survival_cox_all_tops_and_tails.txt", sep = "\t")










########################################################

#### read in the clinical file for metabric
clin2 <- read.table("Clinical_final_METABRIC.txt", sep = "\t",
                   header = TRUE, row.names = 1)
#### Intersect this with the trimmed sample list for METABRIC
list <- intersect(colnames(datExpr2), rownames(clin2))

clin2 <- clin2[list,]
#### Now that you have trimmed the list, set the factors you would like
Pam50_2 <- as.factor(clin2$Pam50_subtype)

clin_bas_meta <- clin2[Pam50_2 == "Basal",]
datExpr2_bas <- datExpr2[,Pam50_2 == "Basal"]

rownames(ME_2A) <- colnames(datExpr2)

list <- intersect(rownames(clin_bas_meta), rownames(ME_2A))
mod_bas_meta <- ME_2A[list,]


iterations = 177  ## number of samples (rows) in the data frame
variables = 13  ## number of variables/genes (cols) in the data frame

output <- matrix(ncol=variables, nrow=iterations)


for(i in 1:ncol(mod_bas_meta))
{
  col = mod_bas_meta[,i]
  brk = unique(quantile(col, probs=c(0, 0.25, 0.75, 1)))
  
  if ( length(brk) < 3)
  { 
    for(j in 1:nrow(mod_bas_meta))
      output[,i] = 'NA'
  }
  else
  {
    output[,i] =  as.factor(cut(col,
                                breaks=brk,
                                labels = FALSE, 
                                include.lowest=TRUE))
  }
}


output <- data.frame(output)
rownames(output) <- rownames(mod_bas_meta)
colnames(output) <- colnames(mod_bas_meta)

bads <- names(which(sapply(output, function(x) any(x == "NA"))))

for(i in 1:length(bads))
{
  output[,bads[i]] <- NULL
}


#write.table(output, "Bins_.txt", sep = "\t")


##################  Then do survival testing (Anova / cox) and order by which is sig. 
output <- read.table("Bins_exp_basal.txt", sep = "\t", row.names =1, header = TRUE)
#check if clinical file is there.  Toggle values
### To do the tops and tails, will need to exclude all the samples with a "level2" 
### ie, the middle group - they are unhelpful. 
### or do the anova one again??

library(survival)

SurvTest <- output

attach(SurvTest)

aov.out <- matrix(ncol = 5, nrow = 13)

## How to use only the top and bottom factors, and exlclude the others
#call <- factor(SurvTest$SOX10, levels = c(1, 3))

for(i in 1:ncol(SurvTest)) {
  
  #level = factor(SurvTest[,i])
  level = factor(SurvTest[,i], levels = c(1,3))
  survival_time = clin_bas_meta$Survival.Time
  survival_event = clin_bas_meta$Survival.Event
  survival = Surv(survival_time,survival_event == 1) ~ level
  cox = coxph(survival, data = SurvTest)
  coeffs <- as.matrix(coef(summary(cox)))  ## coeffs has more info
  #aov = anova(cox)  ### anovas, for eg more than 2 groups
  aov.out[i,] <- coeffs
  #aov.out[i,] = (aov$`Pr(>|Chi|)`)[2]  ### output anovas from above
  
}

#### PUT IN GENE NAMES AS ROWNAMES HERE
rownames(aov.out) <- colnames(SurvTest)
colnames(aov.out) <- c("log Hazard Ratio", "Hazard Ratio", "SE Hazard Ratio", "Z score", "p.value" )

##write.table(aov.out, "Survival_cox_all_tops_and_tails.txt", sep = "\t")


#### none of these are significant - probably the dilution effect of too
## many useless genes at the end of the module


### Survivial curves
library(survival)


mods <- c("MEblue", "MEgreen", "MEyellow")

mods_s <- output[,mods]
attach(mods_s)


survival_time = clin_bas_tcga$OS_Time_nature2012
survival_event = clin_bas_tcga$OS_event_nature2012
level = factor(mods_s$MEblue, levels = c(1,3))


fit.diff = survdiff(Surv(survival_time,survival_event == 1) ~ level) 

chisq2 = signif(1-pchisq(fit.diff$chisq,length(levels(level))-1),3) 
fit1 = survfit(Surv(survival_time,survival_event == 1)~level,
               conf.type="log-log") 
kmcolours <- c("black", "red", "green", "blue")
plot(fit1, conf.int=F,col=kmcolours,xlab="Time to death (days)", 
     ylab="Survival",main=c("All subtype survival by\n Module Expression"), 
     lwd=4,mark.time=TRUE)
legend("bottomleft",legend=levels(factor(mods_s)),
       fill = kmcolours, cex = 1)

###do another anova for signifigance
anova(coxph(Surv(OS_Time_nature2012,OS_event_nature2012 == 1)~factor(Brown), data = datClin))


### add to the end and do a prism one instead


mods <- c("MEblue", "MEgreen", "MEyellow")

mods_m <- output[,mods]

survival_time = clin_bas_meta$Survival.Time
survival_event = clin_bas_meta$Survival.Event

mods_m$survival_time <- survival_time
mods_m$survival_event <- survival_event


write.table(mods_s, "Modules_survival_TCGA.txt", sep = "\t")
write.table(mods_m, "Modules_survival_METABRIC.txt", sep = "\t")


levels(Pam50_2)


fit1 <- survfit(formula=Surv(Survival.Time, Survival.Event) ~ Pam50_2, data=clin2)
plot(fit1, col = c(1, 2, 3, 4, 5), ymin=0.45)
# colors: 1=black; 2=red, 4=blue
summary(fit1)

fit1 <- survfit(formula=Surv(Survival.Time, Survival.Event) ~ Pam50_2, data=clin2)
mdrcox <- coxph(Surv(Survival.Time, Survival.Event) ~ factor(Pam50_2) + factor(LN_Cat) + factor(histological_type),
                data=clin2)
summary(mdrcox)
cox.zph(mdrcox)


ggforest(cox.zph, data = clin2)



fit1 <- survfit(formula=Surv(Survival.Time, Survival.Event) ~ Pam50_2, data=clin2)
mdrcox <- coxph(Surv(Survival.Time, Survival.Event) ~  factor(Pam50_2) + factor(size_CAT) + factor(grade) + factor(Age_Cat),
                data=clin2)
summary(mdrcox)
cox.zph(mdrcox)
## if the global Chi-sq is significant, that suggests that the hazard is not proportionate - that is
# that there is some factor which is contributing a larger proportion of risks
# the factors in this instance are: PAM 50 (p value highly sig) and grade and age

# we can address this by stratefying:

mdrcox <- coxph(Surv(Survival.Time, Survival.Event) ~  strata(Pam50_2) + factor(size_CAT) + factor(grade) + factor(Age_Cat),
                data=clin2)
summary(mdrcox)
cox.zph(mdrcox)
## this stratification has revealed 
plot(cox.zph(mdrcox))




########################

### load in the bins 
bins <- read.table("Dose_mods.txt", sep = "\t", header = TRUE, row.names = 1)

list <- intersect(rownames(bins), rownames(clin2))

clin_bas_meta <- clin2[list,]


## order the dataframes
bins1 <- bins[ order(row.names(bins)), ]
clin_bas_meta1 <- clin_bas_meta[ order(row.names(clin_bas_meta)), ]

bins1$Age <- ifelse(clin_bas_meta1$age_at_diagnosis >= 50, 2, 1)
bins1$TumourSize <- clin_bas_meta1$size_CAT

bins1$Grade <- ifelse(clin_bas_meta1$grade == 2,1, 
                     ifelse(clin_bas_meta1$grade == 3, 2, 0))

#### create the factors

bins1$Bhigh <- ifelse(bins1$blue.per >= 0.33, 1, 0)
bins1$Ghigh <- ifelse(bins1$green.per >= 0.33, 1, 0)
bins1$Yhigh <- ifelse(bins1$yellow.per >= 0.1, 1, 0)

bins1$B_G <- ifelse(bins1$Bhigh == 1 & bins1$Ghigh == 1, 2, 
                   ifelse(bins1$Bhigh == 1 & bins1$Ghigh == 0, 1, 0))
bins1$G_B <- ifelse(bins1$Ghigh == 1 & bins1$Bhigh == 1, 2, 
                    ifelse(bins1$Ghigh == 1 & bins1$Bhigh == 0, 1, 0))

bins1$G_Y <- ifelse(bins1$Ghigh == 1 & bins1$Yhigh == 1, 2, 
                 ifelse(bins1$Ghigh == 1 & bins1$Yhigh == 0, 1, 0))

bins1$B_Y <- ifelse(bins1$Bhigh == 1 & bins1$Yhigh == 1, 2, 
                    ifelse(bins1$Bhigh == 1 & bins1$Yhigh == 0, 1, 0))


write.table(bins1, "Cox_input_META.txt", sep = "\t")

##################  Then do survival testing (Anova / cox) and order by which is sig. 
output <- read.table("Cox_input_META_2.txt", sep = "\t", row.names =1, header = TRUE)
#check if clinical file is there.  Toggle values
### To do the tops and tails, will need to exclude all the samples with a "level2" 
### ie, the middle group - they are unhelpful. 
### or do the anova one again??

library(survival)

SurvTest <- output

attach(SurvTest)

aov.out <- matrix(ncol = 5, nrow = 10)

## How to use only the top and bottom factors, and exclude the others
#call <- factor(SurvTest$SOX10, levels = c(1, 3))

for(i in 1:ncol(SurvTest)) {
  
  #level = factor(SurvTest[,i])
  level = factor(SurvTest[,i], levels = c(1,2))
  survival_time = bins1$Censor.yrs
  survival_event = bins1$Censor
  survival = Surv(survival_time,survival_event == 1) ~ level
  cox = coxph(survival, data = SurvTest)
  coeffs <- as.matrix(coef(summary(cox)))  ## coeffs has more info
  #aov = anova(cox)  ### anovas, for eg more than 2 groups
  aov.out[i,] <- coeffs
  #aov.out[i,] = (aov$`Pr(>|Chi|)`)[2]  ### output anovas from above
  
}

#### PUT IN GENE NAMES AS ROWNAMES HERE
rownames(aov.out) <- colnames(SurvTest)
colnames(aov.out) <- c("log Hazard Ratio", "Hazard Ratio", "SE Hazard Ratio", "Z score", "p.value" )

write.table(aov.out, "Survival_cox_groups_META.txt", sep = "\t")



################



fit1 <- survfit(formula=Surv(Survival.Time, Survival.Event) ~ Pam50_2, data=clin2)
mdrcox <- coxph(Surv(Survival.Time, Survival.Event) ~  factor(Pam50_2) + factor(size_CAT) + factor(grade) + factor(Age_Cat),
                data=clin2)
summary(mdrcox)
cox.zph(mdrcox)
## if the global Chi-sq is significant, that suggests that the hazard is not proportionate - that is
# that there is some factor which is contributing a larger proportion of risks
# the factors in this instance are: PAM 50 (p value highly sig) and grade and age

# we can address this by stratefying:

mdrcox <- coxph(Surv(Survival.Time, Survival.Event) ~  strata(Pam50_2) + factor(size_CAT) + factor(grade) + factor(Age_Cat),
                data=clin2)
summary(mdrcox)
cox.zph(mdrcox)
## this stratification has revealed 
plot(cox.zph(mdrcox))





fit <- survfit(formula=Surv(Censor.yrs, Censor) ~ G_Y, data=bins1)
plot(fit1, col = c(1, 2, 3, 4, 5), ymin=0.45)
# colors: 1=black; 2=red, 4=blue
summary(fit1, times = 20)

ggsurvplot(fit, 
           data = bins1,
           size = 1,                 # change line size
           #palette =
           # c("#E7B800", "#2E9FDF", ),# custom color palettes
           #conf.int = TRUE,          # Add confidence interval
           pval = TRUE,              # Add p-value
           risk.table = TRUE,        # Add risk table
           risk.table.col = "strata",# Risk table color by groups
           legend.labs = c("Other", "G+ B low", "G+ B high"),    # Change legend labels
           risk.table.height = 0.25, # Useful to change when you have multiple groups
           ggtheme = theme_bw()      # Change ggplot2 theme
)



#################
covariates <- c("Age", "TumourSize", "Grade", "Bhigh", "Ghigh", "Yhigh", "B_G", "G_B", "G_Y", "B_Y")
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(Censor.yrs, Censor)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = bins1)})
# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
res <- as.data.frame(res)

write.table(res, "Univariate Coxph Module Dose.txt", sep = "\t")



##### multivariate analysis using a larger number of variables
res.cox <- coxph(Surv(Censor.yrs, Censor) ~ Age + Grade + TumourSize + B_G + G_B + G_Y + B_Y, data =  bins1)#add in a correction for size
test.ph <- cox.zph(res.cox)
test.ph

summary(res.cox)

levels(bins1$B_Y)

## Now just on the CT treated patients
bins2 <- bins1[CT == 1,]

res.cox <- coxph(Surv(Censor.yrs, Censor) ~ Age + Grade + TumourSize + B_G + G_B + G_Y + B_Y, data =  bins2)#add in a correction for size
test.ph <- cox.zph(res.cox)
test.ph
summary(res.cox)


bins1$B_G <- as.factor(bins1$B_G)
bins1$G_B <- as.factor(bins1$G_B)
bins1$G_Y <- as.factor(bins1$G_Y)
bins1$B_Y <- as.factor(bins1$B_Y)

bins1$Age <- as.factor(bins1$Age)
bins1$Grade <- as.factor(bins1$Grade)
bins1$TumourSize <- as.factor(bins1$TumourSize)
levels(bins1$TumourSize)


covariate_names <- c(age="Age at Dx",
                     Grade="Tumour Grade",
                     TumourSize="Tumour Size",
                     B_G="B+ G high",
                     G_B ="G+ B high",
                     G_Y="G+ Y+",
                     B_Y="B+ Y+")

bins1 %>%
  analyse_multivariate(vars(Censor.yrs, Censor),
                       covariates = vars(Age,Grade,TumourSize,B_G,G_B,G_Y,B_Y),
                       covariate_name_dict=covariate_names,
                       reference_level_dict=c(Age = "1", Grade = "0",TumourSize = "1", B_G="1",G_B="1",G_Y="1",B_Y="1"))


write.table(bins1, "Cox_1.txt", sep = "\t")
write.table(bins2, "Cox_2.txt", sep = "\t")
