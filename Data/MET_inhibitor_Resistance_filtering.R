###################################################################################################################
# MET_inhibitor_Resistance_filtering.R
# by g.estevam @ ucsf 
# code to filter resistance and sensitizing mutations 
###################################################################################################################

library(ggplot2)
library(gglorenz)
library(readr)
library(stringr)
library(tidyverse)
library(dplyr)
library(ggridges)
library(tidyquant)
library(bio3d)
library(corrplot)
library(data.table)
source("dms_analysis_utilities.R")
library("RColorBrewer")

library(scales)
library(colorspace)
library(ggpubr)
library(rstatix)


###################################################################################################################
# Notes on scores and significance testing: 

# all scores are normalized to the growth rate of DMSO (no inhibitor selection)

# pval0.5.up tests if variant growth rate is significantly larger than 0.5 (this is the positive tail of the distribution)
# pval0.5 <=  0.1 most confident are resistant
# pval0up tests whether variant growth rate is significantly larger than 0 ( this gives the tails of both distributions )


# "filtered" scores contain variants that pass these rules: 
#    1) 11/12 (3 replicate x 4 data points) counts for the variant == 0
#    2) all counts at time point 0 (before selection) == 0
#    3) mean count for the variant is <= 3


###################################################################################################################

#-------------------------------------- Load Rosace Fitness scores ----------------------------------------------#

### open data file as data frame for all raw inhibitor ALL scores generated from ROSACE

# delta exon14 
met_rosace_norm_scores <- as.data.frame(fread("met_scores_unfiltered.tsv")) # here met means MET delta Ex14
met_rosace_norm_scores  <- met_rosace_norm_scores %>% mutate(position = position + 1058)

# wt 
ex14_rosace_norm_scores <- as.data.frame(fread("ex14_scores_unfiltered.tsv")) # here ex14 means MET+Ex14 
ex14_rosace_norm_scores  <- ex14_rosace_norm_scores  %>% mutate(position = position + 1058)

## DMSO subtraction = detla.score, which is added as a new column 
ex14_rosace_norm_scores <- ex14_rosace_norm_scores %>% group_by(variants) %>% mutate(delta.score = mean - mean[key == "DMSO"]) %>%ungroup()
met_rosace_norm_scores <- met_rosace_norm_scores %>%group_by(variants) %>%mutate(delta.score = mean - mean[key == "DMSO"][1]) %>%ungroup()


#----------------------------------------- Define resistance mutations -------------------------------------#


# Note 1: use pval.up <= 0.1 for resistance mutation
# pval.up > 0.1 means that the variant's score is not significantly LARGER than the threshold.
# Note 2: pval <= 0.1 detects whether the score is significantly DIFFERENT from the threshold.
# So both lower than and greater than the threshold would be included.

# this will generate additional pvalue thresholds 
threshold <- 0.5 # provided arbitrary threshold here

## adds pval 0.75 
# wt
ex14_rosace_norm_scores_filtered <- ex14_rosace_norm_scores %>% rowwise() %>%
  mutate(pval0.5.up = pnorm(threshold, mean = mean, sd = sd, lower.tail = TRUE),
         pval0.5 = min(pnorm(threshold, mean = mean, sd = sd),
                       pnorm(threshold, mean = mean, sd = sd, lower.tail = FALSE))) 


# met del ex14
met_rosace_norm_scores_filtered <- met_rosace_norm_scores %>% rowwise() %>%
  mutate(pval0.5.up = pnorm(threshold, mean = mean, sd = sd, lower.tail = TRUE),
         pval0.5 = min(pnorm(threshold, mean = mean, sd = sd),
                       pnorm(threshold, mean = mean, sd = sd, lower.tail = FALSE))) 


#------------------------------------------ filter resistance mutations ------------------------------------------#


# these are the main dfs for the resistance mutations
# these are the scores that are used for plotting the positional manhattan type plots 
# this also provides the filtering for all the scores plotted in 3D


ex14_rosace_resistant_pval0.5 <- ex14_rosace_norm_scores_filtered %>%
  group_by(variants) %>%
  filter(
    (delta.score >= 0 & pval0.5.up <= 0.1 & mean >= 0 & mean[match("DMSO", key)] <= 0)) %>%
  ungroup()



ex14_rosace_condensed <- data.frame (hgvs = ex14_rosace_norm_scores $variants,
                                           pos = ex14_rosace_norm_scores $position, 
                                           mutation = ex14_rosace_norm_scores $mutation, 
                                           score = ex14_rosace_norm_scores $delta.score,
                                           type = ex14_rosace_norm_scores $type,
                                           inhib= ex14_rosace_norm_scores $key)


ex14_rosace_resis_condensed <- data.frame (hgvs = ex14_rosace_resistant_pval0.5$variants,
                                           pos = ex14_rosace_resistant_pval0.5$position, 
                                           mutation = ex14_rosace_resistant_pval0.5$mutation, 
                                           score = ex14_rosace_resistant_pval0.5$delta.score,
                                           type = ex14_rosace_resistant_pval0.5$type,
                                           inhib= ex14_rosace_resistant_pval0.5$key)


ex14_rosace_resis_wide <- ex14_rosace_resis_condensed %>% pivot_wider(names_from = inhib, values_from = score) 
ex14_rosace_norm_wide <- ex14_rosace_condensed %>% pivot_wider(names_from = inhib, values_from = score) 

### get all scores for individuals ####

# type Is
ex14_rosace_norm_wide_crizo <- data.frame(pos =  ex14_rosace_norm_wide$pos,score = ex14_rosace_norm_wide$Crizo, type =  ex14_rosace_norm_wide$type)
ex14_rosace_norm_wide_tepo <- data.frame(pos =  ex14_rosace_norm_wide$pos,score = ex14_rosace_norm_wide$Tepo, type =  ex14_rosace_norm_wide$type)
ex14_rosace_norm_wide_camp <- data.frame(pos =  ex14_rosace_norm_wide$pos,score = ex14_rosace_norm_wide$Camp, type =  ex14_rosace_norm_wide$type)
ex14_rosace_norm_wide_glu <- data.frame(pos =  ex14_rosace_norm_wide$pos,score = ex14_rosace_norm_wide$Glu, type =  ex14_rosace_norm_wide$type)
#ex14_rosace_norm_wide_savo <- data.frame(pos =  ex14_rosace_norm_wide$pos,score = ex14_rosace_norm_wide$Savo, type =  ex14_rosace_norm_wide$type) # these end up not passing our filters
#ex14_rosace_norm_wide_nvp <- data.frame(pos =  ex14_rosace_norm_wide$pos,score = ex14_rosace_norm_wide$NVP, type =  ex14_rosace_norm_wide$type) # these end up not passing our filters
# type IIs
ex14_rosace_norm_wide_mere <- data.frame(pos =  ex14_rosace_norm_wide$pos, score = ex14_rosace_norm_wide$Mere, type = ex14_rosace_norm_wide$type)
ex14_rosace_norm_wide_cabo <- data.frame(pos =  ex14_rosace_norm_wide$pos, score = ex14_rosace_norm_wide$Cabo, type = ex14_rosace_norm_wide$type)
ex14_rosace_norm_wide_gle <- data.frame(pos =  ex14_rosace_norm_wide$pos, score = ex14_rosace_norm_wide$Gle, type = ex14_rosace_norm_wide$type)
# type I.5
ex14_rosace_norm_wide_A458 = data.frame(pos =  ex14_rosace_norm_wide$pos, score = ex14_rosace_norm_wide$A458, type = ex14_rosace_norm_wide$type)

### get all scores for resistant ###

# type Is
ex14_resis_crizo = data.frame(hgvs = ex14_rosace_resis_wide$hgvs, pos = ex14_rosace_resis_wide$pos, mutation = ex14_rosace_resis_wide$mutation, score = ex14_rosace_resis_wide$Crizo)
ex14_resis_tepo = data.frame(hgvs = ex14_rosace_resis_wide$hgvs, pos = ex14_rosace_resis_wide$pos, mutation = ex14_rosace_resis_wide$mutation, score = ex14_rosace_resis_wide$Tepo)
ex14_resis_camp = data.frame(hgvs = ex14_rosace_resis_wide$hgvs, pos = ex14_rosace_resis_wide$pos, mutation = ex14_rosace_resis_wide$mutation, score = ex14_rosace_resis_wide$Camp)
ex14_resis_glu = data.frame(hgvs = ex14_rosace_resis_wide$hgvs, pos = ex14_rosace_resis_wide$pos, mutation = ex14_rosace_resis_wide$mutation, score = ex14_rosace_resis_wide$Glu)
#ex14_resis_savo = data.frame(hgvs = ex14_rosace_resis_wide$hgvs, pos = ex14_rosace_resis_wide$pos, mutation = ex14_rosace_resis_wide$mutation, score = ex14_rosace_resis_wide$Savo)
#ex14_resis_nvp = data.frame(hgvs = ex14_rosace_resis_wide$hgvs, pos = ex14_rosace_resis_wide$pos, mutation = ex14_rosace_resis_wide$mutation, score = ex14_rosace_resis_wide$NVP)
# type IIs
ex14_resis_mere = data.frame(hgvs = ex14_rosace_resis_wide$hgvs, pos = ex14_rosace_resis_wide$pos, mutation = ex14_rosace_resis_wide$mutation, score = ex14_rosace_resis_wide$Mere)
ex14_resis_cabo= data.frame(hgvs = ex14_rosace_resis_wide$hgvs, pos = ex14_rosace_resis_wide$pos, mutation = ex14_rosace_resis_wide$mutation, score = ex14_rosace_resis_wide$Cabo)
ex14_resis_gle= data.frame(hgvs = ex14_rosace_resis_wide$hgvs, pos = ex14_rosace_resis_wide$pos, mutation = ex14_rosace_resis_wide$mutation, score = ex14_rosace_resis_wide$Gle)
# type I.5
ex14_resis_A458 = data.frame(hgvs = ex14_rosace_resis_wide$hgvs, pos = ex14_rosace_resis_wide$pos, mutation = ex14_rosace_resis_wide$mutation, score = ex14_rosace_resis_wide$A458)



ex14_resis_crizo <- na.omit(ex14_resis_crizo) 
ex14_resis_tepo <- na.omit(ex14_resis_tepo) 
ex14_resis_camp <- na.omit(ex14_resis_camp) 
ex14_resis_glu<- na.omit(ex14_resis_glu) 
#ex14_resis_savo<- na.omit(ex14_resis_savo) 
#ex14_resis_nvp <- na.omit(ex14_resis_nvp)

ex14_resis_mere <- na.omit(ex14_resis_mere) 
ex14_resis_cabo <- na.omit(ex14_resis_cabo) 
ex14_resis_gle <- na.omit(ex14_resis_gle) 

ex14_resis_A458 <- na.omit(ex14_resis_A458) 


#--------------------------------- Met del Ex14 filter resistance mutations ------------------------------------------#


# these are the main dfs for the resistance mutations
# these are the scores that are used for plotting the positional manhattan type plots 
# this also provides the filtering for all the scores plotted in 3D


met_rosace_resistant_pval0.5 <- met_rosace_norm_scores_filtered %>%
  group_by(variants) %>%
  filter(
    (delta.score >= 0 & pval0.5.up <= 0.1 & mean >= 0 & mean[match("DMSO", key)] <= 0)) %>%
  ungroup()




met_rosace_condensed <- data.frame (hgvs = met_rosace_norm_scores $variants,
                                     pos = met_rosace_norm_scores $position, 
                                     mutation = met_rosace_norm_scores $mutation, 
                                     score = met_rosace_norm_scores $delta.score,
                                     type = met_rosace_norm_scores $type,
                                     inhib= met_rosace_norm_scores $key)


met_rosace_resis_condensed <- data.frame (hgvs = met_rosace_resistant_pval0.5$variants,
                                           pos = met_rosace_resistant_pval0.5$position, 
                                           mutation = met_rosace_resistant_pval0.5$mutation, 
                                           score = met_rosace_resistant_pval0.5$delta.score,
                                           type = met_rosace_resistant_pval0.5$type,
                                           inhib= met_rosace_resistant_pval0.5$key)


met_rosace_resis_wide <- met_rosace_resis_condensed %>% pivot_wider(names_from = inhib, values_from = score) 
met_rosace_norm_wide <- met_rosace_condensed %>% pivot_wider(names_from = inhib, values_from = score) 

### get all scores for individuals ####

# type Is
met_rosace_norm_wide_crizo <- data.frame(pos =  met_rosace_norm_wide$pos,score = met_rosace_norm_wide$Crizo, type =  met_rosace_norm_wide$type)
met_rosace_norm_wide_tepo <- data.frame(pos =  met_rosace_norm_wide$pos,score = met_rosace_norm_wide$Tepo, type =  met_rosace_norm_wide$type)
met_rosace_norm_wide_camp <- data.frame(pos =  met_rosace_norm_wide$pos,score = met_rosace_norm_wide$Camp, type =  met_rosace_norm_wide$type)
met_rosace_norm_wide_glu <- data.frame(pos =  met_rosace_norm_wide$pos,score = met_rosace_norm_wide$Glu, type =  met_rosace_norm_wide$type)
met_rosace_norm_wide_savo <- data.frame(pos =  met_rosace_norm_wide$pos,score = met_rosace_norm_wide$Savo, type =  met_rosace_norm_wide$type) # these end up not passing our filters
met_rosace_norm_wide_nvp <- data.frame(pos =  met_rosace_norm_wide$pos,score = met_rosace_norm_wide$NVP, type =  met_rosace_norm_wide$type) # these end up not passing our filters
# type IIs
met_rosace_norm_wide_mere <- data.frame(pos =  met_rosace_norm_wide$pos, score = met_rosace_norm_wide$Mere, type = met_rosace_norm_wide$type)
met_rosace_norm_wide_cabo <- data.frame(pos =  met_rosace_norm_wide$pos, score = met_rosace_norm_wide$Cabo, type = met_rosace_norm_wide$type)
met_rosace_norm_wide_gle <- data.frame(pos =  met_rosace_norm_wide$pos, score = met_rosace_norm_wide$Gle, type = met_rosace_norm_wide$type)
# type I.5
met_rosace_norm_wide_A458 = data.frame(pos =  met_rosace_norm_wide$pos, score = met_rosace_norm_wide$A458, type = met_rosace_norm_wide$type)

### get all scores for resistant ###

# type Is
met_resis_crizo = data.frame(hgvs = met_rosace_resis_wide$hgvs, pos = met_rosace_resis_wide$pos, mutation = met_rosace_resis_wide$mutation, score = met_rosace_resis_wide$Crizo)
met_resis_tepo = data.frame(hgvs = met_rosace_resis_wide$hgvs, pos = met_rosace_resis_wide$pos, mutation = met_rosace_resis_wide$mutation, score = met_rosace_resis_wide$Tepo)
met_resis_camp = data.frame(hgvs = met_rosace_resis_wide$hgvs, pos = met_rosace_resis_wide$pos, mutation = met_rosace_resis_wide$mutation, score = met_rosace_resis_wide$Camp)
met_resis_glu = data.frame(hgvs = met_rosace_resis_wide$hgvs, pos = met_rosace_resis_wide$pos, mutation = met_rosace_resis_wide$mutation, score = met_rosace_resis_wide$Glu)
met_resis_savo = data.frame(hgvs = met_rosace_resis_wide$hgvs, pos = met_rosace_resis_wide$pos, mutation = met_rosace_resis_wide$mutation, score = met_rosace_resis_wide$Savo)
met_resis_nvp = data.frame(hgvs = met_rosace_resis_wide$hgvs, pos = met_rosace_resis_wide$pos, mutation = met_rosace_resis_wide$mutation, score = met_rosace_resis_wide$NVP)

# type IIs
met_resis_mere = data.frame(hgvs = met_rosace_resis_wide$hgvs, pos = met_rosace_resis_wide$pos, mutation = met_rosace_resis_wide$mutation, score = met_rosace_resis_wide$Mere)
met_resis_cabo= data.frame(hgvs = met_rosace_resis_wide$hgvs, pos = met_rosace_resis_wide$pos, mutation = met_rosace_resis_wide$mutation, score = met_rosace_resis_wide$Cabo)
met_resis_gle= data.frame(hgvs = met_rosace_resis_wide$hgvs, pos = met_rosace_resis_wide$pos, mutation = met_rosace_resis_wide$mutation, score = met_rosace_resis_wide$Gle)
# type I.5
met_resis_A458 = data.frame(hgvs = met_rosace_resis_wide$hgvs, pos = met_rosace_resis_wide$pos, mutation = met_rosace_resis_wide$mutation, score = met_rosace_resis_wide$A458)



met_resis_crizo <- na.omit(met_resis_crizo) 
met_resis_tepo <- na.omit(met_resis_tepo) 
met_resis_camp <- na.omit(met_resis_camp) 
met_resis_glu<- na.omit(met_resis_glu) 
met_resis_savo<- na.omit(met_resis_savo) 
met_resis_nvp <- na.omit(met_resis_nvp)

met_resis_mere <- na.omit(met_resis_mere) 
met_resis_cabo <- na.omit(met_resis_cabo) 
met_resis_gle <- na.omit(met_resis_gle) 

met_resis_A458 <- na.omit(met_resis_A458) 



#----------------------------------  plot specific mutations:these are positional plots ------------------------------#

# data points in color are classified resistant mutations 


#------------- Met -----------#
scatterplots <- ggarrange(
ggplot() + 
  geom_point(data = ex14_rosace_norm_wide_crizo, aes(x = pos, y = score), color = "ivory4", alpha=0.2)+
  geom_point(data = ex14_resis_crizo , aes(x = pos, y = score),  color = "#e63466ff") + 
  xlab("Position")+
  ylab("Score")+
  ylim(0.5,2.5)+
  ggtitle("Crizotinib ")+
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")+
  theme_classic(),

ggplot() + 
  geom_point(data = ex14_rosace_norm_wide_tepo, aes(x = pos, y = score), color = "ivory4", alpha=0.2)+
  geom_point(data = ex14_resis_tepo , aes(x = pos, y = score),  color = "#e63466ff") + 
  xlab("Position")+
  ylab("Score")+
  ylim(0.5,2.5)+
  ggtitle("Tepotinib")+
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")+
  theme_classic(),

ggplot() + 
  geom_point(data = ex14_rosace_norm_wide_camp, aes(x = pos, y = score), color = "ivory4", alpha=0.2)+
  geom_point(data = ex14_resis_camp , aes(x = pos, y = score),  color = "#e63466ff") + 
  xlab("Position")+
  ylab("Score")+
  ylim(0.5,2.5)+
  ggtitle("Capmatinib")+
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")+
  theme_classic(),

ggplot() + 
  geom_point(data = ex14_rosace_norm_wide_glu, aes(x = pos, y = score), color = "ivory4", alpha=0.2)+
  geom_point(data = ex14_resis_glu , aes(x = pos, y = score),  color = "#e63466ff") + 
  xlab("Position")+
  ylab("Score")+
  ylim(0.5,2.5)+
  ggtitle("Glumetinib")+
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")+
  theme_classic(),

mere_resis <- ggplot() + 
  geom_point(data = ex14_rosace_norm_wide_mere, aes(x = pos, y = score), color = "ivory4", alpha=0.2)+
  geom_point(data = ex14_resis_mere, aes(x = pos, y = score),  color = "deepskyblue") + 
  xlab("Position")+
  ylab("Score")+
  ylim(0.5,2.5)+
  ggtitle("Merestinib")+
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")+
  theme_classic(),

ggplot() + 
  geom_point(data = ex14_rosace_norm_wide_cabo, aes(x = pos, y = score), color = "ivory4", alpha=0.2)+
  geom_point(data = ex14_resis_cabo, aes(x = pos, y = score),  color = "deepskyblue") + 
  xlab("Position")+
  ylab("Score")+
  ylim(0.5,2.5)+
  ggtitle("Cabozantinib")+
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")+
  theme_classic() ,


 ggplot() + 
  geom_point(data = ex14_rosace_norm_wide_gle, aes(x = pos, y = score), color = "ivory4", alpha=0.2)+
  geom_point(data = ex14_resis_gle, aes(x = pos, y = score),  color = "deepskyblue") + 
  xlab("Position")+
  ylab("Score")+
  ylim(0.5,2.5)+
  ggtitle("Glesatinib")+
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")+
  theme_classic() ,


ggplot() + 
  geom_point(data = ex14_rosace_norm_wide_A458, aes(x = pos, y = score), color = "ivory4", alpha=0.2)+
  geom_point(data = ex14_resis_A458, aes(x = pos, y = score),  color = "#8cbe61ff") + 
  xlab("Position")+
  ylab("Score")+
  ylim(0.5,2.5)+
  ggtitle("AMG-458")+
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")+
  theme_classic(),

nrow =2, ncol = 4
)

plot(scatterplots)



#------------- Met del ex14 -----------#

scatterplots <- ggarrange(
  ggplot() + 
    geom_point(data = met_rosace_norm_wide_crizo, aes(x = pos, y = score), color = "ivory4", alpha=0.2)+
    geom_point(data = met_resis_crizo , aes(x = pos, y = score),  color = "#e63466ff") + 
    xlab("Position")+
    ylab("Score")+
    ylim(0.5,2.5)+
    ggtitle("Crizotinib ")+
    geom_hline(yintercept = 1, linetype = "dashed", color = "black")+
    theme_classic(),
  
  ggplot() + 
    geom_point(data = met_rosace_norm_wide_tepo, aes(x = pos, y = score), color = "ivory4", alpha=0.2)+
    geom_point(data = met_resis_tepo , aes(x = pos, y = score),  color = "#e63466ff") + 
    xlab("Position")+
    ylab("Score")+
    ylim(0.5,2.5)+
    ggtitle("Tepotinib")+
    geom_hline(yintercept = 1, linetype = "dashed", color = "black")+
    theme_classic(),
  
  ggplot() + 
    geom_point(data = met_rosace_norm_wide_camp, aes(x = pos, y = score), color = "ivory4", alpha=0.2)+
    geom_point(data = met_resis_camp , aes(x = pos, y = score),  color = "#e63466ff") + 
    xlab("Position")+
    ylab("Score")+
    ylim(0.5,2.5)+
    ggtitle("Capmatinib")+
    geom_hline(yintercept = 1, linetype = "dashed", color = "black")+
    theme_classic(),
  
  ggplot() + 
    geom_point(data = met_rosace_norm_wide_glu, aes(x = pos, y = score), color = "ivory4", alpha=0.2)+
    geom_point(data = met_resis_glu , aes(x = pos, y = score),  color = "#e63466ff") + 
    xlab("Position")+
    ylab("Score")+
    ylim(0.5,2.5)+
    ggtitle("Glumetinib")+
    geom_hline(yintercept = 1, linetype = "dashed", color = "black")+
    theme_classic(),
  
  ggplot() + 
    geom_point(data = met_rosace_norm_wide_nvp, aes(x = pos, y = score), color = "ivory4", alpha=0.2)+
    geom_point(data = met_resis_nvp , aes(x = pos, y = score),  color = "#f088a5ff") + 
    xlab("Position")+
    ylab("Score")+
    ylim(0.5,2.5)+
    ggtitle("NVP-BVU972")+
    geom_hline(yintercept = 1, linetype = "dashed", color = "black")+
    theme_classic(),
  
   ggplot() + 
    geom_point(data = met_rosace_norm_wide_mere, aes(x = pos, y = score), color = "ivory4", alpha=0.2)+
    geom_point(data = met_resis_mere, aes(x = pos, y = score),  color = "deepskyblue") + 
    xlab("Position")+
    ylab("Score")+
    ylim(0.5,2.5)+
    ggtitle("Merestinib")+
    geom_hline(yintercept = 1, linetype = "dashed", color = "black")+
    theme_classic(),
  
  ggplot() + 
    geom_point(data = met_rosace_norm_wide_cabo, aes(x = pos, y = score), color = "ivory4", alpha=0.2)+
    geom_point(data = met_resis_cabo, aes(x = pos, y = score),  color = "deepskyblue") + 
    xlab("Position")+
    ylab("Score")+
    ylim(0.5,2.5)+
    ggtitle("Cabozantinib")+
    geom_hline(yintercept = 1, linetype = "dashed", color = "black")+
    theme_classic() ,
  
  
  ggplot() + 
    geom_point(data = met_rosace_norm_wide_gle, aes(x = pos, y = score), color = "ivory4", alpha=0.2)+
    geom_point(data = met_resis_gle, aes(x = pos, y = score),  color = "deepskyblue") + 
    xlab("Position")+
    ylab("Score")+
    ylim(0.5,2.5)+
    ggtitle("Glesatinib")+
    geom_hline(yintercept = 1, linetype = "dashed", color = "black")+
    theme_classic() ,
  
  
  ggplot() + 
    geom_point(data = met_rosace_norm_wide_A458, aes(x = pos, y = score), color = "ivory4", alpha=0.2)+
    geom_point(data = met_resis_A458, aes(x = pos, y = score),  color = "#8cbe61ff") + 
    xlab("Position")+
    ylab("Score")+
    ylim(0.5,2.5)+
    ggtitle("AMG-458")+
    geom_hline(yintercept = 1, linetype = "dashed", color = "black")+
    theme_classic(),
  
  nrow =3, ncol = 3
)

plot(scatterplots)


#--------------------------------- map resistance positions for individuals on PDBs--------------------------------------------#


# read pdb files 
crizo_2WGJ <- read.pdb("2WGJ")
tepo_4R1V <- read.pdb("4R1V")
tiv_3RHK <- read.pdb("3RHK")
amg458_5T3Q <- read.pdb("5T3Q")
nvp_3QTI <- read.pdb("3QTI")
mere_4EEV <- read.pdb("4EEV")
savo_6SDE <- read.pdb("6SDE")
cap_docked <- read.pdb("capmatinib_2wgj_docked.pdb")
cabo_docked <- read.pdb("cabozantinib_4eev_docked.pdb")
glu_docked <- read.pdb("glumetinib_2wgj_docked.pdb")
gle_docked <- read.pdb("glesatinib_4eev_docked.pdb")



# function to make all pdbs w/ mapped resistance 
map_pdb <- function(df, pdb, file.name) {
  df2 = df %>% group_by(pos) %>% dplyr::summarise(mean = mean(score, na.rm=TRUE))
  x = map_scores_pdb(pdb, df2, "mean")
  write.pdb(x, file=file.name)
}


map_pdb(ex14_resis_crizo, crizo_2WGJ, "ex14_Crizo_resis.pdb")
map_pdb(ex14_resis_tepo, tepo_4R1V, "ex14_Tepo_resis.pdb")
map_pdb(ex14_resis_camp, cap_docked, "ex14_Camp_resis.pdb")
map_pdb(ex14_resis_glu, glu_docked, "ex14_Glu_resis.pdb")
map_pdb(ex14_resis_mere, mere_4EEV, "ex14_Mere_resis.pdb")
map_pdb(ex14_resis_cabo, cabo_docked, "ex14_Cabo_resis.pdb")
map_pdb(ex14_resis_gle, gle_docked, "ex14_Gle_resis.pdb")
map_pdb(ex14_resis_A458, amg458_5T3Q, "ex14_A458_resis.pdb")




#-------------------------------------  plot distributions  ------------------------------------------#


ex14_rosace_norm_wide_tiv <- data.frame(pos =  ex14_rosace_norm_wide$pos,
                                        score = ex14_rosace_norm_wide$Tiv,
                                        type =  ex14_rosace_norm_wide$type)
ex14_rosace_norm_wide_tiv <- ex14_rosace_norm_wide_tiv %>% filter(type != "synonymous" & type != "nonsense")

ex14_rosace_norm_wide_dmso <- data.frame(pos =  ex14_rosace_norm_wide$pos,
                                         score = ex14_rosace_norm_wide$DMSO,
                                         type =  ex14_rosace_norm_wide$type)
ex14_rosace_norm_wide_dmso <- ex14_rosace_norm_wide_dmso %>% filter(type != "synonymous" & type != "nonsense")


plot_Tiv <- ggplot( ) +
  geom_histogram(aes(x = ex14_rosace_norm_wide_dmso$score),color="grey4",fill="grey", alpha=0.5) +
  geom_histogram(aes(x = ex14_rosace_norm_wide_tiv$score),color="darkorange4",fill="darkorange", alpha=0.5) +
  xlab("Fitness Score") + 
  ylab("Count") + 
  ggtitle("Tivantinib")+
  xlim(-4,3)+
  geom_vline(xintercept = 0, linetype="dashed") +
  theme_classic()+
  theme(panel.background = element_rect(colour = "black", size=1))
plot(plot_Tiv)

plot_crizo <- ggplot( ) +
  geom_histogram(aes(x = ex14_rosace_norm_wide_dmso$score),color="grey4",fill="grey", alpha=0.5) +
  geom_histogram(aes(x = ex14_rosace_norm_wide_crizo$score),color="darkorange4",fill="darkorange", alpha=0.5) +
  xlab("Fitness Score") + 
  ylab("Count") + 
  ggtitle("Tivantinib")+
  xlim(-4,3)+
  geom_vline(xintercept = 0, linetype="dashed") +
  theme_classic()+
  theme(panel.background = element_rect(colour = "black", size=1))
plot(plot_crizo)



#------------------------------ General Fitness for ATP and R-Spine residues -------------------------------------#

ATP_site_resi <- c("1084","1085","1086","1087","1088","1089","1090","1091",
                   "1057","1058","1059","1060","1061","1062","1063","1064", 
                   "1211", "1222","1140","1204","1205","1206","1207","1208",
                   "1209","1210","1110","1110")
Rspine_resi <- c("1131", "1142", "1202","1223")

atp_site_fitness <-  ex14_rosace_norm_scores_filtered %>% filter(position %in%  ATP_site_resi)  %>% filter(type != "synonymous" & type!= "nonsense" & key != "DMSO") 
Rspine_fitness <-  ex14_rosace_norm_scores_filtered %>% filter(position %in%  Rspine_resi)  %>% filter(type != "synonymous" & type!= "nonsense" & key != "DMSO")


plot_ATPsite <- ggplot( ) +
  geom_histogram(aes(x = atp_site_fitness$mean),color="deeppink4",fill="maroon1", alpha=0.6) +
  xlab("Fitness Score") + 
  ylab("Count") + 
  ggtitle("ATP neighboring")+
  xlim(-4,3)+
  geom_vline(xintercept = 0, linetype="dashed") +
  theme_classic()+
  theme(panel.background = element_rect(colour = "black", size=1))
plot(plot_ATPsite)

plot_Rspine <- ggplot( ) +
  geom_histogram(aes(x = Rspine_fitness$mean),color="royalblue4",fill="deepskyblue", alpha=0.7) +
  xlab("Fitness Score") + 
  ylab("Count") + 
  ggtitle("R-spine")+
  xlim(-4,3)+
  geom_vline(xintercept = 0, linetype="dashed") +
  theme_classic()+
  theme(panel.background = element_rect(colour = "black", size=1))
plot(plot_Rspine)


core_residues = ex14_rosace_norm_scores_filtered%>% filter(key == "DMSO" & type != "synonymous" & type != "nonsense" & (pos %in% c(1095,1107,1108,1121,1124,1127,1128,1131,1139,1144,1145,
                                                                                                1146,1153,1155,1156,1161,1168,1176,1180,1181,1182,1183,
                                                                                                1184,1185,1188,1189,1191,1195,1202,1205,1206,1207,1210,
                                                                                                1212,1218,1220,1223,1226,1251,1254,1265,1266,1267,1269,
                                                                                                1270,1271,1272,1273,1275,1275,1280,1292,1312,1316,1320,1321,1328,1333,1337)))

surface_residues = ex14_rosace_norm_scores_filtered%>% filter(key == "DMSO" & type != "synonymous" & type != "nonsense" & !(pos %in% c(1095,1107,1108,1121,1124,1127,1128,1131,1139,1144,1145,
                                                                                                    1146,1153,1155,1156,1161,1168,1176,1180,1181,1182,1183,
                                                                                                    1184,1185,1188,1189,1191,1195,1202,1205,1206,1207,1210,
                                                                                                    1212,1218,1220,1223,1226,1251,1254,1265,1266,1267,1269,
                                                                                                    1270,1271,1272,1273,1275,1275,1280,1292,1312,1316,1320,1321,1328,1333,1337)))

core_surface <- ggplot( ) +
  geom_histogram(aes(x = surface_residues$mean),color="royalblue4",fill="purple4", alpha=0.2) +
  geom_histogram(aes(x = core_residues$mean),color="royalblue4",fill="mediumpurple1", alpha=0.5) +
  xlab("Fitness Score") + 
  ylab("Count") + 
  ggtitle("DMSO")+
  xlim(-4,3)+
  geom_vline(xintercept = 0, linetype="dashed") +
  theme_classic()+
  theme(panel.background = element_rect(colour = "black", size=1))
plot(core_surface)



#---------------------------------- Group inhibitors by type and Sum Resistance Mutations -----------------------------------#

custom_palette <- c("#e63466ff", "#56B4E9","#8cbe61ff")

library(dplyr)
library(tidyr)

# define groups
key_list <- c("DMSO", "A458", "Crizo", "Camp", "Tepo", "Tiv", "Glu", "Savo", "NVP", "Cabo", "Mere", "Gle")
group_list <- c("Control", "I5", "I", "I", "I", "III", "I", "I", "I", "II", "II", "II")


#ex14_rosace_resistant_pval0.5 <- ex14_rosace_norm_scores_filtered %>%


# generate data frames 
ex14_scores_type <- ex14_rosace_resistant_pval0.5 %>%
  mutate(inhib = ifelse(key %in% key_list, group_list[match(key, key_list)], NA)) 

ex14_resis_dt <- as.data.table(ex14_scores_type)
#ex14_type_resis_sum <- ex14_resis_sensitive_dt[, .(count = sum(pval0.5.up <= 0.1, na.rm = TRUE)), by = .(variants, position, mutation, inhib)][, dcast(.SD, ... ~ inhib, value.var = "count", fill = 0)]

ex14_type_resis_sum <- ex14_resis_dt[, .(count = .N), by = .(variants, position, mutation, inhib)][, dcast(.SD, variants + position + mutation ~ inhib, value.var = "count", fill = 0)]

resis_counts <- data.frame(position = ex14_type_resis_sum$position,
                           I = ex14_type_resis_sum$I,
                           II = ex14_type_resis_sum$II,
                           I5 = ex14_type_resis_sum$I5)

position_sum <- resis_counts %>%
  group_by(position) %>%
  summarise_all(list(~sum(., na.rm = TRUE))) %>%
  arrange(position)


all_positions <- data.frame(position = seq(1059, 1345, by = 1))

# Left join your existing data frame with the full positions data frame
position_sum <- all_positions %>%
  left_join(position_sum, by = "position") %>%
  mutate(I = ifelse(is.na(I), 0, I),
         II = ifelse(is.na(II), 0, II),
         I5 = ifelse(is.na(I5), 0, I5))



#-------------------------------- Venn diagram of shared resistance positions  -----------------------------------#

library(ggvenn)
library(RVenn)
library(sf)
library(ggVennDiagram)

venn_df <- data.frame(hgvs = ex14_resis_dt$variants, 
                      pos = ex14_resis_dt$pos, 
                      inhib= ex14_resis_dt$inhib)
venn_df$hgvs <- gsub("p.\\(|\\)|[0-9]", "", venn_df$hgvs)
venn_df <- venn_df %>% 
  mutate(mut = paste0(substr(hgvs, 1, 1), pos, substr(hgvs, 2, nchar(hgvs))))

venn_df_mut <- data.frame(mut= venn_df$mut,
                      inhib = venn_df$inhib)

venn_df_pos <- data.frame(pos= venn_df$pos,
                       inhib = venn_df$inhib)

venn_dt_mut <- as.data.table(venn_df_mut)
venn_dt_mut_sum <- venn_dt_mut[, .(count = pmin(.N, 1)), by = .(mut, inhib)][, dcast(.SD, mut ~ inhib, value.var = "count", fill = 0)]
venn_df2_mut <- data.frame(venn_dt_mut_sum)

venn_dt_pos <- as.data.table(venn_df_pos)
venn_dt_pos_sum <- venn_dt_pos[, .(count = pmin(.N, 1)), by = .(pos, inhib)][, dcast(.SD, pos ~ inhib, value.var = "count", fill = 0)]
venn_df2_pos <- data.frame(venn_dt_pos_sum)


#----- Venn 

# mutation based
sets_list2 <- lapply(venn_df2_mut[, -1], function(x) which(x == 1))
ggVennDiagram(
  sets_list2,
  category.names = colnames(venn_df2_mut)[-1],
  filename = NULL,
  output = TRUE,
  labels = list(pos = unique(venn_df2_mut$mut)), 
  category.weights = list(Mutation = unique(venn_df2_mut$mut))
)

mut_1= venn_df2_mut%>% group_by(mut) %>% filter (I==1 & II==0 & I5==0) 
mut_2= venn_df2_mut%>% group_by(mut) %>% filter (I==0 & II==1 & I5==0) 
mut_1.5= venn_df2_mut%>% group_by(mut) %>% filter (I==0 & II==0 & I5==1) 
mut_1_2= venn_df2_mut%>% group_by(mut) %>% filter (I==1 & II==1 & I5==0)
mut_1_1.5= venn_df2_mut%>% group_by(mut) %>% filter (I==1 & II==0 & I5==1) 
mut_2_1.5= venn_df2_mut%>% group_by(mut) %>% filter (I==0 & II==1& I5==1) 
mut_1_2_1.5= venn_df2_mut%>% group_by(mut) %>% filter (I==1 & II==1 & I5==1)



# position based
sets_list <- lapply(venn_df2_pos[, -1], function(x) which(x == 1))
ggVennDiagram(
  sets_list,
  category.names = colnames(venn_df2_pos)[-1],
  filename = NULL,
  output = TRUE,
  labels = list(pos = unique(venn_df2_pos$pos)), 
  category.weights = list(Position = unique(venn_df2_pos$pos))
)


sum_1= venn_df2_pos%>% group_by(pos) %>% filter (I==1 & II==0 & I5==0) %>% dplyr::summarise(sum = sum(pos, na.rm=TRUE))
sum_2= venn_df2_pos%>% group_by(pos) %>% filter (I==0 & II==1 & I5==0) %>% dplyr::summarise(sum = sum(pos, na.rm=TRUE))
sum_1.5= venn_df2_pos%>% group_by(pos) %>% filter (I==0 & II==0 & I5==1) %>% dplyr::summarise(sum = sum(pos, na.rm=TRUE))
sum_1_2= venn_df2_pos%>% group_by(pos) %>% filter (I==1 & II==1 & I5==0) %>% dplyr::summarise(sum = sum(pos, na.rm=TRUE))
sum_1_1.5= venn_df2_pos%>% group_by(pos) %>% filter (I==1 & II==0 & I5==1) %>%dplyr::summarise(sum = sum(pos, na.rm=TRUE))
sum_2_1.5= venn_df2_pos%>% group_by(pos) %>% filter (I==0 & II==1& I5==1) %>% dplyr::summarise(sum = sum(pos, na.rm=TRUE))
sum_1_2_1.5= venn_df2_pos%>% group_by(pos) %>% filter (I==1 & II==1 & I5==1) %>% dplyr::summarise(sum = sum(pos, na.rm=TRUE))

ex14_ATP_crizo = read.pdb("2WGJ")
a = map_scores_pdb(ex14_ATP_crizo,sum_1, "sum")
b = map_scores_pdb(ex14_ATP_crizo,sum_2, "sum")
c = map_scores_pdb(ex14_ATP_crizo,sum_1.5, "sum")
d = map_scores_pdb(ex14_ATP_crizo,sum_1_2, "sum")
e = map_scores_pdb(ex14_ATP_crizo,sum_1_1.5, "sum")
f = map_scores_pdb(ex14_ATP_crizo,sum_2_1.5, "sum")
g = map_scores_pdb(ex14_2WGJ_crizo,sum_1_2_1.5, "sum")

write.pdb(a, file="sum_1.pdb")
write.pdb(b, file="sum_2.pdb")
write.pdb(c, file="sum_1.5.pdb")
write.pdb(d, file="sum_1_2.pdb")
write.pdb(e, file="sum_1_1.5.pdb")
write.pdb(f, file="sum_2_1.5.pdb")
write.pdb(g, file="sum_1_2_1.5.pdb")


#----- PieDonut 

venn_df2_mut <- venn_df2_mut %>%
  mutate(group = case_when(
    I == 1 & II == 0 & I5 == 0 ~ "I",
    I == 0 & II == 1 & I5 == 0 ~ "II",
    I == 0 & II == 0 & I5 == 1 ~ "I½ ",
    I == 1 & II == 1 & I5 == 0 ~ "I+II",
    I == 1 & II == 0 & I5 == 1 ~ "I+I½",
    I == 0 & II == 1 & I5 == 1 ~ "II+I½",
    I == 1 & II == 1 & I5 == 1 ~ "I+II+I½",
    TRUE ~ "Other Groups"  # Adjust this condition based on your specific grouping criteria
  ))

reshaped_data <- data.frame(group = venn_df2_mut$group, 
                            mut = venn_df2_mut$mut)


PieDonut(reshaped_data, aes(group, mut), title = "Shared Mutations")




# ---------------------------------------- heatmap of resistance  ---------------------------------------------------#

complete_data <- expand_grid(position = 1059:1345, mutation = variant_names)

# Merge the complete data with your existing data
ex14_type_resis_sum_complete <- merge(complete_data, ex14_type_resis_sum, by = c("position", "mutation"), all.x = TRUE)

# Fill missing values with 0
ex14_type_resis_sum_complete[is.na(ex14_type_resis_sum_complete$I), "I"] <- 0
ex14_type_resis_sum_complete[is.na(ex14_type_resis_sum_complete$II), "II"] <- 0
ex14_type_resis_sum_complete[is.na(ex14_type_resis_sum_complete$I5), "I5"] <- 0


order <- c('H', 'K', 'R', 'D', 'E','C','M','N','Q','S','T','A','I','L','V','W','F', 'Y','G','P')
variant_names <- c('H', 'K', 'R', 'D', 'E','C','M','N','Q','S','T','A','I','L','V','W','F','Y','G','P')
met_wt_fasta = "Met_wt2.fasta"

met_wt_fasta_con=file(met_wt_fasta, open="r")
met_wt_sequence = readLines(met_wt_fasta_con)[2]
close(met_wt_fasta_con)

met_wt_1 = str_split(substr(met_wt_sequence, 1, 144), '')[[1]]
met_wt_2 = str_split(substr(met_wt_sequence, 145, 287), '')[[1]]


# Create the Row1 plot with black border
Row1 <- ggplot(data = ex14_type_resis_sum_complete %>% filter(position %in% c(1059:1202)),
               aes(x = position, y = factor(mutation, level = order, labels = variant_names), fill = I)) +
  geom_tile(size = 0.2, color = "lightgray") +  # Add black border
  scale_fill_gradientn(colors = brewer.pal(9, "PuRd")) +
  scale_color_manual(values = c(NA, 'green')) +
  scale_x_continuous(breaks = seq(1059, 1202, by = 5),
                     expand = c(0, 0),
                     sec.axis = sec_axis(
                       trans = ~.,
                       name = "Sequence",
                       breaks = seq(1059, 1202),
                       labels = met_wt_1,
                       guide = derive()
                     )) +
  coord_fixed(ratio = 1) +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgray"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 4),
    axis.text = element_text(size = 5),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
    axis.ticks = element_blank(),
    legend.position = "none"
  ) + 
  labs(y = "Mutation", x = "Position")
plot(Row1)

Row1.5 =  ggplot(data = ex14_type_resis_sum_complete %>% filter(position %in% c(1203:1345)),
                 aes(x = position, y = factor(mutation, level = order, labels = variant_names), fill = I)) +
  geom_tile(size = 0.2, color = "lightgray") +  # Add black border
  scale_fill_gradientn(colors = brewer.pal(9, "PuRd")) +
  scale_color_manual(values = c(NA, 'green')) +
  scale_x_continuous(breaks = seq(1203, 1345, by = 5),
                     expand = c(0, 0),
                     sec.axis = sec_axis(
                       trans = ~.,
                       name = "Sequence",
                       breaks = seq(1203, 1345),
                       labels = met_wt_2,
                       guide = derive()
                     )) +
  coord_fixed(ratio = 1) +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgray"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 4),
    axis.text = element_text(size = 5),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
    axis.ticks = element_blank(),
    legend.position = "none"
  ) + 
  labs(y = "Mutation", x = "Position")
heatmap_plot <- ggarrange(Row1,Row1.5,nrow = 3, ncol = 1)
ggsave("Met_typeI_resis.pdf",height = 11 , width = 8.5, heatmap_plot)
plot(heatmap_plot)



Row2 <- ggplot(data = ex14_type_resis_sum_complete %>% filter(position %in% c(1059:1202)),
               aes(x = position, y = factor(mutation, level = order, labels = variant_names), fill = II)) +
  geom_tile(size = 0.2, color = "lightgray") +  # Add black border
  scale_fill_gradientn(colors = brewer.pal(9, "Blues")) +
  scale_color_manual(values = c(NA, 'green')) +
  scale_x_continuous(breaks = seq(1059, 1202, by = 5),
                     expand = c(0, 0),
                     sec.axis = sec_axis(
                       trans = ~.,
                       name = "Sequence",
                       breaks = seq(1059, 1202),
                       labels = met_wt_1,
                       guide = derive()
                     )) +
  coord_fixed(ratio = 1) +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgray"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 4),
    axis.text = element_text(size = 5),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
    axis.ticks = element_blank(),
    legend.position = "none"
  ) + 
  labs(y = "Mutation", x = "Position")
Row2.5 =  ggplot(data = ex14_type_resis_sum_complete %>% filter(position %in% c(1203:1345)),
                 aes(x = position, y = factor(mutation, level = order, labels = variant_names), fill = II)) +
  geom_tile(size = 0.2, color = "lightgray") +  # Add black border
  scale_fill_gradientn(colors = brewer.pal(9, "Blues")) +
  scale_color_manual(values = c(NA, 'green')) +
  scale_x_continuous(breaks = seq(1203, 1345, by = 5),
                     expand = c(0, 0),
                     sec.axis = sec_axis(
                       trans = ~.,
                       name = "Sequence",
                       breaks = seq(1203, 1345),
                       labels = met_wt_2,
                       guide = derive()
                     )) +
  coord_fixed(ratio = 1) +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgray"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 4),
    axis.text = element_text(size = 5),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
    axis.ticks = element_blank(),
    legend.position = "none"
  ) + 
  labs(y = "Mutation", x = "Position")
heatmap_plot <- ggarrange(Row2,Row2.5,nrow = 3, ncol = 1)
ggsave("Met_typeII_resis.pdf",height = 11 , width = 8.5, heatmap_plot)
plot(heatmap_plot)



Row3 <- ggplot(data = ex14_type_resis_sum_complete %>% filter(position %in% c(1059:1202)),
               aes(x = position, y = factor(mutation, level = order, labels = variant_names), fill = I5)) +
  geom_tile(size = 0.2, color = "lightgray") +  # Add black border
  scale_fill_gradientn(colors = brewer.pal(9, "Greens")) +
  scale_color_manual(values = c(NA, 'green')) +
  scale_x_continuous(breaks = seq(1059, 1202, by = 5),
                     expand = c(0, 0),
                     sec.axis = sec_axis(
                       trans = ~.,
                       name = "Sequence",
                       breaks = seq(1059, 1202),
                       labels = met_wt_1,
                       guide = derive()
                     )) +
  coord_fixed(ratio = 1) +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgray"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 4),
    axis.text = element_text(size = 5),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
    axis.ticks = element_blank(),
    legend.position = "none"
  ) + 
  labs(y = "Mutation", x = "Position")
Row3.5 =  ggplot(data = ex14_type_resis_sum_complete %>% filter(position %in% c(1203:1345)),
                 aes(x = position, y = factor(mutation, level = order, labels = variant_names), fill = I5)) +
  geom_tile(size = 0.2, color = "lightgray") +  # Add black border
  scale_fill_gradientn(colors = brewer.pal(9, "Greens")) +
  scale_color_manual(values = c(NA, 'green')) +
  scale_x_continuous(breaks = seq(1203, 1345, by = 5),
                     expand = c(0, 0),
                     sec.axis = sec_axis(
                       trans = ~.,
                       name = "Sequence",
                       breaks = seq(1203, 1345),
                       labels = met_wt_2,
                       guide = derive()
                     )) +
  coord_fixed(ratio = 1) +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgray"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 4),
    axis.text = element_text(size = 5),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
    axis.ticks = element_blank(),
    legend.position = "none"
  ) + 
  labs(y = "Mutation", x = "Position")
heatmap_plot <- ggarrange(Row3,Row3.5,nrow = 3, ncol = 1)
ggsave("Met_typeI.5_resis.pdf",height = 11 , width = 8.5, heatmap_plot)
plot(heatmap_plot)



#----------------------------------- collapsed heatmaps to highlight hotspots -------------------------------------#

### this is the code block to make sum of sums heatmaps for the inhibitor resistance hotspots in figure 3

met_wt_1 = str_split(substr(met_wt_sequence, 1, 144), '')[[1]]
met_wt_2 = str_split(substr(met_wt_sequence, 145, 287), '')[[1]]


resis_map_1 <- ggplot(data = position_sum  %>%filter(position %in% c(1059:1202)),
                      aes(x = position, y = 1, fill = I)) +
  geom_tile(size = 0.2, color = "lightgray") + 
  scale_fill_gradientn(colors = brewer.pal(9, "PuRd"), limits = c(0, 23)) +
  scale_color_manual(values = c(NA,'green')) +
  scale_x_continuous(breaks = seq(1059,1202, by = 5),
                     expand = c(0,0),
                     sec.axis = sec_axis(
                       trans = ~.,
                       name = "Sequence",
                       breaks = seq(1059,1202),
                       labels = met_wt_1,
                       guide = derive()
                     )) +
  coord_fixed(ratio = 1) +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgray"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 4),
    axis.text = element_text(size = 4),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
    axis.ticks = element_blank(),
    #legend.position="none"
    ) + 
  labs(y = "Mutation", x = "Position")

resis_map_1.1 <- ggplot(data = position_sum  %>%filter(position %in% c(1203:1345)),
                         aes(x = position, y = 1, fill = I)) +
  geom_tile(size = 0.2, color = "lightgray") +  
  scale_fill_gradientn(colors = brewer.pal(9, "PuRd"), limits = c(0, 23)) +
  scale_color_manual(values = c(NA,'green')) +
  scale_x_continuous(breaks = seq(1203,1345, by = 5),
                     expand = c(0,0),
                     sec.axis = sec_axis(
                       trans = ~.,
                       name = "Sequence",
                       breaks = seq(1203,1345),
                       labels = met_wt_2,
                       guide = derive()
                     )) +
  coord_fixed(ratio = 1) +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgray"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 4),
    axis.text = element_text(size = 4),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
    axis.ticks = element_blank(),
    #legend.position="none", 
    legend.key.size = unit(0.4, 'cm')
  ) + 
  labs(y = "Mutation", x = "Position")


resis_map_2 <- ggplot(data = position_sum  %>%filter(position %in% c(1059:1202)),
                      aes(x = position, y = 1, fill = II)) +
  geom_tile(size = 0.2, color = "lightgray") +  
  scale_fill_gradientn(colors = brewer.pal(9, "Blues"), limits = c(0, 13)) +
  scale_color_manual(values = c(NA,'green')) +
  scale_x_continuous(breaks = seq(1059,1202, by = 5),
                     expand = c(0,0),
                     sec.axis = sec_axis(
                       trans = ~.,
                       name = "Sequence",
                       breaks = seq(1059,1202),
                       labels = met_wt_1,
                       guide = derive()
                     )) +
  coord_fixed(ratio = 1) +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgray"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 4),
    axis.text = element_text(size = 4),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
    axis.ticks = element_blank(),
    #legend.position="none"
    ) + 
  labs(y = "Mutation", x = "Position")

resis_map_2.2 <- ggplot(data = position_sum  %>%filter(position %in% c(1203:1345)),
                        aes(x = position, y = 1, fill = II)) +
  geom_tile(size = 0.2, color = "lightgray") + 
  scale_fill_gradientn(colors = brewer.pal(9, "Blues"), limits = c(0, 13)) +
  scale_color_manual(values = c(NA,'green')) +
  scale_x_continuous(breaks = seq(1203,1345, by = 5),
                     expand = c(0,0),
                     sec.axis = sec_axis(
                       trans = ~.,
                       name = "Sequence",
                       breaks = seq(1203,1345),
                       labels = met_wt_2,
                       guide = derive()
                     )) +
  coord_fixed(ratio = 1) +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgray"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 4),
    axis.text = element_text(size = 4),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
    axis.ticks = element_blank(),
    #legend.position="none", 
    legend.key.size = unit(0.4, 'cm')
  ) + 
  labs(y = "Mutation", x = "Position")


resis_map_3 <- ggplot(data = position_sum  %>%filter(position %in% c(1059:1202)),
                      aes(x = position, y = 1, fill = I5)) +
  geom_tile(size = 0.2, color = "lightgray") +  
  scale_fill_gradientn(colors = brewer.pal(9, "Greens"), limits = c(0, 10)) +
  scale_color_manual(values = c(NA,'green')) +
  scale_x_continuous(breaks = seq(1059,1202, by = 5),
                     expand = c(0,0),
                     sec.axis = sec_axis(
                       trans = ~.,
                       name = "Sequence",
                       breaks = seq(1059,1202),
                       labels = met_wt_1,
                       guide = derive()
                     )) +
  coord_fixed(ratio = 1) +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgray"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 4),
    axis.text = element_text(size = 4),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
    axis.ticks = element_blank(),
    #legend.position="none"
  ) + 
  labs(y = "Mutation", x = "Position")

resis_map_3.3<- ggplot(data = position_sum  %>%filter(position %in% c(1203:1345)),
                        aes(x = position, y = 1, fill = I5)) +
  geom_tile(size = 0.2, color = "lightgray") +  
  scale_fill_gradientn(colors = brewer.pal(9, "Greens"), limits = c(0, 10)) +
  scale_color_manual(values = c(NA,'green')) +
  scale_x_continuous(breaks = seq(1203,1345, by = 5),
                     expand = c(0,0),
                     sec.axis = sec_axis(
                       trans = ~.,
                       name = "Sequence",
                       breaks = seq(1203,1345),
                       labels = met_wt_2,
                       guide = derive()
                     )) +
  coord_fixed(ratio = 1) +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgray"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 4),
    axis.text = element_text(size = 4),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
    axis.ticks = element_blank(),
    #legend.position="none", 
    legend.key.size = unit(0.4, 'cm')
  ) + 
  labs(y = "Mutation", x = "Position")


heatmap_plot <- ggarrange(resis_map_1, resis_map_1.1, 
                          resis_map_2, resis_map_2.2, 
                          resis_map_3, resis_map_3.3,
                          nrow = 6, ncol = 1)
plot(heatmap_plot)
ggsave("Type1_2_sinlge_resis_hotspot.pdf",height = 11 , width = 8.5, heatmap_plot)



#----------------------------------------- Project hotspot sums onto structures -------------------------------------#

avg_ex14_1 = data.frame(pos = position_sum$position, mean = position_sum$I)
avg_ex14_2 = data.frame(pos = position_sum$position, mean = position_sum$II)
avg_ex14_1.5 = data.frame(pos = position_sum$position, mean = position_sum$I5)


crizo_2WGJ <- read.pdb("2WGJ") 
x = map_scores_pdb(crizo_2WGJ, avg_ex14_1, "mean")
write.pdb(x, file="resis_hotspot_1.pdb")

mere_4EEV <- read.pdb("4EEV")
x = map_scores_pdb(mere_4EEV, avg_ex14_2, "mean")
write.pdb(x, file="resis_hotspot_2.pdb")

amg458_5T3Q <- read.pdb("5T3Q")
x = map_scores_pdb(amg458_5T3Q , avg_ex14_1.5, "mean")
write.pdb(x, file="resis_hotspot_1.5.pdb")


#################################################################################################################
# ----------------------------------- heatmap of reistance ONLY condensed  -----------------------------------------#
#################################################################################################################


# these are the smaller heatmaps shown in figure 3, that only show the resistance mutations 


complete_data <- expand_grid(position = 1059:1345, mutation = variant_names)

# Merge the complete data with your existing data
ex14_type_resis_sum_complete <- merge(complete_data, ex14_type_resis_sum, by = c("position", "mutation"), all.x = TRUE)

met_type_resis_sum_complete <- merge(complete_data, met_type_resis_sum, by = c("position", "mutation"), all.x = TRUE)



# Fill missing values with 0
ex14_type_resis_sum_complete[is.na(ex14_type_resis_sum_complete$I), "I"] <- 0
ex14_type_resis_sum_complete[is.na(ex14_type_resis_sum_complete$II), "II"] <- 0
ex14_type_resis_sum_complete[is.na(ex14_type_resis_sum_complete$I5), "I5"] <- 0



order <- c('H', 'K', 'R', 'D', 'E','C','M','N','Q','S','T','A','I','L','V','W','F', 'Y','G','P')
variant_names <- c('H', 'K', 'R', 'D', 'E','C','M','N','Q','S','T','A','I','L','V','W','F','Y','G','P')
#met_wt_fasta = "Met_wt2.fasta"

#met_wt_fasta_con=file(met_wt_fasta, open="r")
#met_wt_sequence = readLines(met_wt_fasta_con)[2]
#close(met_wt_fasta_con)

#met_wt_1 = str_split(substr(met_wt_sequence, 1, 287), '')[[1]]


#----- type I -----#

I_data <- data.frame(position = ex14_type_resis_sum_complete$position,
                     mutation = ex14_type_resis_sum_complete$mutation,
                     I = ex14_type_resis_sum_complete$I)

I_data <- I_data %>% filter(I != 0) # remove all zero values 

I_data_expanded <- I_data %>%
  filter(mutation != "Stop") %>%
  complete(position, mutation = factor(variant_names, levels = variant_names), fill = list(I = 0))


I_data_expanded$mutation <- factor(I_data_expanded$mutation, levels = variant_names)


I_data_expanded <- I_data_expanded[order(I_data_expanded$position), ]

Row1 <- ggplot(data = I_data_expanded,
               aes(x = factor(position, levels = unique(I_data_expanded$position)), 
                   y = factor(mutation, level = order, labels = variant_names), 
                   fill = I)) +
  geom_tile(size = 0.7, color = "lightgray") +
  scale_fill_gradientn(colors = brewer.pal(9, "PuRd")) +
  scale_color_manual(values = c(NA, 'green')) +
  coord_fixed(ratio = 1) +
  labs(y = "Mutation", x = "Position")+
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgray"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    axis.text = element_text(size = 7),
    axis.ticks = element_blank(),
    #legend.position="right"
  ) 
plot(Row1)

# ---- type II ----#

II_data <- data.frame(position = ex14_type_resis_sum_complete$position,
                     mutation = ex14_type_resis_sum_complete$mutation,
                     II = ex14_type_resis_sum_complete$II)

II_data <- II_data %>% filter(II != 0) # remove all zero values 

II_data_expanded <- II_data %>%
  filter(mutation != "Stop") %>%
  complete(position, mutation = factor(variant_names, levels = variant_names), fill = list(II = 0))


II_data_expanded$mutation <- factor(II_data_expanded$mutation, levels = variant_names)


II_data_expanded <- II_data_expanded[order(II_data_expanded$position), ]

Row2 <- ggplot(data = II_data_expanded,
               aes(x = factor(position, levels = unique(II_data_expanded$position)), 
                   y = factor(mutation, level = order, labels = variant_names), 
                   fill = II)) +
  geom_tile(size = 0.7, color = "lightgray") +
  scale_fill_gradientn(colors = brewer.pal(9, "Blues")) +
  scale_color_manual(values = c(NA, 'green')) +
  coord_fixed(ratio = 1) +
  labs(y = "Mutation", x = "Position")+
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgray"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    axis.text = element_text(size = 7),
    axis.ticks = element_blank(),
    #legend.position="right"
  ) 
plot(Row2)

#--- type 1.5 ---#

I5_data <- data.frame(position = ex14_type_resis_sum_complete$position,
                      mutation = ex14_type_resis_sum_complete$mutation,
                      I5 = ex14_type_resis_sum_complete$I5)

I5_data <- I5_data %>% filter(I5 != 0) # remove all zero values 

I5_data_expanded <- I5_data %>%
  filter(mutation != "Stop") %>%
  complete(position, mutation = factor(variant_names, levels = variant_names), fill = list(I5 = 0))


I5_data_expanded$mutation <- factor(I5_data_expanded$mutation, levels = variant_names)


I5_data_expanded <- I5_data_expanded[order(I5_data_expanded$position), ]

Row3 <- ggplot(data = I5_data_expanded,
               aes(x = factor(position, levels = unique(I5_data_expanded$position)), 
                   y = factor(mutation, level = order, labels = variant_names), 
                   fill = I5)) +
  geom_tile(size = 0.6, color = "lightgray") +
  scale_fill_gradientn(colors = brewer.pal(9, "Greens")) +
  scale_color_manual(values = c(NA, 'green')) +
  coord_fixed(ratio = 1) +
  labs(y = "Mutation", x = "Position")+
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgray"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    axis.text = element_text(size = 7),
    axis.ticks = element_blank(),
    #legend.position="right"
  ) 
plot(Row3)

library(cowplot)
plot_grid(Row1, Row2, Row3)

heatmap_plot <- ggarrange(Row1,Row2,Row3, nrow = 1, ncol = 3)
plot(heatmap_plot)
ggsave("resistance_hotspot_short_heatmap.pdf",height = 11 , width = 8.5, heatmap_plot)



###############################################################################################
#-------------------------------------- MET del Ex14 -----------------------------------------#
###############################################################################################

# ----------------------- heatmap of resistance MET del Ex14  --------------------------------#

met_rosace_resistant_pval0.5 <- met_rosace_norm_scores_filtered %>%
  group_by(variants) %>%
  filter(
    (delta.score >= 0 & pval0.5.up <= 0.1 & mean >= 0 & mean[match("DMSO", key)] <= 0)) %>%
  ungroup()

# generate data frames 
met_scores_type <- met_rosace_resistant_pval0.5 %>%
  mutate(inhib = ifelse(key %in% key_list, group_list[match(key, key_list)], NA)) 

met_resis_dt <- as.data.table(met_scores_type)
#met_type_resis_sum <- met_resis_sensitive_dt[, .(count = sum(pval0.5.up <= 0.1, na.rm = TRUE)), by = .(variants, position, mutation, inhib)][, dcast(.SD, ... ~ inhib, value.var = "count", fill = 0)]

met_type_resis_sum <- met_resis_dt[, .(count = .N), by = .(variants, position, mutation, inhib)][, dcast(.SD, variants + position + mutation ~ inhib, value.var = "count", fill = 0)]

met_resis_counts <- data.frame(position = met_type_resis_sum$position,
                           I = met_type_resis_sum$I,
                           II = met_type_resis_sum$II,
                           I5 = met_type_resis_sum$I5)

met_position_sum <- met_resis_counts %>%
  group_by(position) %>%
  summarise_all(list(~sum(., na.rm = TRUE))) %>%
  arrange(position)


met_all_positions <- data.frame(position = seq(1059, 1345, by = 1))

# Left join your existing data frame with the full positions data frame
met_position_sum <- met_all_positions %>%
  left_join(met_position_sum, by = "position") %>%
  mutate(I = ifelse(is.na(I), 0, I),
         II = ifelse(is.na(II), 0, II),
         I5 = ifelse(is.na(I5), 0, I5))
met_type_resis_sum_complete <- merge(complete_data, met_type_resis_sum, by = c("position", "mutation"), all.x = TRUE)

met_type_resis_sum_complete[is.na(met_type_resis_sum_complete$I), "I"] <- 0
met_type_resis_sum_complete[is.na(met_type_resis_sum_complete$II), "II"] <- 0
met_type_resis_sum_complete[is.na(met_type_resis_sum_complete$I5), "I5"] <- 0




# Create the Row1 plot with black border
Row1 <- ggplot(data = met_type_resis_sum_complete %>% filter(position %in% c(1059:1202)),
               aes(x = position, y = factor(mutation, level = order, labels = variant_names), fill = I)) +
  geom_tile(size = 0.2, color = "lightgray") +  # Add black border
  scale_fill_gradientn(colors = brewer.pal(9, "PuRd")) +
  scale_color_manual(values = c(NA, 'green')) +
  scale_x_continuous(breaks = seq(1059, 1202, by = 5),
                     expand = c(0, 0),
                     sec.axis = sec_axis(
                       trans = ~.,
                       name = "Sequence",
                       breaks = seq(1059, 1202),
                       labels = met_wt_1,
                       guide = derive()
                     )) +
  coord_fixed(ratio = 1) +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgray"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 4),
    axis.text = element_text(size = 5),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
    axis.ticks = element_blank(),
    #legend.position = "none"
  ) + 
  labs(y = "Mutation", x = "Position")
plot(Row1)

Row1.5 =  ggplot(data = met_type_resis_sum_complete %>% filter(position %in% c(1203:1345)),
                 aes(x = position, y = factor(mutation, level = order, labels = variant_names), fill = I)) +
  geom_tile(size = 0.2, color = "lightgray") +  # Add black border
  scale_fill_gradientn(colors = brewer.pal(9, "PuRd")) +
  scale_color_manual(values = c(NA, 'green')) +
  scale_x_continuous(breaks = seq(1203, 1345, by = 5),
                     expand = c(0, 0),
                     sec.axis = sec_axis(
                       trans = ~.,
                       name = "Sequence",
                       breaks = seq(1203, 1345),
                       labels = met_wt_2,
                       guide = derive()
                     )) +
  coord_fixed(ratio = 1) +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgray"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 4),
    axis.text = element_text(size = 5),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
    axis.ticks = element_blank(),
    #legend.position = "none"
  ) + 
  labs(y = "Mutation", x = "Position")
heatmap_plot <- ggarrange(Row1,Row1.5,nrow = 3, ncol = 1)
ggsave("Metdelex14_typeI_resis.pdf",height = 11 , width = 8.5, heatmap_plot)
plot(heatmap_plot)



Row2 <- ggplot(data = met_type_resis_sum_complete %>% filter(position %in% c(1059:1202)),
               aes(x = position, y = factor(mutation, level = order, labels = variant_names), fill = II)) +
  geom_tile(size = 0.2, color = "lightgray") +  # Add black border
  scale_fill_gradientn(colors = brewer.pal(9, "Blues")) +
  scale_color_manual(values = c(NA, 'green')) +
  scale_x_continuous(breaks = seq(1059, 1202, by = 5),
                     expand = c(0, 0),
                     sec.axis = sec_axis(
                       trans = ~.,
                       name = "Sequence",
                       breaks = seq(1059, 1202),
                       labels = met_wt_1,
                       guide = derive()
                     )) +
  coord_fixed(ratio = 1) +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgray"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 4),
    axis.text = element_text(size = 5),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
    axis.ticks = element_blank(),
    #legend.position = "none"
  ) + 
  labs(y = "Mutation", x = "Position")
Row2.5 =  ggplot(data = met_type_resis_sum_complete %>% filter(position %in% c(1203:1345)),
                 aes(x = position, y = factor(mutation, level = order, labels = variant_names), fill = II)) +
  geom_tile(size = 0.2, color = "lightgray") +  # Add black border
  scale_fill_gradientn(colors = brewer.pal(9, "Blues")) +
  scale_color_manual(values = c(NA, 'green')) +
  scale_x_continuous(breaks = seq(1203, 1345, by = 5),
                     expand = c(0, 0),
                     sec.axis = sec_axis(
                       trans = ~.,
                       name = "Sequence",
                       breaks = seq(1203, 1345),
                       labels = met_wt_2,
                       guide = derive()
                     )) +
  coord_fixed(ratio = 1) +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgray"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 4),
    axis.text = element_text(size = 5),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
    axis.ticks = element_blank(),
    #legend.position = "none"
  ) + 
  labs(y = "Mutation", x = "Position")
heatmap_plot <- ggarrange(Row2,Row2.5,nrow = 3, ncol = 1)
ggsave("Metdelex14_typeII_resis.pdf",height = 11 , width = 8.5, heatmap_plot)
plot(heatmap_plot)



Row3 <- ggplot(data = met_type_resis_sum_complete %>% filter(position %in% c(1059:1202)),
               aes(x = position, y = factor(mutation, level = order, labels = variant_names), fill = I5)) +
  geom_tile(size = 0.2, color = "lightgray") +  # Add black border
  scale_fill_gradientn(colors = brewer.pal(9, "Greens")) +
  scale_color_manual(values = c(NA, 'green')) +
  scale_x_continuous(breaks = seq(1059, 1202, by = 5),
                     expand = c(0, 0),
                     sec.axis = sec_axis(
                       trans = ~.,
                       name = "Sequence",
                       breaks = seq(1059, 1202),
                       labels = met_wt_1,
                       guide = derive()
                     )) +
  coord_fixed(ratio = 1) +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgray"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 4),
    axis.text = element_text(size = 5),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
    axis.ticks = element_blank(),
    #legend.position = "none"
  ) + 
  labs(y = "Mutation", x = "Position")
Row3.5 =  ggplot(data = met_type_resis_sum_complete %>% filter(position %in% c(1203:1345)),
                 aes(x = position, y = factor(mutation, level = order, labels = variant_names), fill = I5)) +
  geom_tile(size = 0.2, color = "lightgray") +  # Add black border
  scale_fill_gradientn(colors = brewer.pal(9, "Greens")) +
  scale_color_manual(values = c(NA, 'green')) +
  scale_x_continuous(breaks = seq(1203, 1345, by = 5),
                     expand = c(0, 0),
                     sec.axis = sec_axis(
                       trans = ~.,
                       name = "Sequence",
                       breaks = seq(1203, 1345),
                       labels = met_wt_2,
                       guide = derive()
                     )) +
  coord_fixed(ratio = 1) +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgray"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 4),
    axis.text = element_text(size = 5),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
    axis.ticks = element_blank(),
    #legend.position = "none"
  ) + 
  labs(y = "Mutation", x = "Position")
heatmap_plot <- ggarrange(Row3,Row3.5,nrow = 3, ncol = 1)
ggsave("Metdelex14_typeI.5_resis.pdf",height = 11 , width = 8.5, heatmap_plot)
plot(heatmap_plot)



#####################################################################################################

#----------------------------- Shared resistance Doughnut chart -------------------------------------#

# generate data frame that keeps inhibitor association but then plot counts into a doughnut chart
# this should provide comparative information on mutations that are shared and bespoke to inhibitor classes 


library(webr)

### build data frame 
resis_freq <- data.frame(inhib = ex14_rosace_resistant_pval0.5$key, 
                               variants = ex14_rosace_resistant_pval0.5$variants)

# Filter out rows with specific values in 'inhib'
type1_resis_freq <- resis_freq %>% filter(!(inhib %in% c("Cabo", "Mere", "Gle", "A458")))
type2_resis_freq <- resis_freq %>% filter(!(inhib %in% c("Crizo","Camp", "Tepo", "Glu","Savo", "NVP","A458")))

##__________________________________________type 1 ________________________________#

variant_counts_1 <- table(type1_resis_freq$variant)
result_df1 <- data.frame(
  variant = names(variant_counts_1),
  count = as.numeric(variant_counts_1),
  occur_multi = variant_counts_1 > 1,
  inhibs = character(length(variant_counts_1))  # Initialize 'inhibs' as character
)

# If 'occur_multi' is TRUE, annotate all associated 'inhibs' for each 'variant'
result_df1$inhibs[result_df1$occur_multi] <- sapply(result_df1$variant[result_df1$occur_multi], function(v) {
  paste(unique(type1_resis_freq$inhib[type1_resis_freq$variant == v]), collapse = ", ")
})

# If 'occur_multi' is FALSE, annotate a single 'inhibs' value for each 'variant'
result_df1$inhibs[!result_df1$occur_multi] <- sapply(result_df1$variant[!result_df1$occur_multi], function(v) {
  unique(type1_resis_freq$inhib[type1_resis_freq$variant == v])
})

# Fill in 'occur_multi' column with "Shared" and "Unique"
result_df1$occur_multi <- ifelse(result_df1$occur_multi, "Shared", "Unique")

result_df1 <- result_df1 %>%
  mutate(variant = paste0(sub(".*([A-Za-z]+)(\\d+)([A-Za-z]+).*", "\\1", variant),
                          as.numeric(sub(".*([A-Za-z]+)(\\d+)([A-Za-z]+).*", "\\2", variant)) + 1058,
                          sub(".*([A-Za-z]+)(\\d+)([A-Za-z]+).*", "\\3", variant)))
# Print the resulting data frame
print(result_df1)


# Filter data for occurrences where occur_multi is "Shared"
result_shared_df <- result_df1[result_df1$occur_multi == "Shared", ]
result_unique_df <- result_df1[result_df1$occur_multi == "Unique", ]


#------------------------------ Pie Doughnut Charts --------------------------------------#

# create the doughnut chart!
PieDonut(result_shared_df, aes(inhibs, variant, count = count))
PieDonut(result_unique_df, aes(inhibs, variant, count = count))


# Create the doughnut chart with counts using ggplot2
ggplot(result_shared_df, aes(x = 1, y = count, fill = inhibs)) +
  geom_bar(stat = "identity", width = 6, color = "white") +
  geom_bar(stat = "identity", width = 1.2, color = "white", fill = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = count), position = position_stack(vjust = 0.5), size = 4, color = "black") +
  geom_text(aes(label = paste(variant, sep = "\n")), position = position_stack(vjust = 0.5, reverse = TRUE), size = 3, color = "black") +
  theme_void() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 14),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  ) +
  labs(title = "Shared Occurrence of Variants with Inhibitors", fill = "Inhibitor") +
  guides(fill = guide_legend(title = "Inhibitor"))


library(tibble)
library(dplyr)
library(ggplot2)


# shared resis for type 1

lvl0 <- tibble(variant = "Parent", value = 0, level = 0, fill = NA)

lvl1 <- result_df1 %>%
  group_by(occur_multi) %>%
  summarise(value = sum(count)) %>%
  ungroup() %>%
  mutate(level = 1, fill = as.character(occur_multi))

lvl2 <- result_df1 %>%
  group_by(variant, occur_multi) %>%
  summarise(value = sum(count)) %>%
  ungroup() %>%
  mutate(level = 2)

result_df_combined <- bind_rows(lvl0, lvl1, lvl2) %>%
  arrange(fill, variant) %>%
  mutate(variant = factor(variant, levels = unique(variant)),
         fill = ifelse(occur_multi == "Shared", "Shared", "Unique")) %>%
  mutate(level = as.factor(level))

ggplot(result_df_combined, aes(x = level, y = value, fill = fill, alpha = level)) +
  geom_bar(data = result_df_combined[result_df_combined$level == 1, ], 
           aes(x = level, y = value, fill = fill, alpha = level), 
           stat = "identity", width = 0.7, color = "gray90", size = 0.25, position = position_stack()) +
  geom_bar(data = result_df_combined[result_df_combined$level == 2, ], 
           aes(x = level, y = value, fill = fill, alpha = level), 
           stat = "identity", width = 1.2, color = "gray90", size = 0.25, position = position_stack()) +
  geom_text(aes(label = paste(variant, "\n", value)), size = 5, position = position_stack(vjust = 0.5)) +  # Include raw count information
  coord_polar(theta = "y") +
  scale_alpha_manual(values = c("0" = 0, "1" = 1, "2" = 0.7), guide = FALSE) +
  scale_x_discrete(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  scale_fill_manual(values = c("Shared" = "deeppink1", "Unique" = "plum2"), na.translate = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_minimal()



## shared resis for type 2

variant_counts_2 <- table(type2_resis_freq$variant)
result_df2 <- data.frame(
  variant = names(variant_counts_2),
  count = as.numeric(variant_counts_2),
  occur_multi = variant_counts_2 > 1,
  inhibs = character(length(variant_counts_2))  # Initialize 'inhibs' as character
)

# If 'occur_multi' is TRUE, annotate all associated 'inhibs' for each 'variant'
result_df2$inhibs[result_df2$occur_multi] <- sapply(result_df2$variant[result_df2$occur_multi], function(v) {
  paste(unique(type2_resis_freq$inhib[type2_resis_freq$variant == v]), collapse = ", ")
})

# If 'occur_multi' is FALSE, annotate a single 'inhibs' value for each 'variant'
result_df2$inhibs[!result_df2$occur_multi] <- sapply(result_df2$variant[!result_df2$occur_multi], function(v) {
  unique(type2_resis_freq$inhib[type2_resis_freq$variant == v])
})

# Fill in 'occur_multi' column with "Shared" and "Unique"
result_df2$occur_multi <- ifelse(result_df2$occur_multi, "Shared", "Unique")

result_df2 <- result_df2 %>%
  mutate(variant = paste0(sub(".*([A-Za-z]+)(\\d+)([A-Za-z]+).*", "\\1", variant),
                          as.numeric(sub(".*([A-Za-z]+)(\\d+)([A-Za-z]+).*", "\\2", variant)) + 1058,
                          sub(".*([A-Za-z]+)(\\d+)([A-Za-z]+).*", "\\3", variant)))
# Print the resulting data frame
print(result_df2)


# Filter data for occurrences where occur_multi is "Shared"
result_shared_df <- result_df2[result_df2$occur_multi == "Shared", ]
result_unique_df <- result_df2[result_df2$occur_multi == "Unique", ]
# create the doughnut chart!
PieDonut(result_shared_df, aes(inhibs, variant, count = count))
PieDonut(result_unique_df, aes(inhibs, variant, count = count))


