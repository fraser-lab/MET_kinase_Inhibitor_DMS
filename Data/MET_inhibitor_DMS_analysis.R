#####################################################################################################
# Code by: Gabriella Estevam @ UCSF
#
# DMS analysis and plots of TPR-MET +/- Exon14 under inhibitor conditions 
# Utilizes "dms_analysis_utilities.R" and "MET_Inhibitor_V1.Rmd" code written by Chris MacDonald @ UCSF 
#
#####################################################################################################


#####################################################################################################
###---------------------------------------- load libraries ------------------------------------------
#####################################################################################################

library(ggplot2)
library(gglorenz)
library(ggpubr)
library(ggridges)
library(readr)
library(stringr)
library(tidyverse)
library(plyr)
library(dplyr)
library(ggridges)
library(ggbraid)
library(tidyquant)
library(bio3d)
source("dms_analysis_utilities.R")
library(scales)
library(colorspace)
library(ggpubr)
library(rstatix)
library(patchwork)
library(cowplot)
library(dendextend)
library("RColorBrewer")
library(gghighlight)
library(stringdist)
library(forcats)

# color palette
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#0072B2","#F0E442","#D55E00", "#CC79A7")

#####################################################################################################
###------------------------------------ define initial data frames ----------------------------------
#####################################################################################################

require(data.table)

### writes this data frame into a permanent data file - this should be a ONE time run, then commented out

#ex14_inhib_scores = met_ex_scores
#met_inhib_scores = met_scores

#ex14_inhib_scores = ex14_inhib_scores %>% filter(mutation_type != "X")
#met_inhib_scores  = met_inhib_scores %>% filter(mutation_type != "X")


##################### Enrich2 scores #####################
### open data file as data frame for all raw inhibitor scores generated from ENRICH2 
ex14_inhib_scores <- read.csv("ex14_inhib_scores.csv")#+Ex14
met_inhib_scores <- read.csv("met_inhib_scores.csv")#-Ex14 



##################### Rosace scores #####################
### open data file as data frame for all raw inhibitor ALL scores generated from ROSACE
met_rosace_scores <- as.data.frame(fread("EX_rosace_effect_all.tsv")) # here met is delta exon14
met_rosace_scores <- met_rosace_scores  %>% mutate(position = position + 1058) 
ex14_rosace_scores <- as.data.frame(fread("WT_rosace_effect_all.tsv")) # here ex14 means MET+Ex14 
ex14_rosace_scores <- ex14_rosace_scores  %>% mutate(position = position + 1058) 


### open data file as data frame for all raw inhibitor RESISTANCE scores generated from ROSACE
met_rosace_resis_scores <- as.data.frame(fread("EX_inhibitor_resistance_count.tsv")) # here met is delta exon14
met_rosace_resis_scores <- met_rosace_resis_scores  %>% mutate(position = position + 1058) 
ex14_rosace_resis_scores <- as.data.frame(fread("WT_inhibitor_resistance_count.tsv")) # here ex14 means MET+Ex14 
ex14_rosace_resis_scores <- ex14_rosace_resis_scores %>% mutate(position = position + 1058) 

# these are scores with statistically determined large differences between drug and DMSO 
# from ROSACE 
met_roscace_DMSO_diff_scores <- as.data.frame(fread("EX14_drug_compare_DMSO_sig_all.tsv")) 
met_roscace_DMSO_diff_scores  <- met_roscace_DMSO_diff_scores %>% mutate(position = position + 1058) 
ex14_roscace_DMSO_diff_scores  <- as.data.frame(fread("WT_drug_compare_DMSO_sig_all.tsv")) 
ex14_roscace_DMSO_diff_scores   <- ex14_roscace_DMSO_diff_scores  %>% mutate(position = position + 1058) 
  




### create a new data frame where the scores for the DMSO control is subtracted from each inhibitor 
# this is a SINGLE data frame with delta score, standard errors, etc.
ex14_delta_scores <- data.frame(hgvs=ex14_inhib_scores$hgvs,
                                pos=ex14_inhib_scores$pos, 
                                variants=ex14_inhib_scores$variants,
                                ex14_DMSO_score= (ex14_inhib_scores$DMSO_score),
                                ex14_DMSO_SE = (ex14_inhib_scores$DMSO_SE),
                                ex14_AMG458_delta = (ex14_inhib_scores$A458_score)-(ex14_inhib_scores$DMSO_score),
                                ex14_AMG458_SE = ex14_inhib_scores$A458_SE,
                                ex14_Crizo_delta = (ex14_inhib_scores$Crizo_score)-(ex14_inhib_scores$DMSO_score),
                                ex14_Crizo_SE = ex14_inhib_scores$Crizo_SE,
                                ex14_Cabo_delta = (ex14_inhib_scores$Cabo_score)-(ex14_inhib_scores$DMSO_score),
                                ex14_Cabo_SE = ex14_inhib_scores$Cabo_SE,
                                ex14_Camp_delta = (ex14_inhib_scores$Camp_score)-(ex14_inhib_scores$DMSO_score),
                                ex14_Camp_SE = ex14_inhib_scores$Camp_SE,
                                ex14_Tepo_delta =(ex14_inhib_scores$Tepo_score)-(ex14_inhib_scores$DMSO_score),
                                ex14_Tepo_SE = ex14_inhib_scores$Tepo_SE,
                                ex14_Tiv_delta = (ex14_inhib_scores$Tiv_score)-(ex14_inhib_scores$DMSO_score),
                                ex14_Tiv_SE = ex14_inhib_scores$Tiv_SE,
                                ex14_Gle_delta = (ex14_inhib_scores$Gle_score)-(ex14_inhib_scores$DMSO_score),
                                ex14_GleSE = ex14_inhib_scores$Gle_SE,
                                ex14_Glu_delta = (ex14_inhib_scores$Glu_score)-(ex14_inhib_scores$DMSO_score),
                                ex14_Glu_SE = ex14_inhib_scores$Glu_SE,
                                ex14_Savo_delta = (ex14_inhib_scores$Savo_score)-(ex14_inhib_scores$DMSO_score),
                                ex14_Savo_SE = ex14_inhib_scores$Savo_SE,
                                ex14_Mere_delta = (ex14_inhib_scores$Mere_score)-(ex14_inhib_scores$DMSO_score),
                                ex14_Mere_SE = ex14_inhib_scores$Mere_SE,
                                ex14_NVP_delta = (ex14_inhib_scores$NVP_score)-(ex14_inhib_scores$DMSO_score),
                                ex14_NVP_SE = ex14_inhib_scores$NVP_SE
                                )



met_delta_scores  <- data.frame(hgvs=met_inhib_scores$hgvs,
                                pos=met_inhib_scores$pos, 
                                variants=met_inhib_scores$variants,
                                met_DMSO_score= (met_inhib_scores$DMSO_score),
                                met_DMSO_SE = (met_inhib_scores$DMSO_SE),
                                met_AMG458_delta = (met_inhib_scores$A458_score)-(met_inhib_scores$DMSO_score),
                                met_AMG458_SE = met_inhib_scores$A458_SE,
                                met_Crizo_delta = (met_inhib_scores$Crizo_score)-(met_inhib_scores$DMSO_score),
                                met_Crizo_SE = met_inhib_scores$Crizo_SE,
                                met_Cabo_delta = (met_inhib_scores$Cabo_score)-(met_inhib_scores$DMSO_score),
                                met_Cabo_SE = met_inhib_scores$Cabo_SE,
                                met_Camp_delta = (met_inhib_scores$Camp_score)-(met_inhib_scores$DMSO_score),
                                met_Camp_SE = met_inhib_scores$Camp_SE,
                                met_Gle_delta = (met_inhib_scores$Gle_score)-(met_inhib_scores$DMSO_score),
                                met_GleSE = met_inhib_scores$Gle_SE,
                                met_Glu_delta = (met_inhib_scores$Glu_score)-(met_inhib_scores$DMSO_score),
                                met_Glu_SE = met_inhib_scores$Glu_SE,
                                met_Mere_delta = (met_inhib_scores$Mere_score)-(met_inhib_scores$DMSO_score),
                                met_Mere_SE = met_inhib_scores$Mere_SE,
                                met_NVP_delta = (met_inhib_scores$NVP_score)-(met_inhib_scores$DMSO_score),
                                met_NVP_SE = met_inhib_scores$NVP_SE,
                                met_Savo_delta = (met_inhib_scores$Savo_score)-(met_inhib_scores$DMSO_score),
                                met_Savo_SE = met_inhib_scores$Savo_SE,
                                met_Tepo_delta =(met_inhib_scores$Tepo_score)-(met_inhib_scores$DMSO_score),
                                met_Tepo_SE = met_inhib_scores$Tepo_SE,
                                met_Tiv_delta = (met_inhib_scores$Tiv_score)-(met_inhib_scores$DMSO_score),
                                met_Tiv_SE = met_inhib_scores$Tiv_SE
                                )

met_DMSO_scores <- data.frame(hgvs=met_inhib_scores$hgvs,
                              pos=met_inhib_scores$pos, 
                              variants=met_inhib_scores$variants,
                              mutation_type=(met_inhib_scores$mutation_type),
                              met_DMSO_score= (met_inhib_scores$DMSO_score),
                              met_DMSO_SE = (met_inhib_scores$DMSO_SE))

ex14_DMSO_scores <- data.frame(hgvs=ex14_inhib_scores$hgvs,
                               pos=ex14_inhib_scores$pos, 
                               mutation_type=(ex14_inhib_scores$mutation_type),
                               variants=ex14_inhib_scores$variants,
                               ex14_DMSO_score= (ex14_inhib_scores$DMSO_score),
                               ex14_DMSO_SE = (ex14_inhib_scores$DMSO_SE))

### same information as before (DMSO control is subtracted from each inhibitor), but this creates a LIST of DATA FRAMES 
### this data structure is much easier to iterate over 
### this data structure allows for analysis of a single structure independently, or in the deltas_df_list, iterate over all data frames
ex14_AMG458_delta <-data.frame(hgvs=ex14_inhib_scores$hgvs,mutation_type=(ex14_inhib_scores$mutation_type),pos=ex14_inhib_scores$pos,variants=ex14_inhib_scores$variants,ex14_AMG458_delta=(ex14_inhib_scores$A458_score)-(ex14_inhib_scores$DMSO_score),ex14_AMG458_scores=(ex14_inhib_scores$A458_score),ex14_AMG458_SE = ex14_inhib_scores$A458_SE,ex14_DMSO_score=(ex14_inhib_scores$DMSO_score),ex14_DMSO_SE = ex14_inhib_scores$DMSO_SE)
ex14_Crizo_delta <-data.frame(hgvs=ex14_inhib_scores$hgvs,mutation_type=(ex14_inhib_scores$mutation_type),pos=ex14_inhib_scores$pos,variants=ex14_inhib_scores$variants,ex14_Crizo_delta=(ex14_inhib_scores$Crizo_score)-(ex14_inhib_scores$DMSO_score),ex14_Crizo_scores=(ex14_inhib_scores$Crizo_score),ex14_Crizo_SE = ex14_inhib_scores$Crizo_SE,ex14_DMSO_score=(ex14_inhib_scores$DMSO_score),ex14_DMSO_SE = ex14_inhib_scores$DMSO_SE)
ex14_Cabo_delta <-data.frame(hgvs=ex14_inhib_scores$hgvs,mutation_type=(ex14_inhib_scores$mutation_type),pos=ex14_inhib_scores$pos,variants=ex14_inhib_scores$variants,ex14_Cabo_delta=(ex14_inhib_scores$A458_score)-(ex14_inhib_scores$DMSO_score),ex14_Cabo_scores=(ex14_inhib_scores$Cabo_score),ex14_Cabo_SE = ex14_inhib_scores$Cabo_SE,ex14_DMSO_score=(ex14_inhib_scores$DMSO_score),ex14_DMSO_SE = ex14_inhib_scores$DMSO_SE)
ex14_Camp_delta <-data.frame(hgvs=ex14_inhib_scores$hgvs,mutation_type=(ex14_inhib_scores$mutation_type),pos=ex14_inhib_scores$pos,variants=ex14_inhib_scores$variants,ex14_Camp_delta=(ex14_inhib_scores$Camp_score)-(ex14_inhib_scores$DMSO_score),ex14_Camp_scores=(ex14_inhib_scores$Camp_score),ex14_Camp_SE = ex14_inhib_scores$Camp_SE,ex14_DMSO_score=(ex14_inhib_scores$DMSO_score),ex14_DMSO_SE = ex14_inhib_scores$DMSO_SE)
ex14_Tepo_delta <-data.frame(hgvs=ex14_inhib_scores$hgvs,mutation_type=(ex14_inhib_scores$mutation_type),pos=ex14_inhib_scores$pos,variants=ex14_inhib_scores$variants,ex14_Tepo_delta=(ex14_inhib_scores$Tepo_score)-(ex14_inhib_scores$DMSO_score),ex14_Tepo_scores=(ex14_inhib_scores$Tepo_score),ex14_Tepo_SE = ex14_inhib_scores$Tepo_SE,ex14_DMSO_score=(ex14_inhib_scores$DMSO_score),ex14_DMSO_SE = ex14_inhib_scores$DMSO_SE)
ex14_Tiv_delta <-data.frame(hgvs=ex14_inhib_scores$hgvs,mutation_type=(ex14_inhib_scores$mutation_type),pos=ex14_inhib_scores$pos,variants=ex14_inhib_scores$variants,ex14_Tiv_delta=(ex14_inhib_scores$Tiv_score)-(ex14_inhib_scores$DMSO_score),ex14_Tiv_scores=(ex14_inhib_scores$Tiv_score),ex14_Tiv_SE = ex14_inhib_scores$Tiv_SE,ex14_DMSO_score=(ex14_inhib_scores$DMSO_score),ex14_DMSO_SE = ex14_inhib_scores$DMSO_SE)
ex14_Gle_delta <-data.frame(hgvs=ex14_inhib_scores$hgvs,mutation_type=(ex14_inhib_scores$mutation_type),pos=ex14_inhib_scores$pos,variants=ex14_inhib_scores$variants,ex14_Gle_delta=(ex14_inhib_scores$Gle_score)-(ex14_inhib_scores$DMSO_score),ex14_Gle_scores=(ex14_inhib_scores$Gle_score),ex14_Gle_SE = ex14_inhib_scores$Gle_SE,ex14_DMSO_score=(ex14_inhib_scores$DMSO_score),ex14_DMSO_SE = ex14_inhib_scores$DMSO_SE)
ex14_Glu_delta <-data.frame(hgvs=ex14_inhib_scores$hgvs,mutation_type=(ex14_inhib_scores$mutation_type),pos=ex14_inhib_scores$pos,variants=ex14_inhib_scores$variants,ex14_Glu_delta=(ex14_inhib_scores$Glu_score)-(ex14_inhib_scores$DMSO_score),ex14_Glu_scores=(ex14_inhib_scores$Glu_score),ex14_Glu_SE = ex14_inhib_scores$Glu_SE,ex14_DMSO_score=(ex14_inhib_scores$DMSO_score),ex14_DMSO_SE = ex14_inhib_scores$DMSO_SE)
ex14_Savo_delta <-data.frame(hgvs=ex14_inhib_scores$hgvs,mutation_type=(ex14_inhib_scores$mutation_type),pos=ex14_inhib_scores$pos,variants=ex14_inhib_scores$variants,ex14_Savo_delta=(ex14_inhib_scores$Savo_score)-(ex14_inhib_scores$DMSO_score),ex14_Savo_scores=(ex14_inhib_scores$Savo_score),ex14_Savo_SE = ex14_inhib_scores$Savo_SE,ex14_DMSO_score=(ex14_inhib_scores$DMSO_score),ex14_DMSO_SE = ex14_inhib_scores$DMSO_SE)
ex14_Mere_delta <-data.frame(hgvs=ex14_inhib_scores$hgvs,mutation_type=(ex14_inhib_scores$mutation_type),pos=ex14_inhib_scores$pos,variants=ex14_inhib_scores$variants,ex14_Mere_delta=(ex14_inhib_scores$Mere_score)-(ex14_inhib_scores$DMSO_score),ex14_Mere_scores=(ex14_inhib_scores$Mere_score),ex14_Mere_SE = ex14_inhib_scores$Mere_SE,ex14_DMSO_score=(ex14_inhib_scores$DMSO_score),ex14_DMSO_SE = ex14_inhib_scores$DMSO_SE)
ex14_NVP_delta <-data.frame(hgvs=ex14_inhib_scores$hgvs,mutation_type=(ex14_inhib_scores$mutation_type), pos=ex14_inhib_scores$pos,variants=ex14_inhib_scores$variants,ex14_NVP_delta=(ex14_inhib_scores$NVP_score)-(ex14_inhib_scores$DMSO_score),ex14_NVP_scores=(ex14_inhib_scores$NVP_score),ex14_NVP_SE = ex14_inhib_scores$NVP_SE,ex14_DMSO_score=ex14_inhib_scores$DMSO_score,ex14_DMSO_SE = ex14_inhib_scores$DMSO_SE)
ex14_deltas_df_list <- list(ex14_AMG458_delta,ex14_Crizo_delta,ex14_Cabo_delta,ex14_Camp_delta,ex14_Tepo_delta,ex14_Tiv_delta,ex14_Gle_delta,ex14_Glu_delta,ex14_Savo_delta,ex14_Mere_delta,ex14_NVP_delta)

met_AMG458_delta <-data.frame(hgvs=met_inhib_scores$hgvs,mutation_type=(met_inhib_scores$mutation_type),pos=met_inhib_scores$pos,variants=met_inhib_scores$variants,met_AMG458_delta=(met_inhib_scores$A458_score)-(met_inhib_scores$DMSO_score),met_AMG458_scores=(met_inhib_scores$A458_score),met_AMG458_SE = met_inhib_scores$A458_SE,met_DMSO_score=(met_inhib_scores$DMSO_score),met_DMSO_SE = met_inhib_scores$DMSO_SE)
met_Crizo_delta <-data.frame(hgvs=met_inhib_scores$hgvs,mutation_type=(met_inhib_scores$mutation_type),pos=met_inhib_scores$pos,variants=met_inhib_scores$variants,met_Crizo_delta=(met_inhib_scores$Crizo_score)-(met_inhib_scores$DMSO_score),met_Crizo_scores=(met_inhib_scores$Crizo_score),met_Crizo_SE = met_inhib_scores$Crizo_SE,met_DMSO_score=(met_inhib_scores$DMSO_score),met_DMSO_SE = met_inhib_scores$DMSO_SE)
met_Cabo_delta <-data.frame(hgvs=met_inhib_scores$hgvs,mutation_type=(met_inhib_scores$mutation_type),pos=met_inhib_scores$pos,variants=met_inhib_scores$variants,met_Cabo_delta=(met_inhib_scores$A458_score)-(met_inhib_scores$DMSO_score),met_Cabo_scores=(met_inhib_scores$Cabo_score),met_Cabo_SE = met_inhib_scores$Cabo_SE,met_DMSO_score=(met_inhib_scores$DMSO_score),met_DMSO_SE = met_inhib_scores$DMSO_SE)
met_Camp_delta <-data.frame(hgvs=met_inhib_scores$hgvs,mutation_type=(met_inhib_scores$mutation_type),pos=met_inhib_scores$pos,variants=met_inhib_scores$variants,met_Camp_delta=(met_inhib_scores$Camp_score)-(met_inhib_scores$DMSO_score),met_Camp_scores=(met_inhib_scores$Camp_score),met_Camp_SE = met_inhib_scores$Camp_SE,met_DMSO_score=(met_inhib_scores$DMSO_score),met_DMSO_SE = met_inhib_scores$DMSO_SE)
met_Tepo_delta <-data.frame(hgvs=met_inhib_scores$hgvs,mutation_type=(met_inhib_scores$mutation_type),pos=met_inhib_scores$pos,variants=met_inhib_scores$variants,met_Tepo_delta=(met_inhib_scores$Tepo_score)-(met_inhib_scores$DMSO_score),met_Tepo_scores=(met_inhib_scores$Tepo_score),met_Tepo_SE = met_inhib_scores$Tepo_SE,met_DMSO_score=(met_inhib_scores$DMSO_score),met_DMSO_SE = met_inhib_scores$DMSO_SE)
met_Tiv_delta <-data.frame(hgvs=met_inhib_scores$hgvs,mutation_type=(met_inhib_scores$mutation_type),pos=met_inhib_scores$pos,variants=met_inhib_scores$variants,met_Tiv_delta=(met_inhib_scores$Tiv_score)-(met_inhib_scores$DMSO_score),met_Tiv_scores=(met_inhib_scores$Tiv_score),met_Tiv_SE = met_inhib_scores$Tiv_SE,met_DMSO_score=(met_inhib_scores$DMSO_score),met_DMSO_SE = met_inhib_scores$DMSO_SE)
met_Gle_delta <-data.frame(hgvs=met_inhib_scores$hgvs,mutation_type=(met_inhib_scores$mutation_type),pos=met_inhib_scores$pos,variants=met_inhib_scores$variants,met_Gle_delta=(met_inhib_scores$Gle_score)-(met_inhib_scores$DMSO_score),met_Gle_scores=(met_inhib_scores$Gle_score),met_Gle_SE = met_inhib_scores$Gle_SE,met_DMSO_score=(met_inhib_scores$DMSO_score),met_DMSO_SE = met_inhib_scores$DMSO_SE)
met_Glu_delta <-data.frame(hgvs=met_inhib_scores$hgvs,mutation_type=(met_inhib_scores$mutation_type),pos=met_inhib_scores$pos,variants=met_inhib_scores$variants,met_Glu_delta=(met_inhib_scores$Glu_score)-(met_inhib_scores$DMSO_score),met_Glu_scores=(met_inhib_scores$Glu_score),met_Glu_SE = met_inhib_scores$Glu_SE,met_DMSO_score=(met_inhib_scores$DMSO_score),met_DMSO_SE = met_inhib_scores$DMSO_SE)
met_Savo_delta <-data.frame(hgvs=met_inhib_scores$hgvs,mutation_type=(met_inhib_scores$mutation_type),pos=met_inhib_scores$pos,variants=met_inhib_scores$variants,met_Savo_delta=(met_inhib_scores$Savo_score)-(met_inhib_scores$DMSO_score),met_Savo_scores=(met_inhib_scores$Savo_score),met_Savo_SE = met_inhib_scores$Savo_SE,met_DMSO_score=(met_inhib_scores$DMSO_score),met_DMSO_SE = met_inhib_scores$DMSO_SE)
met_Mere_delta <-data.frame(hgvs=met_inhib_scores$hgvs,mutation_type=(met_inhib_scores$mutation_type),pos=met_inhib_scores$pos,variants=met_inhib_scores$variants,met_Mere_delta=(met_inhib_scores$Mere_score)-(met_inhib_scores$DMSO_score),met_Mere_scores=(met_inhib_scores$Mere_score),met_Mere_SE = met_inhib_scores$Mere_SE,met_DMSO_score=(met_inhib_scores$DMSO_score),met_DMSO_SE = met_inhib_scores$DMSO_SE)
met_NVP_delta <-data.frame(hgvs=met_inhib_scores$hgvs,mutation_type=(met_inhib_scores$mutation_type), pos=met_inhib_scores$pos,variants=met_inhib_scores$variants,met_NVP_delta=(met_inhib_scores$NVP_score)-(met_inhib_scores$DMSO_score),met_NVP_scores=(met_inhib_scores$NVP_score),met_NVP_SE = met_inhib_scores$NVP_SE,met_DMSO_score=met_inhib_scores$DMSO_score,met_DMSO_SE = met_inhib_scores$DMSO_SE)
met_deltas_df_list <- list(met_AMG458_delta,met_Crizo_delta,met_Cabo_delta,met_Camp_delta,met_Tepo_delta,met_Tiv_delta,met_Gle_delta,met_Glu_delta,met_Savo_delta,met_Mere_delta,met_NVP_delta)


# create data frames with only the scores, NOT subtracted from the DMSO 
# this will the data frames used for plotting of syn, missense, and nonsense distributions 
# this is creates a LIST of DATA FRAMES 
#this data structure is much easier to iterate over 
ex14_DMSO <- data.frame(mutation_type = (ex14_inhib_scores$mutation_type),DMSO = (ex14_inhib_scores$DMSO_score))
ex14_AMG458 <- data.frame(mutation_type = (ex14_inhib_scores$mutation_type),AMG458 = (ex14_inhib_scores$A458_score))
ex14_Crizo <- data.frame(mutation_type = (ex14_inhib_scores$mutation_type),Crizo = (ex14_inhib_scores$Crizo_score))
ex14_Cabo <- data.frame(mutation_type = (ex14_inhib_scores$mutation_type),Cabo = (ex14_inhib_scores$Cabo_score))
ex14_Camp <- data.frame(mutation_type = (ex14_inhib_scores$mutation_type),Camp = (ex14_inhib_scores$Camp_score))
ex14_Gle <- data.frame(mutation_type = (ex14_inhib_scores$mutation_type),Gle= (ex14_inhib_scores$Gle_score))
ex14_Glu <- data.frame(mutation_type = (ex14_inhib_scores$mutation_type),Glu = (ex14_inhib_scores$Glu_score))
ex14_Mere <- data.frame(mutation_type = (ex14_inhib_scores$mutation_type),Mere = (ex14_inhib_scores$Mere_score))
ex14_NVP <- data.frame(mutation_type = (ex14_inhib_scores$mutation_type),NVP = (ex14_inhib_scores$NVP_score))
ex14_Savo <- data.frame(mutation_type = (ex14_inhib_scores$mutation_type),Savo = (ex14_inhib_scores$Savo_score))
ex14_Tepo <- data.frame(mutation_type = (ex14_inhib_scores$mutation_type),Tepo = (ex14_inhib_scores$Tepo_score))
ex14_Tiv <- data.frame(mutation_type = (ex14_inhib_scores$mutation_type),Tiv = (ex14_inhib_scores$Tiv_score))

ex14_inhib_scores_only <- list(ex14_DMSO,ex14_AMG458,ex14_Crizo,ex14_Cabo,ex14_Camp,ex14_Tepo,ex14_Tiv,
                               ex14_Gle,ex14_Glu,ex14_Savo,ex14_Mere,ex14_NVP)


met_DMSO <- data.frame(mutation_type = (met_inhib_scores$mutation_type),DMSO = (met_inhib_scores$DMSO_score))
met_AMG458 <- data.frame(mutation_type = (met_inhib_scores$mutation_type),AMG458 = (met_inhib_scores$A458_score))
met_Crizo <- data.frame(mutation_type = (met_inhib_scores$mutation_type),Crizo = (met_inhib_scores$Crizo_score))
met_Cabo <- data.frame(mutation_type = (met_inhib_scores$mutation_type),Cabo = (met_inhib_scores$Cabo_score))
met_Camp <- data.frame(mutation_type = (met_inhib_scores$mutation_type),Camp = (met_inhib_scores$Camp_score))
met_Gle <- data.frame(mutation_type = (met_inhib_scores$mutation_type),Gle = (met_inhib_scores$Gle_score))
met_Glu <- data.frame(mutation_type = (met_inhib_scores$mutation_type),Glu = (met_inhib_scores$Glu_score))
met_Mere <- data.frame(mutation_type = (met_inhib_scores$mutation_type),Mere = (met_inhib_scores$Mere_score))
met_NVP <- data.frame(mutation_type = (met_inhib_scores$mutation_type),NVP = (met_inhib_scores$NVP_score))
met_Savo <- data.frame(mutation_type = (met_inhib_scores$mutation_type),Savo = (met_inhib_scores$Savo_score))
met_Tepo <- data.frame(mutation_type = (met_inhib_scores$mutation_type),Tepo = (met_inhib_scores$Tepo_score))
met_Tiv <- data.frame(mutation_type = (met_inhib_scores$mutation_type),Tiv = (met_inhib_scores$Tiv_score))
met_inhib_scores_only <- list(met_AMG458,met_Crizo,met_Cabo,met_Camp,met_Tepo,met_Tiv, met_Gle,met_Glu,met_Savo,met_Mere,met_NVP)


# create data frames with only the scores, mutations, hgvs NOT subtracted from the DMSO 
# this will be the data frame for comparing scires within conditions 
# this is creates a LIST of DATA FRAMES 
ex14_AMG458_scores <-data.frame(hgvs=ex14_inhib_scores$hgvs,mutation_type=(ex14_inhib_scores$mutation_type),pos=ex14_inhib_scores$pos,variants=ex14_inhib_scores$variants,ex14_AMG458_scores=(ex14_inhib_scores$A458_score),ex14_AMG458_SE = ex14_inhib_scores$A458_SE)
ex14_Crizo_scores <-data.frame(hgvs=ex14_inhib_scores$hgvs,mutation_type=(ex14_inhib_scores$mutation_type),pos=ex14_inhib_scores$pos,variants=ex14_inhib_scores$variants,ex14_Crizo_scores=(ex14_inhib_scores$Crizo_score),ex14_Crizo_SE = ex14_inhib_scores$Crizo_SE)
ex14_Cabo_scores <-data.frame(hgvs=ex14_inhib_scores$hgvs,mutation_type=(ex14_inhib_scores$mutation_type),pos=ex14_inhib_scores$pos,variants=ex14_inhib_scores$variants,ex14_Cabo_scores=(ex14_inhib_scores$Cabo_score),ex14_Cabo_SE = ex14_inhib_scores$Cabo_SE)
ex14_Camp_scores <-data.frame(hgvs=ex14_inhib_scores$hgvs,mutation_type=(ex14_inhib_scores$mutation_type),pos=ex14_inhib_scores$pos,variants=ex14_inhib_scores$variants,ex14_Camp_scores=(ex14_inhib_scores$Camp_score),ex14_Camp_SE = ex14_inhib_scores$Camp_SE)
ex14_Tepo_scores <-data.frame(hgvs=ex14_inhib_scores$hgvs,mutation_type=(ex14_inhib_scores$mutation_type),pos=ex14_inhib_scores$pos,variants=ex14_inhib_scores$variants,ex14_Tepo_scores=(ex14_inhib_scores$Tepo_score),ex14_Tepo_SE = ex14_inhib_scores$Tepo_SE)
ex14_Tiv_scores <-data.frame(hgvs=ex14_inhib_scores$hgvs,mutation_type=(ex14_inhib_scores$mutation_type),pos=ex14_inhib_scores$pos,variants=ex14_inhib_scores$variants,ex14_Tiv_scores=(ex14_inhib_scores$Tiv_score),ex14_Tiv_SE = ex14_inhib_scores$Tiv_SE)
ex14_Gle_scores <-data.frame(hgvs=ex14_inhib_scores$hgvs,mutation_type=(ex14_inhib_scores$mutation_type),pos=ex14_inhib_scores$pos,variants=ex14_inhib_scores$variants,ex14_Gle_scores=(ex14_inhib_scores$Gle_score),ex14_Gle_SE = ex14_inhib_scores$Gle_SE)
ex14_Glu_scores <-data.frame(hgvs=ex14_inhib_scores$hgvs,mutation_type=(ex14_inhib_scores$mutation_type),pos=ex14_inhib_scores$pos,variants=ex14_inhib_scores$variants,ex14_Glu_scores=(ex14_inhib_scores$Glu_score),ex14_Glu_SE = ex14_inhib_scores$Glu_SE)
ex14_Savo_scores<-data.frame(hgvs=ex14_inhib_scores$hgvs,mutation_type=(ex14_inhib_scores$mutation_type),pos=ex14_inhib_scores$pos,variants=ex14_inhib_scores$variants,ex14_Savo_scores=(ex14_inhib_scores$Savo_score),ex14_Savo_SE = ex14_inhib_scores$Savo_SE)
ex14_Mere_scores <-data.frame(hgvs=ex14_inhib_scores$hgvs,mutation_type=(ex14_inhib_scores$mutation_type),pos=ex14_inhib_scores$pos,variants=ex14_inhib_scores$variants,ex14_Mere_scores=(ex14_inhib_scores$Mere_score),ex14_Mere_SE = ex14_inhib_scores$Mere_SE)
ex14_NVP_scores <-data.frame(hgvs=ex14_inhib_scores$hgvs,mutation_type=(ex14_inhib_scores$mutation_type), pos=ex14_inhib_scores$pos,variants=ex14_inhib_scores$variants,ex14_NVP_scores=(ex14_inhib_scores$NVP_score),ex14_NVP_SE = ex14_inhib_scores$NVP_SE)

ex14_inhib_scores_list <- list(#ex14_DMSO_scores,
                               ex14_AMG458_scores,
                               ex14_Crizo_scores,
                               ex14_Cabo_scores,
                               ex14_Camp_scores,
                               ex14_Tepo_scores,
                               ex14_Tiv_scores,
                               ex14_Gle_scores,
                               ex14_Glu_scores,
                               ex14_Savo_scores,
                               ex14_Mere_scores,
                               ex14_NVP_scores)

ex14_inhib_scores_list  %>% reduce(full_join, by='hgvs')
ex14_inhib_scores_list  %>% reduce(full_join, by='hgvs')

met_AMG458_scores <-data.frame(hgvs=met_inhib_scores$hgvs,mutation_type=(met_inhib_scores$mutation_type),pos=met_inhib_scores$pos,variants=met_inhib_scores$variants,met_AMG458_scores=(met_inhib_scores$A458_score),met_AMG458_SE = met_inhib_scores$A458_SE)
met_Crizo_scores <-data.frame(hgvs=met_inhib_scores$hgvs,mutation_type=(met_inhib_scores$mutation_type),pos=met_inhib_scores$pos,variants=met_inhib_scores$variants,met_Crizo_scores=(met_inhib_scores$Crizo_score),met_Crizo_SE = met_inhib_scores$Crizo_SE)
met_Cabo_scores <-data.frame(hgvs=met_inhib_scores$hgvs,mutation_type=(met_inhib_scores$mutation_type),pos=met_inhib_scores$pos,variants=met_inhib_scores$variants,met_Cabo_scores=(met_inhib_scores$Cabo_score),met_Cabo_SE = met_inhib_scores$Cabo_SE)
met_Camp_scores <-data.frame(hgvs=met_inhib_scores$hgvs,mutation_type=(met_inhib_scores$mutation_type),pos=met_inhib_scores$pos,variants=met_inhib_scores$variants,met_Camp_scores=(met_inhib_scores$Camp_score),met_Camp_SE = met_inhib_scores$Camp_SE)
met_Tepo_scores <-data.frame(hgvs=met_inhib_scores$hgvs,mutation_type=(met_inhib_scores$mutation_type),pos=met_inhib_scores$pos,variants=met_inhib_scores$variants,met_Tepo_scores=(met_inhib_scores$Tepo_score),met_Tepo_SE = met_inhib_scores$Tepo_SE)
met_Tiv_scores <-data.frame(hgvs=met_inhib_scores$hgvs,mutation_type=(met_inhib_scores$mutation_type),pos=met_inhib_scores$pos,variants=met_inhib_scores$variants,met_Tiv_scores=(met_inhib_scores$Tiv_score),met_Tiv_SE = met_inhib_scores$Tiv_SE)
met_Gle_scores <-data.frame(hgvs=met_inhib_scores$hgvs,mutation_type=(met_inhib_scores$mutation_type),pos=met_inhib_scores$pos,variants=met_inhib_scores$variants,met_Gle_scores=(met_inhib_scores$Gle_score),met_Gle_SE = met_inhib_scores$Gle_SE)
met_Glu_scores <-data.frame(hgvs=met_inhib_scores$hgvs,mutation_type=(met_inhib_scores$mutation_type),pos=met_inhib_scores$pos,variants=met_inhib_scores$variants,met_Glu_scores=(met_inhib_scores$Glu_score),met_Glu_SE = met_inhib_scores$Glu_SE)
met_Savo_scores<-data.frame(hgvs=met_inhib_scores$hgvs,mutation_type=(met_inhib_scores$mutation_type),pos=met_inhib_scores$pos,variants=met_inhib_scores$variants,met_Savo_scores=(met_inhib_scores$Savo_score),met_Savo_SE = met_inhib_scores$Savo_SE)
met_Mere_scores <-data.frame(hgvs=met_inhib_scores$hgvs,mutation_type=(met_inhib_scores$mutation_type),pos=met_inhib_scores$pos,variants=met_inhib_scores$variants,met_Mere_scores=(met_inhib_scores$Mere_score),met_Mere_SE = met_inhib_scores$Mere_SE)
met_NVP_scores <-data.frame(hgvs=met_inhib_scores$hgvs,mutation_type=(met_inhib_scores$mutation_type), pos=met_inhib_scores$pos,variants=met_inhib_scores$variants,met_NVP_scores=(met_inhib_scores$NVP_score),met_NVP_SE = met_inhib_scores$NVP_SE)

met_inhib_scores_list <- list(#met_DMSO_scores,
                               met_AMG458_scores,
                               met_Crizo_scores,
                               met_Cabo_scores,
                               met_Camp_scores,
                               met_Tepo_scores,
                               met_Tiv_scores,
                               met_Gle_scores,
                               met_Glu_scores,
                               met_Savo_scores,
                               met_Mere_scores,
                               met_NVP_scores)




#####################################################################################################
###----------------------- ROSACE general score mapping on structures ----------------------------
#####################################################################################################


#--------------------------------- METdEx14 --------------------------------
# add a "pos" columns 
met_rosace_scores$pos = met_rosace_scores$position

# average crizo scores 
crizo_met_rosace_avg <- met_rosace_scores %>% filter(inhibitor == "Crizo") %>% 
  group_by(pos) %>% 
  dplyr::summarise(avg = mean(ROSACE_effects))%>% ungroup()

# average DMSO scores 
DMSO_met_rosace_avg <-  met_rosace_scores %>% filter(inhibitor == "DMSO") %>% 
  group_by(pos) %>% 
  dplyr::summarise(avg = mean(ROSACE_effects))%>% ungroup()

# map scores 
x = map_scores_pdb(met_2WGJ_crizo , crizo_met_rosace_avg, "avg")
write.pdb(x, file="Crizo_avg_met_rosace.pdb")

x = map_scores_pdb(met_2WGJ_crizo , DMSO_met_rosace_avg, "avg")
write.pdb(x, file="DMSO_avg_met_rosace.pdb")


#--------------------------------- WT MET--------------------------------
# add a "pos" columns 
ex14_rosace_scores$pos = ex14_rosace_scores$position

# average crizo scores 
crizo_ex14_rosace_avg <- ex14_rosace_scores %>% filter(inhibitor == "Crizo") %>% 
  group_by(pos) %>% 
  dplyr::summarise(avg = mean(ROSACE_effects))%>% ungroup()

# average DMSO scores 
DMSO_ex14_rosace_avg <-  ex14_rosace_scores %>% filter(inhibitor == "DMSO") %>% 
  group_by(pos) %>% 
  dplyr::summarise(avg = mean(ROSACE_effects))%>% ungroup()

# map scores 
x = map_scores_pdb(ex14_2WGJ_crizo , crizo_ex14_rosace_avg, "avg")
write.pdb(x, file="Crizo_avg_ex14_rosace.pdb")

x = map_scores_pdb(ex14_2WGJ_crizo , DMSO_ex14_rosace_avg, "avg")
write.pdb(x, file="DMSO_avg_ex14_rosace.pdb")



#######################################u#############################################################
###------------------- ROSACE grouped resistance score mapping on structures ------------------------
#####################################################################################################

#--------------------------------- METdEx14 --------------------------------
met_rosace_resis_scores$pos = met_rosace_resis_scores$position
avg_met_rosace_resis = ddply(met_rosace_resis_scores,c("pos"),summarise,mean=mean(I))
met_2WGJ_crizo <- read.pdb("2WGJ")
x = map_scores_pdb(met_2WGJ_crizo, avg_met_rosace_resis, "mean")
write.pdb(x, file="avg_met_rosace_resis_scores_I.pdb")

met_rosace_resis_scores$pos = met_rosace_resis_scores$position
met_2WGJ_crizo <- read.pdb("2WGJ")

avg_met_rosace_I = met_rosace_resis_scores%>% group_by(pos) %>% dplyr::summarise(I_mean = mean(I, na.rm=TRUE))
x = map_scores_pdb(met_2WGJ_crizo , avg_met_rosace_I, "I_mean")
write.pdb(x, file="avg_met_rosace_I.pdb")

avg_met_rosace_II = met_rosace_resis_scores%>% group_by(pos) %>% dplyr::summarise(II_mean = mean(II, na.rm=TRUE))
met_5DG5  <- read.pdb("5DG5")
x = map_scores_pdb(met_5DG5  , avg_met_rosace_II, "II_mean")
write.pdb(x, file="avg_met_rosace_II.pdb")


#--------------------------------- WT MET--------------------------------
ex14_rosace_resis_scores$pos = ex14_rosace_resis_scores$position
avg_ex14_rosace_resis = ddply(ex14_rosace_resis_scores,c("pos"),summarise,mean=mean(I))
ex14_2WGJ_crizo <- read.pdb("2WGJ")
x = map_scores_pdb(ex14_2WGJ_crizo, avg_ex14_rosace_resis, "mean")
write.pdb(x, file="avg_ex14_rosace_resis_scores_I.pdb")

ex14_rosace_resis_scores$pos = ex14_rosace_resis_scores$position
ex14_2WGJ_crizo <- read.pdb("2WGJ")

avg_ex14_rosace_I = ex14_rosace_resis_scores%>% group_by(pos) %>% dplyr::summarise(I_mean = mean(I, na.rm=TRUE))
x = map_scores_pdb(ex14_2WGJ_crizo , avg_ex14_rosace_I, "I_mean")
write.pdb(x, file="avg_ex14_rosace_I.pdb")

avg_ex14_rosace_II = ex14_rosace_resis_scores%>% group_by(pos) %>% dplyr::summarise(II_mean = mean(II, na.rm=TRUE))
ex14_5DG5  <- read.pdb("5DG5")
x = map_scores_pdb(ex14_5DG5  , avg_ex14_rosace_II, "II_mean")
write.pdb(x, file="avg_ex14_rosace_II.pdb")


#######################################u#############################################################
###----------------------------- ROSACE subtraction calculations ------------------------------------
#####################################################################################################

# Mut_drug - WT_drug 
# for ALL conditions
ex14_rosace_scores_DRUG_SUBTRACTION <- ex14_rosace_scores %>% group_by(position, inhibitor) %>% 
  mutate(WT_drug_delta = (ROSACE_effects - ROSACE_effects[match("synonymous", type)])) %>% ungroup()

ex14_rosace_scores_DRUG_SUBTRACTION

# crizo_specific filtering 
ex14_rosace_scores_DRUG_SUBTRACTION$pos = ex14_rosace_scores_DRUG_SUBTRACTION$position
crizo_ex14_rosace_DRUG_diff <- ex14_rosace_scores_DRUG_SUBTRACTION %>% filter(inhibitor == "Crizo")

# map scores onto crizotinb structure 
avg_ex14_rosace_diff_crizo = crizo_ex14_rosace_DRUG_diff%>% group_by(pos) %>% dplyr::summarise(mean = mean(WT_drug_delta, na.rm=TRUE))
x = map_scores_pdb(ex14_2WGJ_crizo,avg_ex14_rosace_diff_crizo, "mean")
write.pdb(x, file="avg_ex14_rosace_diff_crizo.pdb")



# generate pocket heatmap slices 

ex14_rosace_scores_crizo <- ex14_rosace_scores %>% filter(inhibitor == "Crizo")
ex14_rosace_scores_crizo$pos = ex14_rosace_scores_crizo$position


Crizo_Y1230 <- ex14_rosace_scores_crizo %>% filter(pos == "1230")
Crizo_Y1230_ <- (data.frame(variants = Crizo_Y1230$mutation, score = Crizo_Y1230$ROSACE_effects, pos = Crizo_Y1230$pos))

Crizo_G1163 <- ex14_rosace_scores_crizo %>% filter(pos == "1163" )
Crizo_G1163_ <- (data.frame(variants = Crizo_G1163$mutation, score = Crizo_G1163$ROSACE_effects, pos = Crizo_G1163$pos))

Crizo_P1158 <- ex14_rosace_scores_crizo %>% filter(pos == "1158" )
Crizo_P1158_ <- (data.frame(variants = Crizo_P1158$mutation, score = Crizo_P1158$ROSACE_effects, pos = Crizo_P1158$pos))

Crizo_M1160 <-ex14_rosace_scores_crizo %>% filter(pos == "1160" )
Crizo_M1160_ <- (data.frame(variants = Crizo_M1160$mutation, score = Crizo_M1160$ROSACE_effects, pos = Crizo_M1160$pos))

Crizo_M1211 <- ex14_rosace_scores_crizo %>% filter(pos == "1211" )
Crizo_M1211_ <- (data.frame(variants = Crizo_M1211$mutation, score = Crizo_M1211$ROSACE_effects, pos = Crizo_M1211$pos))

Crizo_D1228 <- ex14_rosace_scores_crizo %>% filter(pos == "1228" )
Crizo_D1228_ <- (data.frame(variants = Crizo_D1228$mutation, score = Crizo_D1228$ROSACE_effects, pos = Crizo_D1228$pos))

Crizo_V1092 <- ex14_rosace_scores_crizo %>% filter(pos == "1092" )
Crizo_V1092_ <- (data.frame(variants = Crizo_V1092$mutation, score = Crizo_V1092$ROSACE_effects, pos = Crizo_V1092$pos))

Crizo_A1108 <- ex14_rosace_scores_crizo %>% filter(pos == "1108" )
Crizo_A1108_ <- (data.frame(variants = Crizo_A1108$mutation, score = Crizo_A1108$ROSACE_effects, pos = Crizo_A1108$pos))

Crizo_Y1159 <- ex14_rosace_scores_crizo %>% filter(pos == "1159" )
Crizo_Y1159_ <- (data.frame(variants = Crizo_Y1159$mutation, score = Crizo_Y1159$ROSACE_effects, pos = Crizo_Y1159$pos))



ggplot(Crizo_Y1230_, aes(pos, variants)) +
  geom_tile(aes(fill = score)) +
  geom_text(aes(label = round(score, 1))) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(-12,5))+
  ggtitle("Y1230")

ggplot(Crizo_G1163_, aes(pos, variants)) +
  geom_tile(aes(fill = score)) +
  geom_text(aes(label = round(score, 1))) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(-12,5))+
  ggtitle("G1163")

ggplot(Crizo_P1158_, aes(pos, variants)) +
  geom_tile(aes(fill = score)) +
  geom_text(aes(label = round(score, 1))) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(-12,5))+
  ggtitle("P1158")

ggplot(Crizo_M1160_, aes(pos, variants)) +
  geom_tile(aes(fill = score)) +
  geom_text(aes(label = round(score, 1))) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(-12,5))+
  ggtitle("M1160")

ggplot(Crizo_M1211_, aes(pos, variants)) +
  geom_tile(aes(fill = score)) +
  geom_text(aes(label = round(score, 1))) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(-12,5))+
  ggtitle("M1211")

ggplot(Crizo_V1092_, aes(pos, variants)) +
  geom_tile(aes(fill = score)) +
  geom_text(aes(label = round(score, 1))) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(-12,5))+
  ggtitle("V1092")

ggplot(Crizo_D1228_, aes(pos, variants)) +
  geom_tile(aes(fill = score)) +
  geom_text(aes(label = round(score, 1))) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(-12,5))+
  ggtitle("D1228")


ggplot(Crizo_A1108_, aes(pos, variants)) +
  geom_tile(aes(fill = score)) +
  geom_text(aes(label = round(score, 1))) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(-12,5))+
  ggtitle("A1108")

ggplot(Crizo_Y1159_, aes(pos, variants)) +
  geom_tile(aes(fill = score)) +
  geom_text(aes(label = round(score, 1))) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(-12,5))+
  ggtitle("Y1159")






# generate subtraction heatmap 
order <- c('X','H', 'K', 'R', 'D', 'E',
           'C','M','N','Q','S','T','A','I','L','V','W','F',
           'Y','G','P')

variant_names <- c('Stop','H', 'K', 'R', 'D', 'E',
                   'C','M','N','Q','S','T','A','I','L','V','W','F',
                   'Y','G','P')


met_wt_fasta = "Met_wt.fasta"

met_wt_fasta_con=file(met_wt_fasta, open="r")
met_wt_sequence = readLines(met_wt_fasta_con)[2]
close(met_wt_fasta_con)

met_wt_1 = str_split(substr(met_wt_sequence, 1, 144), '')[[1]]
met_wt_2 = str_split(substr(met_wt_sequence, 145, 287), '')[[1]]

Met_row1 = ggplot(data = ex14_rosace_scores_crizo %>%filter(pos %in% c(1059:1202)),
                  aes(x = pos, y = factor(mutation, level = order, labels = variant_names), fill = ROSACE_effects )) +
  geom_tile(size = 0.2) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(-10,4)) + 
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
    legend.position="none"
  ) +
  labs(y = "Mutation", x = "Position")

Met_row2 = ggplot(data = ex14_rosace_scores_crizo %>% filter(pos %in% c(1203:1345)), 
                  aes(x = pos, y = factor(mutation, level = order, labels = variant_names), fill = ROSACE_effects)) +
  geom_tile(size = 0.2) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(-10,4)) + 
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
    legend.position="none"
  ) +
  labs(y = "Mutation", x = "Position")

Met_DMS = ggarrange(Met_row1,Met_row2,nrow = 2, ncol = 1)

ggsave("Ex14_MET_ROSACE_Crizo_.pdf", height = 7, width = 8.5, Met_DMS)
ggsave("Ex14_MET_ROSACE_Crizo_.png", height = 7, width = 8.5, Met_DMS)


#####################################################################################################
###-------------------- ROSACE DMSO difference  3D mapping and analysis -----------------------------
#####################################################################################################

#--------------------------------- METdEx14 --------------------------------
met_roscace_DMSO_diff_scores$pos= met_roscace_DMSO_diff_scores$position
avg_met_rosace_DMSO_diff <- met_roscace_DMSO_diff_scores %>% group_by(inhibitor) %>%group_by(pos) %>% dplyr::summarise(avg = mean(diff))%>% ungroup()

crizo_met_rosace_DMSO_diff_avg <- met_roscace_DMSO_diff_scores %>% filter(inhibitor == "Crizo") %>% 
  group_by(inhibitor) %>%group_by(pos) %>% 
  dplyr::summarise(avg = mean(diff))%>% ungroup()

met_2WGJ_crizo_DMSO <- read.pdb("2WGJ")
x = map_scores_pdb(met_2WGJ_crizo_DMSO, crizo_met_rosace_DMSO_diff_avg , "avg")
write.pdb(x, file="crizo_met_rosace_DMSO_diff_avg.pdb")


#--------------------------------- WT MET--------------------------------
ex14_roscace_DMSO_diff_scores$pos= ex14_roscace_DMSO_diff_scores$position
avg_ex14_rosace_DMSO_diff <- ex14_roscace_DMSO_diff_scores %>% group_by(inhibitor) %>%group_by(pos) %>% dplyr::summarise(avg = mean(diff))%>% ungroup()

crizo_ex14_rosace_DMSO_diff_avg <- ex14_roscace_DMSO_diff_scores %>% filter(inhibitor == "Crizo") %>% 
  group_by(inhibitor) %>%group_by(pos) %>% 
  dplyr::summarise(avg = mean(diff))%>% ungroup()

ex14_2WGJ_crizo_DMSO <- read.pdb("2WGJ")
x = map_scores_pdb(ex14_2WGJ_crizo_DMSO, crizo_ex14_rosace_DMSO_diff_avg , "avg")
write.pdb(x, file="crizo_ex14_rosace_DMSO_diff_avg.pdb")


#####################################################################################################
###---------------------------- Histograms of all conditions for validation  ------------------------
#####################################################################################################

# function that takes each inhibitor condition for a given data frame and plots:
# synonymous, missense, and nonsense distributions 

# loop over list of data frames for each condition 
# filter mutations based on type and plot distibutions 

# Ex14 screen distributions 
for (i in ex14_inhib_scores_only){
  for (j in names(i[2])){
    syn = (i %>% filter(mutation_type != "X" & mutation_type !="N" & mutation_type != "M"))
    missense = (i %>% filter(mutation_type != "X" & mutation_type !="S" & mutation_type != "N")) 
    nonsense = (i %>% filter(mutation_type != "X" & mutation_type !="S" & mutation_type != "M")) 
    
    plot_syn <- ggplot(syn)+
      geom_histogram(aes(x=syn[[j]]))+
      #geom_density(aes(x=syn[[j]]),bw = "nrd", kernel='gaussian', size=1)+
      xlab("Activity score") + 
      ylab("Count") + 
      xlim(-15,5)+
      theme_classic() +
      ggtitle(paste("Ex14",j,"Syn"))
    plot_missense <- ggplot(missense)+
      geom_histogram(aes(x=missense[[j]]))+
      #geom_density(aes(x=syn[[j]]))+
      xlab("Activity score") + 
      ylab("Count") + 
      xlim(-15,5)+
      theme_classic() +
      ggtitle(paste("Ex14",j,"Missense"))
    plot_nonsense <- ggplot(nonsense)+
      geom_histogram(aes(x=nonsense[[j]]))+
      #geom_density(aes(x=syn[[j]]))+
      xlab("Activity score") + 
      ylab("Count") + 
      xlim(-15,5)+
      theme_classic() +
      ggtitle(paste("Ex14",j,"Stops"))
    plot_density <- ggplot()+
      geom_histogram(aes(x=syn[[j]]),color="gray")+
      geom_histogram(aes(x=missense[[j]]), color="blue")+
      geom_histogram(aes(x=nonsense[[j]]), color="red")+
      #geom_density(aes(x=syn[[j]]),bw = "nrd", kernel='gaussian', size=1,color="gray")+
      #geom_density(aes(x=missense[[j]]),bw = "nrd", kernel='gaussian', size=1, color="blue")+
      #geom_density(aes(x=nonsense[[j]]),bw = "nrd", kernel='gaussian', size=1,color="red")+
      xlab("Activity score") + 
      ylab("Count") + 
      xlim(-15,5)+
      geom_vline(xintercept = 0, linetype="dashed")+
      theme_classic()+
      ggtitle(paste("Ex14",j))
    plot_all <- plot_grid(plot_syn, plot_missense, plot_nonsense, ncol=3, nrow=1)
    print(plot_all)
    print(plot_density)
  }
}


# MET screen distributions 
for (i in met_inhib_scores_only){
  for (j in names(i[2])){
    syn = (i %>% filter(mutation_type != "X" & mutation_type !="N" & mutation_type != "M"))
    missense = (i %>% filter(mutation_type != "X" & mutation_type !="S" & mutation_type != "N")) 
    nonsense = (i %>% filter(mutation_type != "X" & mutation_type !="S" & mutation_type != "M")) 
    
    plot_syn <- ggplot(syn)+
      geom_histogram(aes(x=syn[[j]]))+
      xlab("Activity score") + 
      ylab("Count") + 
      xlim(-15,5)+
      theme_classic() +
      ggtitle(paste("MET",j,"Syn"))
    plot_missense <- ggplot(missense)+
      geom_histogram(aes(x=missense[[j]]))+
      xlab("Activity score") + 
      ylab("Count") + 
      xlim(-15,5)+
      theme_classic() +
      ggtitle(paste("MET",j,"Missense"))
    plot_nonsense <- ggplot(nonsense)+
      geom_histogram(aes(x=nonsense[[j]]))+
      xlab("Activity score") + 
      ylab("Count") + 
      xlim(-15,5)+
      theme_classic() +
      ggtitle(paste("MET",j,"Stops"))
    plot_all <- plot_grid(plot_syn, plot_missense, plot_nonsense, ncol=3, nrow=1)
    print(plot_all)
  }
}


#####################################################################################################
###---------------------------- Distributions of control DMSO condition -----------------------------
#####################################################################################################


### distribution of wt syn mutations in DMSO 

met_DMSO = data.frame(met_DMSO_score = met_inhib_scores$DMSO_score,
                      met_DMSO_SE = met_inhib_scores$DMSO_SE,
                      met_mutation_type = met_inhib_scores$mutation_type,
                      met_pos=met_inhib_scores$pos)

ex14_DMSO = data.frame(ex14_DMSO_score = ex14_inhib_scores$DMSO_score,
                       ex14_DMSO_SE = ex14_inhib_scores$DMSO_SE,
                       ex14_mutation_type = ex14_inhib_scores$mutation_type,
                       ex14_pos=ex14_inhib_scores$pos)

met_wt_DMSO_score= (met_DMSO %>% filter(met_mutation_type != "M" & met_mutation_type != "N" & met_mutation_type != "X"))
ex14_wt_DMSO_score= (ex14_DMSO %>% filter(ex14_mutation_type != "M" & ex14_mutation_type != "N" & ex14_mutation_type != "X"))

# mean all wt scores 
met_wt_DMSO_score_mean = mean(met_wt_DMSO_score$met_DMSO_score, na.rm = T) # -0.3112883
ex14_wt_DMSO_score_mean = mean(ex14_wt_DMSO_score$ex14_DMSO_score, na.rm = T)# -0.425165

# mean of all wt standard error scores
met_wt_DMSO_SE_mean = mean(met_wt_DMSO_score$met_DMSO_SE)
ex14_wt_DMSO_SE_mean = mean(ex14_wt_DMSO_score$ex14_DMSO_SE)

# correlation plot of DMSO +/-Exon14 

#####################################################################################################
###------------------------- DMSO STOPS average score and distributions -------------------------
#####################################################################################################


met_stop_DMSO_score= (met_DMSO %>% filter(met_mutation_type != "M" & met_mutation_type != "S" & met_mutation_type != "X"))
ex14_stop_DMSO_score= (ex14_DMSO %>% filter(ex14_mutation_type != "M" & ex14_mutation_type != "S" & ex14_mutation_type != "X"))

# mean all stop scores 
met_stop_DMSO_score_mean = mean(met_stop_DMSO_score$met_DMSO_score, na.rm = T) # -4.440125
ex14_stop_DMSO_score_mean = mean(ex14_stop_DMSO_score$ex14_DMSO_score, na.rm = T) # -4.166853

# mean of all stop standard error scores
met_stop_DMSO_SE_mean = mean(met_stop_DMSO_score$met_DMSO_SE)
ex14_stop_DMSO_SE_mean = mean(ex14_stop_DMSO_score$ex14_DMSO_SE)



#####################################################################################################
###------------------------- DMSO subtracted scores (done more efficiently) -------------------------
#####################################################################################################

# data frame looking at missense mutation subtracted from the same position in dmso

# 1) Missensee_drug > avg_STOP_DMSO
# 2) Delta = Missense_drug - Missesnse_DMSO 
# 3) Delta > propagation of error


ex14_inhib_scores_long <-
  ex14_inhib_scores %>% pivot_longer(cols = ends_with(c("_score", "_SE")), names_to = c("inhibitor", "SE"), names_sep = "_", values_to = "score") %>% select(!ends_with("_epsilon")) %>% pivot_wider(values_from = "score", names_from = "SE")


ex14_Enrich2_crizo <- ex14_inhib_scores_long  %>% filter(inhibitor == "Crizo")

ex14_DMSO_subtacted <- ex14_inhib_scores_long %>% group_by(hgvs) %>% 
  mutate(DMSO_delta = (score - score[match("DMSO", inhibitor)]), 
         prop_error =  (((SE)^2 + (SE[match("DMSO", inhibitor)])^2)^1/2) ) %>% ungroup()


ex14_DMSO_stop_filtered <- subset(ex14_DMSO_subtacted , ex14_DMSO_subtacted$score > ex14_stop_DMSO_score_mean) 
write.csv(ex14_DMSO_stop_filtered, "ex14_DMSO_stop_filtered.csv", row.names=FALSE)

ex14_DMSO_stop_error_filtered <- subset(ex14_DMSO_subtacted , ex14_DMSO_subtacted$score > ex14_stop_DMSO_score_mean & ex14_DMSO_subtacted$DMSO_delta >= ex14_DMSO_subtacted$prop_error) 
write.csv(ex14_DMSO_stop_error_filtered, "ex14_DMSO_stop_error_filtered.csv", row.names=FALSE)


#Missense_drug > Avg_WT_DRUG?
ex14_missense_drug_avg_WT_drug <- ex14_inhib_scores_long %>% group_by(inhibitor) %>% 
  mutate(score[match("S", inhibitor)])

ex14_WT_Crizo <- subset(ex14_Crizo_scores, mutation_type == "S") 
ex14_Avg_WT_Crizo = mean(ex14_WT_Crizo$ex14_Crizo_scores)
ex14_missenseCrizo_avgWTcrizo <- subset(ex14_Crizo_scores, ex14_Crizo_scores>ex14_Avg_WT_Crizo & mutation_type != "S" )
ex14_missenseCrizo_avgWTDMSO <- subset(ex14_Crizo_scores, ex14_Crizo_scores>ex14_wt_DMSO_score_mean & mutation_type != "S" )

  
  
#####################################################################################################
###---------------------------- Ridge plot of all conditions for validation  ------------------------
#####################################################################################################

ex14_DMSO_2 <- data.frame(mutation_type = (ex14_inhib_scores$mutation_type),score = (ex14_inhib_scores$DMSO_score))
ex14_AMG458_2 <- data.frame(mutation_type = (ex14_inhib_scores$mutation_type),score = (ex14_inhib_scores$A458_score))
ex14_Crizo_2 <- data.frame(mutation_type = (ex14_inhib_scores$mutation_type),score = (ex14_inhib_scores$Crizo_score))
ex14_Cabo_2 <- data.frame(mutation_type = (ex14_inhib_scores$mutation_type),score = (ex14_inhib_scores$Cabo_score))
ex14_Camp_2 <- data.frame(mutation_type = (ex14_inhib_scores$mutation_type),score = (ex14_inhib_scores$Camp_score))
ex14_Gle_2 <- data.frame(mutation_type = (ex14_inhib_scores$mutation_type),score = (ex14_inhib_scores$Gle_score))
ex14_Glu_2 <- data.frame(mutation_type = (ex14_inhib_scores$mutation_type),score = (ex14_inhib_scores$Glu_score))
ex14_Mere_2 <- data.frame(mutation_type = (ex14_inhib_scores$mutation_type),score= (ex14_inhib_scores$Mere_score))
ex14_NVP_2 <- data.frame(mutation_type = (ex14_inhib_scores$mutation_type),score = (ex14_inhib_scores$NVP_score))
ex14_Savo_2 <- data.frame(mutation_type = (ex14_inhib_scores$mutation_type),score = (ex14_inhib_scores$Savo_score))
ex14_Tepo_2 <- data.frame(mutation_type = (ex14_inhib_scores$mutation_type),score = (ex14_inhib_scores$Tepo_score))
ex14_Tiv_2 <- data.frame(mutation_type = (ex14_inhib_scores$mutation_type),score = (ex14_inhib_scores$Tiv_score))


ex14_DMSO_2$inhibitor <- 'DMSO'
ex14_AMG458_2$inhibitor <- 'AMG-458'
ex14_Crizo_2$inhibitor <- 'Crizotinib'
ex14_Cabo_2$inhibitor <- 'Cabozantinib'
ex14_Camp_2$inhibitor <- 'Capmatinib'
ex14_Gle_2$inhibitor <-'Glesatinib'
ex14_Glu_2$inhibitor <- 'Glumetinib'
ex14_Mere_2$inhibitor <- 'Merestinib'
ex14_NVP_2$inhibitor <- 'NVP-BVU972'
ex14_Savo_2$inhibitor <- 'Savolitinib'
ex14_Tepo_2$inhibitor <- 'Tepotinib'
ex14_Tiv_2$inhibitor <- 'Tivantinib'

ridge_ex14 <- rbind(ex14_DMSO_2,
                    ex14_AMG458_2,
                    ex14_Crizo_2,
                    ex14_Cabo_2,
                    ex14_Camp_2,
                    ex14_Gle_2,
                    ex14_Glu_2,
                    ex14_Mere_2,
                    ex14_NVP_2,
                    ex14_Savo_2,
                    ex14_Tepo_2,
                    ex14_Tiv_2)


ridge_ex14 <- ridge_ex14 %>%
  mutate(mutation_type = fct_relevel(mutation_type, levels = 'DMSO', 'Crizotinib','Capmatinib','Tepotinib','Tivantinib','Glumetinib','Savolitinib','NVP-BVU972',
                              'AMG-458','Cabozantinib','Merestinib','Glesatinib'))

Ex14_Functional_score <-  ggplot(ridge_ex14) +
  geom_density_ridges(data =ridge_ex14, aes(x = score, y = inhibitor), alpha=0.6, stat="binline") +
  theme_ridges() + 
  theme()+
  ggtitle("MET+Ex14")
plot(Ex14_Functional_score)

met_DMSO_2 <- data.frame(mutation_type = (met_inhib_scores$mutation_type),score = (met_inhib_scores$DMSO_score))
met_AMG458_2 <- data.frame(mutation_type = (met_inhib_scores$mutation_type),score = (met_inhib_scores$A458_score))
met_Crizo_2 <- data.frame(mutation_type = (met_inhib_scores$mutation_type),score = (met_inhib_scores$Crizo_score))
met_Cabo_2 <- data.frame(mutation _type = (met_inhib_scores$mutation_type),score = (met_inhib_scores$Cabo_score))
met_Camp_2 <- data.frame(mutation_type = (met_inhib_scores$mutation_type),score = (met_inhib_scores$Camp_score))
met_Gle_2 <- data.frame(mutation_type = (met_inhib_scores$mutation_type),score = (met_inhib_scores$Gle_score))
met_Glu_2 <- data.frame(mutation_type = (met_inhib_scores$mutation_type),score = (met_inhib_scores$Glu_score))
met_Mere_2 <- data.frame(mutation_type = (met_inhib_scores$mutation_type),score= (met_inhib_scores$Mere_score))
met_NVP_2 <- data.frame(mutation_type = (met_inhib_scores$mutation_type),score = (met_inhib_scores$NVP_score))
met_Savo_2 <- data.frame(mutation_type = (met_inhib_scores$mutation_type),score = (met_inhib_scores$Savo_score))
met_Tepo_2 <- data.frame(mutation_type = (met_inhib_scores$mutation_type),score = (met_inhib_scores$Tepo_score))
met_Tiv_2 <- data.frame(mutation_type = (met_inhib_scores$mutation_type),score = (met_inhib_scores$Tiv_score))


met_DMSO_2$inhibitor <- 'DMSO'
met_AMG458_2$inhibitor <- 'AMG-458'
met_Crizo_2$inhibitor <- 'Crizotinib'
met_Cabo_2$inhibitor <- 'Cabozantinib'
met_Camp_2$inhibitor <- 'Capmatinib'
met_Gle_2$inhibitor <-'Glesatinib'
met_Glu_2$inhibitor <- 'Glumetinib'
met_Mere_2$inhibitor <- 'Merestinib'
met_NVP_2$inhibitor <- 'NVP-BVU972'
met_Savo_2$inhibitor <- 'Savolitinib'
met_Tepo_2$inhibitor <- 'Tepotinib'
met_Tiv_2$inhibitor <- 'Tivantinib'

ridge_met <- rbind(met_DMSO_2,
                    met_AMG458_2,
                    met_Crizo_2,
                    met_Cabo_2,
                    met_Camp_2,
                    met_Gle_2,
                    met_Glu_2,
                    met_Mere_2,
                    met_NVP_2,
                    met_Savo_2,
                    met_Tepo_2,
                    met_Tiv_2)
ggplot(ridge_met, aes(x = score, y = inhibitor)) +
  geom_density_ridges(alpha=0.6, stat="binline") +
  theme_ridges() + 
  theme()+
  ggtitle("MET-Ex14")


#####################################################################################################
###------------------------------ (Mut_drug - WT_drug) subtraction ----------------------------------
#####################################################################################################

# create data frame for each inhibitor with condition delta and SE condition delta 

# this code finds only the TRUE resistance mutation for a specific inhibitor 
# this is normalized to the synonymous mutations


###------------------------- Crizotinib subtracted resistance mutations --------------------------------
 
ex14_inhib_scores_long <-
  ex14_inhib_scores %>% pivot_longer(cols = ends_with(c("_score", "_SE")), names_to = c("inhibitor", "SE"), names_sep = "_", values_to = "score") %>% select(!ends_with("_epsilon")) %>% pivot_wider(values_from = "score", names_from = "SE")

ex14_inhib_scores_long <- ex14_inhib_scores_long %>% group_by(pos, inhibitor) %>% 
  mutate(delta = (score - score[match("S", mutation_type)]),
         prop_error =  (((SE)^2 + (SE[match("S", mutation_type)])^2)^1/2)) %>% ungroup()

ex14_inhib_mutants <- subset(ex14_inhib_scores_long, ex14_inhib_scores_long$delta >= ex14_inhib_scores_long$prop_error) #positive resistance scores only
ex14_inhib_resis_mutants <- subset(ex14_inhib_scores_long, ex14_inhib_scores_long$delta >= ex14_inhib_scores_long$prop_error 
                                   & ex14_inhib_scores_long$score > 0) #positive resistance scores only


ex14_Crizo_inhib_scores <- ex14_inhib_scores_long  %>% filter(inhibitor == "Crizo" & mutation_type != "N") #crizo syn pos delta scores
ex14_Crizo_inhib_scores_2 <- ex14_inhib_mutants %>% filter(inhibitor == "Crizo" & mutation_type != "N") # syn delta scores filtered by error 
ex14_Crizo_resis <-ex14_inhib_resis_mutants %>% filter(inhibitor == "Crizo" & mutation_type != "N") # error filtered and resis filtered 


## crizotinib crystal structure 
ex14_2WGJ_crizo <- read.pdb("2WGJ")

## this produces a PDB of the avg pos delta score (ie each missense score subtracted by the syn score for every position)
# this is the one I used for mapping
ex14_Crizo_mean_inhib_scores <- ddply(ex14_Crizo_inhib_scores,c("pos"),summarise,mean=mean(delta))
x = map_scores_pdb(ex14_2WGJ_crizo, ex14_Crizo_mean_inhib_scores, "mean")
write.pdb(x, file="ex14_2WGJ_crizo_avg.pdb")


## this produces a PDB of the error propagation filtered score 
## same dataset as above but with error taken into account 
ex14_Crizo_mean_inhib_scores_2 <- ddply(ex14_Crizo_inhib_scores_2,c("pos"),summarise,mean=mean(delta))
x = map_scores_pdb(ex14_2WGJ_crizo, ex14_Crizo_mean_inhib_scores_2, "mean")
write.pdb(x, file="ex14_2WGJ_crizo_avg_error_filtered.pdb")

## this only maps resistance mutations 
ex14_Crizo_mean_inhib_scores_3 <- ddply(ex14_Crizo_resis,c("pos"),summarise,mean=mean(delta))
x = map_scores_pdb(ex14_2WGJ_crizo, ex14_Crizo_mean_inhib_scores_3, "mean")
write.pdb(x, file="ex14_2WGJ_crizo_avg_error_filtered_resistant.pdb")

ggplot()+
  geom_histogram(aes(x=ex14_Crizo_mean_inhib_scores$mean))

Crizo_Y1230 <- ex14_Crizo_inhib_scores %>% filter(pos == "1230")
Crizo_Y1230_ <- (data.frame(variants = Crizo_Y1230$variants, score = Crizo_Y1230$delta, pos = Crizo_Y1230$pos))

Crizo_G1163 <- ex14_Crizo_inhib_scores %>% filter(pos == "1163" )
Crizo_G1163_ <- (data.frame(variants = Crizo_G1163$variants, score = Crizo_G1163$delta, pos = Crizo_G1163$pos))

Crizo_P1158 <- ex14_Crizo_inhib_scores %>% filter(pos == "1158" )
Crizo_P1158_ <- (data.frame(variants = Crizo_P1158$variants, score = Crizo_P1158$delta, pos = Crizo_P1158$pos))

Crizo_M1160 <- ex14_Crizo_inhib_scores %>% filter(pos == "1160" )
Crizo_M1160_ <- (data.frame(variants = Crizo_M1160$variants, score = Crizo_M1160$delta, pos = Crizo_M1160$pos))

Crizo_M1211 <- ex14_Crizo_inhib_scores %>% filter(pos == "1211" )
Crizo_M1211_ <- (data.frame(variants = Crizo_M1211$variants, score = Crizo_M1211$delta, pos = Crizo_M1211$pos))

Crizo_D1228 <- ex14_Crizo_inhib_scores %>% filter(pos == "1228" )
Crizo_D1228_ <- (data.frame(variants = Crizo_D1228$variants, score = Crizo_D1228$delta, pos = Crizo_D1228$pos))

Crizo_V1092 <- ex14_Crizo_inhib_scores %>% filter(pos == "1092" )
Crizo_V1092_ <- (data.frame(variants = Crizo_V1092$variants, score = Crizo_V1092$delta, pos = Crizo_V1092$pos))

Crizo_A1108 <- ex14_Crizo_inhib_scores %>% filter(pos == "1108" )
Crizo_A1108_ <- (data.frame(variants = Crizo_A1108$variants, score = Crizo_A1108$delta, pos = Crizo_A1108$pos))

Crizo_Y1159 <- ex14_Crizo_inhib_scores %>% filter(pos == "1159" )
Crizo_Y1159_ <- (data.frame(variants = Crizo_Y1159$variants, score = Crizo_Y1159$delta, pos = Crizo_Y1159$pos))




ggplot(Crizo_Y1230_, aes(pos, variants)) +
  geom_tile(aes(fill = score)) +
  geom_text(aes(label = round(score, 1))) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(-12,5))+
  ggtitle("Y1230")

ggplot(Crizo_G1163_, aes(pos, variants)) +
  geom_tile(aes(fill = score)) +
  geom_text(aes(label = round(score, 1))) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(-12,5))+
  ggtitle("G1163")

ggplot(Crizo_P1158_, aes(pos, variants)) +
  geom_tile(aes(fill = score)) +
  geom_text(aes(label = round(score, 1))) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(-12,5))+
  ggtitle("P1158")

ggplot(Crizo_M1160_, aes(pos, variants)) +
  geom_tile(aes(fill = score)) +
  geom_text(aes(label = round(score, 1))) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(-12,5))+
  ggtitle("M1160")

ggplot(Crizo_M1211_, aes(pos, variants)) +
  geom_tile(aes(fill = score)) +
  geom_text(aes(label = round(score, 1))) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(-12,5))+
  ggtitle("M1211")
  
ggplot(Crizo_V1092_, aes(pos, variants)) +
  geom_tile(aes(fill = score)) +
  geom_text(aes(label = round(score, 1))) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(-12,5))+
  ggtitle("V1092")

ggplot(Crizo_D1228_, aes(pos, variants)) +
  geom_tile(aes(fill = score)) +
  geom_text(aes(label = round(score, 1))) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(-12,5))+
  ggtitle("D1228")


ggplot(Crizo_A1108_, aes(pos, variants)) +
  geom_tile(aes(fill = score)) +
  geom_text(aes(label = round(score, 1))) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(-12,5))+
  ggtitle("A1108")

ggplot(Crizo_Y1159_, aes(pos, variants)) +
  geom_tile(aes(fill = score)) +
  geom_text(aes(label = round(score, 1))) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(-12,5))+
  ggtitle("Y1159")







###------------------------- DMSO subtracted resistance mutations --------------------------------

## same as above but with DMSO correction (ie N1059A_crizo - N1059A_DMSO)

ex14_Crizo_DMSO_delta_avg <- ddply(ex14_Crizo_delta,c("pos"),summarise,mean=mean(ex14_Crizo_delta))

Crizo_Y1230_DMSO <- data.frame(variants = ex14_Crizo_delta$variants, score = ex14_Crizo_delta$ex14_Crizo_delta, pos = ex14_Crizo_delta$pos)
Crizo_Y1230_DMSO <- Crizo_Y1230_DMSO %>% filter(pos == "1230")

Crizo_G1163_DMSO <- data.frame(variants = ex14_Crizo_delta$variants, score = ex14_Crizo_delta$ex14_Crizo_delta, pos = ex14_Crizo_delta$pos)
Crizo_G1163_DMSO <- Crizo_G1163_DMSO %>% filter(pos == "1163")

Crizo_P1158_DMSO <- data.frame(variants = ex14_Crizo_delta$variants, score = ex14_Crizo_delta$ex14_Crizo_delta, pos = ex14_Crizo_delta$pos)
Crizo_P1158_DMSO <- Crizo_P1158_DMSO %>% filter(pos == "1158")

Crizo_M1160_DMSO <- data.frame(variants = ex14_Crizo_delta$variants, score = ex14_Crizo_delta$ex14_Crizo_delta, pos = ex14_Crizo_delta$pos)
Crizo_M1160_DMSO <- Crizo_M1160_DMSO %>% filter(pos == "1160")

Crizo_M1211_DMSO <- data.frame(variants = ex14_Crizo_delta$variants, score = ex14_Crizo_delta$ex14_Crizo_delta, pos = ex14_Crizo_delta$pos)
Crizo_M1211_DMSO <- Crizo_M1211_DMSO %>% filter(pos == "1211")

ggplot(Crizo_Y1230_DMSO , aes(pos, variants)) +
  geom_tile(aes(fill = score)) +
  geom_text(aes(label = round(score, 1))) +
  scale_fill_continuous_divergingx(palette = 'PuOr', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(-12,4))+
  ggtitle("Y1230")

ggplot(Crizo_G1163_DMSO , aes(pos, variants)) +
  geom_tile(aes(fill = score)) +
  geom_text(aes(label = round(score, 1))) +
  scale_fill_continuous_divergingx(palette = 'PuOr', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(-12,4))+
  ggtitle("G1163")

ggplot(Crizo_P1158_DMSO , aes(pos, variants)) +
  geom_tile(aes(fill = score)) +
  geom_text(aes(label = round(score, 1))) +
  scale_fill_continuous_divergingx(palette = 'PuOr', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(-12,4))+
  ggtitle("P1158")

ggplot(Crizo_M1160_DMSO , aes(pos, variants)) +
  geom_tile(aes(fill = score)) +
  geom_text(aes(label = round(score, 1))) +
  scale_fill_continuous_divergingx(palette = 'PuOr', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(-12,4))+
  ggtitle("M1160")

ggplot(Crizo_M1211_DMSO , aes(pos, variants)) +
  geom_tile(aes(fill = score)) +
  geom_text(aes(label = round(score, 1))) +
  scale_fill_continuous_divergingx(palette = 'PuOr', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(-12,4))+
  ggtitle("M1211")




x = map_scores_pdb(ex14_2WGJ_crizo, ex14_Crizo_DMSO_delta_avg, "mean")
write.pdb(x, file="ex14_2WGJ_crizo_DMSO_sub_avg.pdb")


###------------------------- Drug and dmso subtracted scores--------------------------------

####  Missense_drug - (WT_drug + Missense_dmso)



ex14_inhib_scores_long <-
  ex14_inhib_scores %>% pivot_longer(cols = ends_with(c("_score", "_SE")), names_to = c("inhibitor", "SE"), names_sep = "_", values_to = "score") %>% select(!ends_with("_epsilon")) %>% pivot_wider(values_from = "score", names_from = "SE")

ex14_inhib_scores_long_2 <- ex14_inhib_scores_long %>% group_by(pos, inhibitor) %>% 
  mutate(delta_combined = (score - ((score[match("S", mutation_type)])))) %>% ungroup()

ex14_inhib_mutants <- subset(ex14_inhib_scores_long, ex14_inhib_scores_long$delta >= ex14_inhib_scores_long$prop_error) #positive resistance scores only
ex14_inhib_resis_mutants <- subset(ex14_inhib_scores_long, ex14_inhib_scores_long$delta >= ex14_inhib_scores_long$prop_error 
                                   & ex14_inhib_scores_long$score > 0) #positive resistance scores only


ex14_Crizo_inhib_scores <- ex14_inhib_scores_long  %>% filter(inhibitor == "Crizo" & mutation_type != "N") #crizo syn pos delta scores
ex14_Crizo_inhib_scores_2 <- ex14_inhib_mutants %>% filter(inhibitor == "Crizo" & mutation_type != "N") # syn delta scores filtered by error 
ex14_Crizo_resis <-ex14_inhib_resis_mutants %>% filter(inhibitor == "Crizo" & mutation_type != "N") # error filtered and resis filtered 


#####################################################################################################
###---------------------- (Mut_drug/Mut_DMSO)/ [(WT_drug)/WT_DMSO)] -------------------------
#####################################################################################################

##### normalization method from Bolon Lab UMass

# drug score ratios
ex14_inhib_normalization <- ex14_inhib_scores_long %>% group_by(hgvs) %>% 
  mutate(Mut_drug_Mut_DMSO_ratio = (score/((score[match("DMSO", inhibitor)])))) %>% ungroup() 

ex14_inhib_normalization <- ex14_inhib_normalization %>% group_by(pos, inhibitor)%>%
  mutate(WT_drug_WT_DMSO_ratio = (Mut_drug_Mut_DMSO_ratio[match("S", mutation_type)]), 
         DrugScore = Mut_drug_Mut_DMSO_ratio/WT_drug_WT_DMSO_ratio)%>% ungroup() 

# Ridge plot of all conditions for all conditions to show normalization 

ex14_inhib_normalization_noDMSO <- subset(ex14_inhib_normalization, ex14_inhib_normalization$inhibitor != "DMSO")

Ex14_Drug_Score <- ggplot(ex14_inhib_normalization_noDMSO, aes(x = DrugScore, y = inhibitor)) +
  geom_density_ridges(alpha=0.6, stat="binline") +
  theme_ridges() + 
  xlim(-10,10)+
  theme()+
  ggtitle("MET+Ex14 Drug Score")

plot_grid(Ex14_Functional_score,Ex14_Drug_Score)

#####################################################################################################
###-----------------Correlation Plots and statistical filtering  for +/-Ex14 conditions  ------------
###------------------------------------- inhibitor chemistry ----------------------------------------
#####################################################################################################

### data frame with propagation of error from +/- Exon14 conditions 

# merge the data frames that have been DMSO subtracted
# create a list of merged data frames to iterate over and add a score difference column and propagation of error 
#merged_DMSO<-merge(ex14_DMSO_scores,met_DMSO_scores)
merged_AMG458_delta<-merge(ex14_AMG458_delta,met_AMG458_delta)
merged_Crizo_delta<-merge(ex14_Crizo_delta,met_Crizo_delta)
merged_Cabo_delta<-merge(ex14_Cabo_delta,met_Cabo_delta)
merged_Camp_delta<-merge(ex14_Camp_delta,met_Camp_delta)
merged_Tepo_delta<-merge(ex14_Tepo_delta,met_Tepo_delta)
merged_Tiv_delta<-merge(ex14_Tiv_delta,met_Tiv_delta)
merged_Glu_delta<-merge(ex14_Glu_delta,met_Glu_delta)
merged_Gle_delta<-merge(ex14_Gle_delta,met_Gle_delta)
merged_Mere_delta<-merge(ex14_Mere_delta,met_Mere_delta)
merged_NVP_delta<-merge(ex14_NVP_delta,met_NVP_delta)
merged_Savo_delta<-merge(ex14_Savo_delta,met_Savo_delta)

merged_delta_df_list <- list(merged_AMG458_delta, merged_Crizo_delta, merged_Cabo_delta, merged_Camp_delta,
                             merged_Tepo_delta, merged_Tiv_delta, merged_Gle_delta, merged_Glu_delta,
                             merged_Savo_delta, merged_Mere_delta, merged_NVP_delta)

merged_inhib_name_list1 <- list("AMG458","Crizo","Cabo","Camp","Tepo","Tiv","Gle","Glu","Savo","Mere","NVP")

merged_inhib_name_list2 <- list("AMG458","Crizotinib","Cabozantinib","Capmantinib",
                               "Tepotinib","Tivantinib","Glesatinib","Glumetinib",
                               "Savolinitinb","Merestinib","NVP-BVU972")

merged_inhib_name_list3 <- list("AMG458_filtered","Crizo_filtered","Cabo_filtered","Camp_filtered",
                                "Tepo_filtered","Tiv_filtered","Gle_filtered","Glu_filtered","Savo_filtered",
                                "Mere_filtered","NVP_filtered")

merged_inhib_name_list4 <- list("met_AMG458_filtered","met_Crizo_filtered","met_Cabo_filtered","met_Camp_filtered",
                                "met_Tepo_filtered","met_Tiv_filtered","met_Gle_filtered","met_Glu_filtered","met_Savo_filtered",
                                "met_Mere_filtered","met_NVP_filtered")

merged_inhib_name_list5 <- list("ex14_AMG45_filtered","ex14_Crizo_filtered","ex14_Cabo_filtered","ex14_Camp_filtered",
                                "ex14_Tepo_filtered","ex14_Tiv_filtered","ex14_Gle_filtered","ex14_Glu_filtered","ex14_Savo_filtered",
                                "ex14_Mere_filtered","ex14_NVP_filtered")

merged_inhib_name_list6 <- list("ex14_AMG458_GOF_LOF","ex14_Crizo_GOF_LOF","ex14_Cabo_GOF_LOF","ex14_Camp_GOF_LOF",
                                "ex14_Tepo_GOF_LOF","ex14_Tiv_GOF_LOF","ex14_Gle_GOF_LOF","ex14_Glu_GOF_LOF","ex14_Savo_GOF_LOF",
                                "ex14_Mere_GOF_LOF","ex14_NVP_GOF_LOF")

merged_inhib_name_list7 <- list("ex14_AMG458_GOF","ex14_Crizo_GOF","ex14_Cabo_GOF","ex14_Camp_GOF",
                                "ex14_Tepo_GOF","ex14_Tiv_GOF","ex14_Gle_GOF","ex14_Glu_GOF","ex14_Savo_GOF",
                                "ex14_Mere_GOF","ex14_NVP_GOF")

merged_inhib_name_list8 <- list("met_AMG458_GOF","met_Crizo_GOF","met_Cabo_GOF","met_Camp_GOF",
                                "met_Tepo_GOF","met_Tiv_GOF","met_Gle_GOF","met_Glu_GOF","met_Savo_GOF",
                                "met_Mere_GOF","met_NVP_GOF")




# color palette for figs below 
cpb2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# correlation plots of the DMSO subtracted scores for +/-Ex14 backgrounds in comparison to the correlation of DMSO for both backgrounds 
q = 1
for (i in merged_delta_df_list){
  plot_inhib_corr <- ggplot(data=i, aes(x=i[[10]], y=i[[5]])) + 
    geom_point(shape=21,aes(y=i[[5]], colour = mutation_type))+
    scale_color_manual(values=cpb2)+
    theme_linedraw()+
    xlab("MET-Ex14") + ylab("MET+Ex14")+
    ggtitle(paste(merged_inhib_name_list2[q]))
  q = q+1
  print(plot_inhib_corr)
}

# iterate over list of merged +/-Ex1 dfs and add a score difference column and propagation of error 
a = 1
merged_prop_error_df_list <- list() # this is a list of dfs for each condition, now with score diff and prop error values 
for (i in merged_delta_df_list){
  # met vs ex14 
  i[ , paste0(merged_inhib_name_list1[a],"_met_ex14_score_diff")] <- abs(i[11]- i[6])
  i[ , paste0(merged_inhib_name_list1[a],"_met_ex14_prop_error")] <- ((i[12])^2 +(i[7])^2)^(1/2)
  # met vs DMSO
  i[ , paste0(merged_inhib_name_list1[a],"_met_DMSO_score_diff")] <- abs(i[11]- i[13])
  i[ , paste0(merged_inhib_name_list1[a],"_met_DMSO_prop_error")] <- ((i[12])^2 +(i[14])^2)^(1/2)
  # ex14 vs DMSO 
  i[ , paste0(merged_inhib_name_list1[a],"_ex14_DMSO_score_diff")] <- abs(i[6]- i[8])
  i[ , paste0(merged_inhib_name_list1[a],"_ex14_DMSO_prop_error")] <- ((i[7])^2 +(i[9])^2)^(1/2)

  merged_prop_error_df_list[[length(merged_prop_error_df_list) + 1]] <- i
  print(i)
  a = a+1
}

# creates a new INDEPENDENT data frame for each inhibitor condition that has a high score diff based on prop error 
# each table has it's own name (ie "ex14_Crizo_filtered")

# prop error for met vs. DMSO filtering 
w = 1 
for (i in merged_prop_error_df_list){
  j <- subset(i, i[17] >= i[18] & mutation_type != "S"& mutation_type != "N") # met_DMSO_score_diff >= met_DMSO_prop_error
  k <- subset(j%>% filter(j[11] >= 2*sd(met_wt_DMSO_score$met_DMSO_score, na.rm = T))) # met_delta >= 2 SD of DMSO syn mean 
  l <- k[-c(5:9,15,16,19,20)] # remove other inhibitor columns 
  assign(paste(merged_inhib_name_list4[w], sep = ""),l)
  print(l)
  w=w+1
}

# prop error ex14 vs. DMSO filtering 
z = 1 
for (i in merged_prop_error_df_list){
  j <- subset(i, i[19] >= i[20] & mutation_type != "S"& mutation_type != "N") # ex14_DMSO_score_diff >= ex14_DMSO_prop_error
  k <- subset(j%>% filter(j[6] >= 2*sd(ex14_wt_DMSO_score$ex14_DMSO_score, na.rm = T))) # ex14_delta >= 2 SD of DMSO syn mean 
  l <- k[-c(10:18)] # remove other inhibitor columns 
  assign(paste(merged_inhib_name_list5[z], sep = ""),l) # example : ex14_Crizo_filtered
  print(j)
  z=z+1
}

# prop error for met vs. ex14 filtering 
#e = 1 
#for (i in merged_prop_error_df_list){
#  j <- subset(i, i[15] >= i[16] & mutation_type != "S"& mutation_type != "N") # met_ex_score_diff >= met_ex_prop_error
#  k <- subset(j%>% filter(j[10] >= 2*sd(met_wt_DMSO_score$met_DMSO_score, na.rm = T))) # met scores filtered 
#  l <- subset(j%>% filter(j[5] >= 2*sd(met_wt_DMSO_score$ex14_DMSO_score, na.rm = T))) # ex14 scores filtered 
  #assign(paste(merged_inhib_name_list3[e], sep = ""),j)
#  print(j)
#  e=e+1
#}

#met_test <- subset(Crizo_filtered %>% filter(met_Crizo_delta- met_wt_DMSO_score_mean >= 2*sd(met_wt_DMSO_score$met_DMSO_score, na.rm = T)))
#print(met_test)

#####################################################################################################
###------------------------------- MET+Ex14 GOF/resistance mutations --------------------------------------
#####################################################################################################

# ex14_deltas_df_list contains all the delta scores which can be iterated over to get the prop error
# iterate over list of merged +/-Ex1 dfs and add a score difference column and propagation of error 
b = 1
ex14_prop_error_df_list <- list() # this is a list of dfs for each condition, now with score diff and prop error values 
for (i in ex14_deltas_df_list){
  i[ , paste0(merged_inhib_name_list1[b],"_ex14_DMSO_score_diff")] <- abs(i[6]- i[8]) # score - DMSO score 
  i[ , paste0(merged_inhib_name_list1[b],"_ex14_DMSO_prop_error")] <- ((i[7])^2 +(i[9])^2)^(1/2) # [(inhibitorSE^2) + (dmsoSE^2)]^(1/2) 
  ex14_prop_error_df_list[[length(ex14_prop_error_df_list) + 1]] <- i
  print(i)
  b = b+1
}

# creates a new INDEPENDENT data frame for each inhibitor condition that has a high score diff based on prop error 
# each table has it's own name (ie "ex14_Crizo_GOF)
# prop error ex14 vs. DMSO filtering 
p = 1 
for (i in ex14_prop_error_df_list){
  j <- subset(i, i[10] >= i[11] & mutation_type != "S"& mutation_type != "N") # ex14_DMSO_score_diff >= ex14_DMSO_prop_error
  k <- subset(j%>% filter(j[5] > 0 & (j[6]) >= 2*sd(ex14_wt_DMSO_score$ex14_DMSO_score, na.rm = T))) # ex14_delta >= 2 SD of DMSO syn mean, this only gives you GOF
  #k <- subset(j%>% filter(abs(j[6]) >= 2*sd(ex14_wt_DMSO_score$ex14_DMSO_score, na.rm = T))) # this will give you GOF+LOF
  l <- k[-c(12:22)] # remove other inhibitor columns 
  assign(paste(merged_inhib_name_list7[p], sep = ""),l) # example : ex14_Crizo_GOF_LOF
  print(j)
  p=p+1
}


#####################################################################################################
###------------------------------- METdEx14 GOF/resistance mutations --------------------------------
#####################################################################################################

# met_deltas_df_list contains all the delta scores which can be iterated over to get the prop error
# iterate over list of merged +/-Ex1 dfs and add a score difference column and propagation of error 
d = 1
met_prop_error_df_list <- list() # this is a list of dfs for each condition, now with score diff and prop error values 
for (i in met_deltas_df_list){
  i[ , paste0(merged_inhib_name_list1[d],"_met_DMSO_score_diff")] <- abs(i[6]- i[8]) # score - DMSO score 
  i[ , paste0(merged_inhib_name_list1[d],"_met_DMSO_prop_error")] <- ((i[7])^2 +(i[9])^2)^(1/2) # [(inhibitorSE^2) + (dmsoSE^2)]^(1/2) 
  met_prop_error_df_list[[length(met_prop_error_df_list) + 1]] <- i
  print(i)
  d = d+1
}

# creates a new INDEPENDENT data frame for each inhibitor condition that has a high score diff based on prop error 
# each table has it's own name (ie "met_Crizo_GOF)
# prop error met vs. DMSO filtering 
e = 1 
for (i in met_prop_error_df_list){
  j <- subset(i, i[10] >= i[11] & mutation_type != "S"& mutation_type != "N") # met_DMSO_score_diff >= met_DMSO_prop_error
  k <- subset(j%>% filter(j[5] > 0 & (j[6]) >= 2*sd(met_wt_DMSO_score$met_DMSO_score, na.rm = T))) # met_delta >= 2 SD of DMSO syn mean, this only gives you GOF
  #k <- subset(j%>% filter(abs(j[6]) >= 2*sd(met_wt_DMSO_score$met_DMSO_score, na.rm = T))) # this will give you GOF+LOF
  l <- k[-c(12:22)] # remove other inhibitor columns 
  assign(paste(merged_inhib_name_list8[e], sep = ""),l) # example : met_Crizo_GOF_LOF
  print(j)
  e=e+1
}


#####################################################################################################
###------------------ heatmap of GAIN of Function Met+Ex14 mutations --------------------------------
#####################################################################################################

print_heatmap <- function(df, label, low, high, output_file) {
  
  label <- enquo(label)

  order <- c('X','H', 'K', 'R', 'D', 'E','C','M','N','Q','S','T','A','I','L','V','W','F','Y','G','P')
  variant_names <- c('Stop','H', 'K', 'R', 'D', 'E','C','M','N','Q','S','T','A','I','L','V','W','F','Y','G','P')
  met_wt_1 = str_split(substr(met_wt_sequence, 1, 287), '')[[1]]
  
  row1 = ggplot(data = df %>% filter(pos %in% c(1059:1345)), 
                aes(x = pos, y = factor(variants, level = order, labels = variant_names), fill = !!label )) +
    geom_tile(size = 0.2) +
    scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(low,high)) + 
    scale_color_manual(values = c(NA,'green')) +
    scale_x_continuous(breaks = seq(1059,1345, by = 5),
                       expand = c(0,0),
                       sec.axis = sec_axis(
                         trans = ~.,
                         name = "Sequence",
                         breaks = seq(1059,1345),
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
      legend.position="right") +
    labs(y = "Mutation", x = "Position")
  
  heatmap = ggarrange(row1, nrow = 1, ncol = 1)
  plot(heatmap)
  ggsave(output_file, height = 7, width = 8.5, heatmap)
  
}

#met+ex14
print_heatmap(ex14_AMG458_GOF, ex14_AMG458_delta, 0, 6, "ex14_AMG458_Resistance_heatmap.png")
print_heatmap(ex14_Crizo_GOF, ex14_Crizo_delta, 0, 4, "ex14_Crizo_Resistance_heatmap.png")
print_heatmap(ex14_Camp_GOF, ex14_Camp_delta, 1, 6, "ex14_Camp_Resistance_heatmap.png")
print_heatmap(ex14_Tepo_GOF, ex14_Tepo_delta, 0, 4, "ex14_Tepo_Resistance_heatmap.png")
print_heatmap(ex14_NVP_GOF, ex14_NVP_delta, 0, 10, "ex14_NVP_Resistance_heatmap.png")
print_heatmap(ex14_Glu_GOF, ex14_Glu_delta, 0, 4, "ex14_Glu_Resistance_heatmap.png")
#print_heatmap(ex14_Savo_GOF, ex14_Savo_delta, 0, 4, "ex14_Savo_Resistance_heatmap.png")
print_heatmap(ex14_Cabo_GOF, ex14_Cabo_delta, 0, 6, "ex14_Cabo_Resistance_heatmap.png")
print_heatmap(ex14_Gle_GOF, ex14_Gle_delta, 0, 8, "ex14_Gle_Resistance_heatmap.png")
print_heatmap(ex14_Mere_GOF, ex14_Mere_delta, 1, 11, "ex14_Mere_Resistance_heatmap.png")

#met-ex14
print_heatmap(met_AMG458_GOF, met_AMG458_delta, 0, 5, "met_AMG458_Resistance_heatmap.png")
print_heatmap(met_Crizo_GOF, met_Crizo_delta, 0, 8, "met_Crizo_Resistance_heatmap.png")
print_heatmap(met_Camp_GOF, met_Camp_delta, 0, 7, "met_Camp_Resistance_heatmap.png")
print_heatmap(met_Tepo_GOF, met_Tepo_delta, 0, 6, "met_Tepo_Resistance_heatmap.png")
print_heatmap(met_NVP_GOF, met_NVP_delta, 0, 6, "met_NVP_Resistance_heatmap.png")
print_heatmap(met_Glu_GOF, met_Glu_delta, 1, 4, "met_Glu_Resistance_heatmap.png")
print_heatmap(met_Savo_GOF, met_Savo_delta, 0, 5, "met_Savo_Resistance_heatmap.png")
print_heatmap(met_Cabo_GOF, met_Cabo_delta, 0, 4, "met_Cabo_Resistance_heatmap.png")
print_heatmap(met_Gle_GOF, met_Gle_delta, 0, 9, "met_Gle_Resistance_heatmap.png")
print_heatmap(met_Mere_GOF, met_Mere_delta, 0, 4, "met_Mere_Resistance_heatmap.png")


#####################################################################################################
###-------------------------------MET Resistance Balloon Plots  -------------------------------------
###------------------------------- from DMSO subtracted data  ---------------------------------------
#####################################################################################################


# Type 1 (crizo, camp, tepo, glu, nvp, savo)
# Type 2 (cabo, mere, gle)
# Type 1.5 (AMG)

ex14_Crizo_resis <- data.frame(pos =ex14_Crizo_GOF$pos,Crizo = ex14_Crizo_GOF$ex14_Crizo_delta)
ex14_Cabo_resis <- data.frame(pos =ex14_Cabo_GOF$pos,Cabo = ex14_Cabo_GOF$ex14_Cabo_delta)
ex14_Camp_resis <- data.frame(pos =ex14_Camp_GOF$pos,Camp = ex14_Camp_GOF$ex14_Camp_delta)
#ex14_Savo_resis <- data.frame(pos =ex14_Savo_GOF$pos,Savo = ex14_Savo_GOF$ex14_Savo_delta)
ex14_NVP_resis <- data.frame(pos =ex14_NVP_GOF$pos,NVP = ex14_NVP_GOF$ex14_NVP_delta)
ex14_Tepo_resis <- data.frame(pos =ex14_Tepo_GOF$pos,Tepo = ex14_Tepo_GOF$ex14_Tepo_delta)
ex14_Glu_resis <- data.frame(pos =ex14_Glu_GOF$pos,Glu = ex14_Glu_GOF$ex14_Glu_delta)
ex14_Crizo_resis <- data.frame(pos =ex14_Mere_GOF$pos,Mere = ex14_Mere_GOF$ex14_Mere_delta)
ex14_Gle_resis <- data.frame(pos =ex14_Gle_GOF$pos,Gle = ex14_Gle_GOF$ex14_Gle_delta)
ex14_AMG458_resis <- data.frame(pos =ex14_AMG458_GOF$pos,AMG458 = ex14_AMG458_GOF$ex14_AMG458_delta)
ex14_Mere_resis <- data.frame(pos =ex14_Mere_GOF$pos,Mere = ex14_Mere_GOF$ex14_Mere_delta)


ex14_Crizo_resis2 <- ex14_Crizo_resis %>% group_by(pos) %>% dplyr::summarize(count = n())
ex14_Crizo_resis3 <- ex14_Crizo_resis %>% group_by(pos) %>% dplyr::summarise(avg = mean(Crizo)) 
ex14_Crizo_resis4 <- merge(ex14_Crizo_resis2,ex14_Crizo_resis3)
ex14_Crizo_resis4$inhibitor <- 'Crizotinib'

ex14_Cabo_resis2 <- ex14_Cabo_resis %>% group_by(pos) %>% dplyr::summarize(count = n())
ex14_Cabo_resis3 <- ex14_Cabo_resis %>% group_by(pos) %>% dplyr::summarise(avg = mean(Cabo)) 
ex14_Cabo_resis4 <- merge(ex14_Cabo_resis2,ex14_Cabo_resis3)
ex14_Cabo_resis4$inhibitor <- 'Cabozantinib'

ex14_Camp_resis2 <- ex14_Camp_resis %>% group_by(pos) %>% dplyr::summarize(count = n())
ex14_Camp_resis3 <- ex14_Camp_resis %>% group_by(pos) %>% dplyr::summarise(avg = mean(Camp)) 
ex14_Camp_resis4 <- merge(ex14_Camp_resis2,ex14_Camp_resis3)
ex14_Camp_resis4$inhibitor <- 'Capmatinib'

ex14_Tepo_resis2 <- ex14_Tepo_resis %>% group_by(pos) %>% dplyr::summarize(count = n())
ex14_Tepo_resis3 <- ex14_Tepo_resis %>% group_by(pos) %>% dplyr::summarise(avg = mean(Tepo)) 
ex14_Tepo_resis4 <- merge(ex14_Tepo_resis2,ex14_Tepo_resis3)
ex14_Tepo_resis4$inhibitor <- 'Tepotinib'

ex14_Glu_resis2 <- ex14_Glu_resis %>% group_by(pos) %>% dplyr::summarize(count = n())
ex14_Glu_resis3 <- ex14_Glu_resis %>% group_by(pos) %>% dplyr::summarise(avg = mean(Glu)) 
ex14_Glu_resis4 <- merge(ex14_Glu_resis2,ex14_Glu_resis3)
ex14_Glu_resis4$inhibitor <- 'Glumetinib'

ex14_NVP_resis2 <- ex14_NVP_resis %>% group_by(pos) %>% dplyr::summarize(count = n())
ex14_NVP_resis3 <- ex14_NVP_resis %>% group_by(pos) %>% dplyr::summarise(avg = mean(NVP)) 
ex14_NVP_resis4 <- merge(ex14_NVP_resis2,ex14_NVP_resis3)
ex14_NVP_resis4$inhibitor <- 'NVP-BVU972'

ex14_Mere_resis2 <- ex14_Mere_resis %>% group_by(pos) %>% dplyr::summarize(count = n())
ex14_Mere_resis3 <- ex14_Mere_resis %>% group_by(pos) %>% dplyr::summarise(avg = mean(Mere)) 
ex14_Mere_resis4 <- merge(ex14_Mere_resis2,ex14_Mere_resis3)
ex14_Mere_resis4$inhibitor <- 'Merestinib'

ex14_Gle_resis2 <- ex14_Gle_resis %>% group_by(pos) %>% dplyr::summarize(count = n())
ex14_Gle_resis3 <- ex14_Gle_resis %>% group_by(pos) %>% dplyr::summarise(avg = mean(Gle)) 
ex14_Gle_resis4 <- merge(ex14_Gle_resis2,ex14_Gle_resis3)
ex14_Gle_resis4$inhibitor <- 'Glesatinib'

ex14_AMG458_resis2 <- ex14_AMG458_resis %>% group_by(pos) %>% dplyr::summarize(count = n())
ex14_AMG458_resis3 <- ex14_AMG458_resis %>% group_by(pos) %>% dplyr::summarise(avg = mean(AMG458)) 
ex14_AMG458_resis4 <- merge(ex14_AMG458_resis2,ex14_AMG458_resis3)
ex14_AMG458_resis4$inhibitor <- 'AMG-458'


all_inhibs_ex14 <- rbind(ex14_Crizo_resis4,
                    ex14_Camp_resis4,
                    ex14_Tepo_resis4,
                    ex14_Glu_resis4,
                    #ex14_Savo_resis4, 
                    ex14_NVP_resis4,
                    ex14_Cabo_resis4,
                    ex14_Mere_resis4,
                    ex14_Gle_resis4,
                    ex14_AMG458_resis4)


inhib_ballon_ex14 <- ggballoonplot(all_inhibs_ex14, x = "pos", y = 'inhibitor',
                              size = "count", fill = "avg",
                              ggtheme = theme_bw())+
  scale_x_continuous(limits = c(1059, 1345),breaks = seq(1059,1345, by = 10) )+
  labs(y = "Inhibitor", x = "Position")+
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow') 
plot(inhib_ballon_ex14)



#####################################################################################################
###-------------------------------METdEx14 Resistance Balloon Plots  --------------------------------
#####################################################################################################


# Type 1 (crizo, camp, tepo, glu, nvp, savo)
# Type 2 (cabo, mere, gle)
# Type 1.5 (AMG)

met_Crizo_resis <- data.frame(pos =met_Crizo_GOF$pos,Crizo = met_Crizo_GOF$met_Crizo_delta)
met_Cabo_resis <- data.frame(pos =met_Cabo_GOF$pos,Cabo = met_Cabo_GOF$met_Cabo_delta)
met_Camp_resis <- data.frame(pos =met_Camp_GOF$pos,Camp = met_Camp_GOF$met_Camp_delta)
met_Savo_resis <- data.frame(pos =met_Savo_GOF$pos,Savo = met_Savo_GOF$met_Savo_delta)
met_NVP_resis <- data.frame(pos =met_NVP_GOF$pos,NVP = met_NVP_GOF$met_NVP_delta)
met_Tepo_resis <- data.frame(pos =met_Tepo_GOF$pos,Tepo = met_Tepo_GOF$met_Tepo_delta)
met_Glu_resis <- data.frame(pos =met_Glu_GOF$pos,Glu = met_Glu_GOF$met_Glu_delta)
met_Crizo_resis <- data.frame(pos =met_Mere_GOF$pos,Mere = met_Mere_GOF$met_Mere_delta)
met_Gle_resis <- data.frame(pos =met_Gle_GOF$pos,Gle = met_Gle_GOF$met_Gle_delta)
met_AMG458_resis <- data.frame(pos =met_AMG458_GOF$pos,AMG458 = met_AMG458_GOF$met_AMG458_delta)
met_Mere_resis <- data.frame(pos =met_Mere_GOF$pos,Mere = met_Mere_GOF$met_Mere_delta)


met_Crizo_resis2 <- met_Crizo_resis %>% group_by(pos) %>% dplyr::summarize(count = n())
met_Crizo_resis3 <- met_Crizo_resis %>% group_by(pos) %>% dplyr::summarise(avg = mean(Crizo)) 
met_Crizo_resis4 <- merge(met_Crizo_resis2,met_Crizo_resis3)
met_Crizo_resis4$inhibitor <- 'Crizotinib'

met_Cabo_resis2 <- met_Cabo_resis %>% group_by(pos) %>% dplyr::summarize(count = n())
met_Cabo_resis3 <- met_Cabo_resis %>% group_by(pos) %>% dplyr::summarise(avg = mean(Cabo)) 
met_Cabo_resis4 <- merge(met_Cabo_resis2,met_Cabo_resis3)
met_Cabo_resis4$inhibitor <- 'Cabozantinib'

met_Camp_resis2 <- met_Camp_resis %>% group_by(pos) %>% dplyr::summarize(count = n())
met_Camp_resis3 <- met_Camp_resis %>% group_by(pos) %>% dplyr::summarise(avg = mean(Camp)) 
met_Camp_resis4 <- merge(met_Camp_resis2,met_Camp_resis3)
met_Camp_resis4$inhibitor <- 'Capmatinib'

met_Tepo_resis2 <- met_Tepo_resis %>% group_by(pos) %>% dplyr::summarize(count = n())
met_Tepo_resis3 <- met_Tepo_resis %>% group_by(pos) %>% dplyr::summarise(avg = mean(Tepo)) 
met_Tepo_resis4 <- merge(met_Tepo_resis2,met_Tepo_resis3)
met_Tepo_resis4$inhibitor <- 'Tepotinib'

met_Glu_resis2 <- met_Glu_resis %>% group_by(pos) %>% dplyr::summarize(count = n())
met_Glu_resis3 <- met_Glu_resis %>% group_by(pos) %>% dplyr::summarise(avg = mean(Glu)) 
met_Glu_resis4 <- merge(met_Glu_resis2,met_Glu_resis3)
met_Glu_resis4$inhibitor <- 'Glumetinib'

met_NVP_resis2 <- met_NVP_resis %>% group_by(pos) %>% dplyr::summarize(count = n())
met_NVP_resis3 <- met_NVP_resis %>% group_by(pos) %>% dplyr::summarise(avg = mean(NVP)) 
met_NVP_resis4 <- merge(met_NVP_resis2,met_NVP_resis3)
met_NVP_resis4$inhibitor <- 'NVP-BVU972'

met_Mere_resis2 <- met_Mere_resis %>% group_by(pos) %>% dplyr::summarize(count = n())
met_Mere_resis3 <- met_Mere_resis %>% group_by(pos) %>% dplyr::summarise(avg = mean(Mere)) 
met_Mere_resis4 <- merge(met_Mere_resis2,met_Mere_resis3)
met_Mere_resis4$inhibitor <- 'Merestinib'

met_Gle_resis2 <- met_Gle_resis %>% group_by(pos) %>% dplyr::summarize(count = n())
met_Gle_resis3 <- met_Gle_resis %>% group_by(pos) %>% dplyr::summarise(avg = mean(Gle)) 
met_Gle_resis4 <- merge(met_Gle_resis2,met_Gle_resis3)
met_Gle_resis4$inhibitor <- 'Glesatinib'

met_AMG458_resis2 <- met_AMG458_resis %>% group_by(pos) %>% dplyr::summarize(count = n())
met_AMG458_resis3 <- met_AMG458_resis %>% group_by(pos) %>% dplyr::summarise(avg = mean(AMG458)) 
met_AMG458_resis4 <- merge(met_AMG458_resis2,met_AMG458_resis3)
met_AMG458_resis4$inhibitor <- 'AMG-458'

met_Savo_resis2 <- met_Savo_resis %>% group_by(pos) %>% dplyr::summarize(count = n())
met_Savo_resis3 <- met_Savo_resis %>% group_by(pos) %>% dplyr::summarise(avg = mean(Savo)) 
met_Savo_resis4 <- merge(met_Savo_resis2,met_Savo_resis3)
met_Savo_resis4$inhibitor <- 'Savolitinib'



all_inhibs_met <- rbind(met_Crizo_resis4,
                        met_Camp_resis4,
                        met_Tepo_resis4,
                        met_Glu_resis4,
                        met_Savo_resis4, 
                        met_NVP_resis4,
                        met_Cabo_resis4,
                        met_Mere_resis4,
                        met_Gle_resis4,
                        met_AMG458_resis4,
                        met_Savo_resis4)


inhib_ballon_met <- ggballoonplot(all_inhibs_met, x = "pos", y = 'inhibitor',
                                  size = "count", fill = "avg",
                                  ggtheme = theme_bw())+
  scale_x_continuous(limits = c(1059, 1345),breaks = seq(1059,1345, by = 10) )+
  labs(y = "Inhibitor", x = "Position")+
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow') 
plot(inhib_ballon_met)

#####################################################################################################
###----------------------------- Resistance Heatmap + Clustering  ----------------------------------
###------------------------------- From DMSO subtracted data ---------------------------------------
#####################################################################################################

resis_heatmap_df <- data.frame() 

ex14_Crizo_resis.t <- data.frame(t(ex14_Crizo_resis[-1]))
colnames(ex14_Crizo_resis.t) <-ex14_Crizo_resis[, 1]

ex14_Cabo_resis.t <- data.frame(t(ex14_Cabo_resis[-1]))
colnames(ex14_Cabo_resis.t) <-ex14_Cabo_resis[, 1]

ex14_Mere_resis.t <- data.frame(t(ex14_Mere_resis[-1]))
colnames(ex14_Mere_resis.t) <-ex14_Mere_resis[, 1]

ex14_AMG458_resis.t <- data.frame(t(ex14_AMG458_resis[-1]))
colnames(ex14_AMG458_resis.t) <-ex14_AMG458_resis[, 1]

ex14_NVP_resis.t <- data.frame(t(ex14_NVP_resis[-1]))
colnames(ex14_NVP_resis.t) <-ex14_NVP_resis[, 1]

ex14_Gle_resis.t <- data.frame(t(ex14_Gle_resis[-1]))
colnames(ex14_Gle_resis.t) <-ex14_Gle_resis[, 1]

ex14_Glu_resis.t <- data.frame(t(ex14_Glu_resis[-1]))
colnames(ex14_Glu_resis.t) <-ex14_Glu_resis[, 1]

ex14_Tepo_resis.t <- data.frame(t(ex14_Tepo_resis[-1]))
colnames(ex14_Tepo_resis.t) <-ex14_Tepo_resis[, 1]

ex14_Camp_resis.t <- data.frame(t(ex14_Camp_resis[-1]))
colnames(ex14_Camp_resis.t) <-ex14_Camp_resis[, 1]



m <- rbind.fill(ex14_Crizo_resis.t,
                       ex14_Camp_resis.t,
                       ex14_Tepo_resis.t,
                       ex14_Glu_resis.t,
                       ex14_Gle_resis.t,
                       ex14_AMG458_resis.t,
                       ex14_Mere_resis.t, 
                       ex14_Cabo_resis.t,
                       ex14_NVP_resis.t
                       )



#####################################################################################################
###-------------------MET+Ex14 Common resistance sites for each inhibitor type ----------------------
#####################################################################################################

# lists of positional averages, nonsense and wt syn NOT included 
avg_ex14_AMG458 <- data.frame(ex14_AMG458_delta %>% filter (mutation_type != "N" & mutation_type != "S") %>% group_by(pos) %>% summarise(ex14_AMG458_avg = mean(ex14_AMG458_delta, na.rm=TRUE)))
avg_ex14_Crizo <-  data.frame(ex14_Crizo_delta %>% filter (mutation_type != "N" & mutation_type != "S") %>% group_by(pos) %>% summarise(ex14_Crizo_avg = mean(ex14_Crizo_delta, na.rm=TRUE)))
avg_ex14_Cabo <-  data.frame(ex14_Cabo_delta %>% filter (mutation_type != "N" & mutation_type != "S") %>% group_by(pos) %>% summarise(ex14_Cabo_avg = mean(ex14_Cabo_delta, na.rm=TRUE)))
avg_ex14_Camp <-  data.frame(ex14_Camp_delta %>% filter (mutation_type != "N" & mutation_type != "S") %>% group_by(pos) %>% summarise(ex14_Camp_avg = mean(ex14_Camp_delta, na.rm=TRUE)))
avg_ex14_Tepo <-  data.frame(ex14_Tepo_delta %>% filter (mutation_type != "N" & mutation_type != "S") %>% group_by(pos) %>% summarise(ex14_Tepo_avg = mean(ex14_Tepo_delta, na.rm=TRUE)))
avg_ex14_Gle <-  data.frame(ex14_Gle_delta %>% filter (mutation_type != "N" & mutation_type != "S") %>% group_by(pos) %>% summarise(ex14_Gle_avg = mean(ex14_Gle_delta, na.rm=TRUE)))
avg_ex14_Glu <-  data.frame(ex14_Glu_delta %>% filter (mutation_type != "N" & mutation_type != "S") %>% group_by(pos) %>% summarise(ex14_Glu_avg = mean(ex14_Glu_delta, na.rm=TRUE)))
avg_ex14_Savo <- data.frame( ex14_Savo_delta %>% filter (mutation_type != "N" & mutation_type != "S") %>% group_by(pos) %>% summarise(ex14_Savo_avg = mean(ex14_Savo_delta, na.rm=TRUE)))
avg_ex14_Mere <-  data.frame(ex14_Mere_delta %>% filter (mutation_type != "N" & mutation_type != "S") %>% group_by(pos) %>% summarise(ex14_Mere_avg = mean(ex14_Mere_delta, na.rm=TRUE)))
avg_ex14_NVP <-  data.frame(ex14_NVP_delta %>% filter (mutation_type != "N" & mutation_type != "S") %>% group_by(pos) %>% summarise(ex14_NVP_avg = mean(ex14_NVP_delta, na.rm=TRUE)))


# Type 1 (crizo, camp, tepo, glu, nvp, savo)
# Type 2 (cabo, mere, gle)
# Type 1.5 (AMG)
ex14_avg_inhib_list <- list(avg_ex14_Crizo,  
                            avg_ex14_Camp, 
                            avg_ex14_Tepo,  
                            avg_ex14_Glu,
                            avg_ex14_NVP,
                            avg_ex14_Savo,
                            avg_ex14_Cabo, 
                            avg_ex14_Mere, 
                            avg_ex14_Gle, 
                            avg_ex14_Gle, 
                            avg_ex14_AMG458)
  
ex14_avg_inhib <- data.frame(pos = avg_ex14_Crizo$pos,
                            Crizo = avg_ex14_Crizo$ex14_Crizo_avg,
                            Camp = avg_ex14_Camp$ex14_Camp_avg,
                            Tepo = avg_ex14_Tepo$ex14_Tepo_avg,
                            Glu = avg_ex14_Glu$ex14_Glu_avg,
                            NVP = avg_ex14_NVP$ex14_NVP_avg,
                            Savo = avg_ex14_Savo$ex14_Savo_avg,
                            Cabo = avg_ex14_Cabo$ex14_Cabo_avg,
                            Mere = avg_ex14_Mere$ex14_Mere_avg,
                            Gle = avg_ex14_Gle$ex14_Gle_avg,
                            AMG458 = avg_ex14_AMG458$ex14_AMG458_avg)



# create an alignment w/ averages 
# create an alignment w/ average positional resistance 
ex14_avg_inhib.m <-ex14_avg_inhib%>%
  gather(key="inhib", value ="value",-pos) %>% # convert data to long format
  setNames(c("pos","inhib","value")) %>%   # rename columns
  mutate(pos=factor(pos)) %>%   # convert ypos to factor 
  mutate(value=as.numeric(value)) # convert value to numeric (also converts '-' to NA, gives a warning)


## inhibitor heatmap 
inhib_order <- c('AMG458','Cabo','Mere','Gle','Savo','Tepo','Glu','NVP','Camp','Crizo')

met_wt_1 = str_split(substr(met_wt_sequence, 1, 144), '')[[1]]
met_wt_2 = str_split(substr(met_wt_sequence, 145, 287), '')[[1]]

  
inhib_row1 = ggplot(data = ex14_avg_inhib.m %>%filter(pos %in% c(1059:1202)),
                  aes(x = pos, y = factor(inhib, level = inhib_order), fill =  value )) +
  geom_tile( size = 0.2) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow') + 
  scale_color_manual(values = c(NA,'green')) +
  coord_fixed(ratio = 1) +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgray"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 4),
    axis.text.y = element_text(size = 4),
    axis.text = element_text(size = 4),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
    axis.ticks = element_blank(),
    legend.position="none"
  ) +
  labs(y = "Mutation", x = "Position")

inhib_row2 = ggplot(data = ex14_avg_inhib.m %>% filter(pos %in% c(1203:1345)), 
                  aes(x = pos, y = factor(inhib, level = inhib_order), fill =  value)) +
  geom_tile( size = 0.2) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow') + 
  scale_color_manual(values = c(NA,'green')) +
  coord_fixed(ratio = 1) +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgray"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 4),
    axis.text.y = element_text(size = 4),
    axis.text = element_text(size = 4),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
    axis.ticks = element_blank(),
    legend.position="none"
  ) +
  labs(y = "Mutation", x = "Position")

inhib_heatmap= ggarrange(inhib_row1, inhib_row2, nrow = 2, ncol = 1) 
plot(inhib_heatmap)
ggsave("Ex14_Avg_heatmap.pdf", height = 7, width = 8.5, inhib_heatmap)






#####################################################################################################
###------------------------------------- PCA of METdEx14 --------------------------------------------
#####################################################################################################

#example PCA 
library(factoextra)
data(iris)
res.pca <- prcomp(iris[, -5],  scale = TRUE)
fviz_pca_ind(res.pca, label="none", habillage=iris$Species,
             addEllipses=TRUE, ellipse.level=0.95,
             palette = c("#999999", "#E69F00", "#56B4E9"))


## need to create a data frame that goes: 

# score...inhib....type
#...3.....crizo....1

# met+ex14
ex14_Crizo_pca <- data.frame(pos=ex14_Crizo_filtered$pos, score=ex14_Crizo_filtered$ex14_Crizo_delta, inhibitor='Crizotinib',type = '1')
ex14_Camp_pca <- data.frame(pos=ex14_Camp_filtered$pos, score=ex14_Camp_filtered$ex14_Camp_delta, inhibitor='Capmantinib',type = '1')
ex14_Tepo_pca <- data.frame(pos=ex14_Tepo_filtered$pos, score=ex14_Tepo_filtered$ex14_Tepo_delta, inhibitor='Tepotinib',type = '1')
ex14_Glu_pca <- data.frame(pos=ex14_Glu_filtered$pos, score=ex14_Glu_filtered$ex14_Glu_delta, inhibitor='Glumetinib',type = '1')
ex14_NVP_pca <- data.frame(pos=ex14_NVP_filtered$pos, score=ex14_NVP_filtered$ex14_NVP_delta, inhibitor='NVP',type = '1')
#ex14_Savo_pca <- data.frame(pos=ex14_Savo_filtered$pos, score=ex14_Savo_filtered$ex14_Savo_delta, inhibitor='Savolitinib',type = '1')
ex14_Cabo_pca <- data.frame(pos=ex14_Cabo_filtered$pos, score=ex14_Cabo_filtered$ex14_Cabo_delta, inhibitor='Cabozantinib',type = '2')
ex14_Mere_pca <- data.frame(pos=ex14_Mere_filtered$pos, score=ex14_Mere_filtered$ex14_Mere_delta, inhibitor='Merestinib',type = '2')
ex14_Gle_pca <- data.frame(pos=ex14_Gle_filtered$pos, score=ex14_Gle_filtered$ex14_Gle_delta, inhibitor='Glesatinib',type = '2')
ex14_AMG458_pca <- data.frame(pos=ex14_AMG45_filtered$hgvs, score=ex14_AMG45_filtered$ex14_AMG458_delta, inhibitor='AMG-458',type = '1.5')

ex14_pca <- data.frame(rbind(ex14_Crizo_pca,
                  ex14_Camp_pca,
                  ex14_Tepo_pca,
                  ex14_Glu_pca,
                  ex14_NVP_pca,
                  ex14_Cabo_pca,
                  ex14_Mere_pca,
                  ex14_Gle_pca,
                  ex14_AMG458_pca))

ex14_pca2 <- data.frame(ex14_delta_scores)
ex14_pca3<-data.matrix(ex14_pca2)
ex14_pca3<-na.omit(ex14_pca3)

res.pca2 <- prcomp(ex14_pca3,  scale = TRUE)
fviz_pca_ind(res.pca2,
             addEllipses=TRUE, ellipse.level=0.95,
             palette = c("#999999", "#E69F00", "#56B4E9"))

# met-ex14 
met_Crizo_pca <- data.frame(pos=met_Crizo_filtered$pos, score=met_Crizo_filtered$met_Crizo_delta, inhibitor='Crizotinib',type = '1')
met_Camp_pca <- data.frame(pos=met_Camp_filtered$pos, score=met_Camp_filtered$met_Camp_delta, inhibitor='Capmantinib',type = '1')
met_Tepo_pca <- data.frame(pos=met_Tepo_filtered$pos, score=met_Tepo_filtered$met_Tepo_delta, inhibitor='Tepotinib',type = '1')
met_Glu_pca <- data.frame(pos=met_Glu_filtered$pos, score=met_Glu_filtered$met_Glu_delta, inhibitor='Glumetinib',type = '1')
met_NVP_pca <- data.frame(pos=met_NVP_filtered$pos, score=met_NVP_filtered$met_NVP_delta, inhibitor='NVP',type = '1')
met_Savo_pca <- data.frame(pos=met_Savo_filtered$pos, score=met_Savo_filtered$met_Savo_delta, inhibitor='Savolitinib',type = '1')
met_Cabo_pca <- data.frame(pos=met_Cabo_filtered$pos, score=met_Cabo_filtered$met_Cabo_delta, inhibitor='Cabozantinib',type = '2')
met_Mere_pca <- data.frame(pos=met_Mere_filtered$pos, score=met_Mere_filtered$met_Mere_delta, inhibitor='Merestinib',type = '2')
met_Gle_pca <- data.frame(pos=met_Gle_filtered$pos, score=met_Gle_filtered$met_Gle_delta, inhibitor='Glesatinib',type = '2')
met_AMG458_pca <- data.frame(pos=met_AMG458_filtered$hgvs, score=met_AMG458_filtered$met_AMG458_delta, inhibitor='AMG-458',type = '1.5')

met_pca <- data.frame(rbind(met_Crizo_pca,
                  met_Camp_pca,
                  met_Tepo_pca,
                  met_Glu_pca,
                  met_NVP_pca,
                  met_Cabo_pca,
                  met_Mere_pca,
                  met_Gle_pca,
                  met_Savo_pca,
                  met_AMG458_pca))

met_pca2 <- data.matrix(met_pca)


res.pca3 <- prcomp(met_pca2, scale = TRUE)
fviz_pca_ind(res.pca3, habillage=met_pca$type,
             addEllipses=TRUE, ellipse.level=0.95,
             palette = c("#999999", "#E69F00", "#56B4E9"))





# Type 1 (crizo, camp, tepo, glu, nvp, savo)
# Type 2 (cabo, mere, gle)
# Type 1.5 (AMG)


# prepare data frame for hierarchical clustering 
ex14_avg_inhib.t <- data.frame(t(ex14_avg_inhib[-1]))
colnames(ex14_avg_inhib.t) <- ex14_avg_inhib[, 1]
ex14_avg_inhib.t["inhib_type"] <- c("type1","type1","type1","type1","type1","type1","type2","type2","type2","type1.5") #ID inhibitor types 


ex14_avg_inhib.1 <- data.frame(stack(ex14_avg_inhib[-1]))
ex14_avg_inhib.2 <- data.frame(stack(ex14_avg_inhib[1]))
ex14_avg_inhib.1["pos"] <- c(ex14_avg_inhib.2[1])
ex14_avg_inhib.1<-ex14_avg_inhib.1[,c(3,1,2)]

ex14_avg_inhib.1 <- ex14_avg_inhib.1 %>% mutate(type =
                     case_when(str_detect(ind, "Crizo") ~ "type1", 
                               str_detect(ind, "Camp") ~ "type1", 
                               str_detect(ind, "Tepo") ~ "type1", 
                               str_detect(ind, "Glu") ~ "type1", 
                               str_detect(ind,"NVP") ~ "type1", 
                               str_detect(ind, "Savo") ~ "type1", 
                               str_detect(ind,"Cabo")~ "type2", 
                               str_detect(ind,"Mere") ~ "type2",
                               str_detect(ind,"Gle") ~ "type2",
                               str_detect(ind,"AMG458") ~ "type1.5"
                               )
)

ex14_avg_inhib.transpoaed<- t(ex14_avg_inhib.t)


library(cluster)
autoplot(pam(ex14_avg_inhib[-1],3), frame = TRUE, frame.type = 'norm')

# hierarchical clustering of positions
clusters <- hclust(dist(ex14_avg_inhib.1))
plot(clusters)
# hierarchical clustering of inhibitors
clusters2 <- hclust(dist(ex14_avg_inhib.t))
plot(clusters2)

library(ggfortify)
pca_res <- prcomp(ex14_avg_inhib[-1], scale. = TRUE)
autoplot(pca_res)
autoplot(pca_res, data =ex14_avg_inhib.t, colour = 'inhib_type')






  

#####################################################################################################
###---------------------------- Mapping delta scores on PDBs ----------------------------------------
#####################################################################################################

# map scores for each inhibitor selection to their respective structure 
# there are only experimentally determined structures for 6/11 of the inhibitors used
# for the DMSO control, that info should be similar to the previous DMS, and will be plotted on active and inactive structures

#### DMSO : 3R7O 
map_met_DMSO = met_DMSO_scores %>% filter (mutation_type != "N" & mutation_type != "S")
avg_met_DMSO = map_met_DMSO %>% group_by(pos) %>% summarise(met_DMSO_avg = mean(met_DMSO_score, na.rm=TRUE))
met_3R7O_DMSO <- read.pdb("3R7O")
x = map_scores_pdb(met_3R7O_DMSO, avg_met_DMSO, "met_DMSO_avg")
write.pdb(x, file="met_3R7O_DMSO.pdb")

map_ex14_DMSO = ex14_DMSO_scores %>% filter (mutation_type != "N" & mutation_type != "S")
avg_ex14_DMSO = map_ex14_DMSO %>% group_by(pos) %>% summarise(ex14_DMSO_avg = mean(ex14_DMSO_score, na.rm=TRUE))
ex14_3R7O_DMSO <- read.pdb("3R7O")
x = map_scores_pdb(ex14_3R7O_DMSO, avg_ex14_DMSO, "ex14_DMSO_avg")
write.pdb(x, file="ex14_3R7O_DMSO.pdb")

#### crizotinib: 2WGJ 
# met
map_met_crizo = met_Crizo_delta %>% filter (mutation_type != "N" & mutation_type != "S")
avg_met_crizo = map_met_crizo %>% group_by(pos) %>% summarise(crizo_delta = mean(met_Crizo_delta, na.rm=TRUE))
met_2WGJ_crizo <- read.pdb("2WGJ")
x = map_scores_pdb(met_2WGJ_crizo, avg_met_crizo, "crizo_delta")
write.pdb(x, file="met_2WGJ_crizo_delta.pdb")
# ex14
map_ex14_crizo = ex14_Crizo_delta %>% filter (mutation_type != "N" & mutation_type != "S")
avg_ex14_crizo = ddply(map_ex14_crizo,c("pos"),summarise,mean=mean(ex14_Crizo_delta))
ex14_2WGJ_crizo <- read.pdb("2WGJ")
x = map_scores_pdb(ex14_2WGJ_crizo, avg_ex14_crizo, "ex14_crizo_delta")
write.pdb(x, file="ex14_2WGJ_crizo_delta.pdb")




### Tivantinib: 3RHK 
# met
map_met_Tiv = met_Tiv_delta %>% filter (mutation_type != "N" & mutation_type != "S")
avg_met_Tiv = map_met_Tiv %>% group_by(pos) %>% summarise(Tiv_delta = mean(Tiv_delta, na.rm=TRUE))
met_3RHK_Tiv <- read.pdb("3RHK")
x = map_scores_pdb(met_3RHK_Tiv, avg_met_Tiv, "Tiv_delta")
write.pdb(x, file="met_3RHK_Tiv_delta.pdb")
# ex14
map_ex14_Tiv = ex14_Tiv_delta %>% filter (mutation_type != "N" & mutation_type != "S")
avg_ex14_Tiv = map_ex14_Tiv %>% group_by(pos) %>% summarise(Tiv_delta = mean(Tiv_delta, na.rm=TRUE))
ex14_3RHK_Tiv <- read.pdb("3RHK")
x = map_scores_pdb(ex14_3RHK_Tiv, avg_ex14_Tiv, "Tiv_delta")
write.pdb(x, file="ex14_3RHK_Tiv_delta.pdb")


### AMG458: 5T3Q
# met
map_met_AMG458 = met_AMG458_delta %>% filter (mutation_type != "N" & mutation_type != "S")
avg_met_AMG458 = map_met_AMG458 %>% group_by(pos) %>% summarise(AMG458_delta = mean(AMG458_delta, na.rm=TRUE))
met_5T3Q_AMG458 <- read.pdb("5T3Q")
x = map_scores_pdb(met_5T3Q_AMG458, avg_met_AMG458, "AMG458_delta")
write.pdb(x, file="met_5T3Q_AMG458_delta.pdb")
# ex14
map_ex14_AMG458 = ex14_AMG458_delta %>% filter (mutation_type != "N" & mutation_type != "S")
avg_ex14_AMG458 = map_ex14_AMG458 %>% group_by(pos) %>% summarise(AMG458_delta = mean(AMG458_delta, na.rm=TRUE))
ex14_5T3Q_AMG458 <- read.pdb("5T3Q")
x = map_scores_pdb(ex14_5T3Q_AMG458, avg_ex14_AMG458, "AMG458_delta")
write.pdb(x, file="ex14_5T3Q_AMG458_delta.pdb")


### NVP-BVU972: 3QTI 
# met
map_met_NVP = met_NVP_delta %>% filter (mutation_type != "N" & mutation_type != "S")
avg_met_NVP = map_met_NVP %>% group_by(pos) %>% summarise(NVP_delta = mean(NVP_delta, na.rm=TRUE))
met_3QTI_NVP <- read.pdb("3QTI")
x = map_scores_pdb(met_3QTI_NVP, avg_met_NVP, "NVP_delta")
write.pdb(x, file="met_3QTI_NVP_delta.pdb")
# ex14
map_ex14_NVP = ex14_NVP_delta %>% filter (mutation_type != "N" & mutation_type != "S")
avg_ex14_NVP = map_ex14_NVP %>% group_by(pos) %>% summarise(NVP_delta = mean(NVP_delta, na.rm=TRUE))
ex14_3QTI_NVP <- read.pdb("3QTI")
x = map_scores_pdb(ex14_3QTI_NVP, avg_ex14_NVP, "NVP_delta")
write.pdb(x, file="ex14_3QTI_NVP_delta.pdb")


### Merestinib: 4EEV 
map_met_Mere = met_Mere_delta %>% filter (mutation_type != "N" & mutation_type != "S")
avg_met_Mere = map_met_Mere %>% group_by(pos) %>% summarise(Mere_delta = mean(Mere_delta, na.rm=TRUE))
met_4EEV_Mere <- read.pdb("4EEV")
x = map_scores_pdb(met_4EEV_Mere, avg_met_Mere, "Mere_delta")
write.pdb(x, file="met_4EEV_Mere_delta.pdb")
# ex14
map_ex14_Mere = ex14_Mere_delta %>% filter (mutation_type != "N" & mutation_type != "S")
avg_ex14_Mere = map_ex14_Mere %>% group_by(pos) %>% summarise(Mere_delta = mean(Mere_delta, na.rm=TRUE))
ex14_4EEV_Mere <- read.pdb("4EEV")
x = map_scores_pdb(ex14_4EEV_Mere, avg_ex14_Mere, "Mere_delta")
write.pdb(x, file="ex14_4EEV_Mere_delta.pdb")


### Savolitinib: 6SDE
# met
map_met_Savo = met_Savo_delta %>% filter (mutation_type != "N" & mutation_type != "S")
avg_met_Savo = map_met_Savo %>% group_by(pos) %>% summarise(Savo_delta = mean(Savo_delta, na.rm=TRUE))
met_6SDE_Savo <- read.pdb("6SDE")
x = map_scores_pdb(met_6SDE_Savo, avg_met_Savo, "Savo_delta")
write.pdb(x, file="met_6SDE_Savo_delta.pdb")
# ex14
map_ex14_Savo = ex14_Savo_delta %>% filter (mutation_type != "N" & mutation_type != "S")
avg_ex14_Savo = map_ex14_Savo %>% group_by(pos) %>% summarise(Savo_delta = mean(Savo_delta, na.rm=TRUE))
ex14_6SDE_Savo <- read.pdb("6SDE")
x = map_scores_pdb(ex14_6SDE_Savo, avg_ex14_Savo, "Savo_delta")
write.pdb(x, file="ex14_6SDE_Savo_delta.pdb")


### Cabozantinib (representative Type II inhib): 3F82
# met
map_met_Cabo = met_Cabo_delta %>% filter (mutation_type != "N" & mutation_type != "S")
avg_met_Cabo = map_met_Cabo %>% group_by(pos) %>% summarise(Cabo_delta = mean(met_Cabo_delta, na.rm=TRUE))
met_3F82_Cabo <- read.pdb("3F82")
x = map_scores_pdb(met_3F82_Cabo, avg_met_Cabo, "Cabo_delta")
write.pdb(x, file="met_3F82_Cabo_delta.pdb")
# ex14
map_ex14_Cabo = ex14_Cabo_delta %>% filter (mutation_type != "N" & mutation_type != "S")
avg_ex14_Cabo = map_ex14_Cabo %>% group_by(pos) %>% summarise(ex14_Cabo_delta = mean(ex14_Cabo_delta, na.rm=TRUE))
ex14_3F82_Cabo <- read.pdb("3F82")
x = map_scores_pdb(ex14_3F82_Cabo, avg_ex14_Cabo, "ex14_Cabo_delta")
write.pdb(x, file="ex14_3F82_Cabo_delta.pdb")


###########3



# generate subtraction heatmap 
order <- c('X','H', 'K', 'R', 'D', 'E',
           'C','M','N','Q','S','T','A','I','L','V','W','F',
           'Y','G','P')

variant_names <- c('Stop','H', 'K', 'R', 'D', 'E',
                   'C','M','N','Q','S','T','A','I','L','V','W','F',
                   'Y','G','P')


met_wt_fasta = "Met_wt.fasta"

met_wt_fasta_con=file(met_wt_fasta, open="r")
met_wt_sequence = readLines(met_wt_fasta_con)[2]
close(met_wt_fasta_con)

met_wt_1 = str_split(substr(met_wt_sequence, 1, 144), '')[[1]]
met_wt_2 = str_split(substr(met_wt_sequence, 145, 287), '')[[1]]

Met_row1 = ggplot(data = ex14_Enrich2_crizo %>%filter(pos %in% c(1059:1202)),
                  aes(x = pos, y = factor(variants, level = order, labels = variant_names), fill = score )) +
  geom_tile(size = 0.2) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(-12,4)) + 
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
    legend.position="none"
  ) +
  labs(y = "Mutation", x = "Position")

Met_row2 = ggplot(data = ex14_Enrich2_crizo %>% filter(pos %in% c(1203:1345)), 
                  aes(x = pos, y = factor(variants, level = order, labels = variant_names), fill = score)) +
  geom_tile(size = 0.2) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(-12,4)) + 
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
    legend.position="none"
  ) +
  labs(y = "Mutation", x = "Position")

Met_DMS = ggarrange(Met_row1,Met_row2,nrow = 2, ncol = 1)

ggsave("Ex14_MET_Enrich2_Crizo_.pdf", height = 7, width = 8.5, Met_DMS)
ggsave("Ex14_MET_Enrich2_Crizo_.png", height = 7, width = 8.5, Met_DMS)

