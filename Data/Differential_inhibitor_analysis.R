###################################################################################################################
# Differential_inhibitor_analysis.R 
# by g.estevam @ UCSF 
# this script looks at mutations that have reverse sensitivities for inhibitors

# sequential parings type 1 --> 2 
# Savolitinib, cabozantinb
# Crizotinib, glesatinib 
# Crizotinib, merestinib
# Crizotinib, cabozantinib

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
library(reshape2)
source("dms_analysis_utilities.R")

library(scales)
library(colorspace)
library(ggpubr)
library(rstatix)
library(GGally)
library(lattice)
library(cdata)


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
met_rosace_norm_scores  <- met_rosace_norm_scores   %>% mutate(position = position + 1058)

# wt 
ex14_rosace_norm_scores <- as.data.frame(fread("ex14_scores_unfiltered.tsv")) # here ex14 means MET+Ex14 
ex14_rosace_norm_scores  <- ex14_rosace_norm_scores  %>% mutate(position = position + 1058)

## DMSO subtraction = detla.score, which is added as a new column 
ex14_rosace_norm_scores <- ex14_rosace_norm_scores %>%group_by(variants) %>%mutate(delta.score = mean - mean[key == "DMSO"]) %>%ungroup()


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
# these are the scores that are used for plotting the positional Manhattan type plots 
# this also provides the filtering for all the scores plotted in 3D


ex14_rosace_resistant_pval0.5 <- ex14_rosace_norm_scores_filtered %>%
  group_by(variants) %>%
  filter(
    (delta.score >= 0 & pval0.5.up <= 0.1 & mean >= 0 & mean[match("DMSO", key)] <= 0)) %>%
  ungroup()


DMSO_df<- data.frame(type = ex14_rosace_scores_wide$type, score= ex14_rosace_scores_wide$DMSO)
DMSO_wt <- DMSO_df %>% filter(type == "synonymous")
5*sd(DMSO_wt$score)


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



###-------------------------------  PDBs-------------------------------------------

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


###------------------------------------- clinical positions ---------------------------------------------

clinical_cancer <- c("1250", "1070", "1200", "1094", "1230", "1228", "1235", "1100", 
                     "1092", "1123", "1286", "1162", "1312", "1192", "1096", "1148", 
                     "1352", "1124", "1238", "1325", "1205", "1180", "1022", "1251",
                     "1186", "1287", "1375", "1327", "1063", "1073", "1138", "1199", 
                     "1314", "1279", "1308", "1091", "1140", "1311", "1160", "1351", 
                     "1174", "1153", "1258", "1388", "1319", "1060", "1099", "1170", 
                     "1050", "1121", "1014", "1236")
# documented resistanct
#("1211", "1163", "1155", "1195")




# ---------------------- Crizo Cabozantinib (II)----------------------


# DMSO delta scores plotted against each other 
# filter scores based on 

Crizo_ <- ex14_rosace_norm_scores_filtered %>% filter(key %in% c("Crizo", "DMSO"))
Crizo_resis <- Crizo_ %>% group_by(variants) %>% filter((delta.score >= 0 & pval0.5.up <= 0.1 & mean >= 0 & mean[match("DMSO", key)] <= 0)) %>% ungroup()

Cabo_ <- ex14_rosace_norm_scores_filtered %>% filter(key %in% c("Cabo", "DMSO"))
Cabo_resis <- Cabo_ %>% group_by(variants) %>% filter((delta.score >= 0 & pval0.5.up <= 0.1 & mean >= 0 & mean[match("DMSO", key)] <= 0)) %>% ungroup()


Crizo_df <- data.frame (hgvs = Crizo_$variants,
                        pos = Crizo_$position, 
                        mutation = Crizo_$mutation,
                        Crizo = Crizo_$delta.score,
                        Crizo_pval0.5 = Crizo_$pval0.5,
                        Crizo_pval0.5.up = Crizo_$pval0.5)

Cabo_df <- data.frame (hgvs = Cabo_$variants,
                       pos = Cabo_$position, 
                       mutation = Cabo_$mutation, 
                       Cabo = Cabo_$delta.score,
                       Cabo_pval0.5 = Cabo_$pval0.5,
                       Cabo_pval0.5.up = Cabo_$pval0.5)


Crizo_vs_Cabo <- merge(Crizo_df, Cabo_df, by = "hgvs", all = TRUE)


Crizo_GOF_Cabo_LOF <- Crizo_vs_Cabo %>% filter (Crizo>= 0.75 & Cabo <= 0 ) %>% filter(pos.x %in% Crizo_resis$position)
Cabo_GOF_Crizo_LOF <- Crizo_vs_Cabo %>% filter (Cabo>= 0.75 & Crizo <= 0 ) %>% filter(pos.y %in% Cabo_resis$position)


Crizo_vs_Cabo <- data.frame(pos = ex14_rosace_norm_wide$pos,
                          mutation = ex14_rosace_norm_wide$mutation,
                          Crizo = ex14_rosace_norm_wide$Crizo,
                          Cabo = ex14_rosace_norm_wide$Cabo)


# parse out resistance mutations and clinical significance 
Crizo_GOF_Cabo_LOF <- Crizo_vs_Cabo  %>% filter (Crizo>= 0.75 & Cabo <= 0) %>% filter(pos %in% Crizo_resis$position)
Cabo_GOF_Crizo_LOF <- Crizo_vs_Cabo  %>% filter (Cabo>= 0.75 & Crizo <= 0) %>% filter(pos %in% Cabo_resis$position)


Crizo_GOF_Cabo_LOF <- Crizo_vs_Cabo  %>% filter (Crizo>= 0.75 & Cabo <= 0)
Cabo_GOF_Crizo_LOF <- Crizo_vs_Cabo  %>% filter (Cabo>= 0.75 & Crizo <= 0) 


plot_Crizo_vs_Cabo <- ggplot() +
  geom_point(data = ex14_rosace_norm_wide , aes(x = Cabo, y = Crizo),color = "grey60",alpha=0.3) + 
  geom_point(data = Crizo_GOF_Cabo_LOF, aes(x = Cabo, y = Crizo), color ="deeppink") + 
  geom_point(data = Cabo_GOF_Crizo_LOF, aes(x = Cabo, y = Crizo), color = "deepskyblue") + 
  geom_point(aes(x =-0.95939772, y = 1.8697423), color = "black") + 
  xlab("Cabozantinib") + ylab("Crizotinib") + coord_fixed(ratio = 1) +
  #geom_vline(xintercept=0, color = "red", linetype = "dashed")+
  #geom_hline(yintercept=0, color = "red", linetype = "dashed")+
  theme_classic() 
plot(plot_Crizo_vs_Cabo)




ex14_2WGJ_crizo <- read.pdb("2WGJ")
cabo_docked <- read.pdb("cabozantinib_4eev_docked.pdb")

avg_Crizo_GOF_Cabo_LOF = Crizo_GOF_Cabo_LOF %>% group_by(pos) %>% dplyr::summarise(mean = mean(Crizo, na.rm=TRUE))
x = map_scores_pdb(crizo_2WGJ, avg_Crizo_GOF_Cabo_LOF, "mean")
write.pdb(x, file="avg_Crizo_GOF_cabo_LOF.pdb")

avg_Cabo_GOF_Crizo_LOF = Cabo_GOF_Crizo_LOF %>% group_by(pos) %>% dplyr::summarise(mean = mean(Cabo, na.rm=TRUE))
x = map_scores_pdb(cabo_docked , avg_Cabo_GOF_Crizo_LOF, "mean")
write.pdb(x, file="avg_Cabo_GOF_Crizo_LOF.pdb")


# ---------------------- Crizo Capmatinib(II)----------------------


Crizo_vs_Camp <- data.frame(pos = ex14_rosace_norm_wide$pos, 
                            Crizo = ex14_rosace_norm_wide$Crizo,
                            mutation = ex14_rosace_norm_wide$mutation,
                            Camp = ex14_rosace_norm_wide$Camp )


Crizo_GOF_Camp_LOF <- Crizo_vs_Camp %>% filter (Crizo>= 0.75 & Camp <= 0)
Camp_GOF_Crizo_LOF <- Crizo_vs_Camp  %>% filter (Camp>= 0.75 & Crizo <= 0)

plot_Crizo_vs_Camp <- ggplot() +
  geom_point(data = ex14_rosace_norm_wide , aes(x = Camp, y = Crizo),color = "grey60",alpha=0.3) + 
  geom_point(data = Crizo_GOF_Camp_LOF, aes(x = Camp, y = Crizo), color ="deeppink") + 
  geom_point(data = Camp_GOF_Crizo_LOF, aes(x = Camp, y = Crizo), color = "deepskyblue") + 
  xlab("Capmantinib") + ylab("Crizotinib") + coord_fixed(ratio = 1) +
  #geom_vline(xintercept=0, color = "red", linetype = "dashed")+
  #geom_hline(yintercept=0, color = "red", linetype = "dashed")+
  theme_classic() 
plot(plot_Crizo_vs_Camp)


avg_Crizo_GOF_Camp_LOF = Crizo_GOF_Camp_LOF %>% group_by(pos) %>% dplyr::summarise(mean = mean(Crizo, na.rm=TRUE))
x = map_scores_pdb(crizo_2WGJ, avg_Crizo_GOF_Camp_LOF, "mean")
write.pdb(x, file="avg_Crizo_GOF_Camp_LOF.pdb")

avg_Camp_GOF_Crizo_LOF = Camp_GOF_Crizo_LOF %>% group_by(pos) %>% dplyr::summarise(mean = mean(Camp, na.rm=TRUE))
x = map_scores_pdb(cap_docked, avg_Camp_GOF_Crizo_LOF, "mean")
write.pdb(x, file="avg_Camp_GOF_Crizo_LOF.pdb")



# ---------------------- Capmatinib, (ib) Gle (II) ----------------------

Camp_vs_Gle <- data.frame(pos = ex14_rosace_norm_wide$pos,
                          mutation = ex14_rosace_norm_wide$mutation,
                          Camp = ex14_rosace_norm_wide$Camp,
                          Gle = ex14_rosace_norm_wide$Gle)

Gle_ <- ex14_rosace_norm_scores_filtered %>% filter(key %in% c("Gle", "DMSO"))
Gle_resis <- Gle_ %>% group_by(variants) %>% filter((delta.score >= 0 & pval0.5.up <= 0.1 & mean >= 0 & mean[match("DMSO", key)] <= 0)) %>% ungroup()

Camp_ <- ex14_rosace_norm_scores_filtered %>% filter(key %in% c("Camp", "DMSO"))
Camp_resis <- Camp_ %>% group_by(variants) %>% filter((delta.score >= 0 & pval0.5.up <= 0.1 & mean >= 0 & mean[match("DMSO", key)] <= 0)) %>% ungroup()


# parse out resistance mutations and clinical significance 
Camp_GOF_Gle_LOF <- Camp_vs_Gle %>% filter (Camp>= 0.75 & Gle <= 0)
Gle_GOF_Camp_LOF <- Camp_vs_Gle %>% filter (Gle>= 0.75 & Camp <= 0)

Camp_GOF_Gle_LOF <- Camp_vs_Gle  %>% filter (Camp>= 0.75 & Gle <= 0) %>% filter(pos %in% Camp_resis$position)
Gle_GOF_Camp_LOF <- Camp_vs_Gle  %>% filter (Gle>= 0.75 & Camp<= 0) %>% filter(pos %in% Gle_resis$position)
 

# plots 
plot_Camp_vs_Gle <- ggplot() +
  geom_point(data = ex14_rosace_norm_wide , aes(x = Gle, y = Camp),color = "grey60",alpha=0.3) + 
  geom_point(data = Camp_GOF_Gle_LOF, aes(x = Gle, y = Camp), color ="deeppink") + 
  geom_point(data = Gle_GOF_Camp_LOF, aes(x = Gle, y = Camp), color = "deepskyblue") + 
  xlab("Glesatinib") + ylab("Capmatinib") + coord_fixed(ratio = 1) +
  #geom_vline(xintercept=0, color = "red", linetype = "dashed")+
  #geom_hline(yintercept=0, color = "red", linetype = "dashed")+
  theme_classic() 
plot(plot_Camp_vs_Gle)


#Camptinib GOF mapped structurally 
avg_Camp_GOF_Gle_LOF = Camp_GOF_Gle_LOF %>% group_by(pos) %>% dplyr::summarise(mean = mean(Camp, na.rm=TRUE))
x = map_scores_pdb(cap_docked, avg_Camp_GOF_Gle_LOF, "mean")
write.pdb(x, file="avg_Camp_GOF_Gle_LOF.pdb")

# Glestinib GOF mapped structurally
avg_Gle_GOF_Camp_LOF = Gle_GOF_Camp_LOF %>% group_by(pos) %>% dplyr::summarise(mean = mean(Gle, na.rm=TRUE))
x = map_scores_pdb(gle_docked, avg_Gle_GOF_Camp_LOF, "mean")
write.pdb(x, file="avg_Gle_GOF_Camp_LOF.pdb")


ggarrange(plot_Crizo_vs_Cabo,plot_Camp_vs_Gle, ncol = 2, nrow=1)

# ---------------- Crizotinib (1a), merestinib (II)----------------------

Crizo_vs_Mere <- data.frame(hgvs = ex14_rosace_norm_wide $hgvs,
                            pos = ex14_rosace_norm_wide $pos, 
                            Crizo = ex14_rosace_norm_wide $Crizo, 
                            Mere = ex14_rosace_norm_wide $Mere)


# parse out resistance mutations and clinical significance 
Crizo_GOF_Mere_LOF <- Crizo_vs_Mere %>% filter (Crizo>= 0.7 & Mere <= 0 ) %>% mutate(delta = Crizo - Mere )  %>% mutate(clinical = pos %in% clinical_cancer) 
Mere_GOF_Crizo_LOF <- Crizo_vs_Mere%>% filter (Mere>= 0.7 & Crizo <= 0) %>% mutate(delta =  Mere - Crizo )  %>% mutate(clinical = pos %in% clinical_cancer) 


# plots 
plot_Crizo_vs_Mere <- ggplot() +geom_point(data = ex14_rosace_norm_wide , aes(x = Crizo, y = Mere),alpha=0.3) + 
  xlab("Crizotinib") + ylab("Merestinib") + coord_fixed(ratio = 1) +
  geom_vline(xintercept=0, color = "red", linetype = "dashed")+
  geom_hline(yintercept=0, color = "red", linetype = "dashed")+
  theme_classic() 
plot(plot_Crizo_vs_Mere)


#crizotinib GOF mapped structurally 
avg_Crizo_GOF_Mere_LOF = Crizo_GOF_Mere_LOF %>% group_by(pos) %>% dplyr::summarise(mean = mean(Crizo, na.rm=TRUE))
x = map_scores_pdb(crizo_2WGJ, avg_Crizo_GOF_Mere_LOF, "mean")
write.pdb(x, file="avg_Crizo_GOF_Mere_LOF.pdb")

# merestinib GOF mapped structurally
avg_Mere_GOF_Crizo_LOF = Mere_GOF_Crizo_LOF %>% group_by(pos) %>% dplyr::summarise(mean = mean(Mere, na.rm=TRUE))
x = map_scores_pdb(mere_4EEV, avg_Mere_GOF_Crizo_LOF, "mean")
write.pdb(x, file="avg_Mere_GOF_Crizo_LOF.pdb")



# ---------------- Crizotinib (1a), Glesatinib (II)----------------------

Crizo_vs_Gle <- data.frame(hgvs = ex14_rosace_norm_wide $hgvs,
                            pos = ex14_rosace_norm_wide $pos, 
                            Crizo = ex14_rosace_norm_wide $Crizo, 
                            Gle = ex14_rosace_norm_wide $Gle)


# parse out resistance mutations and clinical significance 
Crizo_GOF_Gle_LOF <- Crizo_vs_Gle %>% filter (Crizo>= 0.7 & Gle <= 0 ) %>% mutate(delta = Crizo - Gle )  %>% mutate(clinical = pos %in% clinical_cancer) 
Gle_GOF_Crizo_LOF <- Crizo_vs_Gle%>% filter (Gle>= 0.7 & Crizo <= 0) %>% mutate(delta =  Gle - Crizo )  %>% mutate(clinical = pos %in% clinical_cancer) 


# plots 
plot_Crizo_vs_Gle <- ggplot() +geom_point(data = ex14_rosace_norm_wide , aes(x = Crizo, y = Gle),alpha=0.3) + 
  xlab("Crizotinib") + ylab("Glestinib") + coord_fixed(ratio = 1) +
  geom_vline(xintercept=0, color = "red", linetype = "dashed")+
  geom_hline(yintercept=0, color = "red", linetype = "dashed")+
  theme_classic() 
plot(plot_Crizo_vs_Gle)


#crizotinib GOF mapped structurally 
avg_Crizo_GOF_Gle_LOF = Crizo_GOF_Gle_LOF %>% group_by(pos) %>% dplyr::summarise(mean = mean(Crizo, na.rm=TRUE))
x = map_scores_pdb(crizo_2WGJ, avg_Crizo_GOF_Gle_LOF, "mean")
write.pdb(x, file="avg_Crizo_GOF_Gle_LOF.pdb")

# Glestinib GOF mapped structurally
avg_Gle_GOF_Crizo_LOF = Gle_GOF_Crizo_LOF %>% group_by(pos) %>% dplyr::summarise(mean = mean(Gle, na.rm=TRUE))
x = map_scores_pdb(Gle_4EEV, avg_Gle_GOF_Crizo_LOF, "mean")
write.pdb(x, file="avg_Gle_GOF_Crizo_LOF.pdb")




# ---------------------- Capmatinib, (ib) Cabozantinib (II) ----------------------

### I like this pairng because they are both targeted to delta ex 

Camp_vs_Cabo <- data.frame(hgvs = ex14_rosace_norm_wide $hgvs,
                            pos = ex14_rosace_norm_wide $pos, 
                            Camp = ex14_rosace_norm_wide $Camp, 
                            Cabo = ex14_rosace_norm_wide $Cabo)


# parse out resistance mutations and clinical significance 
Camp_GOF_Cabo_LOF <- Camp_vs_Cabo %>% filter (Camp>= 0.7 & Cabo <= 0 ) %>% mutate(delta = Camp - Cabo )  %>% mutate(clinical = pos %in% clinical_cancer) 
Cabo_GOF_Camp_LOF <- Camp_vs_Cabo%>% filter (Cabo>= 0.7 & Camp <= 0) %>% mutate(delta =  Cabo - Camp )  %>% mutate(clinical = pos %in% clinical_cancer) 


# plots 
plot_Camp_vs_Cabo <- ggplot() +geom_point(data = ex14_rosace_norm_wide , aes(x = Camp, y = Cabo),alpha=0.3) + 
  xlab("Capmatinib") + ylab("Cabozantinib") + coord_fixed(ratio = 1) +
  geom_vline(xintercept=0, color = "red", linetype = "dashed")+
  geom_hline(yintercept=0, color = "red", linetype = "dashed")+
  theme_classic() 
plot(plot_Camp_vs_Cabo)


#Camptinib GOF mapped structurally 
avg_Camp_GOF_Cabo_LOF = Camp_GOF_Cabo_LOF %>% group_by(pos) %>% dplyr::summarise(mean = mean(Camp, na.rm=TRUE))
x = map_scores_pdb(Camp_2WGJ, avg_Camp_GOF_Cabo_LOF, "mean")
write.pdb(x, file="avg_Camp_GOF_Cabo_LOF.pdb")

# Cabostinib GOF mapped structurally
avg_Cabo_GOF_Camp_LOF = Cabo_GOF_Camp_LOF %>% group_by(pos) %>% dplyr::summarise(mean = mean(Cabo, na.rm=TRUE))
x = map_scores_pdb(Cabo_4EEV, avg_Cabo_GOF_Camp_LOF, "mean")
write.pdb(x, file="avg_Cabo_GOF_Camp_LOF.pdb")





# ----------------------Crizotinib, (Ia) Capmatinib (Ib) -----------------------


Crizo_vs_Camp <- data.frame(hgvs = ex14_rosace_scores_wide$hgvs,
                            pos = ex14_rosace_scores_wide$pos, 
                            Crizo = ex14_rosace_scores_wide$Crizo, 
                            Camp = ex14_rosace_scores_wide$Camp)


# parse out resistance mutations and clinical significance 
Crizo_GOF_Camp_LOF <- Crizo_vs_Camp %>% filter (Crizo>= 0.7 & Camp <= 0 ) %>% mutate(delta = Crizo - Camp )  %>% mutate(clinical = pos %in% clinical_cancer) 
Camp_GOF_Crizo_LOF <- Crizo_vs_Camp%>% filter (Camp>= 0.7 & Crizo <= 0) %>% mutate(delta =  Camp - Crizo )  %>% mutate(clinical = pos %in% clinical_cancer) 


# plots 
plot_Crizo_vs_Camp <- ggplot() +geom_point(data = ex14_rosace_scores_wide, aes(x = Crizo, y = Camp),alpha=0.3) + 
  xlab("Crizotinib") + ylab("Capmatinib") + coord_fixed(ratio = 1) +
  geom_vline(xintercept=0, color = "red", linetype = "dashed")+
  geom_hline(yintercept=0, color = "red", linetype = "dashed")+
  theme_classic() 
plot(plot_Crizo_vs_Camp)


#crizotinib GOF mapped structurally 
avg_Crizo_GOF_Camp_LOF = Crizo_GOF_Camp_LOF %>% group_by(pos) %>% dplyr::summarise(mean = mean(Crizo, na.rm=TRUE))
x = map_scores_pdb(crizo_2WGJ, avg_Crizo_GOF_Camp_LOF, "mean")
write.pdb(x, file="avg_Crizo_GOF_Camp_LOF.pdb")

# Campstinib GOF mapped structurally
avg_Camp_GOF_Crizo_LOF = Camp_GOF_Crizo_LOF %>% group_by(pos) %>% dplyr::summarise(mean = mean(Camp, na.rm=TRUE))
x = map_scores_pdb(Camp_4EEV, avg_Camp_GOF_Crizo_LOF, "mean")
write.pdb(x, file="avg_Camp_GOF_Crizo_LOF.pdb")




# ----------------------Crizotinib, (Ia) Tepotinib (Ib) -----------------------


Crizo_vs_Tepo <- data.frame(hgvs = ex14_rosace_scores_wide$hgvs,
                            pos = ex14_rosace_scores_wide$pos, 
                            Crizo = ex14_rosace_scores_wide$Crizo, 
                            Tepo = ex14_rosace_scores_wide$Tepo)


# parse out resistance mutations and clinical significance 
Crizo_GOF_Tepo_LOF <- Crizo_vs_Tepo %>% filter (Crizo>= 0.7 & Tepo <= 0 ) %>% mutate(delta = Crizo - Tepo )  %>% mutate(clinical = pos %in% clinical_cancer) 
Tepo_GOF_Crizo_LOF <- Crizo_vs_Tepo%>% filter (Tepo>= 0.7 & Crizo <= 0) %>% mutate(delta =  Tepo - Crizo )  %>% mutate(clinical = pos %in% clinical_cancer) 


# plots 
plot_Crizo_vs_Tepo <- ggplot() +geom_point(data = ex14_rosace_scores_wide, aes(x = Crizo, y = Tepo),alpha=0.3) + 
  xlab("Crizotinib") + ylab("Tepostinib") + coord_fixed(ratio = 1) +
  geom_vline(xintercept=0, color = "red", linetype = "dashed")+
  geom_hline(yintercept=0, color = "red", linetype = "dashed")+
  theme_classic() 
plot(plot_Crizo_vs_Tepo)


#crizotinib GOF mapped structurally 
avg_Crizo_GOF_Tepo_LOF = Crizo_GOF_Tepo_LOF %>% group_by(pos) %>% dplyr::summarise(mean = mean(Crizo, na.rm=TRUE))
x = map_scores_pdb(crizo_2WGJ, avg_Crizo_GOF_Tepo_LOF, "mean")
write.pdb(x, file="avg_Crizo_GOF_Tepo_LOF.pdb")

# Tepostinib GOF mapped structurally
avg_Tepo_GOF_Crizo_LOF = Tepo_GOF_Crizo_LOF %>% group_by(pos) %>% dplyr::summarise(mean = mean(Tepo, na.rm=TRUE))
x = map_scores_pdb(Tepo_4EEV, avg_Tepo_GOF_Crizo_LOF, "mean")
write.pdb(x, file="avg_Tepo_GOF_Crizo_LOF.pdb")


# ----------------------Crizotinib, (Ia) Savo(Ib) -----------------------


Crizo_vs_Savo <- data.frame(hgvs = ex14_rosace_scores_wide$hgvs,
                            pos = ex14_rosace_scores_wide$pos, 
                            Crizo = ex14_rosace_scores_wide$Crizo, 
                            Savo = ex14_rosace_scores_wide$Savo)


# parse out resistance mutations and clinical significance 
Crizo_GOF_Savo_LOF <- Crizo_vs_Savo %>% filter (Crizo>= 0.7 & Savo <= 0 ) %>% mutate(delta = Crizo - Savo )  %>% mutate(clinical = pos %in% clinical_cancer) 
Savo_GOF_Crizo_LOF <- Crizo_vs_Savo%>% filter (Savo>= 0.7 & Crizo <= 0) %>% mutate(delta =  Savo - Crizo )  %>% mutate(clinical = pos %in% clinical_cancer) 


# plots 
plot_Crizo_vs_Savo <- ggplot() +geom_point(data = ex14_rosace_scores_wide, aes(x = Crizo, y = Savo),alpha=0.3) + 
  xlab("Crizotinib") + ylab("Savostinib") + coord_fixed(ratio = 1) +
  geom_vline(xintercept=0, color = "red", linetype = "dashed")+
  geom_hline(yintercept=0, color = "red", linetype = "dashed")+
  theme_classic() 
plot(plot_Crizo_vs_Savo)


#Camptinib GOF mapped structurally 
avg_Camp_GOF_Savo_LOF = Camp_GOF_Savo_LOF %>% group_by(pos) %>% dplyr::summarise(mean = mean(Camp, na.rm=TRUE))
x = map_scores_pdb(Camp_2WGJ, avg_Camp_GOF_Savo_LOF, "mean")
write.pdb(x, file="avg_Camp_GOF_Savo_LOF.pdb")

# Savostinib GOF mapped structurally
avg_Savo_GOF_Camp_LOF = Savo_GOF_Camp_LOF %>% group_by(pos) %>% dplyr::summarise(mean = mean(Savo, na.rm=TRUE))
x = map_scores_pdb(Savo_4EEV, avg_Savo_GOF_Camp_LOF, "mean")
write.pdb(x, file="avg_Tepo_GOF_Camp_LOF.pdb")






# ---------------------- Merestinib, Glesatinib (II)----------------------

Mere_vs_Gle <- data.frame(hgvs = ex14_rosace_scores_wide$hgvs,
                          pos = ex14_rosace_scores_wide$pos, 
                          Mere = ex14_rosace_scores_wide$Mere,
                          Gle = ex14_rosace_scores_wide$Gle )

Mere_GOF_Gle_LOF <- Mere_vs_Gle %>% filter (Mere>= 0.7 & Gle <= 0) %>% mutate(delta = Mere - Gle)  %>% mutate(clinical = pos %in% clinical_cancer) 
Gle_GOF_Mere_LOF <- Mere_vs_Gle %>% filter (Gle>= 0.7 & Mere <= 0)%>% mutate(delta = Gle - Mere)  %>% mutate(clinical = pos %in% clinical_cancer) 


plot_Mere_vs_Gle <- ggplot() +geom_point(data = ex14_rosace_scores_wide, aes(x = Mere, y = Gle),alpha=0.3) + 
  xlab("Merestinib") + ylab("Glesatinib") + coord_fixed(ratio = 1) +
  geom_vline(xintercept=0, color = "red", linetype = "dashed")+
  geom_hline(yintercept=0, color = "red", linetype = "dashed")+
  theme_classic() 
plot(plot_Mere_vs_Gle)


ex14_2WGJ_crizo <- read.pdb("2WGJ")
ex14_4EEV_crizo <- read.pdb("4EEV")

avg_Crizo_GOF_Cabo_LOF = Crizo_GOF_Cabo_LOF %>% group_by(pos) %>% dplyr::summarise(mean = mean(Crizo, na.rm=TRUE))
x = map_scores_pdb(ex14_2WGJ_crizo, avg_Crizo_GOF_Cabo_LOF, "mean")
write.pdb(x, file="avg_Crizo_GOF_cabo_LOF.pdb")

avg_Cabo_GOF_Crizo_LOF = Cabo_GOF_Crizo_LOF %>% group_by(pos) %>% dplyr::summarise(mean = mean(Cabo, na.rm=TRUE))
x = map_scores_pdb(ex14_4EEV_crizo, avg_Cabo_GOF_Crizo_LOF, "mean")
write.pdb(x, file="avg_Cabo_GOF_Crizo_LOF.pdb")


###------------------------------------- All inhibitors vs all inhibitors  ---------------------------------------------

# the goal of this code is to see what inhibitor pairs have the largest mutational similarities and differences
# to each plot a pearsons score can be calculated 
# variance can be measured horizontally from the line of best fit to find the strongest outliers

# this plots every condition against all conditions and the correlations for each as well 
cc = ex14_rosace_scores_wide[-(1:4)]
#ggpairs(cc)


# Calculate Pearson correlation coefficients for your data frame cc
cor_matrix <- cor(cc, use = "complete.obs")

# Define the custom order for Var1 and the inverse order for Var2
desired_order_var1 <- c("DMSO", "Tiv", "Camp", "Crizo", "Tepo", "Glu", "Savo", "NVP", "Cabo", "Mere", "Gle", "A458")
desired_order_var2 <- rev(desired_order_var1)

# Create new factor variables to specify the order
cor_matrix <- cor_matrix[desired_order_var1, desired_order_var2]

# Define group breaks
break_indices_var1 <- c(1, 2, 8, 11)
break_indices_var2 <- c(1, 3, 10)

# Convert the correlation matrix to long format
cor_matrix_long <- as.data.frame(as.table(cor_matrix))

# Rename columns for clarity
colnames(cor_matrix_long) <- c("Var1", "Var2", "Correlation")

# Create a scatter plot of all columns against one another
library(ggplot2)

# Function to add white lines as separators between groups
add_group_lines <- function(p, breaks_var1, breaks_var2) {
  for (b in breaks_var1) {
    p <- p + geom_vline(xintercept = b + 0.5, color = "white", size = 2)
    p <- p + geom_hline(yintercept = b + 0.5, color = "white", size = 2)
  }
  for (b in breaks_var2) {
    p <- p + geom_vline(xintercept = b + 0.5, color = "white", size = 2)
    p <- p + geom_hline(yintercept = b + 0.5, color = "white", size = 2)
  }
  return(p)
}

# Create the scatter plot
scatter_plot <- ggplot(cor_matrix_long, aes(Var1, Var2, fill = Correlation)) +
  geom_tile() +
  scale_fill_gradient2(low = "grey", mid = "white", high = "#01b6f0", midpoint = 0.7) +
  labs(title = "Met Pearson Correlation Heatmap") +
  geom_text(aes(label = round(Correlation, 2)), vjust = 1) +  # Display the correlation values +
  scale_x_discrete(labels = function(x) ifelse(as.numeric(x) %in% break_indices_var1, "", x)) +
  scale_y_discrete(labels = function(x) ifelse(as.numeric(x) %in% break_indices_var2, "", x)) +
  theme_minimal()  # Remove the background

# Add white lines as separators
scatter_plot <- add_group_lines(scatter_plot, break_indices_var1, break_indices_var2)

# Print the scatter plot
print(scatter_plot)


#-------- Met del Ex14 -------# 

# this plots every condition against all conditions and the correlations for each as well 
cc = met_rosace_scores_wide[-(1:4)]
#ggpairs(cc)


# Calculate Pearson correlation coefficients for your data frame cc
cor_matrix <- cor(cc, use = "complete.obs")

# Define the custom order for Var1 and the inverse order for Var2
desired_order_var1 <- c("DMSO", "Tiv", "Camp", "Crizo", "Tepo", "Glu", "Savo", "NVP", "Cabo", "Mere", "Gle", "A458")
desired_order_var2 <- rev(desired_order_var1)

# Create new factor variables to specify the order
cor_matrix <- cor_matrix[desired_order_var1, desired_order_var2]

# Define group breaks
break_indices_var1 <- c(1, 2, 8, 11)
break_indices_var2 <- c(1, 3, 10)

# Convert the correlation matrix to long format
cor_matrix_long <- as.data.frame(as.table(cor_matrix))

# Rename columns for clarity
colnames(cor_matrix_long) <- c("Var1", "Var2", "Correlation")

# Create a scatter plot of all columns against one another
library(ggplot2)

# Function to add white lines as separators between groups
add_group_lines <- function(p, breaks_var1, breaks_var2) {
  for (b in breaks_var1) {
    p <- p + geom_vline(xintercept = b + 0.5, color = "white", size = 2)
    p <- p + geom_hline(yintercept = b + 0.5, color = "white", size = 2)
  }
  for (b in breaks_var2) {
    p <- p + geom_vline(xintercept = b + 0.5, color = "white", size = 2)
    p <- p + geom_hline(yintercept = b + 0.5, color = "white", size = 2)
  }
  return(p)
}

# Create the scatter plot
scatter_plot <- ggplot(cor_matrix_long, aes(Var1, Var2, fill = Correlation)) +
  geom_tile() +
  scale_fill_gradient2(low = "grey", mid = "white", high = "#01b6f0", midpoint = 0.85) +
  labs(title = "Met del Ex14 Pearson Correlation Heatmap") +
  geom_text(aes(label = round(Correlation, 2)), vjust = 1) +  # Display the correlation values +
  scale_x_discrete(labels = function(x) ifelse(as.numeric(x) %in% break_indices_var1, "", x)) +
  scale_y_discrete(labels = function(x) ifelse(as.numeric(x) %in% break_indices_var2, "", x)) +
  theme_minimal()  # Remove the background

# Add white lines as separators
scatter_plot <- add_group_lines(scatter_plot, break_indices_var1, break_indices_var2)

# Print the scatter plot
print(scatter_plot)


###------------------------------------- Type 1 vs Type 1---------------------------------------------

type1 <-  subset(ex14_rosace_scores_wide , select=c("Crizo","Camp","Tepo","Glu","NVP","Savo"))
#ggpairs(type1)


###------------------------------------- Type 2 vs Type 2---------------------------------------------

type2 <-  subset(ex14_rosace_scores_wide , select=c("Cabo","Mere","Gle"))
ggpairs(type2)

###------------------------------------- Type 2 vs Type 1---------------------------------------------

# crizotinib
Crizo_vs_Cabo <- ggplot() +geom_point(data = ex14_rosace_norm_wide , aes(x = Cabo, y = Crizo), size = 0.5, alpha=0.2)+theme_classic()
Crizo_vs_Gle <- ggplot() +geom_point(data = ex14_rosace_norm_wide , aes(x = Gle, y = Crizo), size = 0.5, alpha=0.2)+theme_classic()
Crizo_vs_Mere <- ggplot() +geom_point(data = ex14_rosace_norm_wide , aes(x = Mere, y = Crizo), size = 0.5, alpha=0.2)+theme_classic()

# capmatinib
Camp_vs_Cabo <- ggplot() +geom_point(data = ex14_rosace_norm_wide , aes(x = Cabo, y = Camp), size = 0.5, alpha=0.2)+theme_classic()
Camp_vs_Gle <- ggplot() +geom_point(data = ex14_rosace_norm_wide , aes(x = Gle, y = Camp), size = 0.5, alpha=0.2)+theme_classic()
Camp_vs_Mere <- ggplot() +geom_point(data = ex14_rosace_norm_wide , aes(x = Mere, y = Camp),size = 0.5,  alpha=0.2)+theme_classic()

# tepotinib
Tepo_vs_Cabo <- ggplot() +geom_point(data = ex14_rosace_norm_wide , aes(x = Cabo, y = Tepo),size = 0.5,  alpha=0.2)+theme_classic()
Tepo_vs_Gle <- ggplot() +geom_point(data = ex14_rosace_norm_wide , aes(x = Gle, y = Tepo), size = 0.5, alpha=0.2)+theme_classic()
Tepo_vs_Mere <- ggplot() +geom_point(data = ex14_rosace_norm_wide , aes(x = Mere, y = Tepo),size = 0.5,  alpha=0.2)+theme_classic()

# glumetinib
Glu_vs_Cabo <- ggplot() +geom_point(data = ex14_rosace_norm_wide , aes(x = Cabo, y = Glu),size = 0.5,  alpha=0.2)+theme_classic()
Glu_vs_Gle <- ggplot() +geom_point(data = ex14_rosace_norm_wide , aes(x = Gle, y = Glu), size = 0.5, alpha=0.2)+theme_classic()
Glu_vs_Mere <- ggplot() +geom_point(data = ex14_rosace_norm_wide , aes(x = Mere, y = Glu),size = 0.5,  alpha=0.2)+theme_classic()

# savolitinib
Savo_vs_Cabo <- ggplot() +geom_point(data = ex14_rosace_norm_wide , aes(x = Cabo, y = Savo), size = 0.5, alpha=0.2)+theme_classic()
Savo_vs_Gle <- ggplot() +geom_point(data = ex14_rosace_norm_wide , aes(x = Gle, y = Savo), size = 0.5, alpha=0.2)+theme_classic()
Savo_vs_Mere <- ggplot() +geom_point(data = ex14_rosace_norm_wide , aes(x = Mere, y = Savo), size = 0.5, alpha=0.2)+theme_classic()

# NVP
NVP_vs_Cabo <- ggplot() +geom_point(data = ex14_rosace_norm_wide , aes(x = Cabo, y = NVP), size = 0.5, alpha=0.2)+theme_classic()
NVP_vs_Gle <- ggplot() +geom_point(data = ex14_rosace_norm_wide , aes(x = Gle, y = NVP), size = 0.5, alpha=0.2)+theme_classic()
NVP_vs_Mere <- ggplot() +geom_point(data = ex14_rosace_norm_wide , aes(x = Mere, y = NVP), size = 0.5, alpha=0.1)+theme_classic()

ggarrange(Crizo_vs_Cabo, Crizo_vs_Gle, Crizo_vs_Mere,
          Camp_vs_Cabo, Camp_vs_Gle, Camp_vs_Mere,
          Tepo_vs_Cabo, Tepo_vs_Gle, Tepo_vs_Mere,
          Glu_vs_Cabo, Glu_vs_Gle, Glu_vs_Mere,
          Savo_vs_Cabo, Savo_vs_Gle, Savo_vs_Mere,
          NVP_vs_Cabo, NVP_vs_Gle, NVP_vs_Mere,
          nrow=6, ncol=3)


ggarrange(Crizo_vs_Cabo, Camp_vs_Cabo, Tepo_vs_Cabo, Glu_vs_Cabo, Savo_vs_Cabo, NVP_vs_Cabo,
          Crizo_vs_Gle, Camp_vs_Gle, Tepo_vs_Gle, Glu_vs_Gle,Savo_vs_Gle,NVP_vs_Gle,
          Crizo_vs_Mere, Camp_vs_Mere, Tepo_vs_Mere, Glu_vs_Mere, Savo_vs_Mere, NVP_vs_Mere,
          nrow=3, ncol=6)

plot(Crizo_vs_Cabo)


#####------------------- Correlations of DMSO subtracted scores for type I vs type II --------------------#####


numeric_cols <- c("Gle", "Mere", "Cabo", "Camp","NVP", "Savo", "Glu", "Tepo", "Crizo")
cor_matrix <- cor(ex14_rosace_norm_wide[, numeric_cols, drop = FALSE])

# Replace NA with 0 and Inf with a finite value
cor_matrix[is.na(cor_matrix)] <- 0
cor_matrix[is.infinite(cor_matrix)] <- 1e10  # You can use any large finite value

# Load ggplot2 library
library(ggplot2)

# Convert the correlation matrix to a long format for ggplot
library(tidyr)
cor_df <- as.data.frame(as.table(cor_matrix))
colnames(cor_df) <- c("Var1", "Var2", "Correlation")

# Set factor levels to reverse the order of x-axis
cor_df$Var1 <- factor(cor_df$Var1, levels = rev(levels(cor_df$Var1)))

# Plotting the correlation matrix with ggplot2
ggplot(cor_df, aes(Var1, Var2, fill = Correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Correlation, 2)), vjust = 1) +  # Add text annotations
  scale_fill_gradient2(low = "grey", mid = "white", high = "#01b6f0", midpoint = 0.8) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  labs(title = "Correlation Matrix Heatmap")

library(reshape2)
ggheatmap <- ggplot(cor_df, aes(Var2, Var1, fill = Correlation))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "grey", mid = "white", high = "#01b6f0", midpoint = 0.8) +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
plot(ggheatmap)


ggheatmap + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))


