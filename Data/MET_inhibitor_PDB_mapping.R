###################################################################################################################
# MET_inhibitor_PDB_mapping.R
# by g.estevam @ ucsf
# code to generate alll the inhibitor heatmaps from the normalized roscae scores 
###################################################################################################################

library(readr)
library(stringr)
library(tidyverse)
library(dplyr)
library(ggridges)
library(tidyquant)
library(bio3d)
source("dms_analysis_utilities.R")

###################################################################################################################


##-------------------------------------- Load Rosace Fitness scores ----------------------------------------------##
### open data file as data frame for all raw inhibitor ALL scores generated from ROSACE

# delta exon14 
met_rosace_norm_scores_filtered <- as.data.frame(fread("met_scores_filtered.tsv")) # here met means MET delta Ex14
met_rosace_norm_scores_filtered  <- met_rosace_norm_scores_filtered  %>% mutate(pos = position + 1058)

met_rosace_norm_scores <- as.data.frame(fread("met_scores_unfiltered.tsv")) # here met means MET delta Ex14
met_rosace_norm_scores  <- met_rosace_norm_scores   %>% mutate(pos= position + 1058)

# wt 
ex14_rosace_norm_scores_filtered <- as.data.frame(fread("ex14_scores_filtered.tsv")) # here ex14 means MET+Ex14 
ex14_rosace_norm_scores_filtered  <- ex14_rosace_norm_scores_filtered  %>% mutate(pos = position + 1058)

ex14_rosace_norm_score <- as.data.frame(fread("ex14_scores_unfiltered.tsv")) # here ex14 means MET+Ex14 
ex14_rosace_norm_scores  <- ex14_rosace_norm_scores   %>% mutate(pos = position + 1058)



#----------- WT MET+Exon14 --------------#

ex14_rosace_scores_condensed <- data.frame (hgvs = ex14_rosace_norm_scores_filtered$variant,
                                            pos =ex14_rosace_norm_scores_filtered$position, 
                                            mutation = ex14_rosace_norm_scores_filtered$mutation, 
                                            score = ex14_rosace_norm_scores_filtered$mean,
                                            type = ex14_rosace_norm_scores_filtered$type,
                                            inhib= ex14_rosace_norm_scores_filtered$key, 
                                            pval1.5= ex14_rosace_norm_scores_filtered$pval1.5,
                                            pval0=ex14_rosace_norm_scores_filtered$pval0)

ex14_rosace_scores_wide <- ex14_rosace_scores_condensed %>% pivot_wider(names_from = inhib, values_from = score) 

ex14_rosace_scores_crizo <- ex14_rosace_scores_condensed  %>% filter(inhib == "Crizo")
ex14_rosace_scores_tepo <- ex14_rosace_scores_condensed  %>% filter(inhib == "Tepo")
ex14_rosace_scores_tiv <- ex14_rosace_scores_condensed  %>% filter(inhib == "Tiv")
ex14_rosace_scores_amg458 <- ex14_rosace_scores_condensed  %>% filter(inhib == "A458")
ex14_rosace_scores_nvp <- ex14_rosace_scores_condensed  %>% filter(inhib == "NVP")
ex14_rosace_scores_mere <- ex14_rosace_scores_condensed  %>% filter(inhib == "Mere")
ex14_rosace_scores_savo <- ex14_rosace_scores_condensed  %>% filter(inhib == "Savo")
ex14_rosace_scores_cap <- ex14_rosace_scores_condensed  %>% filter(inhib == "Camp")
ex14_rosace_scores_cabo <- ex14_rosace_scores_condensed  %>% filter(inhib == "Cabo")
ex14_rosace_scores_glu <- ex14_rosace_scores_condensed  %>% filter(inhib == "Glu")
ex14_rosace_scores_gle <- ex14_rosace_scores_condensed  %>% filter(inhib == "Gle")


#----------- MET del Exon14 --------------#

met_rosace_scores_condensed <- data.frame (hgvs = met_rosace_scores$variant,
                                            pos =met_rosace_scores$position, 
                                            mutation = met_rosace_scores$mutation, 
                                            score = met_rosace_scores$mean,
                                            type = met_rosace_scores$type,
                                            inhib= met_rosace_scores$key)

met_rosace_scores_wide <- met_rosace_scores_condensed %>% pivot_wider(names_from = inhib, values_from = score) 

met_rosace_scores_crizo <- met_rosace_scores_condensed  %>% filter(inhib == "Crizo")
met_rosace_scores_tepo <- met_rosace_scores_condensed  %>% filter(inhib == "Tepo")
met_rosace_scores_tiv <- met_rosace_scores_condensed  %>% filter(inhib == "Tiv")
met_rosace_scores_amg458 <- met_rosace_scores_condensed  %>% filter(inhib == "A458")
met_rosace_scores_nvp <- met_rosace_scores_condensed  %>% filter(inhib == "NVP")
met_rosace_scores_mere <- met_rosace_scores_condensed  %>% filter(inhib == "Mere")
met_rosace_scores_savo <- met_rosace_scores_condensed  %>% filter(inhib == "Savo")
met_rosace_scores_cap <- met_rosace_scores_condensed  %>% filter(inhib == "Camp")
met_rosace_scores_cabo <- met_rosace_scores_condensed  %>% filter(inhib == "Cabo")
met_rosace_scores_glu <- met_rosace_scores_condensed  %>% filter(inhib == "Glu")
met_rosace_scores_gle <- met_rosace_scores_condensed  %>% filter(inhib == "Gle")



##-------------------------------------- PDBs for mapping  ----------------------------------------------##

# experimetally solved 
# PDB 2WGJ, Crizotinib 
# PDB 4R1V, Tepotinib
# PDB 3RHK, Tivantinib
# PDB 5T3Q, AMG-458
# PDB 3QTI, NVP-BVU972
# PDB 4EEV, Merestinib
# PDB 6SDE, Savolitinib 

# docked  
# cabozantinib 
# capmatinib 
# glumetinib 
# glesatinib

# read pdb files 
crizo_2WGJ <- read.pdb("2WGJ")
tepo_4R1V <- read.pdb("4R1V")
tiv_3RHK <- read.pdb("3RHK")
amg458_5T3Q <- read.pdb("5T3Q")
nvp_3QTI <- read.pdb("3QTI")
mere_4EEV <- read.pdb("4EEV")
savo_6SDE <- read.pdb("6SDE")


######################################################### WT MET+Exon14 ##############################################################

##-------------------------------------- average scores for each condition, position  ----------------------------------------------##

## WT MET+Exon14 ###

# type 1 
avg_ex14_crizo = ex14_rosace_scores_crizo%>% group_by(pos) %>% dplyr::summarise(mean = mean(score, na.rm=TRUE))
avg_ex14_tepo = ex14_rosace_scores_tepo %>% group_by(pos) %>% dplyr::summarise(mean = mean(score, na.rm=TRUE))
avg_ex14_tiv = ex14_rosace_scores_tiv %>% group_by(pos) %>% dplyr::summarise(mean = mean(score, na.rm=TRUE))
avg_ex14_nvp = ex14_rosace_scores_nvp %>% group_by(pos) %>% dplyr::summarise(mean = mean(score, na.rm=TRUE))
avg_ex14_savo = ex14_rosace_scores_savo %>% group_by(pos) %>% dplyr::summarise(mean = mean(score, na.rm=TRUE))
#avg_ex14_glu = ex14_rosace_scores_wide %>% group_by(pos) %>% dplyr::summarise(mean = mean(ROSACE_effects, na.rm=TRUE))
#avg_ex14_cap = ex14_rosace_scores_wide %>% group_by(pos) %>% dplyr::summarise(mean = mean(ROSACE_effects, na.rm=TRUE))

# type 2 
avg_ex14_mere = ex14_rosace_scores_mere %>% group_by(pos) %>% dplyr::summarise(mean = mean(score, na.rm=TRUE))
#avg_ex14_gle = ex14_rosace_scores_wide %>% group_by(pos) %>% dplyr::summarise(mean = mean(ROSACE_effects, na.rm=TRUE))
#avg_ex14_cabo = ex14_rosace_scores_wide %>% group_by(pos) %>% dplyr::summarise(mean = mean(ROSACE_effects, na.rm=TRUE))

# type 1.5
avg_ex14_amg458 = ex14_rosace_scores_amg458 %>% group_by(pos) %>% dplyr::summarise(mean = mean(score, na.rm=TRUE))



##-------------------------------------- map on PDB file  ----------------------------------------------##
a = map_scores_pdb(crizo_2WGJ, avg_ex14_crizo, "mean")
b = map_scores_pdb(tepo_4R1V, avg_ex14_tepo, "mean")
c = map_scores_pdb(tiv_3RHK, avg_ex14_tiv , "mean")
d = map_scores_pdb(amg458_5T3Q, avg_ex14_amg458, "mean")
e = map_scores_pdb(nvp_3QTI, avg_ex14_nvp, "mean")
f = map_scores_pdb(mere_4EEV, avg_ex14_mere , "mean")
g = map_scores_pdb(savo_6SDE ,avg_ex14_savo, "mean")
#h = map_scores_pdb( , avg_ex14_glu, "mean")
#i = map_scores_pdb( , avg_ex14_cap , "mean")
#j = map_scores_pdb( ,avg_ex14_gle, "mean")
#k = map_scores_pdb(, avg_ex14_cabo, "mean")



write.pdb(a, file="ex14_rosace_scores_crizo.pdb")
write.pdb(b, file="ex14_rosace_scores_tepo.pdb")
write.pdb(c, file="ex14_rosace_scores_tiv.pdb")
write.pdb(d, file="ex14_rosace_scores_amg458.pdb")
write.pdb(e, file="ex14_rosace_scores_nvp.pdb")
write.pdb(f, file="ex14_rosace_scores_mere.pdb")
write.pdb(g, file="ex14_rosace_scores_savo.pdb")





########################################################### MET del Exon14 ############################################################

##-------------------------------------- average scores for each condition, position  ----------------------------------------------##

## WT MET del Exon14 ###

# type 1 
avg_met_crizo = met_rosace_scores_crizo%>% group_by(pos) %>% dplyr::summarise(mean = mean(score, na.rm=TRUE))
avg_met_tepo = met_rosace_scores_tepo %>% group_by(pos) %>% dplyr::summarise(mean = mean(score, na.rm=TRUE))
avg_met_tiv = met_rosace_scores_tiv %>% group_by(pos) %>% dplyr::summarise(mean = mean(score, na.rm=TRUE))
avg_met_nvp = met_rosace_scores_nvp %>% group_by(pos) %>% dplyr::summarise(mean = mean(score, na.rm=TRUE))
avg_met_savo = met_rosace_scores_savo %>% group_by(pos) %>% dplyr::summarise(mean = mean(score, na.rm=TRUE))
#avg_met_glu = met_rosace_scores_wide %>% group_by(pos) %>% dplyr::summarise(mean = mean(ROSACE_effects, na.rm=TRUE))
#avg_met_cap = met_rosace_scores_wide %>% group_by(pos) %>% dplyr::summarise(mean = mean(ROSACE_effects, na.rm=TRUE))

# type 2 
avg_met_mere = met_rosace_scores_mere %>% group_by(pos) %>% dplyr::summarise(mean = mean(score, na.rm=TRUE))
#avg_met_gle = met_rosace_scores_wide %>% group_by(pos) %>% dplyr::summarise(mean = mean(ROSACE_effects, na.rm=TRUE))
#avg_met_cabo = met_rosace_scores_wide %>% group_by(pos) %>% dplyr::summarise(mean = mean(ROSACE_effects, na.rm=TRUE))

# type 1.5
avg_met_amg458 = met_rosace_scores_amg458 %>% group_by(pos) %>% dplyr::summarise(mean = mean(score, na.rm=TRUE))



##-------------------------------------- map on PDB file  ----------------------------------------------##

l = map_scores_pdb(crizo_2WGJ, avg_met_crizo, "mean")
m = map_scores_pdb(tepo_4R1V, avg_met_tepo, "mean")
n = map_scores_pdb(tiv_3RHK, avg_met_tiv , "mean")
o = map_scores_pdb(amg458_5T3Q, avg_met_amg458, "mean")
p = map_scores_pdb(nvp_3QTI, avg_met_nvp, "mean")
q = map_scores_pdb(mere_4EEV, avg_met_mere , "mean")
r = map_scores_pdb(savo_6SDE ,avg_met_savo, "mean")
#s = map_scores_pdb( , avg_met_glu, "mean")
#t = map_scores_pdb( , avg_met_cap , "mean")
#u = map_scores_pdb( ,avg_met_gle, "mean")
#v = map_scores_pdb(, avg_met_cabo, "mean")


write.pdb(l, file="met_rosace_scores_crizo.pdb")
write.pdb(m, file="met_rosace_scores_tepo.pdb")
write.pdb(n, file="met_rosace_scores_tiv.pdb")
write.pdb(o, file="met_rosace_scores_amg458.pdb")
write.pdb(p, file="met_rosace_scores_nvp.pdb")
write.pdb(q, file="met_rosace_scores_mere.pdb")
write.pdb(r, file="met_rosace_scores_savo.pdb")



