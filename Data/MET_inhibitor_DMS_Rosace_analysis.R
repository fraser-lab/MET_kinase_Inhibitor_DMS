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
#library(ggbraid)
library(tidyquant)
library(bio3d)
source("dms_analysis_utilities.R")
library(scales)
library(colorspace)
library(rstatix)
library(patchwork)
library(cowplot)
library(dendextend)
library("RColorBrewer")
library(gghighlight)
library(stringdist)
library(forcats)
library(janitor)
library(pheatmap)
library(data.table)
library(factoextra)
library(FactoMineR)

# color palette
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#0072B2","#F0E442","#D55E00", "#CC79A7")

#####################################################################################################
###------------------------------------ define initial data frames ----------------------------------
#####################################################################################################


##-------------------------------------- Rosace scores ----------------------------------------------
### open data file as data frame for all raw inhibitor ALL scores generated from ROSACE

met_rosace_scores <- as.data.frame(fread("WT_rosace_effect_all.tsv")) # here met is delta exon14
met_rosace_scores <- met_rosace_scores  %>% mutate(position = position + 1058) 
met_rosace_scores <- met_rosace_scores  %>% mutate(position = position + 1058)
met_rosace_scores <- met_rosace_scores %>% mutate(pos = position) 

ex14_rosace_scores <- as.data.frame(fread("WT_rosace_effect_all.tsv")) # here ex14 means MET+Ex14 
ex14_rosace_scores <- ex14_rosace_scores  %>% mutate(position = position + 1058)
ex14_rosace_scores <- ex14_rosace_scores %>% mutate(pos = position) 

met_rosace_norm_score <- as.data.frame(fread("met_scores_filtered.tsv")) # here met means MET+met 
met_rosace_norm_score  <- met_rosace_norm_score   %>% mutate(position = position + 1058)
met_rosace_norm_score  <- met_rosace_norm_score  %>% mutate(pos = position) 

ex14_rosace_norm_score <- as.data.frame(fread("ex14_scores_filtered.tsv")) # here ex14 means MET+Ex14 
ex14_rosace_norm_score  <- ex14_rosace_norm_score   %>% mutate(position = position + 1058)
ex14_rosace_norm_score  <- ex14_rosace_norm_score  %>% mutate(pos = position) 


#####################################################################################################
###---------------------------------------- Prelim analysis -----------------------------------------
#####################################################################################################

# type 1 : crizo, cap, tepo , tiv , glu, savo, nvp 
# type 2 : cabo, mere, gle 
# type 1.5: AMG 

# MET WT inhibitor, DMSO correlation plots 

#--------------------------------- Correlation Plots --------------------------------
ex14_rosace_scores_condensed <- data.frame (hgvs = ex14_rosace_scores$variant,
                                            pos =ex14_rosace_scores$position, 
                                            mutation = ex14_rosace_scores$mutation, 
                                            score = ex14_rosace_scores$ROSACE_effects,
                                            type = ex14_rosace_scores$type,
                                            inhib= ex14_rosace_scores$inhibitor
                                            )

# not normalized
ex14_rosace_scores_wide <- ex14_rosace_scores_condensed %>%
  pivot_wider(names_from = inhib, values_from = score) 

#normalizeds
ex14_rosace_GOF_wide <- ex14_rosace_norm_score %>%
  pivot_wider(names_from = key, values_from = mean) 


plot1 <- ggplot() +geom_point(data = ex14_rosace_scores_wide, aes(x = DMSO, y = Crizo), color = "#97c7b8",alpha=0.3) + 
  xlab("DMSO") + ylab("Crizo") + coord_fixed(ratio = 1) +theme_classic() 
#plot(plot1)

plot2 <- ggplot() +geom_point(data = ex14_rosace_scores_wide, aes(x = DMSO, y = Camp), color ="#97c7b8",alpha=0.3) + 
  xlab("DMSO") + ylab("Cap") + coord_fixed(ratio = 1) +theme_classic() 
#plot(plot2) 

plot3 <- ggplot() +geom_point(data = ex14_rosace_scores_wide, aes(x = DMSO, y = Tepo), color ="#97c7b8",alpha=0.3) + 
  xlab("DMSO") + ylab("Tepo") + coord_fixed(ratio = 1) +theme_classic() 
#plot(plot3)

plot4 <- ggplot() +geom_point(data = ex14_rosace_scores_wide, aes(x = DMSO, y = Tiv), color ="#f47c54",alpha=0.3) + 
  xlab("DMSO") + ylab("Tiv") + coord_fixed(ratio = 1) +theme_classic() 
#plot(plot4)

plot5 <- ggplot() +geom_point(data = ex14_rosace_scores_wide, aes(x = DMSO, y = Glu), color ="#97c7b8",alpha=0.3) + 
  xlab("DMSO") + ylab("Glu") + coord_fixed(ratio = 1) +theme_classic() 
#plot(plot5)

plot6 <- ggplot() +geom_point(data = ex14_rosace_scores_wide, aes(x = DMSO, y = Savo),  color ="#97c7b8",alpha=0.3) + 
  xlab("DMSO") + ylab("Savo") + coord_fixed(ratio = 1) +theme_classic() 
#plot(plot6)

plot7 <- ggplot() +geom_point(data = ex14_rosace_scores_wide, aes(x = DMSO, y = NVP),  color ="#97c7b8",alpha=0.3) + 
  xlab("DMSO") + ylab("NVP") + coord_fixed(ratio = 1) +theme_classic() 
#plot(plot7)

plot8 <- ggplot() +geom_point(data = ex14_rosace_scores_wide, aes(x = DMSO, y = Cabo), color ="#0a91a2",alpha=0.3) + 
  xlab("DMSO") + ylab("Cabo") + coord_fixed(ratio = 1) +theme_classic() 
#plot(plot8)

plot9 <- ggplot() +geom_point(data = ex14_rosace_scores_wide, aes(x = DMSO, y = Mere), color ="#0a91a2",alpha=0.3) + 
  xlab("DMSO") + ylab("Mere") + coord_fixed(ratio = 1) +theme_classic() 
#plot(plot9)

plot10 <- ggplot() +geom_point(data = ex14_rosace_scores_wide, aes(x = DMSO, y = Gle), color ="#0a91a2",alpha=0.3) + 
  xlab("DMSO") + ylab("Gle") + coord_fixed(ratio = 1) +theme_classic() 
#plot(plot10)

plot11 <- ggplot() +geom_point(data = ex14_rosace_scores_wide, aes(x = DMSO, y = A458), color ="#f9a953",alpha=0.3) + 
  xlab("DMSO") + ylab("A458") + coord_fixed(ratio = 1) +theme_classic() 
#plot(plot11)

scatterplots = ggarrange(plot1, plot2, plot3, plot5, plot6, plot7, 
                         plot8, plot9, plot10, plot11, plot4,
                         align = "hv",nrow = 2, ncol = 6)
plot(scatterplots)





# Load required libraries
library(ggplot2)
library(dplyr)
library(ggpubr)

# Define your scatter plots

# Function to create scatter plots
create_scatterplot <- function(x_var, y_var, x_label, y_label, color) {
  p <- ggplot(data = ex14_rosace_scores_wide, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_point(color = color, alpha = 0.3) +
    geom_smooth(method = "lm", se = FALSE, color = color) +  # Add regression line
    xlab(x_label) +
    ylab(y_label) +
    coord_fixed(ratio = 1) +
    theme_classic()
  return(p)
}

plot1 <- create_scatterplot("DMSO", "Crizo", "DMSO", "Crizo", "#97c7b8")
plot2 <- create_scatterplot("DMSO", "Camp", "DMSO", "Camp", "#97c7b8")
plot3 <- create_scatterplot("DMSO", "Tepo", "DMSO", "Tepo", "#97c7b8")
plot4 <- create_scatterplot("DMSO", "Tiv", "DMSO", "Tiv", "#f47c54")
plot5 <- create_scatterplot("DMSO", "Glu", "DMSO", "Glu", "#97c7b8")
plot6 <- create_scatterplot("DMSO", "Savo", "DMSO", "Savo", "#97c7b8")
plot7 <- create_scatterplot("DMSO", "NVP", "DMSO", "NVP", "#97c7b8")
plot8 <- create_scatterplot("DMSO", "Cabo", "DMSO", "Cabo", "#0a91a2")
plot9 <- create_scatterplot("DMSO", "Mere", "DMSO", "Mere", "#0a91a2")
plot10 <- create_scatterplot("DMSO", "Gle", "DMSO", "Gle", "#0a91a2")
plot11 <- create_scatterplot("DMSO", "A458", "DMSO", "A458", "#f9a953")

# Calculate Pearson correlation coefficients
cor_coef <- sapply(c("Crizo", "Camp", "Tepo", "Tiv", "Glu", "Savo", "NVP", "Cabo", "Mere", "Gle", "A458"), function(var) {
  cor(ex14_rosace_scores_wide$DMSO, ex14_rosace_scores_wide[[var]], method = "pearson")
})

# Add correlation coefficients to the plots
labels <- lapply(cor_coef, function(cor) {
  as.character(round(cor, 3))
})
plot_labels <- lapply(1:11, function(i) {
  geom_text(data = NULL, aes(x = Inf, y = -Inf, label = labels[i]))
})

# Arrange the scatter plots with labels
scatterplots = ggarrange(plot1, plot2, plot3, plot8, plot9, plot10, 
                         plot5, plot6, plot7, plot11, plot4,
                         labels = plot_labels,
                         align = "hv", nrow = 2, ncol = 6)

plot(scatterplots)


#--------------------------------- Distributions--------------------------------

# not normalized MET+Ex14
rosace_ridge_ex14 <- data.frame (score = ex14_rosace_scores$ROSACE_effects, 
                                 inhibitor = ex14_rosace_scores$inhibitor, 
                                 mutation_type = ex14_rosace_scores$type)

ggplot(rosace_ridge_ex14, aes(x = score, y = inhibitor)) +
  geom_density_ridges(alpha=0.6, stat="binline") +
  theme_ridges() + 
  theme()+
  ggtitle("MET Rosace Not Normalized")

# normalized
rosace_ridge_ex14_norm <- data.frame (score = ex14_rosace_norm_score$mean, 
                                      inhibitor = ex14_rosace_norm_score$key,
                                      mutation_type = ex14_rosace_norm_score$type, 
                                      variant = ex14_rosace_norm_score$variants)

ggplot(rosace_ridge_ex14_norm, aes(x = score, y = inhibitor)) +
  geom_density_ridges(alpha=0.6, stat="binline") +
  theme_ridges() + 
  theme()+
  ggtitle("MET Rosace Normalized to DMSO GR")


# not normalized METdelEx14
rosace_ridge_met <- data.frame (score = met_rosace_scores$ROSACE_effects, 
                                 inhibitor = met_rosace_scores$inhibitor, 
                                 mutation_type = met_rosace_scores$type)

ggplot(rosace_ridge_met, aes(x = score, y = inhibitor)) +
  geom_density_ridges(alpha=0.6, stat="binline") +
  theme_ridges() + 
  theme()+
  ggtitle("METdelEx14 Rosace Not Normalized")

# normalized
rosace_ridge_met_norm <- data.frame (score = met_rosace_norm_score$mean, 
                                      inhibitor = met_rosace_norm_score$key,
                                      mutation_type = met_rosace_norm_score$type)

ggplot(rosace_ridge_met_norm, aes(x = score, y = inhibitor)) +
  geom_density_ridges(alpha=0.6, stat="binline") +
  theme_ridges() + 
  theme()+
  ggtitle("METdelEx14 Rosace Normalized to DMSO GR")





#--------------------------------- Position Plots -------------------------------


plot1 <- ggplot() +geom_point(data = ex14_rosace_scores_wide, aes(x = pos, y = Crizo), alpha=0.3) + 
  xlab("Position") + ylab("Score") + ggtitle("Crizo") +theme_classic() 
plot(plot1)

plot2 <- ggplot() +geom_point(data = ex14_rosace_scores_wide, aes(x = pos, y = Cabo), alpha=0.3) + 
  xlab("Position") + ylab("Score") + ggtitle("Cabo") +theme_classic() 
plot(plot2)

plot3 <- ggplot() +geom_point(data = ex14_rosace_scores_wide, aes(x = pos, y = Camp), alpha=0.3) + 
  xlab("Position") + ylab("Score") + ggtitle("Cap")  +theme_classic() 
plot(plot3) 

plot4 <- ggplot() +geom_point(data = ex14_rosace_scores_wide, aes(x = pos, y = Tepo), alpha=0.3) + 
  xlab("Position") + ylab("Score") + ggtitle("Tepo") +theme_classic() 
plot(plot4)

plot5 <- ggplot() +geom_point(data = ex14_rosace_scores_wide, aes(x = pos, y = Tiv), alpha=0.3) + 
  xlab("Position") + ylab("Score") + ggtitle("Tiv")  +theme_classic() 
plot(plot5)

plot6 <- ggplot() +geom_point(data = ex14_rosace_scores_wide, aes(x = pos, y = Glu), alpha=0.3) + 
  xlab("Position") + ylab("Score") + ggtitle("Glu") +theme_classic() 
plot(plot6)

plot7 <- ggplot() +geom_point(data = ex14_rosace_scores_wide, aes(x = pos, y = Gle), alpha=0.3) + 
  xlab("Position") + ylab("Score") + ggtitle("Gle") +theme_classic() 
plot(plot7)

plot8 <- ggplot() +geom_point(data = ex14_rosace_scores_wide, aes(x = pos, y = Savo), alpha=0.3) + 
  xlab("Position") + ylab("Score") + ggtitle("Savo") +theme_classic() 
plot(plot8)

plot9 <- ggplot() +geom_point(data = ex14_rosace_scores_wide, aes(x = pos, y = NVP), alpha=0.3) + 
  xlab("Position") + ylab("Score") + ggtitle("NVP")  +theme_classic() 
plot(plot9)

plot10 <- ggplot() +geom_point(data = ex14_rosace_scores_wide, aes(x = pos, y = Mere), alpha=0.3) + 
  xlab("Position") + ylab("Score") + ggtitle("Mere") + theme_classic() 
plot(plot10)

plot11 <- ggplot() +geom_point(data = ex14_rosace_scores_wide, aes(x = pos, y = A458), alpha=0.3) + 
  xlab("Position") + ylab("Score") + ggtitle("AMG 458")+theme_classic() 
plot(plot11)

scatterplots_2 = ggarrange(plot1, plot2, plot3,
                         plot10, plot4, plot6,
                         plot7, plot8, plot9,
                         plot11, plot5,align = "hv",
                         nrow = 4, ncol = 3)
plot(scatterplots_2)

#####################################################################################################
###--------------------------------------- GOF mutations --------------------------------------------
#####################################################################################################

ex14_rosace_GOF <- data.frame (hgvs = ex14_rosace_scores$variant,
                               pos =ex14_rosace_scores$position, 
                               mutation = ex14_rosace_scores$mutation, 
                               score = ex14_rosace_scores$ROSACE_effects,
                               GOF = ex14_rosace_scores$ROSACE.GOF,
                               type = ex14_rosace_scores$type,
                               inhib= ex14_rosace_scores$inhibitor)
  
ex14_rosace_GOF <- ex14_rosace_GOF  %>% filter (GOF == "TRUE")

ex14_rosace_GOF_wide <- ex14_rosace_GOF%>%
  pivot_wider(names_from = inhib, values_from = score) 


plot1 <- ggplot() +geom_point(data = ex14_rosace_GOF_wide, aes(x = pos, y = Crizo), alpha=0.3) + 
  xlab("Position") + ylab("Score") + ggtitle("Crizo") +theme_classic() 
plot(plot1)

plot2 <- ggplot() +geom_point(data = ex14_rosace_GOF_wide, aes(x = pos, y = Cabo), alpha=0.3) + 
  xlab("Position") + ylab("Score") + ggtitle("Cabo") +theme_classic() 
plot(plot2)

plot3 <- ggplot() +geom_point(data = ex14_rosace_GOF_wide, aes(x = pos, y = Camp), alpha=0.3) + 
  xlab("Position") + ylab("Score") + ggtitle("Cap")  +theme_classic() 
plot(plot3) 

plot4 <- ggplot() +geom_point(data = ex14_rosace_GOF_wide, aes(x = pos, y = Tepo), alpha=0.3) + 
  xlab("Position") + ylab("Score") + ggtitle("Tepo") +theme_classic() 
plot(plot4)

plot5 <- ggplot() +geom_point(data = ex14_rosace_GOF_wide, aes(x = pos, y = Tiv), alpha=0.3) + 
  xlab("Position") + ylab("Score") + ggtitle("Tiv")  +theme_classic() 
plot(plot5)

plot6 <- ggplot() +geom_point(data = ex14_rosace_GOF_wide, aes(x = pos, y = Glu), alpha=0.3) + 
  xlab("Position") + ylab("Score") + ggtitle("Glu") +theme_classic() 
plot(plot6)

plot7 <- ggplot() +geom_point(data = ex14_rosace_GOF_wide, aes(x = pos, y = Gle), alpha=0.3) + 
  xlab("Position") + ylab("Score") + ggtitle("Gle") +theme_classic() 
plot(plot7)

plot8 <- ggplot() +geom_point(data = ex14_rosace_GOF_wide, aes(x = pos, y = Savo), alpha=0.3) + 
  xlab("Position") + ylab("Score") + ggtitle("Savo") +theme_classic() 
plot(plot8)

plot9 <- ggplot() +geom_point(data = ex14_rosace_GOF_wide, aes(x = pos, y = NVP), alpha=0.3) + 
  xlab("Position") + ylab("Score") + ggtitle("NVP")  +theme_classic() 
plot(plot9)

plot10 <- ggplot() +geom_point(data = ex14_rosace_GOF_wide, aes(x = pos, y = Mere), alpha=0.3) + 
  xlab("Position") + ylab("Score") + ggtitle("Mere") + theme_classic() 
plot(plot10)

plot11 <- ggplot() +geom_point(data = ex14_rosace_GOF_wide, aes(x = pos, y = A458), alpha=0.3) + 
  xlab("Position") + ylab("Score") + ggtitle("AMG 458")+theme_classic() 
plot(plot11)

scatterplots_3 = ggarrange(plot1, plot2, plot3,
                           plot10, plot4, plot6,
                           plot7, plot8, plot9,
                           plot11, plot5, 
                           nrow = 3, ncol = 4)
plot(scatterplots_3)


#####################################################################################################
###------------------------------------- Crizo vs Cabo ---------------------------------------------
#####################################################################################################

plot_Crizo_vs_Cabo <- ggplot() +geom_point(data = ex14_rosace_scores_wide, aes(x = Crizo, y = Cabo),alpha=0.3) + 
  xlab("Crizo") + ylab("Cabo") + coord_fixed(ratio = 1) +
  geom_vline(xintercept=0, color = "red", linetype = "dashed")+
  geom_hline(yintercept=0, color = "red", linetype = "dashed")+
  theme_classic() 
plot(plot_Crizo_vs_Cabo )

plot_Crizo_vs_Mere <- ggplot() +geom_point(data = ex14_rosace_scores_wide, aes(x = Crizo, y = Mere),alpha=0.3) + 
  xlab("Crizo") + ylab("Mere") + coord_fixed(ratio = 1) +
  geom_vline(xintercept=0, color = "red", linetype = "dashed")+
  geom_hline(yintercept=0, color = "red", linetype = "dashed")+
  theme_classic() 
plot(plot_Crizo_vs_Mere )

plot_Crizo_vs_Gle <- ggplot() +geom_point(data = ex14_rosace_scores_wide, aes(x = Crizo, y = Gle),alpha=0.3) + 
  xlab("Crizo") + ylab("Gle") + coord_fixed(ratio = 1) +
  geom_vline(xintercept=0, color = "red", linetype = "dashed")+
  geom_hline(yintercept=0, color = "red", linetype = "dashed")+
  theme_classic() 
plot(plot_Crizo_vs_Gle )


Crizo_vs_Cabo <- data.frame(pos = ex14_rosace_scores_wide$pos, 
                            Crizo = ex14_rosace_scores_wide$Crizo, 
                            Cabo = ex14_rosace_scores_wide$Cabo)

Crizo_GOF_Cabo_LOF <- Crizo_vs_Cabo %>% filter (Crizo>= 1 & Cabo <= 0)
Cabo_GOF_Crizo_LOF <- Crizo_vs_Cabo %>% filter (Cabo>= 1 & Crizo <= 0)

ex14_2WGJ_crizo <- read.pdb("2WGJ")
ex14_4EEV_crizo <- read.pdb("4EEV")

avg_Crizo_GOF_Cabo_LOF = Crizo_GOF_Cabo_LOF %>% group_by(pos) %>% dplyr::summarise(mean = mean(Crizo, na.rm=TRUE))
x = map_scores_pdb(ex14_2WGJ_crizo, avg_Crizo_GOF_Cabo_LOF, "mean")
write.pdb(x, file="avg_Crizo_GOF_cabo_LOF.pdb")

avg_Cabo_GOF_Crizo_LOF = Cabo_GOF_Crizo_LOF %>% group_by(pos) %>% dplyr::summarise(mean = mean(Cabo, na.rm=TRUE))
x = map_scores_pdb(ex14_4EEV_crizo, avg_Cabo_GOF_Crizo_LOF, "mean")
write.pdb(x, file="avg_Cabo_GOF_Crizo_LOF.pdb")


#####################################################################################################
###-------------------------- Crizotinib Sensitive, GOF, Resis -------------------------------------
#####################################################################################################

ex14_rosace_crizo = ex14_rosace_scores %>% filter (inhibitor == "Crizo")
ex14_rosace_crizo_LOF = ex14_rosace_crizo %>% filter (ROSACE.LOF == "TRUE")
ex14_rosace_crizo_GOF = ex14_rosace_crizo %>% filter (ROSACE.GOF == "TRUE" &  type != "synonymous")


ex14_rosace_crizo_dmso <- data.frame(pos = ex14_rosace_scores_wide$pos, 
                                     mutation = ex14_rosace_scores_wide$mutation, 
                                     type = ex14_rosace_scores_wide$type, 
                                     DMSO = ex14_rosace_scores_wide$DMSO, 
                                     Crizo = ex14_rosace_scores_wide$Crizo)


#----------- count number of GOF mutations per position ------------

ex14_crizo_GOF_count <-data.frame ( pos = ex14_rosace_crizo_GOF$position, 
                                    count = stat(ncount))


#----------- map mutations on structure ------------

ex14_2WGJ_crizo <- read.pdb("2WGJ")
avg_ex14_rosace_crizo_GOF = ex14_rosace_crizo_GOF%>% group_by(pos) %>% dplyr::summarise(mean = mean(ROSACE_effects, na.rm=TRUE))
x = map_scores_pdb(ex14_2WGJ_crizo, avg_ex14_rosace_crizo_GOF, "mean")
write.pdb(x, file="avg_ex14_rosace_crizo_GOF.pdb")



#####################################################################################################
###--------------------------- plot all inhibitors by all inhibitors ---------------------------------
#####################################################################################################

# create a plot where every inhibitor is ploted against every other inhihibitor 
score_by_score <- ex14_PCA [-c(2:4)]
pairs(score_by_score )

#####################################################################################################
###----------------------- ROSACE general score mapping on structures ----------------------------
#####################################################################################################


#--------------------------------- METdEx14 --------------------------------
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

avg_ex14_rosace_1.5 = ex14_rosace_resis_scores%>% group_by(pos) %>% dplyr::summarise(I5_mean = mean("I1/2", na.rm=TRUE))
ex14_5DG5  <- read.pdb("5T3Q")
x = map_scores_pdb(ex14_5DG5  , avg_ex14_rosace_II, "II_mean")
write.pdb(x, file="avg_ex14_rosace_I5.pdb")



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

Crizo_Y1093 <- ex14_rosace_scores_crizo %>% filter(pos == "1093" )
Crizo_Y1093_ <- (data.frame(variants = Crizo_Y1159$mutation, score = Crizo_Y1159$ROSACE_effects, pos = Crizo_Y1159$pos))


Crizo_I1084 <- ex14_rosace_scores_crizo %>% filter(pos == "1084" )
Crizo_I1084_ <- (data.frame(variants = Crizo_I1084$mutation, score = Crizo_I1084$ROSACE_effects, pos = Crizo_I1084$pos))


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

ggplot(Crizo_Y1093_, aes(pos, variants)) +
  geom_tile(aes(fill = score)) +
  geom_text(aes(label = round(score, 1))) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(-12,5))+
  ggtitle("Y1093")

ggplot(Crizo_I1084_, aes(pos, variants)) +
  geom_tile(aes(fill = score)) +
  geom_text(aes(label = round(score, 1))) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(-12,5))+
  ggtitle("I1084")


# generate subtraction heatmap 
order <- c('X','H', 'K', 'R', 'D', 'E',
           'C','M','N','Q','S','T','A','I','L','V','W','F',
           'Y','G','P')

variant_names <- c('Stop','H', 'K', 'R', 'D', 'E',
                   'C','M','N','Q','S','T','A','I','L','V','W','F',
                   'Y','G','P')


met_wt_fasta = "Met_wt2.fasta"

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


