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
#$library(ggbraid)
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
library(janitor)
library(pheatmap)
library(data.table)

# color palette
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#0072B2","#F0E442","#D55E00", "#CC79A7")

#####################################################################################################
###------------------------------------ define initial data frames ----------------------------------
#####################################################################################################


##-------------------------------------- Rosace scores ----------------------------------------------
### open data file as data frame for all raw inhibitor ALL scores generated from ROSACE

met_rosace_scores <- as.data.frame(fread("met_scores_unfiltered.tsv")) # here met means MET delta Ex14
met_rosace_scores  <- met_rosace_scores  %>% mutate(position = position + 1058)

met_rosace_norm_scores <- as.data.frame(fread("met_scores_unfiltered.tsv")) # here met means MET delta Ex14
met_rosace_norm_scores  <- met_rosace_norm_scores  %>% mutate(position = position + 1058)

ex14_rosace_scores <-as.data.frame(fread("ex14_scores_unfiltered.tsv")) # here ex14 means MET+Ex14
ex14_rosace_scores  <- ex14_rosace_scores   %>% mutate(pos = position + 1058)

ex14_rosace_norm_scores <- as.data.frame(fread("ex14_scores_unfiltered.tsv")) # here ex14 means MET+Ex14 
ex14_rosace_norm_scores  <- ex14_rosace_norm_scores   %>% mutate(pos = position + 1058)


#####################################################################################################
###---------------------------------------- Prelim analysis -----------------------------------------
#####################################################################################################

# type 1 : crizo, cap, tepo , tiv , glu, savo, nvp 
# type 2 : cabo, mere, gle 
# type 1.5: AMG 

# MET WT inhibitor, DMSO correlation plots 


#--------------------------------- Correlation Plots MET+Ex14 --------------------------------

ex14_rosace_scores_condensed <- data.frame (hgvs = ex14_rosace_scores$variant,
                                            pos =ex14_rosace_scores$position, 
                                            mutation = ex14_rosace_scores$mutation, 
                                            score = ex14_rosace_scores$mean,
                                            type = ex14_rosace_scores$type,
                                            inhib= ex14_rosace_scores$key
)

ex14_rosace_scores_wide <- ex14_rosace_scores_condensed %>%
  pivot_wider(names_from = inhib, values_from = score) 

# Define a custom color palette
#custom_palette <- c("#055b5c","#08979d", "#8474A1","#CCABD8", "grey")
custom_palette <- c("#e63466ff","#f088a5ff", "#56B4E9","#8cbe61ff","#ff9b00ff")
# Filter out non-numeric values and missing data
ex14_rosace_scores_wide_numeric <- ex14_rosace_scores_wide %>%
  filter(
    !is.na(DMSO) & !is.na(Crizo) & !is.na(Camp) & !is.na(Tepo) &
      !is.na(Tiv) & !is.na(Glu) & !is.na(Savo) & !is.na(NVP) &
      !is.na(Cabo) & !is.na(Mere) & !is.na(Gle) & !is.na(A458)
  )

# Create scatterplots with regression lines and Pearson's R values
scatterplots <- ggarrange(
  ggplot(ex14_rosace_scores_wide_numeric, aes(x = DMSO, y = Crizo)) +
    geom_point(color = custom_palette[1], shape=1, size = 0.5) +
    xlab("DMSO") + ylab("Crizotinib") +
    #oord_fixed(ratio = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    stat_cor(method = "pearson",aes(label = ..r.label..), label.x = -3, label.y = 2, label.sep = "\n", label.case = "lower") +
    ylim(-3.5, 3) +
    theme_classic(),
  
  ggplot(ex14_rosace_scores_wide_numeric, aes(x = DMSO, y = Tepo)) +
    geom_point(color = custom_palette[1], shape=1, size = 0.5) +
    xlab("DMSO") + ylab("Tepotinib") +
    # coord_fixed(ratio = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    stat_cor(method = "pearson", aes(label = ..r.label..),label.x = -3, label.y = 2, label.sep = "\n", label.case = "lower") +
    ylim(-3.5, 3) +
    theme_classic(),
  
  ggplot(ex14_rosace_scores_wide_numeric, aes(x = DMSO, y = Camp)) +
    geom_point(color = custom_palette[1], shape=1, size = 0.5) +
    xlab("DMSO") + ylab("Capmatinib") +
  # coord_fixed(ratio = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    stat_cor(method = "pearson", aes(label = ..r.label..),label.x = -3, label.y = 2, label.sep = "\n", label.case = "lower") +
    ylim(-3.5, 3) +
    theme_classic(),
  
  ggplot(ex14_rosace_scores_wide_numeric, aes(x = DMSO, y = Glu)) +
    geom_point(color = custom_palette[1], shape=1, size = 0.5) +
    xlab("DMSO") + ylab("Glumetinib") +
    #oord_fixed(ratio = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    stat_cor(method = "pearson", aes(label = ..r.label..),label.x = -3, label.y = 2, label.sep = "\n", label.case = "lower") +
    ylim(-3.5, 3) +
    theme_classic(),
  
  ggplot(ex14_rosace_scores_wide_numeric, aes(x = DMSO, y = Savo)) +
    geom_point(color = custom_palette[2], shape=1, size = 0.5) +
    xlab("DMSO") + ylab("Savolitinib") +
    #oord_fixed(ratio = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    stat_cor(method = "pearson", aes(label = ..r.label..),label.x = -3, label.y = 2, label.sep = "\n", label.case = "lower") +
    ylim(-3.5, 3) +
    theme_classic(),
  
  ggplot(ex14_rosace_scores_wide_numeric, aes(x = DMSO, y = NVP)) +
    geom_point(color = custom_palette[2], shape=1, size = 0.5) +
    xlab("DMSO") + ylab("NVP-BVU972") +
    #oord_fixed(ratio = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    stat_cor(method = "pearson", aes(label = ..r.label..),label.x = -3, label.y = 2, label.sep = "\n", label.case = "lower") +
    ylim(-3.5, 3) +
    theme_classic(),
  
  ggplot(ex14_rosace_scores_wide_numeric, aes(x = DMSO, y = Cabo)) +
    geom_point(color = custom_palette[3], shape=1, size = 0.5) +
    xlab("DMSO") + ylab("Cabozantinib") +
    #oord_fixed(ratio = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    stat_cor(method = "pearson", aes(label = ..r.label..),label.x = -3, label.y = 2, label.sep = "\n", label.case = "lower") +
    ylim(-3.5, 3) +
    theme_classic(),
  
  ggplot(ex14_rosace_scores_wide_numeric, aes(x = DMSO, y = Mere)) +
    geom_point(color = custom_palette[3], shape=1, size = 0.5) +
    xlab("DMSO") + ylab("Merestinib") +
    #oord_fixed(ratio = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    stat_cor(method = "pearson", aes(label = ..r.label..),label.x = -3, label.y = 2, label.sep = "\n", label.case = "lower") +
    ylim(-3.5, 3) +
    theme_classic(),
  
  ggplot(ex14_rosace_scores_wide_numeric, aes(x = DMSO, y = Gle)) +
    geom_point(color = custom_palette[3], shape=1, size = 0.5) +
    xlab("DMSO") + ylab("Glesatinib") +
    #oord_fixed(ratio = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    stat_cor(method = "pearson", aes(label = ..r.label..),label.x = -3, label.y = 2, label.sep = "\n", label.case = "lower") +
    ylim(-3.5, 3) +
    theme_classic(),
  
  ggplot(ex14_rosace_scores_wide_numeric, aes(x = DMSO, y = A458)) +
    geom_point(color = custom_palette[4], shape=1, size = 0.5) +
    xlab("DMSO") + ylab("AMG-458") +
    #oord_fixed(ratio = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    stat_cor(method = "pearson", aes(label = ..r.label..), label.x = -3, label.y = 2, label.sep = "\n", label.case = "lower") +
    ylim(-3.5, 3) +
    theme_classic(),
  
  ggplot(ex14_rosace_scores_wide_numeric, aes(x = DMSO, y = Tiv)) +
    geom_point(color = custom_palette[5], shape=1, size = 0.5) +
    xlab("DMSO") + ylab("Tivantinib") +
    #oord_fixed(ratio = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    stat_cor(method = "pearson", aes(label = ..r.label..), label.x = -3, label.y = 2, label.sep = "\n", label.case = "lower") +
    ylim(-3.5, 3) +
    theme_classic(),
  
  nrow =4, ncol = 3
)

# Display the scatterplots
print(scatterplots)
# Save the scatterplots as a PDF file
ggsave("scatterplots.png", plot = scatterplots, height = 2, width = 8, device = "png")



#--------------------------------- Correlation Plots MET+Ex14 --------------------------------

met_rosace_scores_condensed <- data.frame (hgvs = met_rosace_scores$variant,
                                            pos =met_rosace_scores$position, 
                                            mutation = met_rosace_scores$mutation, 
                                            score = met_rosace_scores$mean,
                                            type = met_rosace_scores$type,
                                            inhib= met_rosace_scores$key
)

met_rosace_scores_wide <- met_rosace_scores_condensed %>%
  pivot_wider(names_from = inhib, values_from = score) 

# Define a custom color palette
#custom_palette <- c("#055b5c","#08979d", "#8474A1","#CCABD8", "grey")
custom_palette <- c("#e63466ff","#f088a5ff", "#56B4E9","#8cbe61ff","#ff9b00ff")
# Filter out non-numeric values and missing data
met_rosace_scores_wide_numeric <- met_rosace_scores_wide %>%
  filter(
    !is.na(DMSO) & !is.na(Crizo) & !is.na(Camp) & !is.na(Tepo) &
      !is.na(Tiv) & !is.na(Glu) & !is.na(Savo) & !is.na(NVP) &
      !is.na(Cabo) & !is.na(Mere) & !is.na(Gle) & !is.na(A458)
  )

# Create scatterplots with regression lines and Pearson's R values
scatterplots <- ggarrange(
  ggplot(met_rosace_scores_wide_numeric, aes(x = DMSO, y = Crizo)) +
    geom_point(color = custom_palette[1], shape=1, size = 0.5) +
    xlab("DMSO") + ylab("Crizotinib") +
    #oord_fixed(ratio = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    stat_cor(method = "pearson",aes(label = ..r.label..), label.x = -5, label.y = 2, label.sep = "\n") +
    ylim(-5, 3) +
    theme_classic(),
  
  ggplot(met_rosace_scores_wide_numeric, aes(x = DMSO, y = Tepo)) +
    geom_point(color = custom_palette[1], shape=1, size = 0.5) +
    xlab("DMSO") + ylab("Tepotinib") +
    # coord_fixed(ratio = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    stat_cor(method = "pearson", aes(label = ..r.label..),label.x = -5, label.y = 2, label.sep = "\n") +
    ylim(-5, 3) +
    theme_classic(),
  
  ggplot(met_rosace_scores_wide_numeric, aes(x = DMSO, y = Camp)) +
    geom_point(color = custom_palette[1], shape=1, size = 0.5) +
    xlab("DMSO") + ylab("Capmatinib") +
    # coord_fixed(ratio = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    stat_cor(method = "pearson", aes(label = ..r.label..),label.x = -5, label.y = 2, label.sep = "\n") +
    ylim(-5, 3) +
    theme_classic(),
  
  ggplot(met_rosace_scores_wide_numeric, aes(x = DMSO, y = Glu)) +
    geom_point(color = custom_palette[1], shape=1, size = 0.5) +
    xlab("DMSO") + ylab("Glumetinib") +
    #oord_fixed(ratio = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    stat_cor(method = "pearson", aes(label = ..r.label..),label.x = -5, label.y = 2, label.sep = "\n") +
    ylim(-5, 3) +
    theme_classic(),
  
  ggplot(met_rosace_scores_wide_numeric, aes(x = DMSO, y = Savo)) +
    geom_point(color = custom_palette[2], shape=1, size = 0.5) +
    xlab("DMSO") + ylab("Savolitinib") +
    #oord_fixed(ratio = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    stat_cor(method = "pearson", aes(label = ..r.label..),label.x = -5, label.y = 2, label.sep = "\n") +
    ylim(-5, 3) +
    theme_classic(),
  
  ggplot(met_rosace_scores_wide_numeric, aes(x = DMSO, y = NVP)) +
    geom_point(color = custom_palette[2], shape=1, size = 0.5) +
    xlab("DMSO") + ylab("NVP-BVU972") +
    #oord_fixed(ratio = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    stat_cor(method = "pearson", aes(label = ..r.label..),label.x = -5, label.y = 2, label.sep = "\n") +
    ylim(-5, 3) +
    theme_classic(),
  
  ggplot(met_rosace_scores_wide_numeric, aes(x = DMSO, y = Cabo)) +
    geom_point(color = custom_palette[3], shape=1, size = 0.5) +
    xlab("DMSO") + ylab("Cabozantinib") +
    #oord_fixed(ratio = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    stat_cor(method = "pearson", aes(label = ..r.label..),label.x = -5, label.y = 2, label.sep = "\n") +
    ylim(-5, 3) +
    theme_classic(),
  
  ggplot(met_rosace_scores_wide_numeric, aes(x = DMSO, y = Mere)) +
    geom_point(color = custom_palette[3], shape=1, size = 0.5) +
    xlab("DMSO") + ylab("Merestinib") +
    #oord_fixed(ratio = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    stat_cor(method = "pearson", aes(label = ..r.label..),label.x = -5, label.y = 2, label.sep = "\n") +
    ylim(-5, 3) +
    theme_classic(),
  
  ggplot(met_rosace_scores_wide_numeric, aes(x = DMSO, y = Gle)) +
    geom_point(color = custom_palette[3], shape=1, size = 0.5) +
    xlab("DMSO") + ylab("Glesatinib") +
    #oord_fixed(ratio = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    stat_cor(method = "pearson", aes(label = ..r.label..),label.x = -5, label.y = 2, label.sep = "\n") +
    ylim(-5, 3) +
    theme_classic(),
  
  ggplot(met_rosace_scores_wide_numeric, aes(x = DMSO, y = A458)) +
    geom_point(color = custom_palette[4], shape=1, size = 0.5) +
    xlab("DMSO") + ylab("AMG-458") +
    #oord_fixed(ratio = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    stat_cor(method = "pearson", aes(label = ..r.label..), label.x = -5, label.y = 2, label.sep = "\n") +
    ylim(-5, 3) +
    theme_classic(),
  
  ggplot(met_rosace_scores_wide_numeric, aes(x = DMSO, y = Tiv)) +
    geom_point(color = custom_palette[5], shape=1, size = 0.5) +
    xlab("DMSO") + ylab("Tivantinib") +
    #oord_fixed(ratio = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    stat_cor(method = "pearson", aes(label = ..r.label..), label.x = -5, label.y = 2, label.sep = "\n") +
    ylim(-5, 3) +
    theme_classic(),
  
  nrow =4, ncol = 3
)

# Display the scatterplots
print(scatterplots)
# Save the scatterplots as a PDF file
ggsave("metdelex14_scatterplots.png", plot = scatterplots, height = 2, width = 8, device = "png")


#--------------------------------- Distributions--------------------------------


ex14_rosace_norm_scores <- as.data.frame(fread("ex14_scores_filtered.tsv")) # here ex14 means MET+Ex14 
ex14_rosace_norm_scores  <- ex14_rosace_norm_scores   %>% mutate(pos = position + 1058)


# not normalized MET+Ex14
rosace_ridge_ex14 <- data.frame (score = ex14_rosace_scores$ROSACE_effects, 
                                 inhibitor = ex14_rosace_scores$inhibitor, 
                                 mutation_type = ex14_rosace_scores$type)

ggplot(rosace_ridge_ex14, aes(x = score, y = inhibitor,color = mutation_type)) +
  geom_density_ridges(alpha=0.6, stat="binline") +
  theme_ridges() + 
  theme()+
  ggtitle("MET Rosace Not Normalized")


# normalized

custom_palette <- c("#ff9b00ff","#8cbe61ff", "#56B4E9", "#56B4E9", "#56B4E9", "#f088a5ff", "#f088a5ff", "#e63466ff", "#e63466ff","#e63466ff", "#e63466ff", "darkmagenta")
inhibitor_order <- c("Tiv", "A458","Gle","Mere","Cabo","NVP","Savo","Glu","Camp","Tepo","Crizo", "DMSO")

inhibitor_labels <- c("Tivantinib", "AMG-458", "Glesatinib", "Merestinib", "Cabozantinib",
                      "NVP-BVU972", "Savolitinib", "Glumetinib", "Campmatiinb", "Tepotinib",
                      "Crizotinib", "DMSO")


rosace_ridge_ex14_norm <- data.frame(
  score = ex14_rosace_norm_scores$mean,
  inhibitor = factor(ex14_rosace_norm_scores$key, levels = inhibitor_order),
  mutation_type = ex14_rosace_norm_scores$type
)


ggplot(rosace_ridge_ex14_norm, aes(x = score, y = inhibitor, fill = inhibitor, group = inhibitor)) +
  geom_density_ridges(alpha = 0.6, stat = "binline", bins = 20) +
  scale_fill_manual(values = custom_palette) +
  scale_x_continuous(
    breaks = seq(min(rosace_ridge_ex14_norm$score), max(rosace_ridge_ex14_norm$score), by = 1),
    labels = scales::number_format(scale = 1),  # Format the labels if necessary
    expand = c(0.05, 0)  # Adjust the margin on the x-axis
  ) +
  theme_ridges(grid = FALSE) +  # Disable grid lines
  labs(x = "Score", y = "Inhibitor") +
  ggtitle("MET Rosace Normalized to DMSO GR") +
  theme(axis.title.y = element_text(margin = margin(r = 10, unit = "pt")),  # Adjust y-axis title margin
        axis.text.y = element_text(margin = margin(r = 5, unit = "pt")))  # Adjust y-axis text margin



#------ met del ex14 ridge plot -----------#
met_rosace_norm_scores <- as.data.frame(fread("met_scores_unfiltered.tsv")) # here met means METdelEx14
met_rosace_norm_scores  <- met_rosace_norm_scores   %>% mutate(pos = position + 1058)

rosace_ridge_met_norm <- data.frame(
  score = met_rosace_norm_scores$mean,
  inhibitor = factor(met_rosace_norm_scores$key, levels = inhibitor_order),
  mutation_type = met_rosace_norm_scores$type
)


ggplot(rosace_ridge_met_norm, aes(x = score, y = inhibitor, fill = inhibitor, group = inhibitor)) +
  geom_density_ridges(alpha = 0.6, stat = "binline", bins = 20) +
  scale_fill_manual(values = custom_palette) +
  scale_x_continuous(
    breaks = seq(min(rosace_ridge_met_norm$score), max(rosace_ridge_met_norm$score), by = 1),
    labels = scales::number_format(scale = 1),  # Format the labels if necessary
    expand = c(0.05, 0)  # Adjust the margin on the x-axis
  ) +
  theme_ridges(grid = FALSE) +  # Disable grid lines
  labs(x = "Score", y = "Inhibitor") +
  ggtitle("METdelEx14 Rosace Normalized to DMSO GR") +
  theme(axis.title.y = element_text(margin = margin(r = 10, unit = "pt")),  # Adjust y-axis title margin
        axis.text.y = element_text(margin = margin(r = 5, unit = "pt")))  # Adjust y-axis text margin






#--------------------------------- Position Plots -------------------------------


plot1 <- ggplot() +geom_point(data = ex14_rosace_scores_wide, aes(x = pos, y = Crizo), alpha=0.2) + 
  xlab("Position") + ylab("Score") + ggtitle("Crizo") + 
  geom_hline(yintercept=0, linetype="dashed", color = "red", size=0.5)+
  theme_classic() 
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



#------ met del ex14


met_rosace_scores_condensed <- data.frame (hgvs = met_rosace_scores$variant,
                                            pos =met_rosace_scores$position, 
                                            mutation = met_rosace_scores$mutation, 
                                            score = met_rosace_scores$mean,
                                            type = met_rosace_scores$type,
                                            inhib= met_rosace_scores$key
)

met_rosace_scores_wide <- met_rosace_scores_condensed %>%
  pivot_wider(names_from = inhib, values_from = score) 


plot1 <- ggplot() +geom_point(data = met_rosace_scores_wide, aes(x = pos, y = Crizo), alpha=0.2) + 
  xlab("Position") + ylab("Score") + ggtitle("Crizo") + 
  geom_hline(yintercept=0, linetype="dashed", color = "red", size=0.5)+
  theme_classic() 
plot(plot1)

plot2 <- ggplot() +geom_point(data = met_rosace_scores_wide, aes(x = pos, y = Cabo), alpha=0.3) + 
  xlab("Position") + ylab("Score") + ggtitle("Cabo") +theme_classic() 
plot(plot2)

plot3 <- ggplot() +geom_point(data = met_rosace_scores_wide, aes(x = pos, y = Camp), alpha=0.3) + 
  xlab("Position") + ylab("Score") + ggtitle("Cap")  +theme_classic() 
plot(plot3) 

plot4 <- ggplot() +geom_point(data = met_rosace_scores_wide, aes(x = pos, y = Tepo), alpha=0.3) + 
  xlab("Position") + ylab("Score") + ggtitle("Tepo") +theme_classic() 
plot(plot4)

plot5 <- ggplot() +geom_point(data = met_rosace_scores_wide, aes(x = pos, y = Tiv), alpha=0.3) + 
  xlab("Position") + ylab("Score") + ggtitle("Tiv")  +theme_classic() 
plot(plot5)

plot6 <- ggplot() +geom_point(data = met_rosace_scores_wide, aes(x = pos, y = Glu), alpha=0.3) + 
  xlab("Position") + ylab("Score") + ggtitle("Glu") +theme_classic() 
plot(plot6)

plot7 <- ggplot() +geom_point(data = met_rosace_scores_wide, aes(x = pos, y = Gle), alpha=0.3) + 
  xlab("Position") + ylab("Score") + ggtitle("Gle") +theme_classic() 
plot(plot7)

plot8 <- ggplot() +geom_point(data = met_rosace_scores_wide, aes(x = pos, y = Savo), alpha=0.3) + 
  xlab("Position") + ylab("Score") + ggtitle("Savo") +theme_classic() 
plot(plot8)

plot9 <- ggplot() +geom_point(data = met_rosace_scores_wide, aes(x = pos, y = NVP), alpha=0.3) + 
  xlab("Position") + ylab("Score") + ggtitle("NVP")  +theme_classic() 
plot(plot9)

plot10 <- ggplot() +geom_point(data = met_rosace_scores_wide, aes(x = pos, y = Mere), alpha=0.3) + 
  xlab("Position") + ylab("Score") + ggtitle("Mere") + theme_classic() 
plot(plot10)

plot11 <- ggplot() +geom_point(data = met_rosace_scores_wide, aes(x = pos, y = A458), alpha=0.3) + 
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
###-------------------------------------------- PCA -------------------------------------------------
#####################################################################################################


library(factoextra)
library(FactoMineR)

key_list <- c("DMSO","A458","Crizo","Camp","Tepo","Tiv","Glu","Savo","NVP","Cabo","Mere","Gle")
group_list <- c("Control", "I 1/2", "I", "I", "I", "I", "I", "I", "I", "II", "II", "II")
df_inh <- data.frame(key = key_list, group = group_list)


score_PCA <- ex14_rosace_scores %>% 
  dplyr::select(variants, mean, key) %>%
  tidyr::pivot_wider(names_from = key, values_from = mean) %>%
  tidyr::drop_na() # only use complete row, drop NA

variant_label <- score_PCA$variants
inhibitor_label <- colnames(score_PCA)[-1]

group_label <- df_inh$group[as.numeric(sapply(inhibitor_label, function(x) {which(x == df_inh$key)}))]
d <- score_PCA %>% select(-variants) # d is a score matrix
d <- t(d)

rownames(d) <- inhibitor_label
colnames(d) <- variant_label

res.pca <- prcomp(d, center = TRUE, scale = TRUE) # PCA result
fviz_pca_ind(res.pca, col.ind = group_label,
             legend.title = "Type", repel = TRUE,
             axes = c(1, 2)) # axes: change PCA axe. Ex: c(1, 3) 
# individuals 
fviz_pca_ind(res.pca, 
             col.ind = group_label,
             legend.title = "Type", 
             ggtheme = theme_minimal(base_size = 16))


res.pca <- PCA(d[,-5745], scale = TRUE)
res.hcpc <- HCPC(res.pca, graph = FALSE)


# PCA 
fviz_pca_ind(res.pca, 
             col.ind = group_label,
             legend.title = "Type", 
             palette ="jco")+
  theme_classic()+
  ylim(-100,100)


# Dendogram
fviz_dend(res.hcpc, 
          cex = 0.7,                     # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          #labels_track_height = 0.8      # Augment the room for labels
)

# Factor cluster
fviz_cluster(res.hcpc,
             repel = TRUE,            # Avoid label overlapping
             show.clust.cent = TRUE, # Show cluster centers
             palette = "jco",         # Color palette see ?ggpubr::ggpar
             ggtheme = theme_minimal(),
             main = "Factor map")

fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50))


res.pca <- prcomp(ex14_PCA_3[,-5745], scale = TRUE)
ex14_PCA_df <-unclass(res.pca)
ex14_PCA_df <- data.frame(ex14_PCA_df$x) # turns PCA into data frame based on inhibitor and PC's

dt2 <- ex14_PCA_df %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
head(dt2)
pheatmap(t(ex14_PCA_df), 
         border_color = "black",
         #cluster_rows = FALSE,
         cutree_cols = 3,
         #cutree_rows = 3,
         angle_col = 90,             # angle for column labels
         width = 5,                 # plot width in inches
         height = 7, 
         color = colorRampPalette(c("#E7B800","white","dodgerblue"))(100))


#-------------------------------- Positon level PCA -------------------------------------
ex14_rosace_norm_wide <-data.frame( variant = ex14_rosace_norm_score$variant,
                                    mutation =ex14_rosace_norm_score$mutation,
                                    key = ex14_rosace_norm_score$key,
                                    pos = ex14_rosace_norm_score$position,
                                    type =ex14_rosace_norm_score$type,
                                    score =ex14_rosace_norm_score$mean
)

ex14_rosace_norm_wide <- ex14_rosace_norm_wide %>% pivot_wider(names_from = key, values_from = score)


# condensed data frame, only scores, position, and inhibitor 
#ex14_PCA <- ex14_rosace_scores_wide[-c(2:4)]
ex14_PCA <- ex14_rosace_norm_wide[-c(2:4)]

# remove hgvs data to keep numerican
#ex14_PCA_4 <- ex14_PCA %>% remove_rownames %>% column_to_rownames(var="hgvs")
ex14_PCA_4 <- ex14_PCA %>% remove_rownames %>% column_to_rownames(var="variant")
#ex14_PCA_5 <- ex14_PCA_4

ex14_PCA_4 <- t(ex14_PCA_4)
ex14_PCA_4 <- data.frame(ex14_PCA_4)


#res.pca2 <- PCA(ex14_PCA_5, scale = TRUE)
res.pca2 <- prcomp(ex14_PCA_4, scale =FALSE)
res.pca2 <- PCA(ex14_PCA_4, scale = FALSE)
res.hcpc2 <- HCPC(res.pca2, graph = FALSE)


# PCA 
fviz_pca_ind(res.pca2, 
             legend.title = "Type", 
             palette ="jco")+
  theme_classic()


# Dendogram
fviz_dend(res.hcpc2, 
          cex = 0.7,                     # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          #labels_track_height = 0.8      # Augment the room for labels
)

# Factor cluster
fviz_cluster(res.hcpc2,
             repel = TRUE,            # Avoid label overlapping
             show.clust.cent = TRUE, # Show cluster centers
             palette = "jco",         # Color palette see ?ggpubr::ggpar
             ggtheme = theme_minimal(),
             main = "Factor map")


res.pca <- prcomp(ex14_PCA_4, scale = TRUE)
ex14_PCA_df <-unclass(res.pca)
ex14_PCA_df <- data.frame(ex14_PCA_df$x) # turns PCA into data frame based on inhibitor and PC's

dt2 <- ex14_PCA_df %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
head(dt2)

pheatmap(t(ex14_PCA_df), 
         border_color = "black",
         #cluster_rows = FALSE,
         cutree_cols = 3,
         #cutree_rows = 3,
         angle_col = 90,             # angle for column labels
         width = 5,                 # plot width in inches
         height = 7, 
         color = colorRampPalette(c("#E7B800","white","dodgerblue"))(100))

#-------------- Heatmap of each Principle component ---------------------

#subtract Mut_drug from Mut_dmso
ex14_PCA_DMSO_sub <- ex14_PCA_4 %>% mutate_at(1:12, funs(c(first(.), (. - first(.))[-1]))) # (Mut_drug)-(Mut_DMSO)
ex14_PCA_DMSO_sub <- ex14_PCA_DMSO_sub[-1,] #remove DMSO from df for PCA
ex14_PCA_DMSO_sub <-data.frame(ex14_PCA_DMSO_sub)
#ex14_PCA_DMSO_sub <- t(ex14_PCA_DMSO_sub) 

# do PCA 
# scale=TRUE bases the PCA on the correlation matrix and FALSE on the covariance matrix
# scores are already normalized to DMSO, so scaling is FALSE in this case 
#res.pca3 <-prcomp(ex14_PCA_DMSO_sub, center = TRUE, scale = FALSE)
res.pca3 <- PCA(ex14_PCA_4, scale = FALSE)
res.hcpc3 <- HCPC(res.pca3, graph = FALSE)
#ex14_PCs <-unclass(res.pca) # this allows you to look at all the functions 

#plot contributions
mut.contr <- fviz_contrib(res.pca3, choice = "var", axes = 1:2, fill = "lightblue", color = "darkblue",top=45)
print(mut.contr)

# PCA 
fviz_pca_ind(res.pca3, legend.title = "Type", palette ="jco", ylim = c(-100,100))+ theme_classic()
indvidual <- data.frame((res.pca3$ind$coord)) # flip signage to match the DMS 

# plot PCs with correct signage 
indvidual["inhib"] <- rownames(indvidual)

ggplot(indvidual, aes(x=Dim.1, y=Dim.2, label = inhib)) + 
  geom_point()+
  geom_text()+
  xlim(-100,100)+
  ylim(-100,100)+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_vline(xintercept=0, linetype="dashed")+
  ggtitle("Individuals PCA")+
  xlab("PC1 (42.2%")+
  ylab("PC1 (24.4%")+
  theme_classic()

# plots contributions 
fviz_screeplot(res.pca3, addlabels = TRUE, ylim = c(0, 60))

# contribution of each mutation 
ex14_PCs <- data.frame((res.pca3$var$coord)) # loadings ie contributions
#ex14_PCs <- data.frame(res.pca3$rotation) # for prcomp usage to get loadings 
ex14_PCs$pos <- ex14_rosace_norm_wide$pos
ex14_PCs$mutation <- ex14_rosace_norm_wide$mutation


# contribution of inhibitor in PC space 
#ex14_PC_inhib <- data.frame(res.pca$ind$contrib)
#ex14_PC_inhib

#ggplot(ex14_PC_inhib, aes(x=PC1, y=PC2,label=rownames(ex14_PC_inhib)))+
#  geom_point()+geom_text(nudge_x =1.5)+theme_classic()

#fviz_contrib(res.pca, choice = "var", axes = 1, top = 100) # visual of percent contribution for PC1 


#### histogram of dimensions
ggplot(ex14_PCs, aes(x=Dim.1)) + geom_histogram()
ggplot(ex14_PCs, aes(x=Dim.2)) + geom_histogram()
ggplot(ex14_PCs, aes(x=Dim.3)) + geom_histogram()

ggplot()+
  geom_histogram(aes(x=ex14_PCs$Dim.1),color = "skyblue", alpha = 0.3, fill = "skyblue")+
  geom_histogram(aes(x=ex14_PCs$Dim.2), color = "firebrick", alpha = 0.3, fill = "firebrick")+
  geom_histogram(aes(x=ex14_PCs$Dim.3), color = "gray", alpha = 0.3, fill = "gray")+
  xlab("Variable Loadings")+
  ggtitle("PC ditributions")+
  theme_bw()


### heatmap 
met_wt_fasta = "Met_wt.fasta"

met_wt_fasta_con=file(met_wt_fasta, open="r")
met_wt_sequence = readLines(met_wt_fasta_con)[2]
close(met_wt_fasta_con)

order <- c('X','H', 'K', 'R', 'D', 'E',
           'C','M','N','Q','S','T','A','I','L','V','W','F',
           'Y','G','P')

variant_names <- c('Stop','H', 'K', 'R', 'D', 'E',
                   'C','M','N','Q','S','T','A','I','L','V','W','F',
                   'Y','G','P')

met_wt_1 = str_split(substr(met_wt_sequence, 1, 144), '')[[1]]
met_wt_2 = str_split(substr(met_wt_sequence, 145, 287), '')[[1]]

Met_row1 = ggplot(data = ex14_PCs%>%filter(pos %in% c(1059:1202)),
                  aes(x = pos, y = factor(mutation, level = order, labels = variant_names), fill = Dim.1)) +geom_tile() +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow')+
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

Met_row2 = ggplot(data = ex14_PCs %>% filter(pos %in% c(1203:1345)), 
                  aes(x = pos, y = factor(mutation, level = order, labels = variant_names), fill =  Dim.1)) +
  geom_tile() +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow')+
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
    legend.position="right"
  ) +
  labs(y = "Mutation", x = "Position")

PC1_heatmap = ggarrange(Met_row1,Met_row2,nrow = 2, ncol = 1)

plot(PC1_heatmap)

# individuals heatmap 
#ex14_PC_indiv <- data.frame(res.pca3$ind$coord)
pheatmap(t(ex14_PC2), 
         border_color = "black",
         #cluster_rows = FALSE,
         cutree_cols = 4,
         cutree_rows = 3,
         angle_col = 90,             # angle for column labels
         width = 5,                 # plot width in inches
         height = 7, 
         color = colorRampPalette(c("#E7B800","white","dodgerblue"))(100))

#-------------- Project PC space back onto structure 

#PC1
ex14_PC1 <- data.frame(pos=ex14_PCs$pos, PC1=ex14_PCs$Dim.1)
ex14_PC1_avg <- ex14_PC1 %>% group_by(pos) %>% dplyr::summarise(avg = mean(PC1))%>% ungroup()
x = map_scores_pdb(ex14_2WGJ_crizo , ex14_PC1_avg, "avg")
write.pdb(x, file="ex14_PC1.pdb")

#PC2
ex14_PC2 <- data.frame(pos=ex14_PCs$pos, PC2=ex14_PCs$Dim.2)
ex14_PC2_avg <- ex14_PC2 %>% group_by(pos) %>% dplyr::summarise(avg = mean(PC2))%>% ungroup()
x = map_scores_pdb(ex14_2WGJ_crizo , ex14_PC2_avg, "avg")
write.pdb(x, file="ex14_PC2.pdb")

#PC3
ex14_PC3 <- data.frame(pos=ex14_PCs$pos, PC3=ex14_PCs$Dim.3)
ex14_PC3_avg <- ex14_PC3 %>% group_by(pos) %>% dplyr::summarise(avg = mean(Dim.3))%>% ungroup()
x = map_scores_pdb(ex14_2WGJ_crizo , ex14_PC3_avg, "avg")
write.pdb(x, file="ex14_PC3.pdb")


#####################################################################################################
###-------------------------------------------- MCA -------------------------------------------------
#####################################################################################################
library("FactoMineR")
library("factoextra")

res.mca <- MCA(ex14_rosace_scores, graph = FALSE)

fviz_mca_biplot(res.mca, 
                repel = TRUE, # Avoid text overlapping (slow if many point)
                ggtheme = theme_minimal())

fviz_mca_var(res.mca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE, # Avoid text overlapping
             ggtheme = theme_minimal())


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


