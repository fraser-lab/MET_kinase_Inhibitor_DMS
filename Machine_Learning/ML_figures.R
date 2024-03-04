###################################################################################################################
# ML_model_figures.R
# by g.estevam @ UCSF 
# ML model developed by Ashraya and Karson @UCSF
###################################################################################################################


library(ggplot2)
library(readr)
library(stringr)
library(tidyverse)
library(dplyr)
library(ggridges)
library(tidyquant)
library(bio3d)
library(corrplot)
library(data.table)

library(tidyverse)
library(ggExtra)
library(ggpmisc)
library(ggbreak)



inhib_pearsonR_improvements <- read.csv("inhibitor_pearson_correlation_all_models.csv")
#avg_pearsonR_improvements <- read.csv("average_pearson_correlations_all_models.csv")

inhib_pearsonR_improvements2 <- read.csv("inhibitor_pearson_correlation_all_models_no_dmso_a458.csv")
baseline_r <- inhib_pearsonR_improvements %>% filter(model == "esm_baseline")

inhib_pearsonR_improvements2 <- rbind(inhib_pearsonR_improvements2,baseline_r)
inhib_pearsonR_improvements2 <- inhib_pearsonR_improvements2 %>% filter(inhibitor !="A458" & inhibitor != "DMSO")




custom_palette <- c("#ff9b00ff","#8cbe61ff", "#56B4E9", "#56B4E9", "#56B4E9", "#f088a5ff", "#f088a5ff", "#f088a5ff", "#f088a5ff","#e63466ff", "#e63466ff", "darkmagenta")


#---------------------------------- Feature improvements -------------------------------#

# Define the order for the legend
legend_order <- c("esm_baseline", "esm_only", "esm_ddg", "esm_ddg_dist", "esm_ddg_dist_vol","esm_ddg_dist_vol_type", 
                  "esm_ddg_all","esm_ddg_all_dist","esm_ddg_all_dist_vol", "esm_ddg_all_dist_vol_type")
#x_axis_order <- c("DMSO", "Crizo", "Camp", "Tepo", "Glu", "Savo", "NVP", "Cabo", "Mere", "Gle", "A458")
x_axis_order <- c("Crizo", "Camp", "Tepo", "Glu", "Savo", "NVP", "Cabo", "Mere", "Gle")

#filtered_data <- inhib_pearsonR_improvements %>% filter(inhibitor == "DMSO" & inhibitor == "A458")
#filtered_data <- inhib_pearsonR_improvements[!(inhib_pearsonR_improvements$model %in% c("esm_only", "esm_ddg_dist_vol_type")),]


# Grouped bar plot
corr <- ggplot(inhib_pearsonR_improvements2, aes(factor(inhibitor, levels = x_axis_order), correlation, fill = factor(model, levels = legend_order))) +
  geom_bar(colour="black",stat = "identity", position = 'dodge') +
  scale_fill_manual(values = c("grey","#CAF0F8","#ADE8F4", "#90E0EF","#48CAE4","#00B4D8","#0096C7", "#0077B6","#023E8A", "#03045E")) +
  labs(fill = "Model") +  # Optional: Rename the legend,
  geom_hline(yintercept = 0.5, size = 0.9, color = "red", linetype = "dashed")+
  theme_classic()+
  #theme(text = element_text(size=5))+
  theme(axis.text.x=element_text(colour="black"))+
  theme(legend.position="bottom")

print(corr)


#---------------------------------- Scatter plots with marginal dist -------------------------------#
#esm_all_features <- read.csv("ROSACE_vs_best_model_fitness_predictions.csv")
#esm_baseline <- read.csv("ROSACE_vs_esm_baseline_fitness_predictions.csv")

#esm_baseline <- read.csv("model_wise_fitness_predictions.csv")
#esm_baseline <- esm_baseline%>% filter(model == "esm_only")

#esm_all_features <- read.csv("model_wise_fitness_predictions.csv")
#esm_all_features <- esm_all_features%>% filter(model == "esm_ddg_dist_vol")

esm_all_features <- read.csv("model_wise_fitness_predictions_no_dmso_a458.csv")
esm_all_features2 <- esm_all_features%>% filter(inhibitor == "Crizo")
esm_all_features <- esm_all_features%>% filter(model == "esm_ddg_all_dist_vol_type")

esm_baseline <- read.csv("ROSACE_vs_esm_baseline_fitness_predictions.csv")
esm_baseline <- esm_baseline %>% filter(inhibitor != "DMSO" & inhibitor != "A458")

# ESM baseline:

pearson_r <- cor(esm_baseline$ROSACE_fitness, esm_baseline$ESM_baseline)
#mse <- mean((esm_baseline$ROSACE_fitness - esm_baseline$ESM_baseline_fitness)^2)

# Create scatter plot with Pearson R fit and MSE fit
scatter_esm_features <- ggplot(esm_baseline, aes(x = ROSACE_fitness, y = ESM_baseline_fitness)) +
  geom_point(color = "white", fill = "#00B4D8", shape = 21, size = 2, alpha = 0.7) +
  
  geom_point( aes(x=1.4615671, y= -17.858 ), fill = "#03045E")+ # D1228R 
  geom_point( aes(x=1.5474515, y = -14.900), fill = "#03045E")+ # G1163D
  geom_point( aes(x=-0.690337342, y= 0.137), color = "#03045E")+ # M1229I
  geom_point( aes(x=-2.490182, y=-0.606), fill = "#03045E")+ # 1294H
  
  geom_smooth(method = "lm", color = "red", se = FALSE, aes(group = 1, linetype = "Pearson R")) +  # Pearson R fit
  theme(legend.position = "none") +
  xlab("Experimental Fitness") +
  ylab("ESM LLR") +
  theme_classic()+
  theme(legend.position="none")

# Print the Pearson R coefficient and MSE
print(paste("Pearson R coefficient:", round(pearson_r, 2)))
#print(paste("Mean Squared Error (MSE):", round(mse, 2)))

# Plot the scatter plot with Pearson R fit and MSE fit
ggMarginal(scatter_esm_features, type="histogram",col = "black", fill = "#00B4D8",alpha=0.7)


#ESM + all features:

pearson_r <- cor(esm_all_features$ROSACE_fitness, esm_all_features$predicted_fitness) #0.8133086
mse <- mean((esm_all_features$ROSACE_fitness - esm_all_features$predicted_fitness)^2) #0.27

# Create scatter plot with Pearson R fit and MSE fit
scatter_esm_features2 <- ggplot(esm_all_features, aes(x = ROSACE_fitness, y = predicted_fitness)) +
  geom_point(color = "white", fill = "#00B4D8", shape = 21, size = 2, alpha = 0.7) +
  
  geom_point( aes(x=1.4615671, y= -1.157978), color = "#03045E")+ # D1228R, predicted, 170R
  geom_point( aes(x=1.5474515, y = 0.26009680), fill = "#03045E")+ # G1163E, predicted, 105E
  geom_point( aes(x=-0.690337342, y= -0.7449126), color = "#03045E")+ # M1229I, predicted, 171
  geom_point( aes(x=-2.49018193, y=-0.9577531), color = "#03045E")+ # 1294H, predicted, 236
  
  geom_smooth(method = "lm", color = "red", se = FALSE, aes(group = 1, linetype = "Pearson R")) +  # Pearson R fit
  theme(legend.position = "none") +
  xlab("Experimental Fitness") +
  ylab("Predicted Fitness") +
  theme_classic() +
  theme(legend.position = "none")

# Print the Pearson R coefficient and MSE
print(paste("Pearson R coefficient:", round(pearson_r, 2)))
#print(paste("Mean Squared Error (MSE):", round(mse, 2)))

# Plot the scatter plot with Pearson R fit and MSE fit
ggMarginal(scatter_esm_features2, type="histogram", col = "black", fill = "#00B4D8", alpha=0.7)


#--------------------------- Map mutations that gain prediction -------------------------------#

merge_esm <- merge(esm_baseline, esm_all_features) # merge dfs to compare mutations

merge_esm$delta_baseline_rosace <- (merge_esm$ESM_baseline_fitness) - (merge_esm$ROSACE_fitness) # find mutations with the largest delta 
merge_esm$delta_predicted_rosace <- (merge_esm$predicted_fitness) - (merge_esm$ROSACE_fitness) # find mutations with the largest delta 
#merge_esm <- merge_esm %>% filter(model == "esm_ddg_dist_vol")

# mutations with largest delta from ESM to prediction: 
## crizo:
### D1228R, 1228K
### L1195I
# 1294S

poor_well <- merge_esm %>% filter(ROSACE_fitness <=-1 & ESM_baseline_fitness >= -3) 

#---------------------- Step charts showing Feature improvements for residues  ----------------------#

esm_all_features2 <- read.csv("model_wise_fitness_predictions.csv")

# restructure baseline df
esm_baseline$model <- 'esm_baseline'
names(esm_baseline)[names(esm_baseline) == 'ESM_baseline_fitness'] <- 'predicted_fitness'

merge_baseline <- rbind(esm_baseline,esm_all_features2)
merge_baseline <- merge_baseline[!grepl("esm_baseline", merge_baseline$model),]
merge_baseline <- merge_baseline[!(merge_baseline$model %in% c("esm_baseline")),]


D1228R <- merge_baseline %>% filter(mutation == "170R" & inhibitor=="Crizo")
D1228K <- merge_baseline%>% filter(mutation == "170K" & inhibitor=="Crizo")
V1294S <- merge_baseline %>% filter(mutation == "236S" & inhibitor=="Crizo")

V1294H <- merge_baseline%>% filter(mutation == "236H" & inhibitor=="Crizo")

G1163E <- merge_baseline %>% filter(mutation == "105E" & inhibitor=="Crizo")
V1063L <- merge_baseline %>% filter(mutation == "5L" & inhibitor=="Crizo")
N1288<- merge_baseline %>% filter(mutation == "5L" & inhibitor=="Crizo")
M1229I<- merge_baseline %>% filter(mutation == "171I" & inhibitor=="Crizo")

custom_labels <- c("esm_only"= "XGBoost" ,"esm_ddg" = "∆∆∆G", "esm_ddg_dist" = "Dist", "esm_ddg_dist_vol" = "Vol")

D1228R$model <- factor(D1228R$model, levels = c("esm_only", "esm_ddg", "esm_ddg_dist", "esm_ddg_dist_vol"))
D1228R_step <- ggplot(D1228R, aes(x = model, y = predicted_fitness, group = interaction(mutation, inhibitor))) +
  geom_step(size = 0.5, color = "#390099", linetype = "dashed") +
  geom_point(size = 4, color = "#390099") +
  ylab("Predicted fitness") +
  xlab("Model") +
  ggtitle("D1228R_Crizotinib")+
  theme_classic() +
  ylim(-0.5,2)+
  scale_x_discrete(labels = custom_labels)+
  theme(axis.text = element_text(size = 12),  
        axis.title = element_text(size = 14))
plot(D1228R_step)


V1294H$model <- factor(V1294H$model, levels = c("esm_only", "esm_ddg", "esm_ddg_dist", "esm_ddg_dist_vol"))
V1294H_step <- ggplot(V1294H, aes(x = model, y = predicted_fitness, group = interaction(mutation, inhibitor))) +
  geom_step(size = 0.5, color = "#390099", linetype = "dashed") +
  geom_point(size = 4, color = "#390099") +
  geom_hline(yintercept =-2.49018193, size =1,color = "#e63466ff", linetype ="dotted")+ 
  ylab("Predicted fitness") +
  xlab("Model") +
  ggtitle("V1294H_Crizotinib")+
  theme_classic() +
  ylim(-2.5,-2.1)+
  scale_x_discrete(labels = custom_labels)+
  theme(axis.text = element_text(size = 12),  
        axis.title = element_text(size = 14))
print(V1294H_step)


V1294S$model <- factor(V1294S$model, levels = c("esm_only","esm_ddg", "esm_ddg_dist", "esm_ddg_dist_vol"))
V1294S_step <- ggplot(V1294S, aes(x = model, y = predicted_fitness, group = interaction(mutation, inhibitor))) +
  geom_step(size = 0.5, color = "#390099", linetype = "dashed") +
  geom_point(size = 4, color = "#390099") +
  ylab("Predicted fitness") +
  xlab("Model") +
  ggtitle("V1294S_Crizotinib")+
  theme_classic() +
  ylim(-2.35,-2.2)+
  scale_x_discrete(labels = custom_labels)+
  theme(axis.text = element_text(size = 12),  
        axis.title = element_text(size = 14))
print(V1294S_step)


G1163E$model <- factor(G1163E$model, levels = c("esm_only", "esm_ddg", "esm_ddg_dist", "esm_ddg_dist_vol"))
G1163E_step <- ggplot(G1163E, aes(x = model, y = predicted_fitness, group = interaction(mutation, inhibitor))) +
  geom_step(size = 0.5, color = "#390099", linetype = "dashed") +
  geom_point(size = 4, color = "#390099") +
  ylab("Predicted fitness") +
  xlab("Model") +
  ggtitle("G1163E_Crizotinib")+
  theme_classic() +
  # ylim(-0.5,1)+
  scale_x_discrete(labels = custom_labels)+
  theme(axis.text = element_text(size = 12),  
        axis.title = element_text(size = 14))
print(G1163E_step)


M1229I$model <- factor(M1229I$model, levels = c("esm_only", "esm_ddg", "esm_ddg_dist", "esm_ddg_dist_vol"))
M1229I_step <- ggplot(M1229I, aes(x = model, y = predicted_fitness, group = interaction(mutation, inhibitor))) +
  geom_step(size = 0.5, color = "#390099", linetype = "dashed") +
  geom_point(size = 4, color = "#390099") +
  ylab("Predicted fitness") +
  xlab("Model") +
  ggtitle("M1229I_Crizotinib")+
  theme_classic() +
  # ylim(-0.5,1)+
  scale_x_discrete(labels = custom_labels)+
  theme(axis.text = element_text(size = 12),  
        axis.title = element_text(size = 14))
print(M1229I_step)


