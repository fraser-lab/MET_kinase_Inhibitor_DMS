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
library(cowplot)


inhib_pearsonR_improvements <- read.csv("inhibitor_correlations_all_submodels_final.csv")
#inhib_pearsonR_improvements <- read.csv("inhibitor_pearson_correlation_all_models.csv")
#avg_pearsonR_improvements <- read.csv("average_pearson_correlations_all_models.csv")

inhib_pearsonR_improvements2 <- read.csv("inhibitor_correlations_all_submodels_final.csv")
#baseline_r <- inhib_pearsonR_improvements2 %>% filter(model == "esm_baseline")

inhib_pearsonR_improvements2 <- rbind(inhib_pearsonR_improvements2,baseline_r)
inhib_pearsonR_improvements2 <- inhib_pearsonR_improvements2 %>% filter(inhibitor !="A458" & inhibitor != "DMSO")




custom_palette <- c("#ff9b00ff","#8cbe61ff", "#56B4E9", "#56B4E9", "#56B4E9", "#f088a5ff", "#f088a5ff", "#f088a5ff", "#f088a5ff","#e63466ff", "#e63466ff", "darkmagenta")


#---------------------------------- Feature improvements -------------------------------#

# Define the order for the legend
#legend_order <- c("esm_baseline", "esm_only", "esm_ddg", "esm_ddg_dist", "esm_ddg_dist_vol","esm_ddg_dist_vol_type", "esm_ddg_all","esm_ddg_all_dist","esm_ddg_all_dist_vol", "esm_ddg_all_dist_vol_type")
#x_axis_order <- c("DMSO", "Crizo", "Camp", "Tepo", "Glu", "Savo", "NVP", "Cabo", "Mere", "Gle", "A458")

legend_order <- c("esm", 
                  "esm_inhibitor_weight_ligand_rmsd",
                  "esm_dddg_ddg_all",
                  "esm_atp_distance_inhibitor_distance",
                  "esm_crystal_rmsf_residue_rmsd",
                  "esm_dddg_ddg_all_atp_distance_inhibitor_weight_inhibitor_distance_crystal_rmsf_residue_rmsd_ligand_rmsd")

custom_labels <- c(
  'esm' = 'ESM LLR',
  'esm_inhibitor_weight_ligand_rmsd' = 'Inhibitor',
  'esm_dddg_ddg_all' = 'Stability',
  'esm_atp_distance_inhibitor_distance' = 'Distance',
  'esm_crystal_rmsf_residue_rmsd' = 'Conformation',
  'esm_dddg_ddg_all_atp_distance_inhibitor_weight_inhibitor_distance_crystal_rmsf_residue_rmsd_ligand_rmsd' = 'All'
)



x_axis_order <- c("Crizo", "Camp", "Tepo", "Glu", "Savo", "NVP", "Cabo", "Mere", "Gle")

#filtered_data <- inhib_pearsonR_improvements %>% filter(inhibitor == "DMSO" & inhibitor == "A458")
#filtered_data <- inhib_pearsonR_improvements[!(inhib_pearsonR_improvements$model %in% c("esm_only", "esm_ddg_dist_vol_type")),]


# Grouped bar plot
corr <- ggplot(inhib_pearsonR_improvements %>%  filter (model != "esm_baseline" & 
                                                        model != "esm_dddg_ddg_all_atp_distance_inhibitor_distance" &
                                                        model != "esm_dddg_ddg_all_atp_distance_inhibitor_distance_crystal_rmsf_residue_rmsd"), 
               aes(factor(inhibitor, levels = x_axis_order), correlation, fill = factor(model, levels = legend_order))) +
  
  geom_bar(colour="black",stat = "identity", position = 'dodge') +
  scale_fill_manual(values = c("grey","#CAF0F8", "#90E0EF","#00B4D8","#0096C7", "#0077B6","#023E8A", "#03045E")) +
  labs(fill = "model") +  # Optional: Rename the legend,
  geom_hline(yintercept = 0.5, size = 0.9, color = "red", linetype = "dashed")+
  ylab("Correlation")+
  theme_classic()+
  #theme(text = element_text(size=5))+
  theme(axis.text.x=element_text(colour="black"))+
  theme(legend.position="none")

print(corr)

#---------------------------------- cross validation line plot -------------------------------#

cross_val<- read.csv("cross_val_lineplot.csv")

data <- data.frame(
  Features = c("ESM", "Best model"),
  Train.correlation = c(0.2837153, 0.4152881),
  Validation.correlation = c(0.2633262, 0.3072786),
  Test.correlation = c(0.2506842, 0.3233952)
)

data$Features <- factor(data$Features, levels = c("ESM", "Best model"))
data_long <- gather(data, key = "Dataset", value = "Correlation", -Features)

ggplot(data_long, aes(x = Features, y = Correlation, group = Dataset, color = Dataset)) +
  geom_line(size = 1) +
  geom_point() +
  ylim(0,0.5)+
  geom_hline(yintercept = 0.5, size = 0.9, color = "red", linetype = "dashed")+
  scale_color_manual(values = c("#00B4D8","grey", "#0077B6"), 
                     labels = c("Test", "Train", "Validation")) +  # Custom colors
  theme_classic() +
  labs(title = "Cross Validation",
       x = "Features",
       y = "Correlation")

#---------------------------------- Scatter plots with marginal dist -------------------------------#


esm_all_features <- read.csv("model_wise_fitness_predictions_submodels_new (1).csv")
esm_all_features <- esm_all_features%>% filter(model =="esm_dddg_ddg_all_atp_distance_inhibitor_weight_inhibitor_distance_crystal_rmsf_residue_rmsd_ligand_rmsd")

esm_all_features <- esm_all_features %>%
  mutate(new_column = as.numeric(gsub("\\D", "", mutation)) + 1058)


#esm_baseline <- read.csv("ROSACE_vs_esm_baseline_fitness_predictions.csv")
esm_baseline <- read.csv("ROSACE_vs_esm_baseline_fitness_predictions (4).csv")
esm_baseline <- esm_baseline %>% filter(inhibitor != "DMSO" & inhibitor != "A458")

# ESM baseline:

pearson_r <- cor(esm_baseline$ROSACE_fitness, esm_baseline$ESM_baseline)
#mse <- mean((esm_baseline$ROSACE_fitness - esm_baseline$ESM_baseline_fitness)^2)

# Create scatter plot with Pearson R fit and MSE fit
scatter_esm_features <- ggplot(esm_baseline, aes(x = ROSACE_fitness, y = ESM_baseline_fitness)) +

  geom_point(color = "white", fill = "#00B4D8", shape = 21, size = 2, alpha = 0.7) +
  #geom_smooth(method = "lm", color = "red", se = FALSE, aes(group = 1, linetype = "Pearson R")) +  # Pearson R fit
  #geom_abline(intercept =, slope = 1, color = "red", linetype= "dashed", size=1 )+
  
  geom_point( aes(x=-2.5430146, y= -19.395), fill = "#03045E")+ # G1090D, conserved 
  geom_point( aes(x=1.438725589, y= -11.750), fill = "#03045E")+ # Y1230D, inhibitor contact 
  geom_point( aes(x= -0.2681378, y= -7.935), fill = "#03045E")+ # N1167K
 
  
  theme(legend.position = "none") +
  xlab("Experimental Fitness") +
  ylab("ESM") +
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
 # geom_smooth(method = "lm", color = "red", se = FALSE, aes(group = 1, linetype = "Pearson R")) +  # Pearson R fit
  geom_abline(intercept = -0.5, slope = 0.6, color = "red", linetype= "dashed", size=1 )+
  
  geom_point( aes(x=-2.5430146, y= -2.3182776), fill = "#03045E")+ # G1090D, conserved 
  geom_point( aes(x=1.438725589, y= 0.47242300), fill = "#03045E")+ # Y1230D
  geom_point( aes(x= -0.2681378, y= -0.3160836), fill = "#03045E")+ # N1167K
  
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



poor_well <- merge_esm %>% filter(inhibitor == "Crizo")

poor_well <- merge_esm %>% filter(inhibitor == "Crizo" &
                                  delta_predicted_rosace < -0.2
                                  ) 



#---------------------- Step charts showing Feature improvements for residues  ----------------------#

#esm_all_features2 <- read.csv("model_wise_fitness_predictions (1).csv")
merge_baseline<- read.csv("step_plot_df.csv")

G1090D <- merge_baseline %>% filter(mutation == "32D")
Y1230D <- merge_baseline %>% filter(mutation == "172D")
N1167K <- merge_baseline %>% filter(mutation == "109K")


custom_labels <- c(
  'esm' = 'ESM',
  'esm_dddg_ddg_all' = 'Stability',
  'esm_dddg_ddg_all_atp_distance_inhibitor_distance' = 'Dist.',
  'esm_dddg_ddg_all_atp_distance_inhibitor_distance_crystal_rmsf_residue_rmsd' = 'Conf.',
  'esm_dddg_ddg_all_atp_distance_inhibitor_weight_inhibitor_distance_crystal_rmsf_residue_rmsd_ligand_rmsd' = 'All'
)

G1090D_step <- ggplot(G1090D, aes(x = factor(model), y = predicted_fitness, group = interaction(mutation))) +
  geom_step(size = 0.5, color = "#390099", linetype = "dashed") +
  geom_point(size = 2, color = "#390099") +
  geom_hline(yintercept = G1090D$ROSACE_fitness, color = "red")+
  ylab("Predicted fitness") +
  xlab("Added Feature") +
  ggtitle("G1090D_Crizotinib") +
  theme_classic() +
  #ylim(-0.5,2)+
  scale_x_discrete(labels = custom_labels) +  # Use custom labels for all x-axis values
  theme(axis.text = element_text(size = 7),  
        axis.title = element_text(size = 7), 
        title = element_text(size = 8))
plot(G1090D_step)

Y1230D_step <- ggplot(Y1230D, aes(x = factor(model), y = predicted_fitness, group = interaction(mutation))) +
  geom_step(size = 0.5, color = "#390099", linetype = "dashed") +
  geom_point(size = 2, color = "#390099") +
  geom_hline(yintercept = Y1230D$ROSACE_fitness, color = "red")+
  ylab("Predicted fitness") +
  xlab("Added Feature") +
  ggtitle("Y1230D_Crizotinib") +
  theme_classic() +
  #ylim(-0.5,2)+
  scale_x_discrete(labels = custom_labels) +  # Use custom labels for all x-axis values
  theme(axis.text = element_text(size = 7),  
        axis.title = element_text(size = 7), 
        title = element_text(size = 8))
plot(Y1230D_step)

N1167K_step <- ggplot(N1167K, aes(x = factor(model), y = predicted_fitness, group = interaction(mutation))) +
  geom_step(size = 0.5, color = "#390099", linetype = "dashed") +
  geom_point(size = 2, color = "#390099") +
  geom_hline(yintercept = N1167K$ROSACE_fitness, color = "red")+
  ylab("Predicted fitness") +
  xlab("Added Feature") +
  ggtitle("N1167K_Crizotinib") +
  theme_classic() +
  #ylim(-0.5,2)+
  scale_x_discrete(labels = custom_labels) +  # Use custom labels for all x-axis values
  theme(axis.text = element_text(size = 7),  
        axis.title = element_text(size = 7), 
        title = element_text(size = 8))
plot(N1167K_step)


plot_grid(Y1230D_step, N1167K_step, G1090D_step,ncol =1, nrow= 3)


#----------- 

# plots average of distance between the CA atom and every inhibitor atom


dist_diff <- read.csv("predicted_fitness_scores_vs_ca_distance.csv")
lowest_dist <- dist_diff %>% filter (ca_distance < 11 & inhibitor == "Crizo")
lowest_dist$delta_gxboost_best <- (lowest_dist$predicted_fitness_esm_xgboost - lowest_dist$predicted_fitness_best_model)
lowest_dist$delta_rosace_best <- (lowest_dist$ROSACE_fitness_x- lowest_dist$predicted_fitness_best_model)
lowest_dist$abs_delta <- abs(lowest_dist$delta_gxboost_best- lowest_dist$delta_rosace_best)
lowest_dist$delta_best_rosase <- abs(lowest_dist$ROSACE_fitness_x - lowest_dist$predicted_fitness_best_model)


lowest_dist_filtered <- lowest_dist %>% filter (delta_best_rosase > -0.2 & delta_best_rosase < 0.2 )




