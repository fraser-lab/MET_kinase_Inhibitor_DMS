###############################################################################################
# MET_inhibitor_DMS_replicates_corr.R
# code by g.estevam @ ucsf
# this code plots the correlations between conditions and time points for the MET inhibitor DMS 
###############################################################################################

library(dplyr)
library(ggplot2)
library(gridExtra)

########################################################################################################
#---------------------------------------- Replicate validation  ---------------------------------------#

ex14_inhib_scores <- read.csv("ex14_inhib_scores.csv")#+Ex14
met_inhib_scores <- read.csv("met_inhib_scores.csv")#-Ex14 


condition_names <- c("DMSO", "Crizo", "Camp", "Tepo", "Savo", "Glu", "NVP", "Mere", "Cabo", "Gle", "A458", "Tiv")

ex14_plots_list <- list()

# MET replicates
for (condition in condition_names) {
  rep1 <- as.data.frame(fread(paste0("pathname", condition, "_R1_sel/main_identifiers_scores.tsv")))
  rep2 <- as.data.frame(fread(paste0("pathname", condition, "_R2_sel/main_identifiers_scores.tsv")))
  rep3 <- as.data.frame(fread(paste0("pathname", condition, "_R3_sel/main_identifiers_scores.tsv")))
  
  rep1_df <- data.frame(hgvs = rep1$hgvs, R1 = rep1$score)
  rep1_df <- merge(rep1_df, ex14_inhib_scores)
  
  rep2_df <- data.frame(hgvs = rep2$hgvs, R2 = rep2$score)
  rep2_df <- merge(rep2_df, ex14_inhib_scores)
  
  rep3_df <- data.frame(hgvs = rep3$hgvs, R3 = rep3$score)
  rep3_df <- merge(rep3_df, ex14_inhib_scores)
  
  rep1_avg <- rep1_df %>% group_by(pos) %>% dplyr::summarise(R1 = mean(R1)) %>% ungroup
  rep2_avg <- rep2_df %>% group_by(pos) %>% dplyr::summarise(R2 = mean(R2)) %>% ungroup
  rep3_avg <- rep3_df %>% group_by(pos) %>% dplyr::summarise(R3 = mean(R3)) %>% ungroup
  
  rep12_df <- merge(rep1_avg, rep2_avg)
  rep13_df <- merge(rep1_avg, rep3_avg)
  rep23_df <- merge(rep2_avg, rep3_avg)
  
  corr_rep12_df <- cor(rep12_df$R1, rep12_df$R2, method = "pearson")
  corr_rep13_df <- cor(rep13_df$R1, rep13_df$R3, method = "pearson")
  corr_rep23_df <- cor(rep23_df$R2, rep13_df$R3, method = "pearson")
  
  plot_rep1_rep3 <- ggplot(data = rep13_df, aes(x = R1, y = R3)) +
    ggtitle(paste0(condition, " R1 vs R3")) +
    xlab("Fitness score R1") +
    ylab("Fitness score R3") +
    xlim(-10, 3) +
    ylim(-10, 3) +
    theme_classic() +
    geom_point()+
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    annotate("text", x = -8, y = 0, label = paste("r =", round(corr_rep13_df, 2)), size = 2, color = "black")
  
  
  plot_rep2_rep3 <- ggplot(data = rep23_df, aes(x = R2, y = R3)) +
    ggtitle(paste0(condition, " R2 vs R3")) +
    xlab("Fitness score R2") +
    ylab("Fitness score R3") +
    xlim(-10, 3) +
    ylim(-10, 3) +
    theme_classic() +
    geom_point()+
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    annotate("text", x = -8, y = 0, label = paste("r =", round(corr_rep23_df, 2)), size = 2, color = "black")
  
  # Store each plot in the list
  ex14_plots_list[[paste0(condition, "_R1_R3")]] <- plot_rep1_rep3
  ex14_plots_list[[paste0(condition, "_R2_R3")]] <- plot_rep2_rep3
}

pdf("MET_replicates_plot.pdf", width = 8, height = 10.5)
grid.arrange(grobs = ex14_plots_list, ncol = 4)
dev.off()



met_del_ex14_plots_list <- list()

# MET del Ex14 replicates
for (condition in condition_names) {
  rep1 <- as.data.frame(fread(paste0("pathname", condition, "_R1_sel/main_identifiers_scores.tsv")))
  rep2 <- as.data.frame(fread(paste0("pathname", condition, "_R2_sel/main_identifiers_scores.tsv")))
  rep3 <- as.data.frame(fread(paste0("pathname", condition, "_R3_sel/main_identifiers_scores.tsv")))
  
  rep1_df <- data.frame(hgvs = rep1$hgvs, R1 = rep1$score)
  rep1_df <- merge(rep1_df, met_inhib_scores)
  
  rep2_df <- data.frame(hgvs = rep2$hgvs, R2 = rep2$score)
  rep2_df <- merge(rep2_df, met_inhib_scores)
  
  rep3_df <- data.frame(hgvs = rep3$hgvs, R3 = rep3$score)
  rep3_df <- merge(rep3_df, met_inhib_scores)
  
  rep1_avg <- rep1_df %>% group_by(pos) %>% dplyr::summarise(R1 = mean(R1)) %>% ungroup
  rep2_avg <- rep2_df %>% group_by(pos) %>% dplyr::summarise(R2 = mean(R2)) %>% ungroup
  rep3_avg <- rep3_df %>% group_by(pos) %>% dplyr::summarise(R3 = mean(R3)) %>% ungroup
  
  rep12_df <- merge(rep1_avg, rep2_avg)
  rep13_df <- merge(rep1_avg, rep3_avg)
  rep23_df <- merge(rep2_avg, rep3_avg)
  
  corr_rep12_df <- cor(rep12_df$R1, rep12_df$R2, method = "pearson")
  corr_rep13_df <- cor(rep13_df$R1, rep13_df$R3, method = "pearson")
  corr_rep23_df <- cor(rep23_df$R2, rep13_df$R3, method = "pearson")
  
  plot_rep1_rep3 <- ggplot(data = rep13_df, aes(x = R1, y = R3)) +
    ggtitle(paste0(condition, " R1 vs R3")) +
    xlab("Fitness score R1") +
    ylab("Fitness score R3") +
    xlim(-10, 3) +
    ylim(-10, 3) +
    theme_classic() +
    geom_point()+
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    annotate("text", x = -8, y = 0, label = paste("r =", round(corr_rep13_df, 2)), size = 2, color = "black")
  
  
  plot_rep2_rep3 <- ggplot(data = rep23_df, aes(x = R2, y = R3)) +
    ggtitle(paste0(condition, " R2 vs R3")) +
    xlab("Fitness score R2") +
    ylab("Fitness score R3") +
    xlim(-10, 3) +
    ylim(-10, 3) +
    theme_classic() +
    geom_point()+
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    annotate("text", x = -8, y = 0, label = paste("r =", round(corr_rep23_df, 2)), size = 2, color = "black")
  
  # Store each plot in the list
  met_del_ex14_plots_list[[paste0(condition, "_R1_R2")]] <- plot_rep1_rep3
  met_del_ex14_plots_list[[paste0(condition, "_R2_R3")]] <- plot_rep2_rep3
}

pdf("MET_delEx14_replicates_plot.pdf", width = 8, height = 10.5)
grid.arrange(grobs = met_del_ex14_plots_list, ncol = 4)
dev.off()


