# MET_inhibitor_heatmaps.R
# code to generate alll the inhibitor heatmaps from the normalized roscae scores 
# by g.estevam @ ucsf


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
source("dms_analysis_utilities.R")
library("data.table")
library(scales)
library(colorspace)
library(ggpubr)
library(rstatix)


##-------------------------------------- Load Rosace Fitness scores ----------------------------------------------##
### open data file as data frame for all raw inhibitor ALL scores generated from ROSACE

# wt 
ex14_rosace_norm_scores_filtered <- as.data.frame(fread("ex14_scores_filtered.tsv")) # here ex14 means MET+Ex14 
ex14_rosace_norm_scores_filtered  <- ex14_rosace_norm_scores_filtered  %>% mutate(position = position + 1058)
ex14_rosace_norm_scores_filtered  <- ex14_rosace_norm_scores_filtered %>% mutate(pos = position) 
ex14_rosace_norm_scores_filtered <- ex14_rosace_norm_scores_filtered %>% group_by(variants) %>% mutate(mean = (mean - mean[match("DMSO", key)])) %>% ungroup()

# wt 
ex14_rosace_norm_scores <- as.data.frame(fread("ex14_scores_unfiltered.tsv")) # here ex14 means MET+Ex14 
ex14_rosace_norm_scores  <- ex14_rosace_norm_scores  %>% mutate(position = position + 1058)
ex14_rosace_norm_scores <- ex14_rosace_norm_scores %>% mutate(pos = position) 
ex14_rosace_norm_scores_filtered <- ex14_rosace_norm_scores %>%group_by(variants) %>%mutate(mean = mean - mean[key == "DMSO"]) %>%ungroup()


# delta exon14 
met_rosace_norm_scores_filtered <- as.data.frame(fread("met_scores_filtered.tsv")) # here met means MET+met 
met_rosace_norm_scores_filtered  <- met_rosace_norm_scores_filtered  %>% mutate(position = position + 1058)
met_rosace_norm_scores_filtered  <- met_rosace_norm_scores_filtered %>% mutate(pos = position) 
met_rosace_norm_scores_filtered <- met_rosace_norm_scores_filtered %>% group_by(variants) %>%mutate(mean = mean - mean[key == "DMSO"][1]) %>%ungroup()

# delta exon14 
met_rosace_norm_scores <- as.data.frame(fread("met_scores_unfiltered.tsv")) # here met means MET+met 
met_rosace_norm_scores  <- met_rosace_norm_scores  %>% mutate(position = position + 1058)
met_rosace_norm_scores <- met_rosace_norm_scores %>% mutate(pos = position) 
met_rosace_norm_scores_filtered <- met_rosace_norm_scores_filtered %>% group_by(variants) %>%mutate(mean = mean - mean[key == "DMSO"][1]) %>%ungroup()

##-------------------------------------- make df wide----------------------------------------------##
ex14_rosace_scores_condensed <- data.frame (hgvs = ex14_rosace_norm_scores_filtered$variants,
                                            pos =ex14_rosace_norm_scores_filtered$position, 
                                            mutation = ex14_rosace_norm_scores_filtered$mutation, 
                                            score = ex14_rosace_norm_scores_filtered$mean,
                                            type = ex14_rosace_norm_scores_filtered$type,
                                            inhib= ex14_rosace_norm_scores_filtered$key, 
                                            pval1.5= ex14_rosace_norm_scores_filtered$pval1.5,
                                            pval0=ex14_rosace_norm_scores_filtered$pval0, 
                                            mean = ex14_rosace_norm_scores$mean)

ex14_rosace_scores_wide <- ex14_rosace_scores_condensed %>% pivot_wider(names_from = inhib, values_from = score) 

ex14_rosace_scores_crizo <- ex14_rosace_scores_condensed  %>% filter(inhib == "Crizo")


met_rosace_scores_condensed <- data.frame (hgvs = met_rosace_norm_scores_filtered$variants,
                                            pos =met_rosace_norm_scores_filtered$position, 
                                            mutation = met_rosace_norm_scores_filtered$mutation, 
                                            score = met_rosace_norm_scores_filtered$mean,
                                            type = met_rosace_norm_scores_filtered$type,
                                            inhib= met_rosace_norm_scores_filtered$key, 
                                            pval1.5= met_rosace_norm_scores_filtered$pval1.5,
                                            pval0=met_rosace_norm_scores_filtered$pval0, 
                                            mean = met_rosace_norm_scores$mean)

met_rosace_scores_wide <- met_rosace_scores_condensed %>% pivot_wider(names_from = inhib, values_from = score) 

met_rosace_scores_crizo <- met_rosace_scores_condensed  %>% filter(inhib == "Crizo")



# Assuming met_rosace_norm_scores_filtered is your primary data frame
# Check the number of rows in each vector
n_rows_primary <- nrow(met_rosace_norm_scores_filtered)

# Assuming met_rosace_norm_scores is the secondary data frame
# Check the number of rows in the secondary data frame
n_rows_secondary <- nrow(met_rosace_norm_scores)

# Identify the minimum number of rows between the two data frames
min_rows <- min(n_rows_primary, n_rows_secondary)

# Select the first min_rows rows from each vector in met_rosace_norm_scores_filtered
met_rosace_scores_condensed <- data.frame(
  hgvs = met_rosace_norm_scores_filtered$variants[1:min_rows],
  pos = met_rosace_norm_scores_filtered$position[1:min_rows],
  mutation = met_rosace_norm_scores_filtered$mutation[1:min_rows],
  score = met_rosace_norm_scores_filtered$mean[1:min_rows],
  type = met_rosace_norm_scores_filtered$type[1:min_rows],
  inhib = met_rosace_norm_scores_filtered$key[1:min_rows],
  pval1.5 = met_rosace_norm_scores_filtered$pval1.5[1:min_rows],
  pval0 = met_rosace_norm_scores_filtered$pval0[1:min_rows],
  mean = met_rosace_norm_scores$mean[1:min_rows]
)

# Print the number of rows in the resulting data frame
cat("Number of rows in met_rosace_scores_condensed: ", nrow(met_rosace_scores_condensed), "\n")


##-------------------------------------- MET Heatmap fucntion ----------------------------------------------##

# labels and order for the heatmaps 
order <- c('X','H', 'K', 'R', 'D', 'E','C','M','N','Q','S','T','A','I','L','V','W','F', 'Y','G','P')
variant_names <- c('Stop','H', 'K', 'R', 'D', 'E','C','M','N','Q','S','T','A','I','L','V','W','F','Y','G','P')
met_wt_fasta = "Met_wt2.fasta"

met_wt_fasta_con=file(met_wt_fasta, open="r")
met_wt_sequence = readLines(met_wt_fasta_con)[2]
close(met_wt_fasta_con) 

ex14_rosace_scores_crizo $is.wt = ex14_rosace_scores_crizo$type == "synonymous"

met_wt_1 = str_split(substr(met_wt_sequence, 1, 144), '')[[1]]
met_wt_2 = str_split(substr(met_wt_sequence, 145, 287), '')[[1]]

# define the heatmap layout and rows 

Row1 = ggplot(data = ex14_rosace_scores_crizo %>%filter(pos %in% c(1059:1202)),
                  aes(x = pos, y = factor(mutation, level = order, labels = variant_names), fill = score )) +
  geom_tile(aes(color = as.factor(is.wt)), size = 0.2) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(-3,3)) + 
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
    legend.position="none") + 
  labs(y = "Mutation", x = "Position")

Row2 = ggplot(data = ex14_rosace_scores_crizo %>% filter(pos %in% c(1203:1345)), 
                  aes(x = pos, y = factor(mutation, level = order, labels = variant_names), fill = score)) +
  geom_tile(aes(color = as.factor(is.wt)), size = 0.2) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(-3,3)) + 
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
    ) + 
  labs(y = "Mutation", x = "Position")


heatmap_plot = ggarrange(Row1, Row2, nrow = 2, ncol = 1)
ggsave("ex14_rosace_norm_crizo_delta_heatmap.pdf", height = 11, width = 8.5, heatmap_plot)



##-------------------------------------- Print heatmaps ----------------------------------------------##

# plots and saves all inhibitor heatmaps


met_wt_fasta = "Met_wt2.fasta"
met_wt_fasta_con = file(met_wt_fasta, open = "r")
met_wt_sequence = readLines(met_wt_fasta_con)[2]
close(met_wt_fasta_con)
met_wt_1 = str_split(substr(met_wt_sequence, 1, 144), '')[[1]]
met_wt_2 = str_split(substr(met_wt_sequence, 145, 287), '')[[1]]


print_heatmap <- function(df, label, low, high, output_file, met_wt_1, met_wt_2) {
  
  label <- enquo(label)
  
  order <- c('X','H', 'K', 'R', 'D', 'E',
             'C','M','N','Q','S','T','A','I','L','V','W','F',
             'Y','G','P')
  
  variant_names <- c('Stop','H', 'K', 'R', 'D', 'E',
                     'C','M','N','Q','S','T','A','I','L','V','W','F',
                     'Y','G','P')
  
  met_wt_1 = str_split(substr(met_wt_sequence, 1, 144), '')[[1]]
  met_wt_2 = str_split(substr(met_wt_sequence, 145, 287), '')[[1]]
  
  df$is.wt = df$type == "synonymous"

  
  
  row1 = ggplot(data = df %>% filter(pos %in% c(1059:1202)),
                aes(x = pos, y = factor(mutation, level = order, labels = variant_names), fill = !!label )) +
    geom_tile(aes(color = as.factor(is.wt)), size = 0.2) +
    scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9, rev = FALSE, na.value = 'lightyellow', limits = c(low, high)) + 
    scale_color_manual(values = c(NA,'green')) +
    scale_x_continuous(breaks = seq(1059,1202, by = 5),
                       expand = c(0, 0),
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
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
      axis.text.y = element_text(size = 4),
      axis.text = element_text(size =5),
      axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
      axis.ticks = element_blank(),
      legend.position = "none"
    ) +
    labs(y = "Mutation", x = "Position")
  
  row2 = ggplot(data = df %>% filter(pos %in% c(1203:1345)),
                aes(x = pos, y = factor(mutation, level = order, labels = variant_names), fill = !!label )) +
    geom_tile(aes(color = as.factor(is.wt)), size = 0.2) +
    scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9, rev = FALSE, na.value = 'lightyellow', limits = c(low, high)) + 
    scale_color_manual(values = c(NA,'green')) +
    scale_x_continuous(breaks = seq(1203,1345, by = 5),
                       expand = c(0, 0),
                       sec.axis = sec_axis(
                         trans = ~.,
                         breaks = seq(1203,1345),
                         labels = met_wt_2,
                         guide = derive()
                       )) +
    coord_fixed(ratio = 1) +
    theme(
      panel.background = element_rect(fill = "white",
                                      colour = "white",
                                      size = 0.5, linetype = "solid"),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
      axis.text.y = element_text(size = 4),
      axis.text = element_text(size = 5),
      axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
      axis.ticks = element_blank(),
      legend.position = "none"
    ) +
    labs(y = "Mutation", x = "Position")
  
  
  heatmap = ggarrange(row1,row2, ncol = 1, nrow = 2)
  
  ggsave(output_file, height = 4, width = 7.5, heatmap)
  
}

ex14_rosace_scores_dmso <- ex14_rosace_scores_condensed   %>% filter(inhib == "DMSO")
ex14_rosace_scores_crizo <- ex14_rosace_scores_condensed  %>% filter(inhib == "Crizo")
ex14_rosace_scores_tepo <- ex14_rosace_scores_condensed  %>% filter(inhib == "Tepo")
ex14_rosace_scores_cap <- ex14_rosace_scores_condensed  %>% filter(inhib == "Camp")
ex14_rosace_scores_glu<- ex14_rosace_scores_condensed  %>% filter(inhib == "Glu")
ex14_rosace_scores_nvp <- ex14_rosace_scores_condensed  %>% filter(inhib == "NVP")
ex14_rosace_scores_savo <- ex14_rosace_scores_condensed  %>% filter(inhib == "Savo")
ex14_rosace_scores_cabo <- ex14_rosace_scores_condensed  %>% filter(inhib == "Cabo")
ex14_rosace_scores_mere <- ex14_rosace_scores_condensed  %>% filter(inhib == "Mere")
ex14_rosace_scores_gle <- ex14_rosace_scores_condensed  %>% filter(inhib == "Gle")
ex14_rosace_scores_tiv <- ex14_rosace_scores_condensed  %>% filter(inhib == "Tiv")
ex14_rosace_scores_A458 <- ex14_rosace_scores_condensed  %>% filter(inhib == "A458")

print_heatmap(ex14_rosace_scores_dmso, mean, -3, 3, "Ex14_DMSO_heatmap.pdf", met_wt_1, met_wt_2)
print_heatmap(ex14_rosace_scores_A458, score,-3, 3,  "Ex14_A458_heatmap.pdf", met_wt_1, met_wt_2)
print_heatmap(ex14_rosace_scores_cap, score, -3, 3,  "Ex14_Cap_heatmap.pdf", met_wt_1, met_wt_2)
print_heatmap(ex14_rosace_scores_cabo, score, -3, 3,  "Ex14_Cabo_heatmap.pdf", met_wt_1, met_wt_2)
print_heatmap(ex14_rosace_scores_crizo, score, -3, 3,  "Ex14_Crizo_heatmap.pdf", met_wt_1, met_wt_2)
print_heatmap(ex14_rosace_scores_gle, score, -3, 3,  "Ex14_Gle_heatmap.pdf", met_wt_1, met_wt_2)
print_heatmap(ex14_rosace_scores_glu, score, -3, 3,  "Ex14_Glu_heatmap.pdf", met_wt_1, met_wt_2)
print_heatmap(ex14_rosace_scores_mere, score, -3, 3,  "Ex14_Mere_heatmap.pdf", met_wt_1, met_wt_2)
print_heatmap(ex14_rosace_scores_nvp, score, -3, 3,  "Ex14_NVP_heatmap.pdf", met_wt_1, met_wt_2)
print_heatmap(ex14_rosace_scores_savo, score, -3, 3,  "Ex14_Savo_heatmap.pdf", met_wt_1, met_wt_2)
print_heatmap(ex14_rosace_scores_tepo, score, -3, 3,  "Ex14_Tepo_heatmap.pdf", met_wt_1, met_wt_2)
print_heatmap(ex14_rosace_scores_tiv, score, -3, 3,  "Ex14_Tiv_heatmap.pdf", met_wt_1, met_wt_2)


met_rosace_scores_dmso <- met_rosace_scores_condensed   %>% filter(inhib == "DMSO")
met_rosace_scores_crizo <- met_rosace_scores_condensed  %>% filter(inhib == "Crizo")
met_rosace_scores_tepo <- met_rosace_scores_condensed  %>% filter(inhib == "Tepo")
met_rosace_scores_cap <- met_rosace_scores_condensed  %>% filter(inhib == "Camp")
met_rosace_scores_glu<- met_rosace_scores_condensed  %>% filter(inhib == "Glu")
met_rosace_scores_nvp <- met_rosace_scores_condensed  %>% filter(inhib == "NVP")
met_rosace_scores_savo <- met_rosace_scores_condensed  %>% filter(inhib == "Savo")
met_rosace_scores_cabo <- met_rosace_scores_condensed  %>% filter(inhib == "Cabo")
met_rosace_scores_mere <- met_rosace_scores_condensed  %>% filter(inhib == "Mere")
met_rosace_scores_gle <- met_rosace_scores_condensed  %>% filter(inhib == "Gle")
met_rosace_scores_tiv <- met_rosace_scores_condensed  %>% filter(inhib == "Tiv")
met_rosace_scores_A458 <- met_rosace_scores_condensed  %>% filter(inhib == "A458")

print_heatmap(met_rosace_scores_dmso, mean, -5, 3, "met_DMSO_heatmap.pdf", met_wt_1, met_wt_3)
print_heatmap(met_rosace_scores_A458, score,-5, 3,  "met_A458_heatmap.pdf", met_wt_1, met_wt_3)
print_heatmap(met_rosace_scores_cap, score, -5, 3,  "met_Cap_heatmap.pdf", met_wt_1, met_wt_3)
print_heatmap(met_rosace_scores_cabo, score, -5, 3,  "met_Cabo_heatmap.pdf", met_wt_1, met_wt_3)
print_heatmap(met_rosace_scores_crizo, score, -5, 3,  "met_Crizo_heatmap.pdf", met_wt_1, met_wt_3)
print_heatmap(met_rosace_scores_gle, score, -5, 3,  "met_Gle_heatmap.pdf", met_wt_1, met_wt_3)
print_heatmap(met_rosace_scores_glu, score, -5, 3,  "met_Glu_heatmap.pdf", met_wt_1, met_wt_3)
print_heatmap(met_rosace_scores_mere, score, -5, 3,  "met_Mere_heatmap.pdf", met_wt_1, met_wt_3)
print_heatmap(met_rosace_scores_nvp, score, -5, 3,  "met_NVP_heatmap.pdf", met_wt_1, met_wt_3)
print_heatmap(met_rosace_scores_savo, score, -5, 3,  "met_Savo_heatmap.pdf", met_wt_1, met_wt_3)
print_heatmap(met_rosace_scores_tepo, score, -5, 3,  "met_Tepo_heatmap.pdf", met_wt_1, met_wt_3)
print_heatmap(met_rosace_scores_tiv, score, -5, 3,  "met_Tiv_heatmap.pdf", met_wt_1, met_wt_3)




#-----------------------------------------------------------------------------------------------------##
pocket = ex14_rosace_scores_crizo %>%filter(pos %in% c(1108,1092,1163,1158,1159,1211,1260,1228,1230))

pocket_contacts = ggplot(data = pocket, aes(x = pos, y = factor(mutation, level = order, labels = variant_names), fill = score )) +
  geom_tile(aes(color = as.factor(is.wt)), size = 0.2) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(-3,2)) + 
  scale_color_manual(values = c(NA,'green')) +
  scale_x_discrete(limits=c(1108,1092,1163,1158,1159,1211,1260,1228,1230),expand = c(0,0)) +
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
    axis.ticks = element_blank(),
    legend.position="none") + 
  labs(y = "Mutation", x = "Position")

heatmap_pocket = ggarrange(pocket_contacts, nrow = 1, ncol = 1)
ggsave("ex14_pocket_heatmap.pdf", height = 11, width = 8.5, heatmap_pocket)
plot(heatmap_pocket)
