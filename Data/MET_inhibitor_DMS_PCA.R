##---------------------------------------------------------------------------------------------------
# MET_inhibitor_DMS_PCA.R
# by: g.estevam @ Fraser Lab UCSF 
# the goal of this code is to reduce the dimensionality of the MET DMS inhibitor fitnesses with PCA 
# in this code I analyze the PCA for all conditions (ungrouped) and the PCA for inhibitor groups 
# overall the aim is to understand inhibitor driven effects and pattens at a stuctural level 
##---------------------------------------------------------------------------------------------------

library(ggplot2)
library(gglorenz)
library(ggpubr)
library(ggridges)
library(readr)
library(stringr)
library(tidyverse)
library(dplyr)
library(ggridges)
library(ggbraid)
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
library(GGally)


cbp1 <- c("#000000", "#57a1b0","#ae4f74","#544487")

###-------------------------------  initialize data tables -------------------------------------------

ex14_rosace_scores <- as.data.frame(fread("ex14_scores_filtered.tsv")) # normalized scores, unfiltered 
ex14_rosace_scores <- ex14_rosace_scores  %>% mutate(position = position + 1058)
ex14_rosace_scores <- ex14_rosace_scores %>% mutate(pos = position) 
#ex14_rosace_scores  <- ex14_rosace_scores  %>% filter (type == "missense") # removes WT suyn and stops 

ex14_rosace_norm_score <- as.data.frame(fread("ex14_scores_filtered.tsv")) # normalized and filtered scores
ex14_rosace_norm_score  <- ex14_rosace_norm_score   %>% mutate(position = position + 1058)
ex14_rosace_norm_score  <- ex14_rosace_norm_score  %>% mutate(pos = position) 
#ex14_rosace_norm_score  <- ex14_rosace_norm_score  %>% filter (type == "missense") # removes WT suyn and stops 

###------------------------------- PCA on all inhibitors  -------------------------------------------


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
d <- d

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
          palette = cbp1 ,               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          # Rectangle color
          #labels_track_height = 0.8      # Augment the room for labels
)

# Factor cluster
fviz_cluster(res.hcpc,
             repel = TRUE,            # Avoid label overlapping
             palette = cbp1 ,         # Color palette see ?ggpubr::ggpar
             ggtheme = theme_classic(),
             main = "Factor map")

fviz_screeplot(res.pca, addlabels = TRUE, 
               barfill = "#58a1b0",
               barcolor = "black",
               ylim = c(0, 65), 
               ggtheme = theme_classic()) 


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

ex14_PCA <-data.frame( variant = ex14_rosace_scores$variant,
                       mutation =ex14_rosace_scores$mutation,
                       inhibitor = ex14_rosace_scores$key,
                       pos = ex14_rosace_scores$position,
                       type =ex14_rosace_scores$type,
                       score = ex14_rosace_scores$mean)

ex14_PCA <- ex14_PCA %>% filter(type != "nonsense")

ex14_PCA <- ex14_PCA %>% pivot_wider(names_from = inhibitor, values_from = score)
ex14_PCA_4 <- ex14_PCA [-c(2:4)]


# remove hgvs data to keep numerical
ex14_PCA_4 <- ex14_PCA_4 %>% remove_rownames %>% column_to_rownames(var="variant")
#ex14_PCA_5 <- ex14_PCA_4

#ex14_PCA_4 <- t(ex14_PCA_4)
ex14_PCA_4 <- data.frame(ex14_PCA_4)


#not necessary for the normalized data set
#subtract Mut_drug from Mut_dmso
#ex14_PCA_DMSO_sub <- ex14_PCA_4 %>% mutate_at(1:12, funs(c(first(.), (. - first(.))[-1]))) # (Mut_drug)-(Mut_DMSO)
#ex14_PCA_DMSO_sub <- ex14_PCA_DMSO_sub[-1,] #remove DMSO from df for PCA
#ex14_PCA_DMSO_sub <-data.frame(ex14_PCA_DMSO_sub)
#ex14_PCA_DMSO_sub<- ex14_PCA_DMSO_sub[!(row.names(ex14_PCA_DMSO_sub) %in% c("Tiv")),] # remove tivantinib from var
#ex14_PCA_DMSO_sub <- t(ex14_PCA_DMSO_sub) 

# do PCA 
# scale=TRUE bases the PCA on the correlation matrix and FALSE on the covariance matrix
res.pca3 <-PCA(ex14_PCA_4, scale = TRUE) # this is PCA for DMSO sutracted scores 

# plot loadings
fviz_pca_ind(res.pca3, col.ind = group_label,legend.title = "Type", palette =cbp1)+
  theme_classic()+
  ylim(-100,100)

# plot scree 
fviz_screeplot(res.pca3, addlabels = TRUE, 
               barfill = "#58a1b0",
               barcolor = "black",
               ylim = c(0, 100), 
               ggtheme = theme_classic()) 

# plots contributions 
fviz_screeplot(res.pca3, addlabels = TRUE, ylim = c(0, 100))

# contribution of each mutation 
ex14_PCs <- data.frame((res.pca3$ind$coord)) # for prcomp usage to get loadings for base R
ex14_PCs$pos <- ex14_PCA$pos
ex14_PCs$mutation <- ex14_PCA$mutation
ex14_PCs$type <- ex14_PCA$type


#### histogram of dimensions
ggplot(ex14_PCs, aes(x=Dim.1)) + geom_histogram()
ggplot(ex14_PCs, aes(x=Dim.2)) + geom_histogram()
ggplot(ex14_PCs, aes(x=Dim.3)) + geom_histogram()

ggplot()+
  geom_histogram(aes(x=ex14_PCs$Dim.1),color = "black", alpha = 0.8,fill = "#57a1b0")+
  geom_histogram(aes(x=ex14_PCs$Dim.2), color = "black", alpha = 0.8,fill = "#ae4f74")+
  geom_histogram(aes(x=ex14_PCs$Dim.3), color = "black", alpha = 0.8, fill = "gray")+
  xlab("Variable Loadings")+
  ggtitle("PC ditributions")+
  theme_classic()


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


ex14_PCs$is.wt = ex14_PCs$type == "synonymous"
met_wt_3 = str_split(substr(met_wt_sequence, 1, 287), '')[[1]]

Met_row1 = ggplot(data = ex14_PCs %>%filter(pos %in% c(1059:1345)),
                  aes(x = pos, y = factor(mutation, level = order, labels = variant_names), fill = Dim.1)) +
  geom_tile(aes(color = as.factor(is.wt)), size = 0.2) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow') + 
  scale_color_manual(values = c(NA,'green')) +
  scale_x_continuous(breaks = seq(1058,1345, by = 5),
                     expand = c(0,0),
                     sec.axis = sec_axis(
                       trans = ~.,
                       name = "Sequence",
                       breaks = seq(1059,1345),
                       labels = met_wt_3,
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
    axis.text.y = element_text(angle = 90, hjust = 1, size = 4),
    axis.text = element_text(size = 4),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
    axis.ticks = element_blank(),
    legend.position="right"
  ) +
  
  labs(y = "Mutation", x = "Position")
Met_row2 = ggplot(data = ex14_PCs %>%filter(pos %in% c(1059:1345)),
                  aes(x = pos, y = factor(mutation, level = order, labels = variant_names), fill =  Dim.2)) +
  geom_tile(aes(color = as.factor(is.wt)), size = 0.2) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow') + 
  scale_color_manual(values = c(NA,'green')) +
  scale_x_continuous(breaks = seq(1058,1345, by = 5),
                     expand = c(0,0),
                     sec.axis = sec_axis(
                       trans = ~.,
                       name = "Sequence",
                       breaks = seq(1059,1345),
                       labels = met_wt_3,
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
    axis.text.y = element_text(angle = 90, hjust = 1, size = 4),
    axis.text = element_text(size = 4),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
    axis.ticks = element_blank(),
    legend.position="right"
  ) +
  labs(y = "Mutation", x = "Position")
Met_row3 = ggplot(data = ex14_PCs %>%filter(pos %in% c(1059:1345)),
                  aes(x = pos, y = factor(mutation, level = order, labels = variant_names), fill = Dim.3)) +
  geom_tile(aes(color = as.factor(is.wt)), size = 0.2) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow') + 
  scale_color_manual(values = c(NA,'green')) +
  scale_x_continuous(breaks = seq(1058,1345, by = 5),
                     expand = c(0,0),
                     sec.axis = sec_axis(
                       trans = ~.,
                       name = "Sequence",
                       breaks = seq(1059,1345),
                       labels = met_wt_3,
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
    axis.text.y = element_text(angle = 90, hjust = 1, size = 4),
    axis.text = element_text(size = 4),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
    axis.ticks = element_blank(),
    legend.position="right"
  ) +
  labs(y = "Mutation", x = "Position")

PC1_heatmap = ggarrange(Met_row1,ncol = 1)
PC2_heatmap = ggarrange(Met_row2,ncol = 1)
PC3_heatmap = ggarrange(Met_row3, ncol = 1)
plot(PC1_heatmap)
plot(PC2_heatmap)
plot(PC3_heatmap)

ggsave("PC1_heatmap.pdf", height = 3, width = 15, PC1_heatmap)
ggsave("PC2_heatmap.pdf", height = 3, width = 15, PC2_heatmap)
ggsave("PC3_heatmap.pdf", height = 3, width = 15, PC3_heatmap)

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
ex14_2G15_inactive <- read.pdb("2G15") 
ex14_2WGJ_crizo <- read.pdb("2WGJ") 
ex14_4EEV_mere <- read.pdb("4EEV")

#PC1
ex14_PC1 <- data.frame(pos=ex14_PCs$pos, PC1=ex14_PCs$PC1)
ex14_PC1_avg <- ex14_PC1 %>% group_by(pos) %>% dplyr::summarise(avg = mean(PC1))%>% ungroup()
x = map_scores_pdb(ex14_2G15_inactive  , ex14_PC1_avg, "avg")
write.pdb(x, file="ex14_PC1.pdb")

#PC2
ex14_PC2 <- data.frame(pos=ex14_PCs$pos, PC2=ex14_PCs$PC2)
ex14_PC2_avg <- ex14_PC2 %>% group_by(pos) %>% dplyr::summarise(avg = mean(PC2))%>% ungroup()
x = map_scores_pdb(ex14_2WGJ_crizo , ex14_PC2_avg, "avg")
write.pdb(x, file="ex14_PC2.pdb")

#PC3
ex14_PC3 <- data.frame(pos=ex14_PCs$pos, PC3=ex14_PCs$PC3)
ex14_PC3_avg <- ex14_PC3 %>% group_by(pos) %>% dplyr::summarise(avg = mean(PC3))%>% ungroup()
x = map_scores_pdb(ex14_2WGJ_crizo , ex14_PC3_avg, "avg")
write.pdb(x, file="ex14_PC3.pdb")


#####################################################################################################
###---------------------------------- -- PCA grouped------------------------------------------
#####################################################################################################

# scores grouped by inhibitor

#-------------------------------- type I 
ex14_PCA_type_I <-  subset(ex14_PCA_4 , select=c("Crizo","Camp","Tepo","Glu","NVP","Savo"))
#ex14_PCA_type_I <- na.omit(ex14_PCA_type_I)
ex14_type_I.pca <- PCA(ex14_PCA_type_I, scale = FALSE)

#ex14_PCA_type_I  <- t(ex14_PCA_type_I)
#ex14_type_I.pca <- PCA(ex14_PCA_type_I*-1, scale = TRUE)

#scree plot
fviz_screeplot(ex14_type_I.pca, addlabels = TRUE, ylim = c(0, 100))
fviz_pca_ind(ex14_type_I.pca)

fviz_screeplot(ex14_type_I.pca, addlabels = TRUE, 
               barfill = "#58a1b0",
               barcolor = "black",
               ylim = c(0, 100), 
               ggtheme = theme_classic()) 


ex14_PCA_type_I_ex14_PCs <- data.frame((ex14_type_I.pca$ind$coord)) # lodgings ie contributions
ex14_PCA_type_I_ex14_PCs$pos <- ex14_PCA$pos
ex14_PCA_type_I_ex14_PCs$mutation <- ex14_PCA$mutation
ex14_PCA_type_I_ex14_PCs$type<- ex14_PCA$type


ex14_PCA_type_I_ex14_PCs$is.wt = ex14_PCA_type_I_ex14_PCs$type == "synonymous"
met_wt_3 = str_split(substr(met_wt_sequence, 1, 287), '')[[1]]

order <- c('X','H', 'K', 'R', 'D', 'E',
           'C','M','N','Q','S','T','A','I','L','V','W','F',
           'Y','G','P')

variant_names <- c('Stop','H', 'K', 'R', 'D', 'E',
                   'C','M','N','Q','S','T','A','I','L','V','W','F',
                   'Y','G','P')

met_wt_3 = str_split(substr(met_wt_sequence, 1, 287), '')[[1]]

Met_row1 = ggplot(data = ex14_PCA_type_I_ex14_PCs%>%filter(pos %in% c(1059:1345)),
                  aes(x = pos, y = factor(mutation, level = order, labels = variant_names), fill = Dim.1)) +
  geom_tile(aes(color = as.factor(is.wt)), size = 0.2) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow') + 
  scale_color_manual(values = c(NA,'green')) +
  scale_x_continuous(breaks = seq(1058,1345, by = 5),
                     expand = c(0,0),
                     sec.axis = sec_axis(
                       trans = ~.,
                       name = "Sequence",
                       breaks = seq(1059,1345),
                       labels = met_wt_3,
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
    axis.text.y = element_text(angle = 90, hjust = 1, size = 4),
    axis.text = element_text(size = 4),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
    axis.ticks = element_blank(),
    legend.position="right"
  ) +
  labs(y = "Mutation", x = "Position")

typeI_heatmap_PC1= ggarrange(Met_row1, ncol = 1)
plot(typeI_heatmap_PC1) # plot heatmap
ggsave("typeI_heatmap_PC1.pdf", height = 7 , width = 8.5, typeI_heatmap_PC1)



# map scores onto 3d stucture 
ex14_I_PC1 <- data.frame(pos=ex14_PCA_type_I_ex14_PCs$pos, PC1=ex14_PCA_type_I_ex14_PCs$Dim.1)
ex14_I_PC1_avg <- ex14_I_PC1 %>% group_by(pos) %>% dplyr::summarise(avg = mean(PC1))%>% ungroup()
x = map_scores_pdb(ex14_2WGJ_crizo , ex14_I_PC1_avg , "avg")
write.pdb(x, file="ex14_I_PC1.pdb")


#-------------------------------- type II

ex14_PCA_type_II <-  subset(ex14_PCA_4 , select=c("Cabo","Mere","Gle"))
ex14_type_II.pca <- PCA(ex14_PCA_type_II, scale = FALSE)

#scree plot
fviz_screeplot(ex14_type_II.pca, addlabels = TRUE, ylim = c(0, 100))
fviz_pca_ind(ex14_type_II.pca)

fviz_screeplot(ex14_type_II.pca, addlabels = TRUE, 
               barfill = "#58a1b0",
               barcolor = "black",
               ylim = c(0, 100), 
               ggtheme = theme_classic()) 


ex14_PCA_type_II_ex14_PCs <- data.frame((ex14_type_II.pca$ind$coord)) # lodgings ie contributions
ex14_PCA_type_II_ex14_PCs$pos <- ex14_PCA$pos
ex14_PCA_type_II_ex14_PCs$mutation <- ex14_PCA$mutation
ex14_PCA_type_II_ex14_PCs$type<- ex14_PCA$type


ex14_PCA_type_II_ex14_PCs$is.wt = ex14_PCA_type_II_ex14_PCs$type == "synonymous"
met_wt_3 = str_split(substr(met_wt_sequence, 1, 287), '')[[1]]


order <- c('H', 'K', 'R', 'D', 'E',
           'C','M','N','Q','S','T','A','I','L','V','W','F',
           'Y','G','P')

variant_names <- c('H', 'K', 'R', 'D', 'E',
                   'C','M','N','Q','S','T','A','I','L','V','W','F',
                   'Y','G','P')

met_wt_3 = str_split(substr(met_wt_sequence, 1, 287), '')[[1]]

Met_row2 = ggplot(data = ex14_PCA_type_II_ex14_PCs%>%filter(pos %in% c(1059:1345)),
                  aes(x = pos, y = factor(mutation, level = order, labels = variant_names), fill = Dim.1)) +
  geom_tile(aes(color = as.factor(is.wt)), size = 0.2) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow') + 
  scale_color_manual(values = c(NA,'green')) +
  scale_x_continuous(breaks = seq(1058,1345, by = 5),
                     expand = c(0,0),
                     sec.axis = sec_axis(
                       trans = ~.,
                       name = "Sequence",
                       breaks = seq(1059,1345),
                       labels = met_wt_3,
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
    axis.text.y = element_text(angle = 90, hjust = 1, size = 4),
    axis.text = element_text(size = 4),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
    axis.ticks = element_blank(),
    legend.position="right"
  ) +
  labs(y = "Mutation", x = "Position")

typeII_heatmap_PC1= ggarrange(Met_row2, ncol = 1)
plot(typeII_heatmap_PC1) # plot heatmap
ggsave("typeII_heatmap_PC1.pdf", height = 7 , width = 8.5, typeII_heatmap_PC1)

# map scores onto 3d stucture 
ex14_II_PC1 <- data.frame(pos=ex14_PCA_type_II_ex14_PCs$pos, PC1=ex14_PCA_type_II_ex14_PCs$Dim.1)
ex14_II_PC1_avg <- ex14_II_PC1 %>% group_by(pos) %>% dplyr::summarise(avg = mean(PC1))%>% ungroup()
ex14_4EEV_mere <- read.pdb("4EEV")
x = map_scores_pdb(ex14_4EEV_mere , ex14_II_PC1_avg , "avg")
write.pdb(x, file="ex14_II_PC.pdb")

#---------------------------------- PCA Difference Map-------------------------------------------##

diff_PC1_df <- data.frame(pos = ex14_PCA_type_I_ex14_PCs$pos, 
                          mutation = ex14_PCA_type_I_ex14_PCs$mutation,
                          type = ex14_PCA_type_I_ex14_PCs$type,
                          I_PC1= ex14_PCA_type_I_ex14_PCs$Dim.1,
                          II_PC1=ex14_PCA_type_II_ex14_PCs$Dim.1)

diff_PC1_df$diff <- (diff_PC1_df$I_PC1) - (diff_PC1_df$II_PC1)


diff_PC1_df$is.wt = diff_PC1_df$type == "synonymous"
met_wt_3 = str_split(substr(met_wt_sequence, 1, 287), '')[[1]]

Met_row3 = ggplot(data = diff_PC1_df%>%filter(pos %in% c(1059:1345)),
                  aes(x = pos, y = factor(mutation, level = order, labels = variant_names), fill = diff)) +
  geom_tile(aes(color = as.factor(is.wt)), size = 0.2) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow') + 
  scale_color_manual(values = c(NA,'green')) +
  scale_x_continuous(breaks = seq(1058,1345, by = 5),
                     expand = c(0,0),
                     sec.axis = sec_axis(
                       trans = ~.,
                       name = "Sequence",
                       breaks = seq(1059,1345),
                       labels = met_wt_3,
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
    axis.text.y = element_text(angle = 90, hjust = 1, size = 4),
    axis.text = element_text(size = 4),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
    axis.ticks = element_blank(),
    legend.position="right"
  ) +
  labs(y = "Mutation", x = "Position")

diff_map= ggarrange(Met_row3, ncol = 1)
plot(diff_map) # plot heatmap
ggsave("diff_map_PC.pdf", height = 7 , width = 8.5, diff_map)

# map scores onto 3d stucture 
diff_PCA <- data.frame(pos=diff_PC1_df$pos, PC1=diff_PC1_df$diff)
diff_PCA_avg <- diff_PCA %>% group_by(pos) %>% dplyr::summarise(avg = mean(PC1))%>% ungroup()
ex14_2WGJ_crizo <- read.pdb("2WGJ")
x = map_scores_pdb(ex14_2WGJ_crizo , diff_PCA_avg , "avg")
write.pdb(x, file="diff_PCA.pdb")



#---------------------------------- PCA correlations ---------------------------------------------

# here the idea is to plot the PCA data against the raw data for each condition to see where the corlelations are 
library(dplyr)
library(tidyr)

PCA_correlation <- t(ex14_PCA_4)
PCA_correlation_merge <- merge(PCA_correlation, ex14_PCs, by=0, all=TRUE) #merges raw values with pcs 

PC1_correlations <- PCA_correlation_merge[-(15:28)]
PC2_correlations <- PCA_correlation_merge[-(16:28)]
PC3_correlations <- PCA_correlation_merge[-(17:28)]
  
df.m<- melt(PC1_correlations[-1],"PC1")
df.n<- melt(PC2_correlations[-1],"PC2")
df.o<- melt(PC3_correlations[-1],"PC3")

PC1_correlation_plot = ggplot() + ggscatter(df.m, x="PC1", y="value", add = "reg.line", color="variable") + stat_cor(method = "pearson") +facet_wrap(~ variable, ncol = 2)
PC2_correlation_plot = ggplot() + ggscatter(df.n, x="PC2", y="value", add = "reg.line", color="variable") + stat_cor(method = "pearson") +facet_wrap(~ variable, ncol = 2)
PC3_correlation_plot = ggplot() + ggscatter(df.o, x="PC3", y="value", add = "reg.line", color="variable") + stat_cor(method = "pearson") +facet_wrap(~ variable, ncol = 2)

ggsave("PC1_correlation_plot.pdf", height = 7, width = 8, PC1_correlation_plot)
ggsave("PC2_correlation_plot.pdf", height = 7, width = 8, PC2_correlation_plot)
ggsave("PC3_correlation_plot.pdf", height = 7, width = 8, PC3_correlation_plot)


# corelation heatmap
ggpair(df, columns = 1: ncol(df), title = NULL,
       upper = list(continuous = "cor"),
       lower = list(continuous = "smooth"),
       mapping = NULL)

