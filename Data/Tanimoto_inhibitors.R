###################################################################################################################
# Tanimoto_inhibitors.R 
# by g.estevam @ UCSF 
# this script finds common inhibitor chemistry 
# acts as a basis for inhibitor comparisons 
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

library(scales)
library(colorspace)
library(ggpubr)
library(rstatix)
library(scales)
library(stringdist)
#library(circlize)
library(RecordLinkage)
library(ChemmineR)
library(ChemmineOB)


smiles_list <- list(
  Crizotinib = "CC(C1=C(C=CC(=C1Cl)F)Cl)OC2=C(N=CC(=C2)C3=CN(N=C3)C4CCNCC4)N",
  Cabozantinib = "COC1=CC2=C(C=CN=C2C=C1OC)OC3=CC=C(C=C3)NC(=O)C4(CC4)C(=O)NC5=CC=C(C=C5)F",
  Capmatinib = "CNC(=O)C1=C(C=C(C=C1)C2=NN3C(=CN=C3N=C2)CC4=CC5=C(C=C4)N=CC=C5)F",
  Tepotinib = "CN1CCC(CC1)COC2=CN=C(N=C2)C3=CC=CC(=C3)CN4C(=O)C=CC(=N4)C5=CC=CC(=C5)C",
  Tivantinib = "C1CC2=C3C(=CC=C2)C(=CN3C1)C4C(C(=O)NC4=O)C5=CNC6=CC=CC=C65",
  AMG458 = "CC1=C(C(=O)N(N1CC(C)(C)O)C2=CC=CC=C2)C(=O)NC3=NC=C(C=C3)OC4=C5C=CC(=CC5=NC=C4)OC",
  Glumetinib = "CN1C=C(C=N1)C2=CN3C(=NC=C3S(=O)(=O)N4C5=C(C=N4)N=CC(=C5)C6=CN(N=C6)C)C=C2",
  NVPBVU972 = "CN1C=C(C=N1)C2=NN3C(=NC=C3CC4=CC5=C(C=C4)N=CC=C5)C=C2",
  Merestinib = "CC1=CC=C(C(=O)N1C2=CC=C(C=C2)F)C(=O)NC3=CC(=C(C=C3)OC4=C(C=C5C(=C4)C=NN5C)C6=CNN=C6)F",
  Glesatinib = "CN1C=C(N=C1)C2=CC3=NC=CC(=C3S2)OC4=C(C=C(C=C4)NC(=S)NC(=O)CC5=CC=CC=C5)F",
  Savolitinib = "CC(C1=CN2C=CN=C2C=C1)N3C4=NC(=CN=C4N=N3)C5=CN(N=C5)"
  #Erlotinib = 	"COCCOC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC=CC(=C3)C#C)OCCOC.Cl" , 
  #Osimertinib = "CN1C=C(C2=CC=CC=C21)C3=NC(=NC=C3)NC4=C(C=C(C(=C4)NC(=O)C=C)N(C)CCN(C)C)OC", 
  
)

compound_names <- c(
  "Crizotinib",
  "Cabozantinib",
  "Capmatinib",
  "Tepotinib",
  "Tivantinib",
  "AMG458",
  "Glumetinib",
  "NVPBVU972",
  "Merestinib",
  "Glesatinib",
  "Savolitinib"
)

n <- length(smiles_list)
distance_matrix <- matrix(0, n, n)

for (i in 1:(n - 1)) {
  for (j in (i + 1):n) {
    levenshtein_similarity <- levenshteinSim(smiles_list[[i]], smiles_list[[j]])
    tanimoto_distance <- 1 - levenshtein_similarity
    distance_matrix[i, j] <- tanimoto_distance
    distance_matrix[j, i] <- tanimoto_distance
  }
}


print(distance_matrix)



smiles_list <- c(
  "CC(C1=C(C=CC(=C1Cl)F)Cl)OC2=C(N=CC(=C2)C3=CN(N=C3)C4CCNCC4)N",
  "COC1=CC2=C(C=CN=C2C=C1OC)OC3=CC=C(C=C3)NC(=O)C4(CC4)C(=O)NC5=CC=C(C=C5)F",
  "CNC(=O)C1=C(C=C(C=C1)C2=NN3C(=CN=C3N=C2)CC4=CC5=C(C=C4)N=CC=C5)F",
  "CN1CCC(CC1)COC2=CN=C(N=C2)C3=CC=CC(=C3)CN4C(=O)C=CC(=N4)C5=CC=CC(=C5)C",
  "C1CC2=C3C(=CC=C2)C(=CN3C1)C4C(C(=O)NC4=O)C5=CNC6=CC=CC=C65",
  "CC1=C(C(=O)N(N1CC(C)(C)O)C2=CC=CC=C2)C(=O)NC3=NC=C(C=C3)OC4=C5C=CC(=CC5=NC=C4)OC",
  "CN1C=C(C=N1)C2=CN3C(=NC=C3S(=O)(=O)N4C5=C(C=N4)N=CC(=C5)C6=CN(N=C6)C)C=C2",
  "CN1C=C(C=N1)C2=NN3C(=NC=C3CC4=CC5=C(C=C4)N=CC=C5)C=C2",
  "CC1=CC=C(C(=O)N1C2=CC=C(C=C2)F)C(=O)NC3=CC(=C(C=C3)OC4=C(C=C5C(=C4)C=NN5C)C6=CNN=C6)F",
  "CN1C=C(N=C1)C2=CC3=NC=CC(=C3S2)OC4=C(C=C(C=C4)NC(=S)NC(=O)CC5=CC=CC=C5)F",
  "CC(C1=CN2C=CN=C2C=C1)N3C4=NC(=CN=C4N=N3)C5=CN(N=C5)"
)


molecules <- lapply(smiles_list, function(smiles) parse.smiles(smiles))
                    
similarity_matrix <- Tanimoto(molecules)

print(similarity_matrix)





#----------------------- heatmap -------------------------------# 

desired_order <- c("Crizotinib", "Capmatinib", "Tepotinib", "Glumetinib", "Savolitinib", "NVPBVU972", "Tivantinib", "Cabozantinib", "Merestinib", "Glesatinib", "AMG458")


distance_matrix_ordered <- distance_matrix[match(desired_order, rownames(distance_matrix)), match(desired_order, colnames(distance_matrix))]

heatmap_data <- expand.grid(x = names(smiles_list), y = names(smiles_list))
heatmap_data$Distance <- as.vector(distance_matrix)


heatmap_data_ordered <- heatmap_data[heatmap_data$x %in% desired_order & heatmap_data$y %in% desired_order, ]

heatmap_data_ordered$x <- factor(heatmap_data_ordered$x, levels = desired_order)
heatmap_data_ordered$y <- factor(heatmap_data_ordered$y, levels = desired_order)

ggplot(heatmap_data_ordered, aes(x = x, y = y, fill = Distance, label = round(Distance, 2))) +
  geom_tile() +
  scale_fill_distiller(palette = "Blues", type = "seq", direction = 1, na.value = "lightyellow", name = "Distance") +
  geom_text(color = "black", size = 3) +  # set text color to black
  labs(title = "Tanimoto Distance Heatmap",
       x = "Compounds", y = "Compounds") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # rotate x-axis labels
  scale_color_distiller(palette = "Blues", type = "seq") +  # color text labels
  guides(fill = guide_colorbar(title = "Distance"), color = FALSE)  # hide color legend

######

# Tanimoto distances for all pairs of SMILES 
n <- length(smiles_list)
distance_matrix <- matrix(0, n, n)

for (i in 1:(n - 1)) {
  for (j in (i + 1):n) {
    mol1 <- smilesString(smiles_list[[i]])
    mol2 <- smilesString(smiles_list[[j]])
    tanimoto_distance <- Tanimoto(mol1, mol2)
    distance_matrix[i, j] <- tanimoto_distance
    distance_matrix[j, i] <- tanimoto_distance
  }
}

print(distance_matrix)

                    
if (!requireNamespace("webchem", quietly = TRUE)) {
  install.packages("webchem")
}

library(webchem)
get_2d_structure <- function(smiles, name) {
  csid <- webchem::get_csid_structure(smiles, fileformat = "png")
  plot_path <- paste0(name, "_structure.png")
  file.copy(csid, plot_path)
  print(plot_path)
}


for (name in names(smiles_list)) {
  get_2d_structure(smiles_list[[name]], name)
}
