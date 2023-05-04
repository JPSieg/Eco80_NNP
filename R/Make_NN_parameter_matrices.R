library(MeltR)
library(tidyverse)

####1M NaCl####

####Load in helix energy data####

list.files("Tables/SI_Table_1_fitting_statistics")

df = read.csv("Tables/SI_Table_1_fitting_statistics/SI_Table_1_data.csv") %>%
  filter(Condition == "1 M NaCl") %>%
  filter(Method == "2 Global fit")

####Create NN parameter matrix####

list.df.param.count = {}

for (i in 1:length(df$Helix)){
  list.df.param.count[[i]] = Helix.energy(df$Sequence.F[i],
                                          df$Sequence.R[i],
                                          F.Q = TRUE, Double.label = 0)
}

df.param = bind_rows(list.df.param.count)

df.param$dG = df$G
df.param$dG.error = df$SE.G
df.param$Helix = df$Helix
write.csv(df.param, "Tables/SI_Table_2_established_NN_comparison/SI_Table_2_NN_matrix_and_data_1MNaCl.csv", row.names = F)

####Eco80####

####Load in helix energy data####

list.files("Tables/SI_Table_1_fitting_statistics")

df = read.csv("Tables/SI_Table_1_fitting_statistics/SI_Table_1_data.csv") %>%
  filter(Condition == "Eco80") %>%
  filter(Method == "2 Global fit")

####Create NN parameter matrix####

list.df.param.count = {}

for (i in 1:length(df$Helix)){
  list.df.param.count[[i]] = Helix.energy(df$Sequence.F[i],
                                          df$Sequence.R[i],
                                          F.Q = TRUE, Double.label = 0)
}

df.param = bind_rows(list.df.param.count)

df.param$dG = df$G
df.param$dG.error = df$SE.G
df.param$Helix = df$Helix

write.csv(df.param, "Tables/SI_Table_2_established_NN_comparison/SI_Table_2_NN_matrix_and_data_Eco80.csv", row.names = F)
