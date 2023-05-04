library(tidyverse)
library(ggbeeswarm)
library(viridis)
library(cowplot)
library(NatParksPalettes)

list.files("Tables/Table_3_Prediction_accuracy")

####E coli RNA set accuracy####

df = read.csv("Figures/Figure_2_corrections_and_predictions/Predict_Ecoli_RNA_structure/results/001_Accuracy.csv")

df$Program = factor(df$Program,
                    levels = c("UV_1MNaCl_data_tables",
                               "p7kcalmol_penalty_data_tables",
                               "1kcalmol_penalty_data_tables",
                               "1p7kcalmol_data_tables",
                               "FDBI_1MNaCl_data_tables",
                               "FBDI_Eco80_data_tables",
                               "mxfold2"),
                    labels = c("UV 1MNaCl",
                               "0.7 kcal/mol penalty",
                               "1.0 kcal/mol penalty",
                               "1.7 kcal/mol penalty",
                               "FDBI 1MNaCl",
                               "FBDI Eco80",
                               "mxfold2"))

df.accuracy = df %>%
  group_by(Program) %>%
  summarise(Mean.sens = mean(PPV),
            Median.sens = median(PPV),
            Mean.PPV = mean(PPV),
            Median.PPV = median(PPV),
            Mean.Accuracy = mean(Accuracy),
            Median.Accuracy = median(Accuracy)) %>%
  filter(Program %in% c("UV 1MNaCl", "FDBI 1MNaCl", "1.0 kcal/mol penalty", "FBDI Eco80"))

write.csv(df.accuracy, "Tables/Table_3_Prediction_accuracy/Ecoli_set_accuracy.csv", row.names = F)

####Ecoli Accuracy changes####

v.tables = unique(df$Program)[-which(unique(df$Program) %in% c("UV 1MNaCl"))]

l.df.change = {}

for (i in 1:length(v.tables)){
  print(i)
  df.UV.1M = df %>% filter(Program == "UV 1MNaCl") %>% arrange(RNA)
  df.i = df %>% filter(Program == v.tables[i]) %>% arrange(RNA)
  Change = df.i$Accuracy - df.UV.1M$Accuracy
  N = c(length(which(Change > 0)),
        length(which(Change < 0)),
        length(which(Change == 0)))
  Effect = c("Improved", "Got worse", "No change")
  Data.table = df.i$Program[1]
  l.df.change[[i]] = data.frame(Data.table, Effect, N)
}

df.change = bind_rows(l.df.change)

df.change$Data.table = factor(df.change$Data.table,
                              levels = c("0.7 kcal/mol penalty",
                                         "1.0 kcal/mol penalty",
                                         "1.7 kcal/mol penalty",
                                         "FDBI 1MNaCl",
                                         "FBDI Eco80",
                                         "mxfold2"))
df.change = df.change %>%
  filter(Data.table %in% c("UV 1MNaCl", "FDBI 1MNaCl", "1.0 kcal/mol penalty", "FBDI Eco80")) 

write.csv(df.change, "Tables/Table_3_Prediction_accuracy/Ecoli_set_change.csv", row.names = F)

####mxfold2 RNA set accuracy####

df = read.csv("Figures/Figure_2_corrections_and_predictions/Predict_RNA_structure_3/results/001_Accuracy.csv")

df$Program = factor(df$Program,
                    levels = c("UV_1MNaCl_data_tables",
                               "p7kcalmol_penalty_data_tables",
                               "1kcalmol_penalty_data_tables",
                               "1p7kcalmol_data_tables",
                               "FDBI_1MNaCl_data_tables",
                               "FBDI_Eco80_data_tables",
                               "mxfold2"),
                    labels = c("UV 1MNaCl",
                               "0.7 kcal/mol penalty",
                               "1.0 kcal/mol penalty",
                               "1.7 kcal/mol penalty",
                               "FDBI 1MNaCl",
                               "FBDI Eco80",
                               "mxfold2"))

df.accuracy = df %>%
  group_by(Program) %>%
  summarise(Mean.sens = mean(PPV),
            Median.sens = median(PPV),
            Mean.PPV = mean(PPV),
            Median.PPV = median(PPV),
            Mean.Accuracy = mean(Accuracy),
            Median.Accuracy = median(Accuracy)) %>%
  filter(Program %in% c("UV 1MNaCl", "FDBI 1MNaCl", "1.0 kcal/mol penalty", "FBDI Eco80"))

write.csv(df.accuracy, "Tables/Table_3_Prediction_accuracy/mxfold2_set_accuracy.csv", row.names = F)

####mxfold2 Accuracy changes####

v.tables = unique(df$Program)[-which(unique(df$Program) %in% c("UV 1MNaCl"))]

l.df.change = {}

for (i in 1:length(v.tables)){
  print(i)
  df.UV.1M = df %>% filter(Program == "UV 1MNaCl") %>% arrange(RNA)
  df.i = df %>% filter(Program == v.tables[i]) %>% arrange(RNA)
  Change = df.i$Accuracy - df.UV.1M$Accuracy
  N = c(length(which(Change > 0)),
        length(which(Change < 0)),
        length(which(Change == 0)))
  Effect = c("Improved", "Got worse", "No change")
  Data.table = df.i$Program[1]
  l.df.change[[i]] = data.frame(Data.table, Effect, N)
}

df.change = bind_rows(l.df.change)

df.change$Data.table = factor(df.change$Data.table,
                              levels = c("0.7 kcal/mol penalty",
                                         "1.0 kcal/mol penalty",
                                         "1.7 kcal/mol penalty",
                                         "FDBI 1MNaCl",
                                         "FBDI Eco80",
                                         "mxfold2"))
df.change = df.change %>%
  filter(Data.table %in% c("UV 1MNaCl", "FDBI 1MNaCl", "1.0 kcal/mol penalty", "FBDI Eco80")) 

write.csv(df.change, "Tables/Table_3_Prediction_accuracy/mxfold2_set_change.csv", row.names = F)

####mxfold2 RNA substructures####

df = read.csv("Figures/Figure_2_corrections_and_predictions/Predict_RNA_structure_3/results/002_Substructures.csv")

df$Program = factor(df$Program,
                    levels = c("UV_1MNaCl_data_tables",
                               "p7kcalmol_penalty_data_tables",
                               "1kcalmol_penalty_data_tables",
                               "1p7kcalmol_data_tables",
                               "FDBI_1MNaCl_data_tables",
                               "FBDI_Eco80_data_tables",
                               "mxfold2"),
                    labels = c("UV 1MNaCl",
                               "0.7 kcal/mol penalty",
                               "1.0 kcal/mol penalty",
                               "1.7 kcal/mol penalty",
                               "FDBI 1MNaCl",
                               "FBDI Eco80",
                               "mxfold2"))

df.sub = df %>%
  group_by(Program) %>%
  summarise(Mean.substructures = mean(Structures),
            Median.substructures = median(Structures),
            Total.substructures = sum(Structures))  %>%
  filter(Program %in% c("UV 1MNaCl", "FDBI 1MNaCl", "1.0 kcal/mol penalty", "FBDI Eco80"))

write.csv(df.sub, "Tables/Table_3_Prediction_accuracy/mxfold2_set_substructures.csv", row.names = F)

####mxfold2 features####

df = read.csv("Figures/Figure_2_corrections_and_predictions/Predict_RNA_structure_3/results/003_Structural_features.csv")

df$Program = factor(df$Program,
                    levels = c("UV_1MNaCl_data_tables",
                               "p7kcalmol_penalty_data_tables",
                               "1kcalmol_penalty_data_tables",
                               "1p7kcalmol_data_tables",
                               "FDBI_1MNaCl_data_tables",
                               "FBDI_Eco80_data_tables",
                               "mxfold2"),
                    labels = c("UV 1MNaCl",
                               "0.7 kcal/mol penalty",
                               "1.0 kcal/mol penalty",
                               "1.7 kcal/mol penalty",
                               "FDBI 1MNaCl",
                               "FBDI Eco80",
                               "mxfold2"))


df$pPaired = 100*df$Paired.n/df$Length
df$pAU.pairs = 100*df$AU.pairs/df$Pairs
df$pGC.pairs = 100*df$GC.pairs/df$Pairs
df$pGU.pairs = 100*df$GU.pairs/df$Pairs


df.features = df %>%
  group_by(Program) %>%
  summarise(Mean.pPaired = mean(pPaired),
            Median.pPaired = median(pPaired),
            Mean.pAU.pairs = mean(pAU.pairs),
            Median.pAU.pairs = median(pAU.pairs),
            Mean.pGC.pairs = mean(pGC.pairs),
            Median.pGC.pairs = median(pGC.pairs),
            Mean.pGU.pairs = mean(pGU.pairs),
            Median.pGU.pairs = median(pGU.pairs))  %>%
  filter(Program %in% c("UV 1MNaCl", "FDBI 1MNaCl", "1.0 kcal/mol penalty", "FBDI Eco80"))

write.csv(df.features, "Tables/Table_3_Prediction_accuracy/mxfold2_set_features.csv", row.names = F)
