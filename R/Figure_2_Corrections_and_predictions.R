library(tidyverse)
library(ggbeeswarm)
library(viridis)
library(cowplot)
library(NatParksPalettes)

list.files()

####PA Input RNA####

list.files("Figures/Figure_2_corrections_and_predictions/Predict_Ecoli_RNA_structure/resources")

df.RNA = read.delim("Figures/Figure_2_corrections_and_predictions/Predict_Ecoli_RNA_structure/resources/Ecoli_RNA_names.txt")

df.n = data.frame(table(df.RNA$Class)) %>%
  arrange(Freq)

df.n$Var1 = factor(df.n$Var1,
                   levels = unique(df.n$Var1))

A = ggplot(df.n, aes(x="", y = Freq, fill = Var1)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  theme_classic() +
  ggtitle(paste("N = ", nrow(df.RNA))) +
  geom_text(aes(label = paste0(Freq)),
            position = position_stack(vjust=0.5),
            color = "white") +
  scale_fill_manual(values = viridis(8), name = "RNA type") +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.key.size = unit(0.05, 'in'))

####Folding results####

list.files("Figures/Figure_2_corrections_and_predictions/Predict_Ecoli_RNA_structure/results")

df = read.csv("Figures/Figure_2_corrections_and_predictions/Predict_Ecoli_RNA_structure/results/001_Accuracy.csv")

####PB Accuracy results####

df$Program = factor(df$Program,
                    levels = c("UV_1MNaCl_data_tables",
                               "FDBI_1MNaCl_data_tables",
                               "p7kcalmol_penalty_data_tables",
                               "1kcalmol_penalty_data_tables",
                               "1p7kcalmol_data_tables",
                               "FBDI_Eco80_data_tables",
                               "mxfold2"),
                    labels = c("UV 1MNaCl",
                               "FDBI 1MNaCl",
                               "0.7 kcal/mol penalty",
                               "1.0 kcal/mol penalty",
                               "1.7 kcal/mol penalty",
                               "FBDI Eco80",
                               "mxfold2"))
B = ggplot(df %>% filter(Program %in% c("UV 1MNaCl", "FDBI 1MNaCl", "1.0 kcal/mol penalty", "FBDI Eco80")), aes(x = Program, y = Accuracy)) +
  geom_line(mapping = aes(group = RNA), color = "grey", alpha = 0.5) +
  geom_boxplot() +
  geom_beeswarm(size = 0.5) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black",
                                   angle = 45,
                                   hjust = 1,
                                   size = 8),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_blank(),
        strip.text = element_text(color = "black"),
        legend.title = element_blank(),
        legend.position = "none")

####PC Accuracy changes####

v.tables = unique(df$Program)[-which(unique(df$Program) == "UV 1MNaCl")]

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
                                         "FDBI 1MNaCl",
                                         "1.0 kcal/mol penalty",
                                         "1.7 kcal/mol penalty",
                                         "FBDI Eco80",
                                         "mxfold2"))
unique(df.change$Effect)
df.change$Effect = factor(df.change$Effect,
                          levels = c("Improved", "No change", "Got worse"))

Bryce = natparks.pals("BryceCanyon", 9)

C = ggplot(df.change %>% filter(Data.table %in% c("UV 1MNaCl", "FDBI 1MNaCl", "1.0 kcal/mol penalty", "FBDI Eco80")), aes(x = Data.table, y = N, fill = Effect)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_manual(values = Bryce[c(5,6,1)]) +
  scale_y_continuous(limits = c(0, 115), breaks = c(0, 25, 50, 75, 100)) +
  ylab("Count") +
  theme(axis.text.x = element_text(color = "black",
                                   angle = 45,
                                   hjust = 1,
                                   size = 8),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_blank(),
        strip.text = element_text(color = "black"),
        legend.title = element_blank())

####PD Input RNA ML set####

list.files("Figures/Figure_2_corrections_and_predictions/Predict_RNA_structure_3/resources")

df.RNA = read.csv("Figures/Figure_2_corrections_and_predictions/Predict_RNA_structure_3/resources/Stucture_classifications.csv")%>%
  group_by(RNA) %>%
  summarise(Freq = sum(Count)) %>%
  arrange(Freq)
df.RNA$RNA = factor(df.RNA$RNA,
                    levels = unique(df.RNA$RNA))

D = ggplot(df.RNA, aes(x="", y = Freq, fill = RNA)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  theme_classic() +
  ggtitle(paste("N =", sum(df.RNA$Freq))) +
  guides(fill=guide_legend(ncol=2)) +
  geom_text(aes(label = paste0(Freq)),
            position = position_stack(vjust=0.5),
            color = "white") +
  scale_fill_manual(values = viridis(18), name = "RNA type") +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.key.size = unit(0.05, 'in'))

####Folding results####

list.files("Figures/Figure_2_corrections_and_predictions/Predict_RNA_structure_3")

df = read.csv("Figures/Figure_2_corrections_and_predictions/Predict_RNA_structure_3/results/001_Accuracy.csv")

####PE Accuracy results ML set####

df$Program = factor(df$Program,
                    levels = c("UV_1MNaCl_data_tables",
                               "FDBI_1MNaCl_data_tables",
                               "p7kcalmol_penalty_data_tables",
                               "1kcalmol_penalty_data_tables",
                               "1p7kcalmol_data_tables",
                               "FBDI_Eco80_data_tables",
                               "mxfold2"),
                    labels = c("UV 1MNaCl",
                               "FDBI 1MNaCl",
                               "0.7 kcal/mol penalty",
                               "1.0 kcal/mol penalty",
                               "1.7 kcal/mol penalty",
                               "FBDI Eco80",
                               "mxfold2"))
E = ggplot(df %>% filter(Program %in% c("UV 1MNaCl", "FDBI 1MNaCl", "1.0 kcal/mol penalty", "FBDI Eco80")), aes(x = Program, y = Accuracy)) +
  geom_violin() +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black",
                                   angle = 45,
                                   hjust = 1,
                                   size = 8),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_blank(),
        strip.text = element_text(color = "black"),
        legend.title = element_blank(),
        legend.position = "none")

####PF Accuracy changes####

v.tables = unique(df$Program)[-which(unique(df$Program) == "UV 1MNaCl")]

l.df.change = {}

for (i in 1:length(v.tables)){
  print(i)
  df.i = df %>% filter(Program == v.tables[i]) %>% arrange(RNA) 
  df.UV.1M = df %>% filter(Program == "UV 1MNaCl") %>% filter(RNA %in% df.i$RNA) %>% arrange(RNA)
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
unique(df.change$Effect)
df.change$Effect = factor(df.change$Effect,
                          levels = c("Improved", "No change", "Got worse"))

Bryce = natparks.pals("BryceCanyon", 9)

f = ggplot(df.change %>% filter(Data.table %in% c("UV 1MNaCl", "FDBI 1MNaCl", "1.0 kcal/mol penalty", "FBDI Eco80")), aes(x = Data.table, y = N, fill = Effect)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_manual(values = Bryce[c(5,6,1)]) +
  ylab("Count") +
  theme(axis.text.x = element_text(color = "black",
                                   angle = 45,
                                   hjust = 1,
                                   size = 8),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_blank(),
        strip.text = element_text(color = "black"),
        legend.title = element_blank())

####Read populations####

list.files("Figures/Figure_2_corrections_and_predictions/Predict_RNA_structure_3/results")

df = read.csv("Figures/Figure_2_corrections_and_predictions/Predict_RNA_structure_3/results/002_Substructures.csv")

df$Program = factor(df$Program,
                    levels = c("UV_1MNaCl_data_tables",
                               "FDBI_1MNaCl_data_tables",
                               "p7kcalmol_penalty_data_tables",
                               "1kcalmol_penalty_data_tables",
                               "1p7kcalmol_data_tables",
                               "FBDI_Eco80_data_tables",
                               "mxfold2"),
                    labels = c("UV 1MNaCl",
                               "FDBI 1MNaCl",
                               "0.7 kcal/mol penalty",
                               "1.0 kcal/mol penalty",
                               "1.7 kcal/mol penalty",
                               "FBDI Eco80",
                               "mxfold2"))

Lines = unique(df$RNA)[floor(seq(1, length(unique(df$RNA)), length.out = 200))]

####PG Substructures####

G = ggplot(df %>% filter(Program %in% c("UV 1MNaCl", "FDBI 1MNaCl", "1.0 kcal/mol penalty", "FBDI Eco80")), aes(x = Program, y = Structures)) +
  geom_boxplot() +
  theme_classic() +
  ylab("Number of structures\nin ensemble") +
  theme(axis.text.x = element_text(color = "black",
                                   angle = 45,
                                   hjust = 1,
                                   size = 8),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_blank(),
        strip.text = element_text(color = "black"),
        legend.title = element_blank())

####Read in structure features####

list.files("Figures/Figure_2_corrections_and_predictions/Predict_RNA_structure_3/results")

df = read.csv("Figures/Figure_2_corrections_and_predictions/Predict_RNA_structure_3/results/003_Structural_features.csv")

head(df)

####PH% paired####

df$Percent.paired = df$Paired.n/df$Length

df$Program = factor(df$Program,
                    levels = c("UV_1MNaCl_data_tables",
                               "FDBI_1MNaCl_data_tables",
                               "p7kcalmol_penalty_data_tables",
                               "1kcalmol_penalty_data_tables",
                               "1p7kcalmol_data_tables",
                               "FBDI_Eco80_data_tables",
                               "mxfold2"),
                    labels = c("UV 1MNaCl",
                               "FDBI 1MNaCl",
                               "0.7 kcal/mol penalty",
                               "1.0 kcal/mol penalty",
                               "1.7 kcal/mol penalty",
                               "FBDI Eco80",
                               "mxfold2"))

H = ggplot(df %>% filter(Program %in% c("UV 1MNaCl", "FDBI 1MNaCl", "1.0 kcal/mol penalty", "FBDI Eco80")), aes(x = Program, y = 100*Percent.paired)) +
  geom_violin() +
  geom_boxplot() +
  theme_classic() +
  ylab("Percent of\npaired nucleotides") +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100), limits = c(0, 100)) +
  theme(axis.text.x = element_text(color = "black",
                                   angle = 45,
                                   hjust = 1,
                                   size = 8),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_blank(),
        strip.text = element_text(color = "black"),
        legend.title = element_blank())

####I% AU pairs####

df$Percet.AU.pairs = df$AU.pairs/df$Pairs

I = ggplot(df %>% filter(Program %in% c("UV 1MNaCl", "FDBI 1MNaCl", "1.0 kcal/mol penalty", "FBDI Eco80")),
           aes(x = Program, y = 100*Percet.AU.pairs)) +
  geom_violin() +
  geom_boxplot() +
  theme_classic() +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  ylab("Percent of AU pairs") +
  theme(axis.text.x = element_text(color = "black",
                                   angle = 45,
                                   hjust = 1,
                                   size = 8),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_blank(),
        strip.text = element_text(color = "black"),
        legend.title = element_blank())

####J% GC pairs####

df$Percet.GC.pairs = df$GC.pairs/df$Pairs

J = ggplot(df %>% filter(Program %in% c("UV 1MNaCl", "FDBI 1MNaCl", "1.0 kcal/mol penalty", "FBDI Eco80")),
           aes(x = Program, y = 100*Percet.GC.pairs)) +
  geom_violin() +
  geom_boxplot() +
  theme_classic() +
  ylab("Percent of GC pairs") +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  theme(axis.text.x = element_text(color = "black",
                                   angle = 45,
                                   hjust = 1,
                                   size = 8),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_blank(),
        strip.text = element_text(color = "black"),
        legend.title = element_blank())

####K% GU pairs####

df$Percet.GU.pairs = df$GU.pairs/df$Pairs

K = ggplot(df %>% filter(Program %in% c("UV 1MNaCl", "FDBI 1MNaCl", "1.0 kcal/mol penalty", "FBDI Eco80")),
           aes(x = Program, y = 100*Percet.GU.pairs)) +
  geom_violin() +
  geom_boxplot() +
  theme_classic() +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  ylab("Percent GU pairs") +
  theme(axis.text.x = element_text(color = "black",
                                   angle = 45,
                                   hjust = 1,
                                   size = 8),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_blank(),
        strip.text = element_text(color = "black"),
        legend.title = element_blank())


####Compile plots####

AB = plot_grid(A, D, nrow = 1, rel_widths = c(0.6, 1), labels = c("A   (E. coli set)", "B    (ML set)"))
CDEF = plot_grid(B, E, G, H, nrow = 1, rel_widths = c(1, 1, 1, 1), labels = c("C", "D", "E", "F"))
GHI = plot_grid(I, J, K, nrow = 1, labels = c("G", "H", "I"))
P = plot_grid(AB, CDEF, GHI, ncol = 1)

ggsave("Figures/Figure_2_corrections_and_predictions/Figure_2_predict_static.svg", P, 
       scale = 2, bg = "white", width = 5, height = 4)
