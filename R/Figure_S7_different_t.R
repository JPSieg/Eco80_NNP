library(tidyverse)
library(viridis)
library(ggpubr)
library(cowplot)

#####Make a plot####

list.files("Figures")

df = read.csv("Figures/Figure_4_transcriptome_wide_trend/009_Reactivity_Pss.csv")

head(df)
#df$Mean = df$Mean/7
unique(df$Program)

df$Program = factor(df$Program,
                    levels = c("UV_1MNaCl_data_tables",
                               "FDBI_1MNaCl_data_tables",
                               "1kcalmol_penalty_data_tables",
                               "FBDI_Eco80_data_tables"),
                    labels = c("UV 1MNaCl",
                               "FDBI 1MNaCl",
                               "1.0 kcal/mol penalty",
                               "FBDI Eco80"))

df.l = data.frame("X" = 10^seq(-5, 0,length.out = 100))
df.l$Y = 2*7/(1 + 10^(1 - df.l$X))


P1 = ggplot() +
  facet_wrap(~Program) +
  geom_point(data = df, mapping = aes(x = Pss, y = Mean), alpha = 0.05) +
  geom_line(data = df.l, mapping = aes(x = X, y = Y), color = "red") +
  annotate("text", x = 0.1, y = 0.9, label = "T(R = 7)", color = "red") +
  theme_classic() +
  scale_y_continuous() +
  scale_x_continuous(trans = "log10", ) +
  theme(axis.text = element_text(color = "black")) +
  ylab("Reactivity")

####Score data####

score = function(a = 7, df.i = df){
  Threshold = 2*a/(1 + 10^(1 - df.i$Pss))
  df.i$Wrong = 0
  df.i$Wrong[which(df.i$Mean > Threshold)] = 1
  df.c = df.i %>%
    group_by(Program) %>%
    summarise(Count = sum(Wrong))
  df.c$a = a
  df.c$Percent.change = df.c$Count/df.c$Count[which(df.c$Program == "UV 1MNaCl")]
  output = df.c
}


df.a.7 = score(7)
df.a.1 = score(1)
df.a.7
df.a.1

a = seq(0.5, 10, length.out = 20)
l.df.c = lapply(a, score)
df.c = bind_rows(l.df.c)



P2 = ggplot() +
  geom_point(data = df.c,
             mapping = aes(x = a, y = Count, color = Program, shape = Program)) +
  geom_line(data = df.c,
            mapping = aes(x = a, y = Count, color = Program)) +
  theme_classic() +
  scale_color_manual(values = viridis(7, option = "A")) +
  scale_y_continuous(trans = "log10") +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.25, 0.25),
        legend.title = element_blank(),
        legend.background = element_blank()) +
  xlab("Reactivity threshold to be considered wrong (R)") +
  ylab("Number of wrongly predicted\ndouble stranded nucleotides")

P3 = ggplot() +
  geom_point(data = df.c, mapping = aes(x = a, y = Percent.change, color = Program, shape = Program)) +
  geom_line(data = df.c, mapping = aes(x = a, y = Percent.change, color = Program)) +
  theme_classic() +
  scale_y_continuous() +
  scale_color_manual(values = viridis(7, option = "A")) +
  theme(axis.text = element_text(color = "black"),
        legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_blank()) +
  xlab("Reactivity threshold to be considered wrong (R)") +
  ylab("Fold change in wrongly predicted\ndouble stranded nucleotides")

P23 = plot_grid(P2, P3, ncol = 1, labels = c("B", "C"))

P = plot_grid(P1, P23, nrow = 1, labels = c("A", ""), rel_widths = c(1.5, 1))

list.files("Figures/")
ggsave("Figures/Figure_S7_different_t/SI_Figure_x_different_t.png", bg = "white", dpi = 400, width = 4.5, height = 3, scale = 2.5, units = "in")

print(df.a.7)
