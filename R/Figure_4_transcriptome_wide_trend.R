library(tidyverse)
library(viridis)
library(ggpubr)
library(cowplot)

#####Make a plot####

df = read.csv("Figures/Figure_4_transcriptome_wide_trend/009_Reactivity_Pss.csv")

#df$Mean = df$Mean/7

df$Program = factor(df$Program,
                    levels = c("UV_1MNaCl_data_tables",
                               "FDBI_1MNaCl_data_tables",
                               "1kcalmol_penalty_data_tables",
                               "FBDI_Eco80_data_tables"),
                    labels = c("UV 1MNaCl",
                               "FDBI 1MNaCl",
                               "1.0 kcal/mol penalty",
                               "FBDI Eco80"))

df.l = data.frame("X" = 10^seq(-5, 0,length.out = 100),
                "Y" =  1.0 + 7*10^seq(-5, 0,length.out = 100)) %>% filter(Y <= 7)

P1 = ggplot() +
  facet_wrap(~Program) +
  geom_point(data = df, mapping = aes(x = Pss, y = Mean), alpha = 0.05) +
  geom_line(data = df.l, mapping = aes(x = X, y = Y), color = "red") +
  theme_classic() +
  scale_y_continuous(limits = c(0, 7)) +
  scale_x_continuous(trans = "log10", breaks = c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1.0),
                     labels = c("", "0.0001", "", "0.01", "", 1.0)) +
  theme(axis.text = element_text(color = "black")) +
  ylab("Reactivity")

score = function(a = 1, df.i = df){
  Threshold = 7*df.i$Pss + a
  df.i$Wrong = 0
  df.i$Wrong[which(df.i$Mean > Threshold)] = 1
  df.c = df.i %>%
    group_by(Program) %>%
    summarise(Count = sum(Wrong))
  df.c$a = a
  df.c$Percent.change = df.c$Count/df.c$Count[which(df.c$Program == "UV 1MNaCl")]
  output = df.c
}

a = seq(0, 4, length.out = 25)
l.df.c = lapply(a, score)
df.c = bind_rows(l.df.c)


P2 = ggplot(df.c, aes(x = a, y = Count, color = Program, shape = Program)) +
  geom_point() +
  geom_line() +
  theme_classic() +
  scale_color_manual(values = viridis(6, option = "A")) +
  scale_y_continuous(trans = "log10") +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.75, 0.85),
        legend.title = element_blank(),
        legend.background = element_blank()) +
  xlab("Threshold to be considered wrong (a)") +
  ylab("Incorrectly predicted\ndouble stranded nucleotides")

P3 = ggplot(df.c, aes(x = a, y = 100*Percent.change, color = Program, shape = Program)) +
  geom_point() +
  geom_line() +
  theme_classic() +
  scale_y_continuous() +
  scale_color_manual(values = viridis(6, option = "A")) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.28, 0.6),
        legend.title = element_blank(),
        legend.background = element_blank()) +
  xlab("Threshold to be considered wrong (a)") +
  ylab("%Decrease in incorrectly\npredicted double stranded nucleotides")

PBC =  plot_grid(P2, P3, ncol = 1, labels = c("B", "C"), align = "v")
P = plot_grid(P1, PBC, nrow = 1, labels = c("A", ""), rel_widths = c(1.4, 1))

ggsave("Figures/Figure_4_transcriptome_wide_trend/Figure_4_transcriptome_wide_trends.png",
       P,
       width = 5,
       height = 3,
       units = "in",
       scale = 2.0,
       bg = "white",
       dpi = 600)

 df.result = score(1)

print(df.result)
