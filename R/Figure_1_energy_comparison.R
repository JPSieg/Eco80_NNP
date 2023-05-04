library(tidyverse)
library(svd)
library(MeltR)
library(ggpubr)
library(ggrepel)
library(NatParksPalettes)
library(cowplot)

####Generate color pallets####

Bryce = natparks.pals("BryceCanyon", 9)

str(Bryce)

####Figure 1A####

list.files("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/")

df.Eco80 = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/js7005_Helix_V_Eco80.csv")

fit.Eco80 = meltR.F(df.Eco80,
                 Kd_range = c(10, 400),
                 Kd_error_quantile = 1)

PA.gfit = ggplot(fit.Eco80$df.globalfit, aes(x = B, color = Temperature, group = Reading)) +
  geom_point(mapping = aes(y = Emission)) +
  geom_line(mapping = aes(y = Model)) +
  scale_color_viridis_c(option = "C", name = "Temp. (\u00B0C)") +
  scale_x_continuous() +
  theme_classic() +
  xlab("[BHQ1] (nM)") +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.9, 0.7))

####Figure 1B ####

list.files("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/")

df.1MNaCl = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/js6027_Helix_V.csv")

fit.1M = meltR.F(df.1MNaCl,
                 Kd_range = c(10, 200),
                 Kd_error_quantile = 0.5)

df.Eco80 = fit.Eco80$df.vanthoff
df.1M = fit.1M$df.vanthoff

df.Eco80$Condition = "Eco80"
df.1M$Condition = "1 M NaCl"

df = bind_rows(df.Eco80, df.1M)

PB.VHplot = ggplot(df, aes(x = invT, color = Condition, group = Condition)) +
  geom_pointrange(mapping = aes(y = lnK, ymin = lnK - SE.lnK, ymax = lnK + SE.lnK)) +
  geom_line(mapping = aes(y = Model)) +
  annotate("text", x = 0.00309, y = -16.2, label = "Less\nStable",
           size = 3) +
  geom_segment(aes(x = 0.0031005, y = -16.8, xend = 0.0031005, yend = -16),
               arrow = arrow(length = unit(0.25, "cm")),
               color = "black") +
  theme_classic() +
  scale_x_continuous(limits = c(0.00305, 0.00320), breaks = c(0.00305, 0.00310, 0.00315, 0.00320)) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.75, 0.75)) +
  scale_color_manual(values = Bryce[c(2, 5)]) +
  ylab(expression(ln*"["*" "*K[D]*" "*(M)*" "*"]")) +
  xlab("1/Temperature (K)")

####Figure 1C dG comparison####

list.files("Tables/SI_Table_1_fitting_statistics")

df = read.csv("Tables/SI_Table_1_fitting_statistics/SI_Table_1_data.csv")

df.Eco80 = df %>% 
  filter(Condition == "Eco80") %>%
  filter(Method == "2 Global fit") %>%
  arrange(Helix)
nrow(df.Eco80)
df.1M = df %>% 
  filter(Condition == "1 M NaCl") %>%
  filter(Helix %in% df.Eco80$Helix) %>%
  filter(Method == "2 Global fit") %>%
  arrange(Helix)
nrow(df.1M)

Helix = df.Eco80$Helix
f = df.Eco80$Sequence.F
r = df.Eco80$Sequence.R
dG.Eco80 = df.Eco80$G
dG.error.Eco80 = df.Eco80$SE.G
dG.1M = df.1M$G
dG.error.1M = df.1M$SE.G
ddG = dG.Eco80 - dG.1M
ddG.error = sqrt((dG.error.Eco80)^2 + (dG.error.1M)^2)

df = data.frame(Helix, f, r, dG.Eco80, dG.error.Eco80, dG.1M, dG.error.1M, ddG, ddG.error)

mean.ddG = mean(df$ddG)

SS.res = sum((df$dG.Eco80 - (df$dG.1M + mean.ddG))^2)
SS.tot = sum((df$dG.Eco80 - mean(df$dG.Eco80))^2)

R2.plus.mean.ddG = 1 - SS.res/SS.tot

PC.dG_comparison = ggplot(df, aes(x = dG.1M, xmin = dG.1M - dG.error.1M, xmax = dG.1M + dG.error.1M,
               y = dG.Eco80, ymin = dG.Eco80 - dG.error.Eco80, ymax = dG.Eco80 + dG.error.Eco80)) +
  geom_abline(slope = 1, intercept = 0, color = "black", linewidth = 1.25) +
  geom_abline(slope = 1, intercept = mean.ddG, color = Bryce[c(5)], linewidth = 1.25) +
  geom_pointrange(color = Bryce[c(5)]) +
  geom_errorbarh(color = Bryce[c(5)]) +
  annotate("text", x = -14.5, y = -10,
           label = paste("Less stable\nin Eco80"),
           color = "black",
           size = 3) +
  annotate("text", x = -10.5, y = -13,
           label = paste("More stable\nin Eco80"),
           color = "black",
           size = 3) +
  annotate("text", x = -11.5, y = -14.5,
           label = paste("y = x +", round(mean.ddG, digits = 2), "kcal/mol\nR2 =", round(R2.plus.mean.ddG, digits = 2)),
           color = Bryce[5],
           size = 3) +
  theme_classic() +
  coord_fixed(xlim = c(-16, -8), ylim = c(-16, -8)) +
  scale_x_continuous(breaks = c(-16, -14, -12, -10, -8, -15, -13, -11, -9), labels = c("-16", "-14", "-12", "-10", "-8", "", "", "", "")) +
  scale_y_continuous(breaks = c(-16, -14, -12, -10, -8, -15, -13, -11, -9), labels = c("-16", "-14", "-12", "-10", "-8", "", "", "", "")) +
  theme(axis.text = element_text(color = "black")) +
  xlab(expression("1 M NaCl \u0394G\u00B037 (kcal/mol)")) +
  ylab(expression("Eco80 \u0394G\u00B037 (kcal/mol)"))

#####Figure 1D ddG versus AU content####

head(df)

Length = nchar(df$f)
AU.content = 100*(lengths(regmatches(df$f, gregexpr("A", df$f))) + lengths(regmatches(df$f, gregexpr("U", df$f))))/nchar(df$f)

df$Length = Length
df$AU.content = AU.content

fit.AU = lm(ddG~AU.content, df, weights = 1/ddG.error)

SS.res = sum((df$ddG - predict(fit.AU))^2)
SS.tot = sum((df$ddG - mean(df$ddG))^2)

R2.AU.content = 1 - SS.res/SS.tot

PD.dG.v.AU = ggplot(df, aes(x = AU.content, y = ddG, ymin = ddG - ddG.error, ymax = ddG + ddG.error)) +
  annotate("text", x = 70, y = -0.25,
           label = paste("y = ", round(coef(fit.AU)[2], digits = 2),"x +", round(coef(fit.AU)[1], digits = 2), "\nR2 =", round(R2.AU.content, digits = 2)),
           color = Bryce[5],
           size = 3) +
  geom_abline(slope = coef(fit.AU)[2], intercept = coef(fit.AU)[1], color = Bryce[c(5)], linewidth = 1.25) +
  geom_pointrange(color = Bryce[c(5)]) +
  theme_classic() + theme(axis.text = element_text(color = "black")) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
  xlab("%AU content") +
  ylab("\u0394\u0394G\u00B037 (kcal/mol)")

#####Figure 1E ddG versus length####

fit.length = lm(ddG ~ Length, df, weights = 1/ddG.error)

SS.res = sum((df$ddG - predict(fit.length))^2)
SS.tot = sum((df$ddG - mean(df$ddG))^2)

R2.Length = 1 - SS.res/SS.tot

PE.dG.v.L = ggplot(df, aes(x = Length, y = ddG, ymin = ddG - ddG.error, ymax = ddG + ddG.error)) +
  annotate("text", x = 6.5, y = 1.5,
           label = paste("y = ", round(coef(fit.length)[2], digits = 2),"x +", round(coef(fit.length)[1], digits = 2), "\nR2 =", round(R2.Length, digits = 2)),
           color = Bryce[5],
           size = 3) +
  geom_abline(slope = coef(fit.length)[2], intercept = coef(fit.length)[1], color = Bryce[c(5)], size = 1.25) +
  geom_pointrange(color = Bryce[c(5)]) +
  scale_x_continuous(breaks = c(6, 7, 8)) +
  theme_classic()  +
  theme(axis.text = element_text(color = "black")) +
  xlab("Helix length") +
  ylab("\u0394\u0394G\u00B037 (kcal/mol)")

####Figure 1F RSS by model####

#read in NN models

list.files("Tables/")

df.NN.Sieg = read.csv("Tables/Table_2_NN_parameters/NN_parameters.csv") %>%
  select(Parameter, Estimate,  Std..Error, Condition)
colnames(df.NN.Sieg) = c("Parameter", "Value", "error", "Condition")
df.NN = df.NN.Sieg

#Read in 1M NaCl model

df = read.csv("Tables/SI_Table_2_pvalues_and_established_NN_comparison/SI_Table_2_NN_matrix_and_data_1MNaCl.csv")

#filter data and subtract out fluorophore quencher contributions#

df$dG[which(df$FAMC.GBHQ1 == 1)] = df$dG[which(df$FAMC.GBHQ1 == 1)]  + 3.67

df$dG[which(df$FAMC.GBHQ1 != 1)] = df$dG[which(df$FAMC.GBHQ1 != 1)]  + 4.38

Conditions = unique(df.NN$Condition)

l.df.1M = {}

for (i in 1:length(Conditions)){
  df.NN.i = df.NN %>% filter(Condition == Conditions[i])
  S = as.matrix(df %>% select(Initiation, AA.UU, AU.AU, UA.UA,
                              CU.AG,CA.UG, GU.AC, GA.UC, CG.CG,
                              GG.CC, GC.GC, Term.AU))
  
  
  G.NN = matrix(data = df.NN.i$Value, nrow = 12, ncol = 1)
  dG.known = S %*% G.NN
  
  l.df.1M[[i]] = df %>% select(Helix, dG, dG.error)
  l.df.1M[[i]]$Parameter = Conditions[i]
  l.df.1M[[i]]$dG.known = dG.known
}

df.1M = bind_rows(l.df.1M)
df.1M$Condition = "1 M NaCl"

#Read in Eco80 model

list.files("Tables/SI_Table_2_pvalues_and_established_NN_comparison")

df = read.csv("Tables/SI_Table_2_pvalues_and_established_NN_comparison/SI_Table_2_NN_matrix_and_data_Eco80.csv")

#filter data and subtract out fluorophore quencher contributions#

df$dG[which(df$FAMC.GBHQ1 == 1)] = df$dG[which(df$FAMC.GBHQ1 == 1)]  + 3.67

df$dG[which(df$FAMC.GBHQ1 != 1)] = df$dG[which(df$FAMC.GBHQ1 != 1)]  + 4.38

Conditions = unique(df.NN$Condition)

l.df.Eco80 = {}

for (i in 1:length(Conditions)){
  print(i)
  df.NN.i = df.NN %>% filter(Condition == Conditions[i])
  S = as.matrix(df %>% select(Initiation, AA.UU, AU.AU, UA.UA,
                              CU.AG,CA.UG, GU.AC, GA.UC, CG.CG,
                              GG.CC, GC.GC, Term.AU))
  
  
  G.NN = matrix(data = df.NN.i$Value, nrow = 12, ncol = 1)
  dG.known = S %*% G.NN
  
  l.df.Eco80[[i]] = df %>% select(Helix, dG, dG.error)
  l.df.Eco80[[i]]$Parameter = Conditions[i]
  l.df.Eco80[[i]]$dG.known = dG.known
}

df.Eco80 = bind_rows(l.df.Eco80)
df.Eco80$Condition = "Eco80"


df = bind_rows(df.1M, df.Eco80)
df$Error = df$dG - df$dG.known

library(ggbeeswarm)

df$Parameter = factor(df$Parameter,
                      levels = c("1 M NaCl", "Eco80", "Turner 1989 1 M NaCl", "Adams 2019 1M NaCl 20% PEG 200","Ghosh 2023 100 mM NaCl PEG 200"),
                      labels = c("FDBI 1 M NaCl", "FDBI Eco80", "UV 1 M NaCl", "UV 1M NaCl 20% PEG 200","UV 100 mM NaCl PEG 200"))

df.mean.1M = df %>% filter(Condition == "1 M NaCl") %>% group_by(Parameter) %>% summarise(Mean =  mean(Error), SD = sd(Error))
df.mean.1M$Condition = "1 M NaCl"
df.mean.Eco80 = df %>% filter(Condition == "Eco80") %>% group_by(Parameter) %>% summarise(Mean =  mean(Error), SD = sd(Error))
df.mean.Eco80$Condition = "Eco80"
df.mean = bind_rows(df.mean.1M, df.mean.Eco80)
head(df.mean)

df.mean$Labels = paste(round(df.mean$Mean, digits = 1), "\n(", round(df.mean$SD, digits = 1), ")", sep = "")

df.mean = df.mean %>% filter(Condition == "Eco80") %>% filter(Parameter %in% c("FDBI 1 M NaCl", "FDBI Eco80", "UV 1 M NaCl"))

df = df %>% filter(Condition == "Eco80") %>% filter(Parameter %in% c("FDBI 1 M NaCl", "FDBI Eco80", "UV 1 M NaCl"))

PF.residuals = ggplot() +
  geom_text(data = df.mean, mapping = aes(x = Parameter, y = 3.0, label = Labels), size = 3) +
  geom_hline(yintercept = 0) +
  geom_bar(data = df.mean,
           mapping = aes(x = Parameter, y = Mean, fill = Parameter),
           stat = "identity") + 
  geom_errorbar(data = df.mean,
           mapping = aes(x = Parameter, ymin = Mean - SD, ymax = Mean + SD)) +
  geom_beeswarm(data = df,
                mapping = aes(x = Parameter, y = Error), size = 1) +
  facet_wrap(~Condition, ncol = 2) +
  theme_classic() +
  scale_fill_manual(values = Bryce[c(1, 5, 2)]) +
  scale_y_continuous(limits = c(-1, 3.5), breaks = c(-1, 0, 1, 2, 3)) +
  theme(axis.text.x = element_text(color = "black",
                                   angle = 45,
                                   hjust = 1,
                                   size = 8),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_blank(),
        strip.text = element_text(color = "black"),
        legend.position = "none") +
  ylab("NN model residuals\n(kcal/mol)")


####Figure 1G####

list.files()

df.Turner = read.csv("Tables/SI_Table_2_pvalues_and_established_NN_comparison/PublishedNN_parameters.csv") %>%
  filter(Condition == "Turner 1989 1 M NaCl")

df.NN = bind_rows(df.NN, df.Turner)

df.NN$Parameter = factor(df.NN$Parameter,
                         levels = unique(df.NN$Parameter))

df.NN$Condition = factor(df.NN$Condition,
                         levels = c("1 M NaCl", "Eco80", "Turner 1989 1 M NaCl"),
                         labels = c("FDBI 1 M NaCl", "FDBI Eco80", "Turner 1989 1 M NaCl"))

PG.NN = ggplot(df.NN, aes(x = Condition, y = Value, ymin = Value - error, ymax = Value + error, group = Condition, fill = Condition)) +
  geom_bar(stat = "identity") +
  geom_errorbar() +
  facet_wrap(~Parameter, nrow = 1) +
  theme_classic() +
  scale_fill_manual(values = Bryce[c(1, 5, 2)]) +
  theme(axis.text.x = element_text(color = "black",
                                   angle = 45,
                                   hjust = 1,
                                   size = 8),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_blank(),
        strip.text = element_text(color = "black"),
        legend.title = element_blank(),
        legend.position = "none") +
  ylab("\u0394G\u00B037 (kcal/mol)")


####Make final figure####

P.top = plot_grid(PA.gfit, PB.VHplot, PC.dG_comparison, nrow = 1, labels = c("A","B","C"), label_size = 14)
P.middle = plot_grid(PD.dG.v.AU, PE.dG.v.L, PF.residuals, nrow = 1, labels = c("D","E","F"), label_size = 14)
P.bottom = PG.NN

P = plot_grid(P.top, P.middle, P.bottom, ncol = 1, labels = c("","","G"), label_size = 14, rel_heights = c(1, 1, 1.2))

ggsave("Figures/Figure_1_Nearest_neighbor_parameters/Figure_1_energy_comparison.svg",
       scale = 2, bg = "white", width = 5.4, height = 4.5)

