library(tidyverse)
library(svd)
library(MeltR)
library(ggpubr)
library(ggrepel)
library(cowplot)

####dG Eco80####

####Load in helix energy data####

list.files("Tables/SI_Table_2_pvalues_and_established_NN_comparison/")

df = read.csv("Tables/SI_Table_2_pvalues_and_established_NN_comparison/SI_Table_2_NN_matrix_and_data_Eco80.csv")
df.pub = read.csv("Tables/SI_Table_2_pvalues_and_established_NN_comparison/NN_parameters.csv") %>% filter(Condition == "Eco80")

####Generate possible F.Q contributions to NN model####

tib.F.Q = expand_grid(seq((3.67 - 2), (3.67 + 2), length.out = 10), seq((4.38 - 2), (4.38 + 2), length.out = 10))
colnames(tib.F.Q) = c("FAMC.GBHQ1", "FAMU.ABHQ1")

####Function that calculates NN parameters####

NN.calc = function(FAMC.GBHQ1 =  tib.F.Q$FAMC.GBHQ1[1], FAMU.ABHQ1 = tib.F.Q$FAMU.ABHQ1[1], df){
  df$dG[which(df$FAMC.GBHQ1 == 1)] = df$dG[which(df$FAMC.GBHQ1 == 1)]  + FAMC.GBHQ1
  
  df$dG[which(df$FAMC.GBHQ1 != 1)] = df$dG[which(df$FAMC.GBHQ1 != 1)]  + FAMU.ABHQ1
  
  fit = lm(dG ~ AA.UU + AU.AU + UA.UA + CU.AG + CA.UG + GU.AC + GA.UC +CG.CG + GG.CC + GC.GC + Term.AU,
           df,
           weights = 1/df$dG.error)
  
  df.Eco80 = data.frame(coef(summary(fit)))
  df.Eco80 = bind_cols(df.Eco80, df.pub)
  
  df.Eco80$Parameter = row.names(df.Eco80)
  df.Eco80$Parameter[1] = "Initiation"
  df.Eco80$Condition = "Eco80"
  
  R = cor(df.Eco80$Estimate...1[-1], df.pub$Estimate[-1])
  l.output = list(df.Eco80, R)
}


####Calculate correlation####

R2 = c()
l.out = {}

for (i in 1:nrow(tib.F.Q)){
  print(i)
  l = NN.calc(tib.F.Q$FAMC.GBHQ1[i], tib.F.Q$FAMU.ABHQ1[i], df)
  l.out[[i]] = l[[1]]
  R2[i] = l[[2]]
}

tib.F.Q$R2 = R2
tib.F.Q$i = 1:nrow(tib.F.Q)

#####PA correlation plot####

PA = ggplot(l.out[[34]], aes(x = Estimate...5, xmin = Estimate...5 - Std..Error...6, xmax = Estimate...5 + Std..Error...6,
                        y = Estimate...1, ymin = Estimate...1 - Std..Error...2, ymax = Estimate...5 + Std..Error...2)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_pointrange() +
  geom_errorbarh() +
  annotate("text", x = -2.5, y = 1,
           label = paste("R2 =", round(tib.F.Q$R2[34], digits = 2)),
           size = 3) +
  coord_fixed(ylim = c(-3, 3), xlim = c(-3, 3)) +
  theme_classic() +
  xlab("NNP from Table 2 (kcal/mol)") +
  ylab("NNP with FAM/BHQ1 error (kcal/mol)") +
  theme(axis.text = element_text(color = "black")) 

####PB Plot heat map####

circleFun <- function(center = c(3.67,4.38), diameter = 2.5, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

dat = circleFun()

PB = ggplot() +
  geom_tile(data = tib.F.Q, mapping = aes(x = FAMC.GBHQ1, y = FAMU.ABHQ1, fill = R2)) +
  geom_text(data = tib.F.Q, mapping = aes(x = FAMC.GBHQ1, y = FAMU.ABHQ1, label = round(R2, digits = 2)), size = 3) +
  geom_path(data = dat, mapping = aes(x, y)) +
  coord_fixed() +
  theme_classic() +
  scale_fill_viridis_b() +
  theme(axis.text = element_text(color = "black")) +
  xlab(expression("5'-FAMC/GBHQ1-3' \u0394G\u00B037 (kcal/mol)")) +
  ylab(expression("5'-FAMC/GBHQ1-3' \u0394G\u00B037 (kcal/mol)"))

####Make plot####

P = plot_grid(PA, PB, rel_widths = c(1, 1.2), labels = c("A", "B"))

list.files("Figures/SI_Figure_2_Initiation_parameter_error")

ggsave("Figures/SI_Figure_2_Initiation_parameter_error/SI_Figure_2_Initiaion_parameter_error_analysis.svg", P, 
       units = "in", width = 4.5, height = 2,bg = "white", scale = 2)
