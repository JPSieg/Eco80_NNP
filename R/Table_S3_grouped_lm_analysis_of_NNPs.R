library(tidyverse)
library(svd)
library(MeltR)
library(ggpubr)
library(ggrepel)

####dG 1 M NaCl####

####Load in helix energy data####

list.files("Tables/SI_Table_2_pvalues_and_established_NN_comparison")

df = read.csv("Tables/SI_Table_2_pvalues_and_established_NN_comparison/SI_Table_2_NN_matrix_and_data_1MNaCl.csv")

####filter data and subtract out fluorophore quencher contributions####

df$dG[which(df$FAMC.GBHQ1 == 1)] = df$dG[which(df$FAMC.GBHQ1 == 1)]  + 3.67

df$dG[which(df$FAMC.GBHQ1 != 1)] = df$dG[which(df$FAMC.GBHQ1 != 1)]  + 4.38

####plot comparison to known####

S = as.matrix(df %>% select(Initiation, AA.UU, AU.AU, UA.UA,
                            CU.AG,CA.UG, GU.AC, GA.UC, CG.CG,
                            GG.CC, GC.GC, Term.AU))


G.NN = matrix(data = c(4.09, -0.9, -1.1, -1.3, -2.1, -2.1, -2.2, -2.1, -2.4, -3.3, -3.4, 0.45), nrow = 12, ncol = 1)
dG.known = S %*% G.NN
df$dG.known = dG.known

ggplot(df, aes(x = dG.known, y = dG, ymin = dG - dG.error, ymax = dG + dG.error, label = Helix)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_pointrange() +
  geom_text_repel() +
  coord_fixed() +
  theme_classic()

####Remove 3 outliers####

####Determine 1 M NN using linear regression####

colnames(df)

fit = lm(dG ~ AA.UU + AU.AU + UA.UA + CU.AG + CA.UG + GU.AC + GA.UC +CG.CG + GG.CC + GC.GC + Term.AU,
         df,
         weights = 1/df$dG.error)
summary(fit)


df.1M = data.frame(coef(summary(fit)))

df.1M$Parameter = row.names(df.1M)
df.1M$Parameter[1] = "Initiation"
df.1M$Condition = "1 M NaCl"

####dG Eco80####

####Load in helix energy data####

list.files("Tables/SI_Table_2_pvalues_and_established_NN_comparison/")

df = read.csv("Tables/SI_Table_2_pvalues_and_established_NN_comparison/SI_Table_2_NN_matrix_and_data_Eco80.csv")

####filter data and subtract out fluorophore quencher contributions####

df$dG[which(df$FAMC.GBHQ1 == 1)] = df$dG[which(df$FAMC.GBHQ1 == 1)]  + 3.67

df$dG[which(df$FAMC.GBHQ1 != 1)] = df$dG[which(df$FAMC.GBHQ1 != 1)]  + 4.38

fit = lm(dG ~ AA.UU + AU.AU + UA.UA + CU.AG + CA.UG + GU.AC + GA.UC +CG.CG + GG.CC + GC.GC + Term.AU,
         df,
         weights = 1/df$dG.error)

summary(fit)

df.Eco80 = data.frame(coef(summary(fit)))

df.Eco80$Parameter = row.names(df.Eco80)
df.Eco80$Parameter[1] = "Initiation"
df.Eco80$Condition = "Eco80"

df.result = bind_rows(df.1M, df.Eco80)

####Calculate significance####

df.Eco80 = read.csv("Tables/SI_Table_2_pvalues_and_established_NN_comparison/SI_Table_2_NN_matrix_and_data_Eco80.csv")
df.Eco80$Condition = "Eco80"
df.1M = read.csv("Tables/SI_Table_2_pvalues_and_established_NN_comparison/SI_Table_2_NN_matrix_and_data_1MNaCl.csv")
df.1M$Condition = "1 M NaCl"

df = bind_rows(df.1M, df.Eco80)

df$dG[which(df$FAMC.GBHQ1 == 1)] = df$dG[which(df$FAMC.GBHQ1 == 1)]  + 3.67

df$dG[which(df$FAMC.GBHQ1 != 1)] = df$dG[which(df$FAMC.GBHQ1 != 1)]  + 4.38

l.df.stats = {}

#Model assuming no grouping based on condition

lm.No.group = lm(dG ~ AA.UU + AU.AU + UA.UA + CU.AG + CA.UG + GU.AC + GA.UC +CG.CG + GG.CC + GC.GC + Term.AU,
     df,
     weights = 1/df$dG.error)

#Model assuming grouping of every variable including the initiation (intercept)

lm.group = lm(dG ~ AA.UU*Condition + AU.AU*Condition + UA.UA*Condition + CU.AG*Condition + CA.UG*Condition + GU.AC*Condition + GA.UC*Condition +CG.CG*Condition + GG.CC*Condition + GC.GC*Condition + Term.AU*Condition - 1,
              df,
              weights = 1/df$dG.error)

summary(lm.group)


A = anova(lm.group, lm.No.group)
l.df.stats[[1]] = data.frame("Model.1" = A$Res.Df[1],
                             "Model.2" = A$Res.Df[2],
                             "F.value" = A$F[2],
                             "p-value" = A$`Pr(>F)`[2])
l.df.stats[[1]]$Parameter = "Grouped model to no group"

#AA.UU 2

lm.NN = lm(dG ~ AA.UU + AU.AU*Condition + UA.UA*Condition + CU.AG*Condition + CA.UG*Condition + GU.AC*Condition + GA.UC*Condition +CG.CG*Condition + GG.CC*Condition + GC.GC*Condition + Term.AU*Condition - 1,
              df,
              weights = 1/df$dG.error)


A = anova(lm.group, lm.NN)
l.df.stats[[2]] = data.frame("Model.1" = A$Res.Df[1],
                             "Model.2" = A$Res.Df[2],
                             "F.value" = A$F[2],
                             "p-value" = A$`Pr(>F)`[2])
l.df.stats[[2]]$Parameter = "Grouped model to no group for AA.UU"

#AU.AU 3

lm.NN = lm(dG ~ AA.UU*Condition + AU.AU + UA.UA*Condition + CU.AG*Condition + CA.UG*Condition + GU.AC*Condition + GA.UC*Condition +CG.CG*Condition + GG.CC*Condition + GC.GC*Condition + Term.AU*Condition - 1,
           df,
           weights = 1/df$dG.error)


A = anova(lm.group, lm.NN)
l.df.stats[[3]] = data.frame("Model.1" = A$Res.Df[1],
                             "Model.2" = A$Res.Df[2],
                             "F.value" = A$F[2],
                             "p-value" = A$`Pr(>F)`[2])
l.df.stats[[3]]$Parameter = "Grouped model to no group for AU.AU"

#UA.UA 4

lm.NN = lm(dG ~ AA.UU*Condition + AU.AU*Condition + UA.UA + CU.AG*Condition + CA.UG*Condition + GU.AC*Condition + GA.UC*Condition +CG.CG*Condition + GG.CC*Condition + GC.GC*Condition + Term.AU*Condition - 1,
           df,
           weights = 1/df$dG.error)


A = anova(lm.group, lm.NN)
l.df.stats[[4]] = data.frame("Model.1" = A$Res.Df[1],
                             "Model.2" = A$Res.Df[2],
                             "F.value" = A$F[2],
                             "p-value" = A$`Pr(>F)`[2])
l.df.stats[[4]]$Parameter = "Grouped model to no group for UA.UA"

#CU.AG 5

lm.NN = lm(dG ~ AA.UU*Condition + AU.AU*Condition + UA.UA*Condition + CU.AG + CA.UG*Condition + GU.AC*Condition + GA.UC*Condition +CG.CG*Condition + GG.CC*Condition + GC.GC*Condition + Term.AU*Condition - 1,
           df,
           weights = 1/df$dG.error)


A = anova(lm.group, lm.NN)
l.df.stats[[5]] = data.frame("Model.1" = A$Res.Df[1],
                             "Model.2" = A$Res.Df[2],
                             "F.value" = A$F[2],
                             "p-value" = A$`Pr(>F)`[2])
l.df.stats[[5]]$Parameter = "Grouped model to no group for CU.AG"

#CA.UG 6

lm.NN = lm(dG ~ AA.UU*Condition + AU.AU*Condition + UA.UA*Condition + CU.AG*Condition + CA.UG + GU.AC*Condition + GA.UC*Condition +CG.CG*Condition + GG.CC*Condition + GC.GC*Condition + Term.AU*Condition - 1,
           df,
           weights = 1/df$dG.error)


A = anova(lm.group, lm.NN)
l.df.stats[[6]] = data.frame("Model.1" = A$Res.Df[1],
                             "Model.2" = A$Res.Df[2],
                             "F.value" = A$F[2],
                             "p-value" = A$`Pr(>F)`[2])
l.df.stats[[6]]$Parameter = "Grouped model to no group for CA.UG"

#GU.AC 7

lm.NN = lm(dG ~ AA.UU*Condition + AU.AU*Condition + UA.UA*Condition + CU.AG*Condition + CA.UG*Condition + GU.AC + GA.UC*Condition +CG.CG*Condition + GG.CC*Condition + GC.GC*Condition + Term.AU*Condition - 1,
           df,
           weights = 1/df$dG.error)

A = anova(lm.group, lm.NN)
l.df.stats[[7]] = data.frame("Model.1" = A$Res.Df[1],
                             "Model.2" = A$Res.Df[2],
                             "F.value" = A$F[2],
                             "p-value" = A$`Pr(>F)`[2])
l.df.stats[[7]]$Parameter = "Grouped model to no group for GU.AC"

#GA.UC 8

lm.NN = lm(dG ~ AA.UU*Condition + AU.AU*Condition + UA.UA*Condition + CU.AG*Condition + CA.UG*Condition + GU.AC*Condition + GA.UC +CG.CG*Condition + GG.CC*Condition + GC.GC*Condition + Term.AU*Condition - 1,
           df,
           weights = 1/df$dG.error)

A = anova(lm.group, lm.NN)
l.df.stats[[8]] = data.frame("Model.1" = A$Res.Df[1],
                             "Model.2" = A$Res.Df[2],
                             "F.value" = A$F[2],
                             "p-value" = A$`Pr(>F)`[2])
l.df.stats[[8]]$Parameter = "Grouped model to no group for GA.UC"

#CG.CG 9

lm.NN = lm(dG ~ AA.UU*Condition + AU.AU*Condition + UA.UA*Condition + CU.AG*Condition + CA.UG*Condition + GU.AC*Condition + GA.UC*Condition + CG.CG + GG.CC*Condition + GC.GC*Condition + Term.AU*Condition - 1,
           df,
           weights = 1/df$dG.error)

A = anova(lm.group, lm.NN)
l.df.stats[[9]] = data.frame("Model.1" = A$Res.Df[1],
                             "Model.2" = A$Res.Df[2],
                             "F.value" = A$F[2],
                             "p-value" = A$`Pr(>F)`[2])
l.df.stats[[9]]$Parameter = "Grouped model to no group for CG.CG"

#GG.CC 10

lm.NN = lm(dG ~ AA.UU*Condition + AU.AU*Condition + UA.UA*Condition + CU.AG*Condition + CA.UG*Condition + GU.AC*Condition + GA.UC*Condition + CG.CG*Condition + GG.CC + GC.GC*Condition + Term.AU*Condition - 1,
           df,
           weights = 1/df$dG.error)

A = anova(lm.group, lm.NN)
l.df.stats[[10]] = data.frame("Model.1" = A$Res.Df[1],
                             "Model.2" = A$Res.Df[2],
                             "F.value" = A$F[2],
                             "p-value" = A$`Pr(>F)`[2])
l.df.stats[[10]]$Parameter = "Grouped model to no group for GG.CC"

#GC.GC 11

lm.NN = lm(dG ~ AA.UU*Condition + AU.AU*Condition + UA.UA*Condition + CU.AG*Condition + CA.UG*Condition + GU.AC*Condition + GA.UC*Condition + CG.CG*Condition + GG.CC*Condition + GC.GC + Term.AU*Condition - 1,
           df,
           weights = 1/df$dG.error)

A = anova(lm.group, lm.NN)
l.df.stats[[11]] = data.frame("Model.1" = A$Res.Df[1],
                              "Model.2" = A$Res.Df[2],
                              "F.value" = A$F[2],
                              "p-value" = A$`Pr(>F)`[2])
l.df.stats[[11]]$Parameter = "Grouped model to no group for GC.GC"

#Term.AU 12

lm.NN = lm(dG ~ AA.UU*Condition + AU.AU*Condition + UA.UA*Condition + CU.AG*Condition + CA.UG*Condition + GU.AC*Condition + GA.UC*Condition + CG.CG*Condition + GG.CC*Condition + GC.GC*Condition + Term.AU - 1,
           df,
           weights = 1/df$dG.error)

A = anova(lm.group, lm.NN)
l.df.stats[[12]] = data.frame("Model.1" = A$Res.Df[1],
                              "Model.2" = A$Res.Df[2],
                              "F.value" = A$F[2],
                              "p-value" = A$`Pr(>F)`[2])
l.df.stats[[12]]$Parameter = "Grouped model to no group for Term.AU"

#AU.AU CU.AG GU.AC not grouped

lm.NN = lm(dG ~ AA.UU + AU.AU*Condition + UA.UA + CU.AG*Condition + CA.UG + GU.AC*Condition + GA.UC + CG.CG + GG.CC + GC.GC + Term.AU - 1,
           df,
           weights = 1/df$dG.error)

A = anova(lm.group, lm.NN)
l.df.stats[[13]] = data.frame("Model.1" = A$Res.Df[1],
                              "Model.2" = A$Res.Df[2],
                              "F.value" = A$F[2],
                              "p-value" = A$`Pr(>F)`[2])
l.df.stats[[13]]$Parameter = "Grouped model to AU.AU CU.AG GU.AC only grouped"


#AU.AU CU.AG GU.AConly grouped

anova(lm.group, lm.NN)

lm.NN = lm(dG ~ AA.UU + AU.AU*Condition + UA.UA + CU.AG*Condition + CA.UG + GU.AC*Condition + GA.UC + CG.CG + GG.CC + GC.GC + Term.AU - 1,
           df,
           weights = 1/df$dG.error)

A = anova(lm.No.group, lm.NN)
l.df.stats[[14]] = data.frame("Model.1" = A$Res.Df[1],
                              "Model.2" = A$Res.Df[2],
                              "F.value" = A$F[2],
                              "p-value" = A$`Pr(>F)`[2])
l.df.stats[[14]]$Parameter = "Ungrouped model to AU.AU CU.AG GU.AC only grouped"

#AU.AU CU.AG GU.AC not grouped

df.sig = bind_rows(l.df.stats)

write.csv(df.sig,
          "Tables/SI_Table_2_pvalues_and_established_NN_comparison/Sig_for_parameters.csv",
          row.names = FALSE)

####Write SI table 2####

list.files("Tables/SI_Table_2_pvalues_and_established_NN_comparison")

write.csv(df.result,
          "Tables/SI_Table_2_pvalues_and_established_NN_comparison/NN_parameters.csv",
          row.names = FALSE)

