library(tidyverse)
library(svd)
library(MeltR)
library(ggpubr)
library(ggrepel)

####Read in known parameters####

list.files("Tables/Table_2_NN_parameters")

df.Xia = read.csv("Tables/Table_2_NN_parameters/Xia_NN_params.csv")

####Load in helix energy data####

list.files("Tables/SI_Table_2_established_NN_comparison/")

df = read.csv("Tables/SI_Table_2_established_NN_comparison/SI_Table_2_NN_matrix_and_data_1MNaCl.csv")

sd = mean(df$dG.percent.error)/100

####Change labels####

df$Helix = factor(df$Helix,
                  levels = levels(factor(df$Helix)),
                  labels = 1:length(levels(factor(df$Helix))))

####filter data and subtract out fluorophore quencher contributions####

df$dG[which(df$FAMC.GBHQ1 == 1)] = df$dG[which(df$FAMC.GBHQ1 == 1)]  + 3.67

df$dG[which(df$FAMC.GBHQ1 != 1)] = df$dG[which(df$FAMC.GBHQ1 != 1)]  + 4.38

####Determine 1 M NN using linear regression####

colnames(df)

fit = lm(dG ~ AA.UU + AU.AU + UA.UA + CU.AG + CA.UG + GU.AC + GA.UC +CG.CG + GG.CC + GC.GC + Term.AU,
         df,
         weights = 1/df$SE.G)

df.lm = data.frame(coef(summary(fit)))

df.lm$Known = df.Xia$Xia

####Determine 1 M NN using SVD####

S = as.matrix(df %>% select(Initiation, AA.UU, AU.AU, UA.UA,
                         CU.AG,CA.UG, GU.AC, GA.UC, CG.CG,
                         GG.CC, GC.GC, Term.AU))
G = df$dG
s = df$SE.G
Gs = G/s
Ss = S/s

Ss.svd = La.svd(Ss)

U = Ss.svd$u
w = diag(Ss.svd$d)
V.t = Ss.svd$vt

w.under.1 = (1/w)
w.under.1[which(w.under.1 == Inf)] = 0
G.NN = t(V.t) %*% w.under.1 %*% t(U) %*% Gs

sigma.G.NN = c()

V = t(V.t)

nrow(V)

for (j in 1:nrow(V)){
  sigma.G.NN[j] = sqrt(sum((V[j,]/Ss.svd$d[j])^2))
}

df.svd = data.frame(G.NN)

colnames(df.svd) = "Estimate"

df.svd$Std..Error = sigma.G.NN

df.svd$Parameter = colnames(S)

df.svd$Known = df.Xia$Xia.

####Model in data####

df$dG.known = df$dG.known + 3.933333333

df$Modeled.error = abs(df$dG.known)*rnorm(nrow(df), 0, sd)

df$Modeled = df$dG.known + df$Modeled.error

df$Modeled.error = abs(0.05*df$Modeled)

####Determine 1 M NN using linear regression on modeled data####

colnames(df)

fit = lm(Modeled ~ AA.UU + AU.AU + UA.UA + CU.AG + CA.UG + GU.AC + GA.UC +CG.CG + GG.CC + GC.GC + Term.AU,
         df,
         weights = 1/df$Modeled.error)

df.lm.model = data.frame(coef(summary(fit)))

df.lm.model$Known = df.Xia$Xia.

####Determine 1 M NN using SVD####

S = as.matrix(df %>% select(Initiation, AA.UU, AU.AU, UA.UA,
                            CU.AG,CA.UG, GU.AC, GA.UC, CG.CG,
                            GG.CC, GC.GC, Term.AU))
G = df$Modeled
s = df$Modeled.error
Gs = G/s
Ss = S/s

Ss.svd = La.svd(Ss)

U = Ss.svd$u
w = diag(Ss.svd$d)
V.t = Ss.svd$vt

w.under.1 = (1/w)
w.under.1[which(w.under.1 == Inf)] = 0
G.NN = t(V.t) %*% w.under.1 %*% t(U) %*% Gs

sigma.G.NN = c()

V = t(V.t)

nrow(V)

for (j in 1:nrow(V)){
  sigma.G.NN[j] = sqrt(sum((V[j,]/Ss.svd$d[j])^2))
}

df.svd.modeled = data.frame(G.NN)

colnames(df.svd.modeled) = "Estimate"

df.svd.modeled$Std..Error = sigma.G.NN

df.svd.modeled$Parameter = colnames(S)

df.svd.modeled$Known = df.Xia$Xia.

####Consolidate data####

head(df.lm)
head(df.lm.model)
head(df.svd)
head(df.svd.modeled)

df.lm$Parameter = df.svd$Parameter
df.lm.model$Parameter = df.svd$Parameter

df.lm$Data = "Fluorescence"
df.lm.model$Data = "Modeled"
df.svd$Data = "Fluorescence"
df.svd.modeled$Data = "Modeled"

df.lm$Method = "lm"
df.lm.model$Method = "lm"
df.svd$Method = "svd"
df.svd.modeled$Method = "svd"

df.result = bind_rows(df.lm %>% select(Parameter, Known, Estimate, Std..Error, Data, Method),
          df.lm.model %>% select(Parameter, Known, Estimate, Std..Error, Data, Method),
          df.svd %>% select(Parameter, Known, Estimate, Std..Error, Data, Method),
          df.svd.modeled %>% select(Parameter, Known, Estimate, Std..Error, Data, Method))

list.files("Tables/Table_2_NN_parameters")

write.csv(df.result,
          "Tables/Table_2_NN_parameters/NN_parameters.csv",
          row.names = FALSE)

