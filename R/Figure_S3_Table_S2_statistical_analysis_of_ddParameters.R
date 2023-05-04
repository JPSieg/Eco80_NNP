library(tidyverse)
library(cowplot)

list.files("Tables/SI_Table_1_fitting_statistics")

df = read.csv("Tables/SI_Table_1_fitting_statistics/SI_Table_1_data.csv")
  
####dH####

head(df)

df.dH = data.frame("Helix" = df %>%
                     filter(Method == "2 Global fit") %>%
                     filter(Condition == "1 M NaCl") %>%
                     arrange(Helix) %>%
                     select(Helix),
                   "dH.1MNaCl" = df %>%
                     filter(Method == "2 Global fit") %>%
                     filter(Condition == "1 M NaCl") %>%
                     arrange(Helix) %>%
                     select(H),
                   "SE.dH.1MNaCl" = df %>%
                     filter(Method == "2 Global fit") %>%
                     filter(Condition == "1 M NaCl") %>%
                     arrange(Helix) %>%
                     select(SE.H),
                   "dH.Eco80" = df %>%
                     filter(Method == "2 Global fit") %>%
                     filter(Condition == "Eco80") %>%
                     arrange(Helix) %>%
                     select(H),
                   "SE.dH.Eco80" = df %>%
                     filter(Method == "2 Global fit") %>%
                     filter(Condition == "Eco80") %>%
                     arrange(Helix) %>%
                     select(SE.H))

df.dH$ddH = df.dH$H.1 - df.dH$H
df.dH$SE.ddH = sqrt(df.dH$SE.H^2 + df.dH$SE.H.1^2)

Top.dH = ggplot(df.dH, aes(x = H, y = ddH)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  theme_classic() +
  scale_x_continuous(limits = c(-90, -0)) +
  theme(axis.text = element_text(color = "black"),
        axis.title.x = element_blank()) +
  ylab("\u0394\u0394H\u00B0")

Bottom.dH = ggplot(df.dH, aes(x = H, xmin = H - SE.H, xmax = H + SE.H,
                  y = H.1, ymin = H.1 - SE.H.1, ymax = H.1 + SE.H.1)) +
  geom_abline(slope = 1, intercept = 0) +
  coord_fixed(xlim = c(-90, -0), ylim = c(-90, -0)) +
  geom_pointrange() +
  geom_errorbarh() +
  theme_classic() +
  theme(axis.text = element_text(color = "black")) +
  xlab("1 M NaCl \u0394H\u00B0 (kcal/mol)") +
  ylab("Eco80 \u0394H\u00B0 (kcal/mol)")

PA.dH = plot_grid(Top.dH, Bottom.dH, ncol = 1, rel_heights = c(1, 3.5), align = "v")

dd.obs = abs(mean(df.dH$ddH))

model.dH = function(error){
  x = runif(n=25, min=-80, max=-20)
  x1 = x + rnorm(25, 0, error)
  x2 = x + rnorm(25, 0, error)
  ddG.test = abs(mean(x1 -  x2))
  output = data.frame(error, ddG.test)
}

df.e = df %>% 
  filter(Method == "2 Global fit") %>%
  select(SE.H)
SE = mean(df.e$SE.H)
SE.2x = 2*mean(df.e$SE.H)


l.df.SE = lapply(rep(SE, 1000000), model.dH)
l.df.SE.2x = lapply(rep(SE.2x, 1000000), model.dH)

df.test = bind_rows(l.df.SE, l.df.SE.2x)

head(df.test)

df.test$error = as.character(df.test$error)
df.test$error = factor(df.test$error,
                 levels = c("6.9286524757472", "13.8573049514944"),
                 labels = c("SE", "2xSE"))

dH.test = function(Error){
  df.SE = df.test %>% filter(error == Error)
  Parameter = "dH"
  SE = Error
  Mean.exp = mean(df.SE$ddG.test)
  Q95 = quantile(df.SE$ddG.test, 0.95)
  dd.obs = dd.obs
  p = nrow(df.SE %>% filter(ddG.test >= dd.obs))/nrow(df.SE)
  output = data.frame(Parameter, SE, Mean.exp, Q95, dd.obs, p)
}

l.df.table = lapply(levels(df.test$error), dH.test)

df.table.dH = bind_rows(l.df.table)

P.D.dH.hist = ggplot(df.test, aes(x = ddG.test, fill = error)) +
  geom_histogram(data = df.test %>% filter(error == "SE"), mapping = aes(x = ddG.test, fill = error)) +
  geom_vline(xintercept = dd.obs) + 
  theme_classic() +
  scale_fill_manual(values = c("grey")) +
  theme(axis.text = element_text(color = "black"),
        legend.position = "none") +
  xlab("Modeled mean |\u0394\u0394H\u00B0| (kcal/mol)")

####dS####

head(df)

df.dS = data.frame("Helix" = df %>%
                     filter(Method == "2 Global fit") %>%
                     filter(Condition == "1 M NaCl") %>%
                     arrange(Helix) %>%
                     select(Helix),
                   "dS.1MNaCl" = df %>%
                     filter(Method == "2 Global fit") %>%
                     filter(Condition == "1 M NaCl") %>%
                     arrange(Helix) %>%
                     select(S),
                   "SE.dS.1MNaCl" = df %>%
                     filter(Method == "2 Global fit") %>%
                     filter(Condition == "1 M NaCl") %>%
                     arrange(Helix) %>%
                     select(SE.S),
                   "dS.Eco80" = df %>%
                     filter(Method == "2 Global fit") %>%
                     filter(Condition == "Eco80") %>%
                     arrange(Helix) %>%
                     select(S),
                   "SE.dS.Eco80" = df %>%
                     filter(Method == "2 Global fit") %>%
                     filter(Condition == "Eco80") %>%
                     arrange(Helix) %>%
                     select(SE.S))

df.dS$ddS = df.dS$S.1 - df.dS$S
df.dS$SE.ddS = sqrt(df.dS$SE.S^2 + df.dS$SE.S.1^2)

Top.dS = ggplot(df.dS, aes(x = S, y = ddS)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  theme_classic() +
  scale_x_continuous(limits = c(-220, 0)) + 
  theme(axis.text = element_text(color = "black"),
        axis.title.x = element_blank()) +
  ylab("\u0394\u0394S\u00B0")

Bottom.dS = ggplot(df.dS, aes(x = S, xmin = S - SE.S, xmax = S + SE.S,
                              y = S.1, ymin = S.1 - SE.S.1, ymax = S.1 + SE.S.1)) +
  geom_abline(slope = 1, intercept = 0) +
  coord_fixed(xlim = c(-220, 0), ylim = c(-220, 0)) +
  geom_pointrange() +
  geom_errorbarh() +
  theme_classic() +
  theme(axis.text = element_text(color = "black")) +
  xlab("1 M NaCl \u0394S\u00B0 (cal/mol/K)") +
  ylab("Eco80 \u0394S\u00B0 (cal/mol/K)")

PA.dS = plot_grid(Top.dS, Bottom.dS, ncol = 1, rel_heights = c(1, 3.5), align = "v")

dd.obs = abs(mean(df.dS$ddS))

model.dS = function(error){
  x = runif(n=25, min=-220, max=-20)
  x1 = x + rnorm(25, 0, error)
  x2 = x + rnorm(25, 0, error)
  ddG.test = abs(mean(x1 -  x2))
  output = data.frame(error, ddG.test)
}

df.e = df %>% 
  filter(Method == "2 Global fit") %>%
  select(SE.S)
SE = mean(df.e$SE.S)
SE.2x = 2*mean(df.e$SE.S)


l.df.SE = lapply(rep(SE, 1000000), model.dS)
l.df.SE.2x = lapply(rep(SE.2x, 1000000), model.dS)

df.test = bind_rows(l.df.SE, l.df.SE.2x)

head(df.test)

df.test$error = as.character(df.test$error)
df.test$error = factor(df.test$error,
                       levels = c("21.7299032493372", "43.4598064986745"),
                       labels = c("SE", "2xSE"))

dS.test = function(Error){
  df.SE = df.test %>% filter(error == Error)
  Parameter = "dS"
  SE = Error
  Mean.exp = mean(df.SE$ddG.test)
  Q95 = quantile(df.SE$ddG.test, 0.95)
  dd.obs = dd.obs
  p = nrow(df.SE %>% filter(ddG.test >= dd.obs))/nrow(df.SE)
  output = data.frame(Parameter, SE, Mean.exp, Q95, dd.obs, p)
}

l.df.table = lapply(levels(df.test$error), dS.test)

df.table.dS = bind_rows(l.df.table)

P.D.dS.hist = ggplot(df.test, aes(x = ddG.test, fill = error)) +
  geom_histogram(data = df.test %>% filter(error == "SE"), mapping = aes(x = ddG.test, fill = error)) +
  geom_vline(xintercept = dd.obs) + 
  theme_classic() +
  scale_fill_manual(values = c("grey")) +
  theme(axis.text = element_text(color = "black"),
        legend.position = "none") +
  xlab("Modeled mean |\u0394\u0394S\u00B0| (kcal/mol)") 

####dG####

head(df)

df.dG = data.frame("Helix" = df %>%
                     filter(Method == "2 Global fit") %>%
                     filter(Condition == "1 M NaCl") %>%
                     arrange(Helix) %>%
                     select(Helix),
                   "dG.1MNaCl" = df %>%
                     filter(Method == "2 Global fit") %>%
                     filter(Condition == "1 M NaCl") %>%
                     arrange(Helix) %>%
                     select(G),
                   "SE.dG.1MNaCl" = df %>%
                     filter(Method == "2 Global fit") %>%
                     filter(Condition == "1 M NaCl") %>%
                     arrange(Helix) %>%
                     select(SE.G),
                   "dG.Eco80" = df %>%
                     filter(Method == "2 Global fit") %>%
                     filter(Condition == "Eco80") %>%
                     arrange(Helix) %>%
                     select(G),
                   "SE.dG.Eco80" = df %>%
                     filter(Method == "2 Global fit") %>%
                     filter(Condition == "Eco80") %>%
                     arrange(Helix) %>%
                     select(SE.G))

df.dG$ddG = df.dG$G.1 - df.dG$G
df.dG$SE.ddG = sqrt(df.dG$SE.G^2 + df.dG$SE.G.1^2)

Top.dG = ggplot(df.dG, aes(x = G, y = ddG)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        axis.title.x = element_blank()) +
  scale_x_continuous(limits = c(-16, -8)) +
  ylab("\u0394\u0394G\u00B037")

Bottom.dG = ggplot(df.dG, aes(x = G, xmin = G - SE.G, xmax = G + SE.G,
                              y = G.1, ymin = G.1 - SE.G.1, ymax = G.1 + SE.G.1)) +
  geom_abline(slope = 1, intercept = 0) +
  coord_fixed(xlim = c(-16, -8), ylim = c(-16, -8)) +
  geom_pointrange() +
  geom_errorbarh() +
  theme_classic() +
  theme(axis.text = element_text(color = "black")) +
  xlab("1 M NaCl \u0394G\u00B037 (kcal/mol)") +
  ylab("Eco80 \u0394G\u00B037 (kcal/mol)")

PA.dG = plot_grid(Top.dG, Bottom.dG, ncol = 1, rel_heights = c(1, 3.5), align = "v")

dd.obs = abs(mean(df.dG$ddG))

model.dG = function(error){
  x = runif(n=25, min=-16, max=-10)
  x1 = x + rnorm(25, 0, error)
  x2 = x + rnorm(25, 0, error)
  ddG.test = abs(mean(x1 -  x2))
  output = data.frame(error, ddG.test)
}

df.e = df %>% 
  filter(Method == "2 Global fit") %>%
  select(SE.G)
SE = mean(df.e$SE.G)
SE.2x = 2*mean(df.e$SE.G)


l.df.SE = lapply(rep(SE, 1000000), model.dG)
l.df.SE.2x = lapply(rep(SE.2x, 1000000), model.dG)

df.test = bind_rows(l.df.SE, l.df.SE.2x)

head(df.test)
unique(df.test$error)
df.test$error = as.character(df.test$error)
df.test$error = factor(df.test$error,
                       levels = c("0.209380932573337", "0.418761865146675"),
                       labels = c("SE", "2xSE"))

dG.test = function(Error){
  df.SE = df.test %>% filter(error == Error)
  Parameter = "dG"
  SE = Error
  Mean.exp = mean(df.SE$ddG.test)
  Q95 = quantile(df.SE$ddG.test, 0.95)
  dd.obs = dd.obs
  p = nrow(df.SE %>% filter(ddG.test >= dd.obs))/nrow(df.SE)
  output = data.frame(Parameter, SE, Mean.exp, Q95, dd.obs, p)
}

l.df.table = lapply(levels(df.test$error), dG.test)

df.table.dG = bind_rows(l.df.table)

P.D.dG.hist = ggplot() +
  geom_histogram(data = df.test %>% filter(error == "SE"), mapping = aes(x = ddG.test, fill = error)) +
  geom_vline(xintercept = dd.obs) + 
  theme_classic() +
  scale_fill_manual(values = c("grey")) +
  theme(axis.text = element_text(color = "black"),
        legend.position = "none") +
  xlab("Modeled mean |\u0394\u0394G\u00B037| (kcal/mol)")

####Compile plots####

ABC = plot_grid(PA.dH, PA.dS, PA.dG, nrow = 1, labels = c("A", "B", "C"))
DEF = plot_grid(P.D.dH.hist, P.D.dS.hist, P.D.dG.hist, nrow = 1, labels = c("D", "E", "F"))
P = plot_grid(ABC, DEF, ncol = 1, rel_heights = c(1.3, 1))

ggsave("SI_Figure_X_statistical_analysis_of_parameters.svg", P, width = 4, height = 3.2, units = "in", scale = 2.5, bg = "white")
df.final = bind_rows(df.table.dH, df.table.dS, df.table.dG)
write.csv(df.final, "SI_table_X_statistical_significance_of_ddP.csv", row.names = F)
 