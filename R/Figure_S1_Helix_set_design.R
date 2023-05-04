library(tidyverse)
library(MeltR)
library(viridis)

####Count para,eters####

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

head(df.param)

df.param$Sequence.F = df$Sequence.F
df.param$Sequence.R = df$Sequence.R

####Tidy format NN Parameter matrix####

M = as.matrix(df.param %>% select(!Sequence.F) %>% select(!Sequence.R))

l.df.param = {}

for (i in 1:nrow(df.param)){
  Sequence.F = df.param$Sequence.F[i]
  Sequence.R = df.param$Sequence.R[i]
  Parameter = colnames(M)
  Count = as.vector(M[i,])
  l.df.param[[i]] = data.frame(Sequence.F, Sequence.R, Parameter, Count)
}

df = bind_rows(l.df.param)

####Make Plot####

strreverse <- function(x){ strsplit(x, NULL) %>% lapply(rev) %>% sapply(paste, collapse="") }

df$Helix = paste("5'-", df$Sequence.F, "-3'\n3'-", strreverse(df$Sequence.R), "-5'", sep = "")


df = df %>%
  filter(Count != 0) %>%
  filter(!Parameter %in% c("Initiation", "FAMC.GBHQ1", "FAMU.ABHQ1"))

unique(df$Parameter)

df$Parameter = factor(df$Parameter,
       levels = c("AA.UU", "AU.AU", "UA.UA", "CU.AG", "CA.UG", "GU.AC", "GA.UC", "CG.CG", "GG.CC", "GC.GC", "Term.AU"))

####PA####

PA = ggplot(df,
       aes(x = Parameter, y = Count, fill = Helix)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_manual(values = viridis(length(unique(df$Helix)))) +
  theme(axis.text.x = element_text(color = "black",
                                   angle = 45,
                                   hjust = 1,
                                   size = 8),
        axis.text = element_text(color = "black"),
        strip.text = element_text(color = "black"),
        legend.title = element_blank())

####Load in helix energy data####

df = read.csv("Tables/SI_Table_2_pvalues_and_established_NN_comparison/SI_Table_2_NN_matrix_and_data_1MNaCl.csv")
df.known = read.csv("Tables/SI_Table_2_pvalues_and_established_NN_comparison/PublishedNN_parameters.csv") %>% filter(Condition == "Turner 1989 1 M NaCl")

####Calculate dG####

S = as.matrix(df %>% select(Initiation, AA.UU, AU.AU, UA.UA,
                            CU.AG,CA.UG, GU.AC, GA.UC, CG.CG,
                            GG.CC, GC.GC, Term.AU))

G.NN = matrix(data = c(4.09, -0.9, -1.1, -1.3, -2.1, -2.1, -2.2, -2.1, -2.4, -3.3, -3.4, 0.45), nrow = 12, ncol = 1)
dG.known = S %*% G.NN
df$dG = dG.known + rnorm(length(dG.known), 0, mean(df$dG.error))
df$dG.error = abs(rnorm(length(dG.known), 0, mean(df$dG.error)))

####Calculate NNPs on modeled data####

fit = lm(dG ~ AA.UU + AU.AU + UA.UA + CU.AG + CA.UG + GU.AC + GA.UC +CG.CG + GG.CC + GC.GC + Term.AU,
         df,
         weights = 1/df$dG.error)
summary(fit)


df.1M = data.frame(coef(summary(fit)))

df.1M$Parameter = row.names(df.1M)
df.1M$Parameter[1] = "Initiation"
df.1M$Condition = "1 M NaCl"
df.1M$dG.known = df.known$Value

PB = ggplot(df.1M, aes(x = dG.known,
                        y = Estimate, ymin = Estimate - Std..Error, ymax = Estimate + Std..Error)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_pointrange() +
  annotate("text", x = -2.5, y = 1,
           label = paste("R2 =", round(cor(df.1M$dG.known, df.1M$Estimate), digits = 2)),
           size = 3) +
  coord_fixed(ylim = c(-4, 5), xlim = c(-4, 5)) +
  theme_classic() +
  xlab("NNP Turner 1989 (kcal/mol)") +
  ylab("NNP calculated with\nmodeled data  (kcal/mol)") +
  theme(axis.text = element_text(color = "black"))

####Make Final plot####

P = plot_grid(PA, PB, rel_widths = c(1.5, 1), labels = c("A", "B"))

ggsave("Figures/SI_Figure_1_NN_representation/SI_Figure_1_NN_representation.svg", 
       units = "in", width = 4.5, height = 2,bg = "white", scale = 2)

