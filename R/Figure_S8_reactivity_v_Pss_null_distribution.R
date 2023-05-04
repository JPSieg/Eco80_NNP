library(tidyverse)
library(viridis)
library(ggpubr)
library(cowplot)

#####Make a plot####

list.files("Figures")

df = read.csv("Figures/Figure_4_transcriptome_wide_trend/009_Reactivity_Pss.csv")

df$Program = factor(df$Program,
                    levels = c("UV_1MNaCl_data_tables",
                               "FDBI_1MNaCl_data_tables",
                               "1kcalmol_penalty_data_tables",
                               "FBDI_Eco80_data_tables"),
                    labels = c("UV 1MNaCl",
                               "FDBI 1MNaCl",
                               "1.0 kcal/mol penalty",
                               "FBDI Eco80"))


score = function(a = 1, df.i = df.M){
  Threshold = 7*df.i$Pss + a
  df.i$Wrong = 0
  df.i$Wrong[which(df.i$Mean > Threshold)] = 1
  df.c = df.i %>%
    group_by(Program) %>%
    summarise(Count = sum(Wrong))
  df.c$a = a
  output = df.c
}


####Model data based on randomly assigning Pss###

Model = function(P = "UV 1MNaCl"){
  df.M = df %>% filter(Program == P)
  df.M$Pss = sample(df.M$Pss)
  df.M$Random = T
  a = seq(0.25, 4, length.out = 50)
  l.df.c = {}
  for(i in 1:length(a)){
    l.df.c[[i]] = score(a[i], df.M)
  }
  df.c = bind_rows(l.df.c)
}

Model.null = function(Program = "UV 1MNaCl"){
  print("")
  print(Program)
  l.df.M = {}
  pb = txtProgressBar(min = 0,      # Minimum value of the progress bar
                      max = 1000)
  for (i in 1:1000){
    setTxtProgressBar(pb, i)
    l.df.M[[i]] = Model(Program)
  }
  
  df.M = bind_rows(l.df.M) 
  a = unique(df.M$a)
  
  Count = c()
  lower.Count = c()
  upper.Count = c()
  for (i in 1:length(a)){
    df.i = df.M %>% filter(a == a[i])
    Count[i] = mean(df.i$Count, na.rm = T)
    lower.Count[i] = quantile(df.i$Count, 0.025, na.rm = T)
    upper.Count[i] = quantile(df.i$Count, 0.975, na.rm = T)
  }
  df.M = data.frame(a,
                    Count,
                    lower.Count,
                    upper.Count)
  df.M$Program = Program
  output = df.M
}

l.df.m = lapply(unique(df$Program), Model.null)

df.m = bind_rows(l.df.m)

####Calculate for real data####

a = seq(0.25, 4, length.out = 50)
l.df.c = {}
for (i in 1:length(a)){
  l.df.c[[i]] = score(a[i], df)
}
df.c = bind_rows(l.df.c)

df.m$Data = "Raw data"
df.c$Data = "Scrambled data"


P = ggplot() +
  facet_wrap(~Program) +
  geom_point(data = df.c,
             mapping = aes(x = a, y = Count, color = Program)) +
  geom_pointrange(data = df.m,
             mapping = aes(x = a, y = Count, ymin = lower.Count, ymax = upper.Count), color = "red", size = 0.1) +
  geom_line(data = df.c,
            mapping = aes(x = a, y = Count, color = Program)) +
  geom_line(data = df.m,
                  mapping = aes(x = a, y = Count), color = "red") +
  theme_classic() +
  scale_color_manual(values = viridis(7, option = "A")) +
  scale_y_continuous(trans = "log10") +
  theme(axis.text = element_text(color = "black"),
        legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_blank()) +
  xlab("Reactivity threshold to be considered wrong (a)") +
  ylab("Number of wrongly predicted\ndouble stranded nucleotides")

list.files("Figures/SI_Figure_X_reactivity_v_Pss_null_distribution")

ggsave("Figures/SI_Figure_X_reactivity_v_Pss_null_distribution/SI_Figure_X_reactivity_v_Pss_null_distribution.svg", P, scale = 1.5)
