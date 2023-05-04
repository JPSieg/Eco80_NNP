library(tidyverse)
library(MeltR)
library(ggrepel)

####Read in index####

list.files("Tables/SI_Table_1_fitting_statistics")

df.index = read.csv("Tables/SI_Table_1_fitting_statistics/Data_index.csv")

####Make a list of the fits####

l.fits = {}

####Helix A 1 M NaCl####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/01_js4022_Helix_A.csv")

ggplot(df %>% filter(Reading == 60),
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[1]] = meltR.F(df,
                      Kd_range = c(20, 300),
                      Kd_error_quantile = 0.8,
                      Save_results = "all",
                      file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                      file_prefix = "Helix_A_1M_NaCl")

####Helix B 1 M NaCl####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/02_js4022_Helix_B.csv")

ggplot(df %>% filter(Reading == 60),
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[2]] = meltR.F(df,
                      Save_results = "all",
                      Kd_range = c(20, 300),
                      Kd_error_quantile = 1,
                      file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                      file_prefix = "Helix_B_1M_NaCl")

####Helix C 1 M NaCl####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/03_js4022_Helix_C.csv")

ggplot(df %>% filter(Reading == 10),
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[3]] = meltR.F(df,
                      Kd_range = c(10, 150),
                      Kd_error_quantile = 0.6,
                      Save_results = "all",
                      file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                      file_prefix = "Helix_C_1M_NaCl")

####Helix D 1 M NaCl####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/04_js4022_Helix_D.csv")

ggplot(df %>% filter(Reading == 1),
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[4]] = meltR.F(df,
                      Kd_range = c(50, 500),
                      Kd_error_quantile = 0.8,
                      low_K = 0.01,
                      Save_results = "all",
                      file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                      file_prefix = "Helix_D_1M_NaCl")

####Helix E 1 M NaCl####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/05_js4022_Helix_E.csv")

ggplot(df %>% filter(Reading == 1),
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[5]] = meltR.F(df,
                      Kd_range = c(20, 150),
                      Kd_error_quantile = 0.8,
                      Save_results = "all",
                      file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                      file_prefix = "Helix_E_1M_NaCl")

####Helix F 1 M NaCl####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/06_js5003_Helix_F.csv")

ggplot(df %>% filter(Reading == 1),
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[6]] = meltR.F(df,
                      Kd_range = c(10, 100),
                      Kd_error_quantile = 0.8,
                      Save_results = "all",
                      file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                      file_prefix = "Helix_F_1M_NaCl")

####Helix G 1 M NaCl####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/07_js5003_Helix_G.csv")

ggplot(df %>% filter(Reading == 40),
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[7]] = meltR.F(df,
                      Save_results = "all",
                      Kd_range = c(1, 1000),
                      #low_K = 0.01,
                      file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                      file_prefix = "Helix_G_1M_NaCl")

####Helix H 1 M NaCl####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/08_js5003_Helix_H.csv")

ggplot(df %>% filter(Reading == 1) ,
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[8]] = meltR.F(df %>% filter(!Well %in% c("C11", "C12")),
                      Save_results = "all",
                      Kd_range = c(1, 1000),
                      #low_K = 0.01,
                      file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                      file_prefix = "Helix_H_1M_NaCl")

####Helix I 1 M NaCl####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/09_js5003_Helix_I.csv")

ggplot(df %>% filter(Reading == 1) ,
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[9]] = meltR.F(df,
                      Save_results = "all",
                      Kd_range = c(20, 260),
                      Kd_error_quantile = 0.8,
                      file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                      file_prefix = "Helix_I_1M_NaCl")

####Helix J 1 M NaCl####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/10_js5003_Helix_J.csv")

ggplot(df %>% filter(Reading == 1) ,
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[10]] = meltR.F(df,
                      Save_results = "all",
                      Kd_range = c(20, 260),
                      Kd_error_quantile = 1,
                      file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                      file_prefix = "Helix_J_1M_NaCl")

####Helix L 1 M NaCl####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/mh_1034_Helix_L.csv")

ggplot(df %>% filter(Reading == 1) ,
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[11]] = meltR.F(df,
                       Save_results = "all",
                       Kd_range = c(30, 500),
                       Kd_error_quantile = 1,
                       file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                       file_prefix = "Helix_L_1M_NaCl")

####Helix M 1 M NaCl####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/mh_1044_Helix_M.csv")

ggplot(df %>% filter(Reading == 1) ,
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[12]] = meltR.F(df,
                       Save_results = "all",
                       Kd_range = c(8, 150),
                       Kd_error_quantile = 0.7,
                       file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                       file_prefix = "Helix_M_1M_NaCl")

####Helix N 1 M NaCl####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/mh_1034_Helix_N.csv")

ggplot(df %>% filter(Reading == 1) ,
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[13]] = meltR.F(df,
                       Kd_range = c(30, 300),
                       Kd_error_quantile = 1,
                       Save_results = "all",
                       file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                       file_prefix = "Helix_N_1M_NaCl")

####Helix 0 1 M NaCl####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/mh_1034_Helix_O.csv")

ggplot(df %>% filter(Reading == 1) ,
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[14]] = meltR.F(df %>% filter(!Well %in% c("C10", "C11", "C12")),
                       Save_results = "all",
                       file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                       file_prefix = "Helix_O_1M_NaCl")

####Helix P 1 M NaCl####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/mh_1044_Helix_P.csv")

ggplot(df %>% filter(Reading == 1) ,
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[15]] = meltR.F(df,
                       Kd_range = c(10, 100),
                       Kd_error_quantile = 0.6,
                       Save_results = "all",
                       file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                       file_prefix = "Helix_P_1M_NaCl")

####Helix U 1 M NaCl####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/js6027_Helix_U.csv")

ggplot(df %>% filter(Reading == 1) ,
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[16]] = meltR.F(df,
                       Kd_range = c(20, 100),
                       Kd_error_quantile = 0.3,
                       Save_results = "all",
                       file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                       file_prefix = "Helix_U_1M_NaCl")

####Helix V 1 M NaCl####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/js6027_Helix_V.csv")

ggplot(df %>% filter(Reading == 1) ,
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[17]] = meltR.F(df,
                       Kd_range = c(10, 200),
                       Kd_error_quantile = 0.5,
                       Save_results = "all",
                       file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                       file_prefix = "Helix_V_1M_NaCl")

####Helix X 1 M NaCl####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/js6027_Helix_X.csv")

ggplot(df %>% filter(Reading == 1) ,
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[18]] = meltR.F(df,
                       Kd_range = c(10, 300),
                       Kd_error_quantile = 0.4,
                       Save_results = "all",
                       file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                       file_prefix = "Helix_X_1M_NaCl")

####Helix Y 1 M NaCl####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/mh1089_Helix_Y_1MNaCl.csv")

ggplot(df %>% filter(Reading == 1) ,
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[19]] = meltR.F(df,
                       Kd_range = c(15, 150),
                       Kd_error_quantile = 0.4,
                       Save_results = "all",
                       file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                       file_prefix = "Helix_Y_1M_NaCl")

####Helix Z 1 M NaCl####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/mh1089_Helix_Z_1MNaCl.csv")

ggplot(df %>% filter(Reading == 1) ,
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[20]] = meltR.F(df,
                       Kd_range = c(15, 250),
                       Kd_error_quantile = 1,
                       Save_results = "all",
                       file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                       file_prefix = "Helix_Z_1M_NaCl")

####Helix AA 1 M NaCl####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/mh1087_Helix_AA_1MNaCl.csv")

ggplot(df %>% filter(Reading == 1) ,
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[21]] = meltR.F(df,
                       Kd_range = c(50, 350),
                       Kd_error_quantile = 1,
                       Save_results = "all",
                       file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                       file_prefix = "Helix_AA_1M_NaCl")

####Helix CC 1 M NaCl####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/mh1071_Helix_CC_1MNaCl.csv")

ggplot(df %>% filter(Reading == 1) ,
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[22]] = meltR.F(df,
                       Kd_range = c(20, 200),
                       Kd_error_quantile = 1,
                       Save_results = "all",
                       file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                       file_prefix = "Helix_CC_1M_NaCl")

####Helix DD 1 M NaCl####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/mh1073_Helix_DD_1MNaCl.csv")

ggplot(df %>% filter(Reading == 1) ,
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[23]] = meltR.F(df,
                       Kd_range = c(20, 200),
                       Kd_error_quantile = 1,
                       Save_results = "all",
                       file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                       file_prefix = "Helix_DD_1M_NaCl")

####Helix EE 1 M NaCl####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/mh1073_Helix_EE_1MNaCl.csv")

ggplot(df %>% filter(Reading == 1) ,
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[24]] = meltR.F(df,
                       Kd_range = c(15, 120),
                       Kd_error_quantile = 1,
                       Save_results = "all",
                       file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                       file_prefix = "Helix_EE_1M_NaCl")

####Helix FF 1 M NaCl####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/mh1077_Helix_FF_1MNaCl.csv")

ggplot(df %>% filter(Reading == 1) ,
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[25]] = meltR.F(df,
                       Kd_range = c(20, 150),
                       Kd_error_quantile = 1,
                       Save_results = "all",
                       file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                       file_prefix = "Helix_FF_1M_NaCl")

####Helix A Eco80####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/mh_1063_Helix_A.csv")

ggplot(df %>% filter(Reading == 1) ,
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[26]] = meltR.F(df,
                       Kd_range = c(20, 300),
                       Kd_error_quantile = 0.65,
                       Save_results = "all",
                       file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                       file_prefix = "Helix_A_Eco80")

####Helix B Eco80####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/js7009_Helix_B_Eco80.csv")

ggplot(df %>% filter(Reading == 1) ,
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[27]] = meltR.F(df,
                       Kd_range = c(20, 300),
                       Kd_error_quantile = 1,
                       Save_results = "all",
                       file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                       file_prefix = "Helix_B_Eco80")

####Helix B Eco80####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/js7010_Helix_C_Eco80.csv")

ggplot(df %>% filter(Reading == 1) ,
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[28]] = meltR.F(df,
                       Kd_range = c(20, 400),
                       Kd_error_quantile = 1,
                       Save_results = "all",
                       file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                       file_prefix = "Helix_C_Eco80")

####Helix D Eco80####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/mh_1063_Helix_D.csv")

ggplot(df %>% filter(Reading == 1) ,
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[29]] = meltR.F(df,
                       Kd_range = c(50, 180),
                       Kd_error_quantile = 0.8,
                       Save_results = "all",
                       file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                       file_prefix = "Helix_D_Eco80")

####Helix E Eco80####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/mh_1057_Helix_E.csv")

ggplot(df %>% filter(Reading == 1) ,
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[30]] = meltR.F(df,
                       Kd_range = c(50, 150),
                       Kd_error_quantile = 0.8,
                       Save_results = "all",
                       file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                       file_prefix = "Helix_E_Eco80")

####Helix F Eco80####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/mh1018_Helix_F.csv")

ggplot(df %>% filter(Reading == 1) ,
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[31]] = meltR.F(df,
                       Kd_range = c(15, 100),
                       Kd_error_quantile = 0.9,
                       Save_results = "all",
                       file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                       file_prefix = "Helix_F_Eco80")

####Helix G Eco80####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/js6007_Helix_G.csv")

ggplot(df %>% filter(Reading == 1) ,
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[32]] = meltR.F(df,
                       Kd_range = c(50, 200),
                       Kd_error_quantile = 1,
                       Save_results = "all",
                       file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                       file_prefix = "Helix_G_Eco80")

####Helix H Eco80####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/mh1020_Helix_H.csv")

ggplot(df %>% filter(Reading == 1) ,
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[33]] = meltR.F(df,
                       Kd_range = c(30, 130),
                       Kd_error_quantile = 0.9,
                       Save_results = "all",
                       file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                       file_prefix = "Helix_H_Eco80")

####Helix I Eco80####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/js6007_Helix_I.csv")

ggplot(df %>% filter(Reading == 1) ,
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[34]] = meltR.F(df,
                       Kd_range = c(20, 450),
                       Kd_error_quantile = 1,
                       Save_results = "all",
                       file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                       file_prefix = "Helix_I_Eco80")

####Helix J Eco80####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/js6007_Helix_J.csv")

ggplot(df %>% filter(Reading == 1) ,
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[35]] = meltR.F(df,
                       Kd_range = c(30, 400),
                       Kd_error_quantile = 1,
                       Save_results = "all",
                       file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                       file_prefix = "Helix_J_Eco80")

####Helix L Eco80####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/mh_1060_Helix_L.csv")

ggplot(df %>% filter(Reading == 1) ,
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[36]] = meltR.F(df,
                       Kd_range = c(20, 200),
                       Kd_error_quantile = 1,
                       Save_results = "all",
                       file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                       file_prefix = "Helix_L_Eco80")

####Helix M Eco80####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/mh_1060_Helix_M.csv")

ggplot(df %>% filter(Reading == 1) ,
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[37]] = meltR.F(df,
                       Kd_range = c(10, 300),
                       Kd_error_quantile = 0.5,
                       Save_results = "all",
                       file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                       file_prefix = "Helix_M_Eco80")

####Helix N Eco80####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/js7007_Helix_N_Eco80.csv")

ggplot(df %>% filter(Reading == 1) ,
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[38]] = meltR.F(df,
                       Kd_range = c(20, 300),
                       Kd_error_quantile = 1,
                       Save_results = "all",
                       file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                       file_prefix = "Helix_N_Eco80")

####Helix O Eco80####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/js7006_Helix_O_Eco80.csv")

ggplot(df %>% filter(Reading == 1) ,
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[39]] = meltR.F(df,
                       Kd_range = c(20, 300),
                       Kd_error_quantile = 1,
                       Save_results = "all",
                       file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                       file_prefix = "Helix_O_Eco80")

####Helix P Eco80####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/js7006_Helix_P_Eco80.csv")

ggplot(df %>% filter(Reading == 1) ,
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[40]] = meltR.F(df,
                       Kd_range = c(10, 400),
                       Kd_error_quantile = 1,
                       Save_results = "all",
                       file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                       file_prefix = "Helix_P_Eco80")

####Helix U Eco80####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/js7005_Helix_U_Eco80.csv")

ggplot(df %>% filter(Reading == 1) ,
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[41]] = meltR.F(df,
                       Kd_range = c(50, 400),
                       Kd_error_quantile = 1,
                       Save_results = "all",
                       file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                       file_prefix = "Helix_U_Eco80")

####Helix V Eco80####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/js7005_Helix_V_Eco80.csv")

ggplot(df %>% filter(Reading == 1) ,
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[42]] = meltR.F(df,
                       Kd_range = c(10, 400),
                       Kd_error_quantile = 1,
                       Save_results = "all",
                       file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                       file_prefix = "Helix_V_Eco80")

####Helix X Eco80####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/js7003_Helix_X_Eco80.csv")

ggplot(df %>% filter(Reading == 1) ,
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[43]] = meltR.F(df,
                       Kd_range = c(10, 400),
                       Kd_error_quantile = 1,
                       Save_results = "all",
                       file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                       file_prefix = "Helix_X_Eco80")

####Helix Y Eco80####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/mh1093_Helix_Y_Eco80.csv")

ggplot(df %>% filter(Reading == 1) ,
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[44]] = meltR.F(df,
                       Kd_range = c(10, 300),
                       Kd_error_quantile = 1,
                       Save_results = "all",
                       file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                       file_prefix = "Helix_Y_Eco80")

####Helix Z Eco80####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/mh1093_Helix_Z_Eco80.csv")

ggplot(df %>% filter(Reading == 1) ,
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[45]] = meltR.F(df,
                       Kd_range = c(10, 400),
                       Kd_error_quantile = 1,
                       Save_results = "all",
                       file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                       file_prefix = "Helix_Z_Eco80")

####Helix AA Eco80####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/mh1083_Helix_AA_Eco80.csv")

ggplot(df %>% filter(Reading == 1) ,
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[46]] = meltR.F(df,
                       Kd_range = c(25, 400),
                       Kd_error_quantile = 1,
                       Save_results = "all",
                       file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                       file_prefix = "Helix_AA_Eco80")


####Helix CC Eco80####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/mh1081_Helix_CC_Eco80.csv")

ggplot(df %>% filter(Reading == 1) ,
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[47]] = meltR.F(df,
                       Kd_range = c(25, 400),
                       Kd_error_quantile = 1,
                       Save_results = "all",
                       file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                       file_prefix = "Helix_CC_Eco80")

####Helix DD Eco80####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/mh1083_Helix_DD_Eco80.csv")

ggplot(df %>% filter(Reading == 1) ,
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[48]] = meltR.F(df,
                       Kd_range = c(25, 400),
                       Kd_error_quantile = 1,
                       Save_results = "all",
                       file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                       file_prefix = "Helix_DD_Eco80")

####Helix EE Eco80####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/js7001_Helix_EE_Eco80.csv")

ggplot(df %>% filter(Reading == 1) ,
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[49]] = meltR.F(df,
                       Kd_range = c(25, 400),
                       Kd_error_quantile = 1,
                       Save_results = "all",
                       file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                       file_prefix = "Helix_EE_Eco80")

####Helix FF Eco80####

df.index

df = read.csv("Tables/SI_Table_1_fitting_statistics/Fluorescence_data/js7001_Helix_FF_Eco80.csv")

ggplot(df %>% filter(Reading == 1) ,
       aes(x = B, y = Emission, label = Well)) +
  geom_point() +
  geom_text_repel()

l.fits[[50]] = meltR.F(df,
                       Kd_range = c(25, 400),
                       Kd_error_quantile = 1,
                       Save_results = "all",
                       file_path = "Tables/SI_Table_1_fitting_statistics/Fit_data",
                       file_prefix = "Helix_FF_Eco80")



####Consolidate Results####

unique(df.index$file)

df.index = df.index %>% filter(file != "")

l.df.result = {}

for (i in 1:nrow(df.index)){
  experiment = i
  df = l.fits[[df.index$List.index[i]]]$VantHoff
  df$Helix =  df.index$Helix[experiment]
  df$Sequence.F = df.index$Sequence.F[experiment]
  df$Sequence.R = df.index$Sequence.R[experiment]
  df$Condition = df.index$Condition[experiment]
  l.df.result[[i]] = df
}

df.result = bind_rows(l.df.result)

#View(df.result)

####Calculate known dG####

dG.known = c()

for (i in 1:length(df.result$Helix)){
  dG.known[i] = Helix.energy(df.result$Sequence.F[i],
                             df.result$Sequence.R[i],
                             output = "energy",
                             F.Q = TRUE, Double.label = 0)
}

df.result$dG.known = dG.known

####Write.results####

colnames(df.result)

write.csv(df.result %>% select(Helix, Sequence.F, Sequence.R, Condition, Method,
                        H, SE.H, S, SE.S, G, SE.G,
                        K_error, R, Kd.opt, dG.known),
          "Tables/SI_Table_1_fitting_statistics/SI_Table_1_data.csv",
          row.names = FALSE)
