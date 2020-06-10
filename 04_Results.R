#########################################################
## Checking results of bombus integrated models
## Author: Francesca
## Date created: 2020-05-12
## Date modified: 
########################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(sparta)
library(TrendSummaries)

## Assess convergence across all chains ----

## Run on JASMIN ##
# 
# source("./Functions/combine_chains_modified.R")
# 
# 
# combine_chains(data_dir = "Outputs/Occ_model",
#                output_dir = "Outputs/Occ_model",
#                target = 16670)
# 
# combine_chains(data_dir = "Outputs/Integrated",
#                output_dir = "Outputs/Integrated", model_type = "integrated",
#                target = 16670)

#####################

## Extract model results and plot ----

# load("./Data/BWARS_forJAGS.rdata")
# load("./Data/BeeWalks_forJAGS.rdata")
# 
# bwars <- read.csv("./Data/BWARS_cleaned.csv",
#                   header = TRUE, stringsAsFactors = FALSE)
# 
# beewalks <- read.csv("./Data/BeeWalks_cleaned.csv",
#                           header = TRUE, stringsAsFactors = FALSE)
# 
# 
# species_list <- read.csv("./Data/common_species_bombus.csv", 
#                          header = TRUE, stringsAsFactors = FALSE)
# 
# source("./Functions/extract_metrics.R")
# 
# 
# model_outputs <- lapply(species_list$species, extract_metrics, 
#                         int_input_path = "./Outputs/Integrated/", 
#                         occmod_input_path = "./Outputs/Occ_model/", it = 16670, ep = 100000,
#                         df_str = beewalks, df_unstr = bwars, 
#                         formatted_df_str = beewalks_formatted, 
#                         formatted_df_unstr = BWARS_formatted,
#                         species_col_str = "latin", species_col_unstr = "CONCEPT",
#                         gridref_col_str = "GRIDREF_1KM_PREC", 
#                         gridref_col_unstr = "GRIDREF_1KM_PREC")
# 
# 
# model_outputs_df <- do.call(rbind, model_outputs)
# 
# str(model_outputs_df)
# 
# write.csv(model_outputs_df, "./Outputs/bombus_model_metrics.csv", row.names = FALSE)
# 
# model_outputs_df <- read.csv("./Outputs/bombus_model_metrics.csv",
#                              header = TRUE, stringsAsFactors = FALSE)
# 
# str(model_outputs_df)

# summary stats from data

bwars <- read.csv("./Data/BWARS_cleaned.csv", header = TRUE, stringsAsFactors = FALSE)

beewalks <- read.csv("./Data/BeeWalks_cleaned.csv", header = TRUE, stringsAsFactors = FALSE)


dataDiagnostics(taxa = bwars$CONCEPT, site = bwars$GRIDREF_1KM_PREC, time_period = bwars$YEAR)

dataDiagnostics(taxa = beewalks$latin, site = beewalks$GRIDREF_1KM_PREC, time_period = beewalks$Year)

load("/data/input-data/BWARS_forJAGS.rdata")
load("/data/input-data/BeeWalks_forJAGS.rdata")

# how many records went into the models

dim(BWARS_formatted$spp_vis)

dim(beewalks_formatted$spp_vis)

# how many records of the species

table(BWARS_formatted$spp_vis$'bombus lapidarius')

table(beewalks_formatted$spp_vis$'bombus lapidarius')

# how many sites

length(unique(BWARS_formatted$occDetdata$site))

length(unique(beewalks_formatted$occDetdata$site))

## Look at bombus hypnorum

int_model <- read.csv("/data/outputs/Integrated/bombus hypnorum_it33400_ep2e+05.csv",
                      header = TRUE, row.names = 1, stringsAsFactors = FALSE)

int_model <- int_model %>%
  mutate(parameter = rownames(int_model),
         rhat_threshold = case_when(PSRF > 1.1 ~ "Bad",
                                    PSRF <= 1.1 ~ "Good"),
         model = "integrated",
         precision = 1/SD^2) %>%
  rename(quant_025 = "X2.5.", quant_25 = "X25.", quant_50 = "X50.", 
         quant_75 = "X75.", quant_97.5 = "X97.5.")

occ_model <- read.csv("/data/outputs/Occ_model/bombus hypnorumit33400_ep2e+05.csv",
                      header = TRUE, row.names = 1, stringsAsFactors = FALSE)

occ_model <- occ_model %>%
  mutate(parameter = rownames(occ_model),
         rhat_threshold = case_when(PSRF > 1.1 ~ "Bad",
                                    PSRF <= 1.1 ~ "Good"),
         model = "occupancy",
         precision = 1/SD^2) %>%
  rename(quant_025 = "X2.5.", quant_25 = "X25.", quant_50 = "X50.", 
         quant_75 = "X75.", quant_97.5 = "X97.5.")


int_occ_plot <- ggplot(data = int_model %>% slice(28:34), 
                       aes(x = 2010:2016, y = Mean)) +
  geom_point(size = 4, aes(color = rhat_threshold)) +
  geom_ribbon(aes(ymin = quant_025, ymax = quant_97.5), alpha = 0.2) +
  geom_line(size = 1, col = "blue") +
  geom_text(aes(label = round(precision, digits = 2)), nudge_y = 0.05) +
  scale_color_manual(name = 'Rhat', values = c('Bad' = 'red','Good' = 'blue')) +
  ylab("Occupancy") +
  xlab("Year") +
  scale_y_continuous(limits = c(0.5, 1)) +
  ggtitle("Integrated model for Bombus hypnorum") + 
  theme_bw() +
  theme(plot.title = element_text(lineheight = .8, face = "bold"),
        legend.position = 'bottom')

occ_plot <- ggplot(data = occ_model %>% slice(21:27), 
                   aes(x = 2010:2016, y = Mean)) +
  geom_point(size = 4, aes(color = rhat_threshold)) +
  geom_ribbon(aes(ymin = quant_025, ymax = quant_97.5), alpha = 0.2) +
  geom_text(aes(label = round(precision, digits = 2)), nudge_y = 0.05) +
  geom_line(size = 1, col = "blue") +
  scale_color_manual(name = 'Rhat', values = c('Bad' = 'red','Good' = 'blue')) +
  ylab("Occupancy") +
  xlab("Year") +
  scale_y_continuous(limits = c(0.5, 1)) +
  ggtitle("Occupancy model for Bombus hypnorum") + 
  theme_bw() +
  theme(plot.title = element_text(lineheight = .8, face = "bold"),
        legend.position = 'bottom')

int_occ_plot + occ_plot


int_a_plot <- ggplot(data = int_model %>% slice(1:7), 
                     aes(x = 2010:2016, y = Mean)) +
  geom_point(size = 4, aes(color = rhat_threshold)) +
  geom_ribbon(aes(ymin = quant_025, ymax = quant_97.5), alpha = 0.2) +
  geom_line(size = 1, col = "blue") +
  geom_text(aes(label = round(precision, digits = 2)), nudge_y = 0.2) +
  scale_color_manual(name = 'Rhat', values = c('Bad' = 'red','Good' = 'blue')) +
  ylab("a") +
  xlab("Year") +
  # scale_y_continuous(limits = c(0.5, 1)) +
  ggtitle("Integrated model for Bombus hypnorum") + 
  theme_bw() +
  theme(plot.title = element_text(lineheight = .8, face = "bold"),
        legend.position = 'bottom')

occ_a_plot <- ggplot(data = occ_model %>% slice(1:7), 
                     aes(x = 2010:2016, y = Mean)) +
  geom_point(size = 4, aes(color = rhat_threshold)) +
  geom_ribbon(aes(ymin = quant_025, ymax = quant_97.5), alpha = 0.2) +
  geom_line(size = 1, col = "blue") +
  geom_text(aes(label = round(precision, digits = 2)), nudge_y = 0.2) +
  scale_color_manual(name = 'Rhat', values = c('Bad' = 'red','Good' = 'blue')) +
  ylab("a") +
  xlab("Year") +
  # scale_y_continuous(limits = c(0.5, 1)) +
  ggtitle("Occupancy model for Bombus hypnorum") + 
  theme_bw() +
  theme(plot.title = element_text(lineheight = .8, face = "bold"),
        legend.position = 'bottom')

int_a_plot + occ_a_plot


models_output <- rbind(occ_model, int_model) %>%
  pivot_wider(names_from = model, 
              values_from = c(Mean, SD, Naive.SE, 
                              Time.series.SE, quant_025,
                              quant_25, quant_50, quant_75,
                              quant_97.5, PSRF, PSRF..97.5..quantile,
                              precision, rhat_threshold)) 


ggplot(data = models_output, aes(x = precision_occupancy, y = precision_integrated)) +
  geom_point() +
  geom_smooth(method = "lm", colour = "black") +
  geom_abline(intercept = 0, slope = 1,linetype = "dashed") +
  coord_equal() +
  ylab("Parameters precision (integrated model)") +
  xlab("Parameters precision (occupancy model)") +
  ggtitle("Bombus hypnorum") + 
  theme_bw()


## Look at bombus hortorum

int_model <- read.csv("/data/outputs/Integrated/bombus hortorumit33400_ep2e+05.csv",
                      header = TRUE, row.names = 1, stringsAsFactors = FALSE)

int_model <- int_model %>%
  mutate(parameter = rownames(int_model),
         rhat_threshold = case_when(PSRF > 1.1 ~ "Bad",
                                    PSRF <= 1.1 ~ "Good"),
         model = "integrated",
         precision = 1/SD^2) %>%
  rename(quant_025 = "X2.5.", quant_25 = "X25.", quant_50 = "X50.", 
         quant_75 = "X75.", quant_97.5 = "X97.5.")

occ_model <- read.csv("/data/outputs/Occ_model/bombus hortorumit33400_ep2e+05.csv",
                      header = TRUE, row.names = 1, stringsAsFactors = FALSE)

occ_model <- occ_model %>%
  mutate(parameter = rownames(occ_model),
         rhat_threshold = case_when(PSRF > 1.1 ~ "Bad",
                                    PSRF <= 1.1 ~ "Good"),
         model = "occupancy",
         precision = 1/SD^2) %>%
  rename(quant_025 = "X2.5.", quant_25 = "X25.", quant_50 = "X50.", 
         quant_75 = "X75.", quant_97.5 = "X97.5.")


int_occ_plot <- ggplot(data = int_model %>% slice(28:34), 
                       aes(x = 2010:2016, y = Mean)) +
  geom_point(size = 4, aes(color = rhat_threshold)) +
  geom_ribbon(aes(ymin = quant_025, ymax = quant_97.5), alpha = 0.2) +
  geom_line(size = 1, col = "blue") +
  geom_text(aes(label = round(precision, digits = 2)), nudge_y = 0.05) +
  scale_color_manual(name = 'Rhat', values = c('Bad' = 'red','Good' = 'blue')) +
  ylab("Occupancy") +
  xlab("Year") +
  scale_y_continuous(limits = c(0.5, 1)) +
  ggtitle("Integrated model for Bombus hortorum") + 
  theme_bw() +
  theme(plot.title = element_text(lineheight = .8, face = "bold"),
        legend.position = 'bottom')

occ_plot <- ggplot(data = occ_model %>% slice(21:27), 
                   aes(x = 2010:2016, y = Mean)) +
  geom_point(size = 4, aes(color = rhat_threshold)) +
  geom_ribbon(aes(ymin = quant_025, ymax = quant_97.5), alpha = 0.2) +
  geom_text(aes(label = round(precision, digits = 2)), nudge_y = 0.05) +
  geom_line(size = 1, col = "blue") +
  scale_color_manual(name = 'Rhat', values = c('Bad' = 'red','Good' = 'blue')) +
  ylab("Occupancy") +
  xlab("Year") +
  scale_y_continuous(limits = c(0.5, 1)) +
  ggtitle("Occupancy model for Bombus hortorum") + 
  theme_bw() +
  theme(plot.title = element_text(lineheight = .8, face = "bold"),
        legend.position = 'bottom')

int_occ_plot + occ_plot


int_a_plot <- ggplot(data = int_model %>% slice(1:7), 
                     aes(x = 2010:2016, y = Mean)) +
  geom_point(size = 4, aes(color = rhat_threshold)) +
  geom_ribbon(aes(ymin = quant_025, ymax = quant_97.5), alpha = 0.2) +
  geom_line(size = 1, col = "blue") +
  geom_text(aes(label = round(precision, digits = 2)), nudge_y = 0.2) +
  scale_color_manual(name = 'Rhat', values = c('Bad' = 'red','Good' = 'blue')) +
  ylab("a") +
  xlab("Year") +
  # scale_y_continuous(limits = c(0.5, 1)) +
  ggtitle("Integrated model for Bombus hortorum") + 
  theme_bw() +
  theme(plot.title = element_text(lineheight = .8, face = "bold"),
        legend.position = 'bottom')

occ_a_plot <- ggplot(data = occ_model %>% slice(1:7), 
                     aes(x = 2010:2016, y = Mean)) +
  geom_point(size = 4, aes(color = rhat_threshold)) +
  geom_ribbon(aes(ymin = quant_025, ymax = quant_97.5), alpha = 0.2) +
  geom_line(size = 1, col = "blue") +
  geom_text(aes(label = round(precision, digits = 2)), nudge_y = 0.2) +
  scale_color_manual(name = 'Rhat', values = c('Bad' = 'red','Good' = 'blue')) +
  ylab("a") +
  xlab("Year") +
  # scale_y_continuous(limits = c(0.5, 1)) +
  ggtitle("Occupancy model for Bombus hortorum") + 
  theme_bw() +
  theme(plot.title = element_text(lineheight = .8, face = "bold"),
        legend.position = 'bottom')

int_a_plot + occ_a_plot


models_output <- rbind(occ_model, int_model) %>%
  pivot_wider(names_from = model, 
              values_from = c(Mean, SD, Naive.SE, 
                              Time.series.SE, quant_025,
                              quant_25, quant_50, quant_75,
                              quant_97.5, PSRF, PSRF..97.5..quantile,
                              precision, rhat_threshold)) 


ggplot(data = models_output, aes(x = precision_occupancy, y = precision_integrated)) +
  geom_point() +
  geom_smooth(method = "lm", colour = "black") +
  geom_abline(intercept = 0, slope = 1,linetype = "dashed") +
  coord_equal() +
  ylab("Parameters precision (integrated model)") +
  xlab("Parameters precision (occupancy model)") +
  ggtitle("Bombus hortorum") +
  theme_bw()


## Look at bombus lapidarius

int_model <- read.csv("/data/outputs/Integrated/bombus lapidariusit33400_ep2e+05.csv",
                      header = TRUE, row.names = 1, stringsAsFactors = FALSE)

int_model <- int_model %>%
  mutate(parameter = rownames(int_model),
         rhat_threshold = case_when(PSRF > 1.1 ~ "Bad",
                                    PSRF <= 1.1 ~ "Good"),
         model = "integrated",
         precision = 1/SD^2) %>%
  rename(quant_025 = "X2.5.", quant_25 = "X25.", quant_50 = "X50.", 
         quant_75 = "X75.", quant_97.5 = "X97.5.")

occ_model <- read.csv("/data/outputs/Occ_model/bombus lapidariusit33400_ep2e+05.csv",
                      header = TRUE, row.names = 1, stringsAsFactors = FALSE)

occ_model <- occ_model %>%
  mutate(parameter = rownames(occ_model),
         rhat_threshold = case_when(PSRF > 1.1 ~ "Bad",
                                    PSRF <= 1.1 ~ "Good"),
         model = "occupancy",
         precision = 1/SD^2) %>%
  rename(quant_025 = "X2.5.", quant_25 = "X25.", quant_50 = "X50.", 
         quant_75 = "X75.", quant_97.5 = "X97.5.")


int_occ_plot <- ggplot(data = int_model %>% slice(28:34), 
                       aes(x = 2010:2016, y = Mean)) +
  geom_point(size = 4, aes(color = rhat_threshold)) +
  geom_ribbon(aes(ymin = quant_025, ymax = quant_97.5), alpha = 0.2) +
  geom_line(size = 1, col = "blue") +
  geom_text(aes(label = round(precision, digits = 2)), nudge_y = 0.05) +
  scale_color_manual(name = 'Rhat', values = c('Bad' = 'red','Good' = 'blue')) +
  ylab("Occupancy") +
  xlab("Year") +
  scale_y_continuous(limits = c(0.5, 1)) +
  ggtitle("Integrated model for Bombus lapidarius") + 
  theme_bw() +
  theme(plot.title = element_text(lineheight = .8, face = "bold"),
        legend.position = 'bottom')

occ_plot <- ggplot(data = occ_model %>% slice(21:27), 
                   aes(x = 2010:2016, y = Mean)) +
  geom_point(size = 4, aes(color = rhat_threshold)) +
  geom_ribbon(aes(ymin = quant_025, ymax = quant_97.5), alpha = 0.2) +
  geom_text(aes(label = round(precision, digits = 2)), nudge_y = 0.05) +
  geom_line(size = 1, col = "blue") +
  scale_color_manual(name = 'Rhat', values = c('Bad' = 'red','Good' = 'blue')) +
  ylab("Occupancy") +
  xlab("Year") +
  scale_y_continuous(limits = c(0.5, 1)) +
  ggtitle("Occupancy model for Bombus lapidarius") + 
  theme_bw() +
  theme(plot.title = element_text(lineheight = .8, face = "bold"),
        legend.position = 'bottom')

int_occ_plot + occ_plot


int_a_plot <- ggplot(data = int_model %>% slice(1:7), 
                     aes(x = 2010:2016, y = Mean)) +
  geom_point(size = 4, aes(color = rhat_threshold)) +
  geom_ribbon(aes(ymin = quant_025, ymax = quant_97.5), alpha = 0.2) +
  geom_line(size = 1, col = "blue") +
  geom_text(aes(label = round(precision, digits = 2)), nudge_y = 1) +
  scale_color_manual(name = 'Rhat', values = c('Bad' = 'red','Good' = 'blue')) +
  ylab("a") +
  xlab("Year") +
  scale_y_continuous(limits = c(0, 10)) +
  ggtitle("Integrated model for Bombus lapidarius") + 
  theme_bw() +
  theme(plot.title = element_text(lineheight = .8, face = "bold"),
        legend.position = 'bottom')

occ_a_plot <- ggplot(data = occ_model %>% slice(1:7), 
                     aes(x = 2010:2016, y = Mean)) +
  geom_point(size = 4, aes(color = rhat_threshold)) +
  geom_ribbon(aes(ymin = quant_025, ymax = quant_97.5), alpha = 0.2) +
  geom_line(size = 1, col = "blue") +
  geom_text(aes(label = round(precision, digits = 2)), nudge_y = 1) +
  scale_color_manual(name = 'Rhat', values = c('Bad' = 'red','Good' = 'blue')) +
  ylab("a") +
  xlab("Year") +
  scale_y_continuous(limits = c(0, 10)) +
  ggtitle("Occupancy model for Bombus lapidarius") + 
  theme_bw() +
  theme(plot.title = element_text(lineheight = .8, face = "bold"),
        legend.position = 'bottom')

int_a_plot + occ_a_plot


models_output <- rbind(occ_model, int_model) %>%
  pivot_wider(names_from = model, 
              values_from = c(Mean, SD, Naive.SE, 
                              Time.series.SE, quant_025,
                              quant_25, quant_50, quant_75,
                              quant_97.5, PSRF, PSRF..97.5..quantile,
                              precision, rhat_threshold)) 

models_output$parameter_grouped <- sub("\\[.*", "", models_output$parameter)


ggplot(data = models_output, aes(x = precision_occupancy, y = precision_integrated)) +
  geom_point(aes(color = parameter_grouped, fill = parameter_grouped),
             position=position_jitter(h=100, w=100),
             shape = 21, alpha = 0.5, size = 3) +
  scale_fill_viridis_d() +
  scale_colour_viridis_d() +
  geom_smooth(method = "lm", colour = "black") +
  geom_abline(intercept = 0, slope = 1,linetype = "dashed") +
  coord_equal() +
  ylab("Parameters precision (integrated model)") +
  xlab("Parameters precision (occupancy model)") +
  ggtitle("Bombus lapidarius") +
  theme_bw()



## Look at bombus sylvestris

int_model <- read.csv("/data/outputs/Integrated/bombus sylvestrisit33400_ep2e+05.csv",
                      header = TRUE, row.names = 1, stringsAsFactors = FALSE)

int_model <- int_model %>%
  mutate(parameter = rownames(int_model),
         rhat_threshold = case_when(PSRF > 1.1 ~ "Bad",
                                    PSRF <= 1.1 ~ "Good"),
         model = "integrated",
         precision = 1/SD^2) %>%
  rename(quant_025 = "X2.5.", quant_25 = "X25.", quant_50 = "X50.", 
         quant_75 = "X75.", quant_97.5 = "X97.5.")

occ_model <- read.csv("/data/outputs/Occ_model/bombus sylvestrisit33400_ep2e+05.csv",
                      header = TRUE, row.names = 1, stringsAsFactors = FALSE)

occ_model <- occ_model %>%
  mutate(parameter = rownames(occ_model),
         rhat_threshold = case_when(PSRF > 1.1 ~ "Bad",
                                    PSRF <= 1.1 ~ "Good"),
         model = "occupancy",
         precision = 1/SD^2) %>%
  rename(quant_025 = "X2.5.", quant_25 = "X25.", quant_50 = "X50.", 
         quant_75 = "X75.", quant_97.5 = "X97.5.")


int_occ_plot <- ggplot(data = int_model %>% slice(28:34), 
                       aes(x = 2010:2016, y = Mean)) +
  geom_point(size = 4, aes(color = rhat_threshold)) +
  geom_ribbon(aes(ymin = quant_025, ymax = quant_97.5), alpha = 0.2) +
  geom_line(size = 1, col = "blue") +
  geom_text(aes(label = round(precision, digits = 2)), nudge_y = 0.05) +
  scale_color_manual(name = 'Rhat', values = c('Bad' = 'red','Good' = 'blue')) +
  ylab("Occupancy") +
  xlab("Year") +
  scale_y_continuous(limits = c(0, 1)) +
  ggtitle("Integrated model for Bombus sylvestris") + 
  theme_bw() +
  theme(plot.title = element_text(lineheight = .8, face = "bold"),
        legend.position = 'bottom')

occ_plot <- ggplot(data = occ_model %>% slice(21:27), 
                   aes(x = 2010:2016, y = Mean)) +
  geom_point(size = 4, aes(color = rhat_threshold)) +
  geom_ribbon(aes(ymin = quant_025, ymax = quant_97.5), alpha = 0.2) +
  geom_text(aes(label = round(precision, digits = 2)), nudge_y = 0.05) +
  geom_line(size = 1, col = "blue") +
  scale_color_manual(name = 'Rhat', values = c('Bad' = 'red','Good' = 'blue')) +
  ylab("Occupancy") +
  xlab("Year") +
  scale_y_continuous(limits = c(0, 1)) +
  ggtitle("Occupancy model for Bombus sylvestris") + 
  theme_bw() +
  theme(plot.title = element_text(lineheight = .8, face = "bold"),
        legend.position = 'bottom')

int_occ_plot + occ_plot


int_a_plot <- ggplot(data = int_model %>% slice(1:7), 
                     aes(x = 2010:2016, y = Mean)) +
  geom_point(size = 4, aes(color = rhat_threshold)) +
  geom_ribbon(aes(ymin = quant_025, ymax = quant_97.5), alpha = 0.2) +
  geom_line(size = 1, col = "blue") +
  geom_text(aes(label = round(precision, digits = 2)), nudge_y = 1) +
  scale_color_manual(name = 'Rhat', values = c('Bad' = 'red','Good' = 'blue')) +
  ylab("a") +
  xlab("Year") +
  scale_y_continuous(limits = c(-3, 3)) +
  ggtitle("Integrated model for Bombus sylvestris") + 
  theme_bw() +
  theme(plot.title = element_text(lineheight = .8, face = "bold"),
        legend.position = 'bottom')

occ_a_plot <- ggplot(data = occ_model %>% slice(1:7), 
                     aes(x = 2010:2016, y = Mean)) +
  geom_point(size = 4, aes(color = rhat_threshold)) +
  geom_ribbon(aes(ymin = quant_025, ymax = quant_97.5), alpha = 0.2) +
  geom_line(size = 1, col = "blue") +
  geom_text(aes(label = round(precision, digits = 2)), nudge_y = 1) +
  scale_color_manual(name = 'Rhat', values = c('Bad' = 'red','Good' = 'blue')) +
  ylab("a") +
  xlab("Year") +
  scale_y_continuous(limits = c(-3, 3)) +
  ggtitle("Occupancy model for Bombus sylvestris") + 
  theme_bw() +
  theme(plot.title = element_text(lineheight = .8, face = "bold"),
        legend.position = 'bottom')

int_a_plot + occ_a_plot


models_output <- rbind(occ_model, int_model) %>%
  pivot_wider(names_from = model, 
              values_from = c(Mean, SD, Naive.SE, 
                              Time.series.SE, quant_025,
                              quant_25, quant_50, quant_75,
                              quant_97.5, PSRF, PSRF..97.5..quantile,
                              precision, rhat_threshold)) 

models_output$parameter_grouped <- sub("\\[.*", "", models_output$parameter)



ggplot(data = models_output, aes(x = precision_occupancy, y = precision_integrated)) +
  geom_point(aes(color = parameter_grouped, fill = parameter_grouped),
             position=position_jitter(h=1, w=1),
             shape = 21, alpha = 0.5, size = 3) +
  scale_fill_viridis_d() +
  scale_colour_viridis_d() +
  geom_smooth(method = "lm", colour = "black") +
  geom_abline(intercept = 0, slope = 1,linetype = "dashed") +
  # coord_equal() +
  ylab("Parameters precision (integrated model)") +
  xlab("Parameters precision (occupancy model)") +
  ggtitle("Bombus sylvestris") +
  theme_bw()



