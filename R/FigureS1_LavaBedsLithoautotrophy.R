##### load packages and data #####
pacman::p_load(tidyverse, plotly, broom, ggdendro, RColorBrewer, ggpubr)
irms_stds <- read_csv("../data/irms_standards.csv")

##### IRMS tank values conversion #####

# split irms standards into known VPDB d13C values and measured tank values
irms_vpdb <- irms_stds %>% filter(ID == "VPDB") %>%
  select(1:23, -c(6:8)) %>%
  pivot_longer(cols = 6:ncol(.), names_to = "peak", values_to = "d13C_vpdb") %>%
  select(Batch, peak, d13C_vpdb)

irms_tank <- irms_stds %>%
  filter(!ID == "VPDB") %>%
  select(1:23, -c(6:8)) %>%
  pivot_longer(cols = 6:ncol(.), names_to = "peak", values_to = "d13C_tank") %>%
  group_by(Batch, peak) %>%
  mutate(avg_d13C_tank = mean(d13C_tank),
         sd_d13C_tank = sd(d13C_tank))

# predict d13C values (VPDB) from tank values from respective linear models

# merge residuals from standard linear models with i45 data
irms_i45 <- irms_stds %>%
  filter(!ID == "VPDB") %>%
  select(Analysis, ID, Batch, 27:ncol(.)) %>%
  pivot_longer(cols = 4:ncol(.), names_to = "peak", values_to = "i45_intensity") %>%
  mutate(peak = str_remove(peak, "_i45"))

irms_joined <- left_join(irms_tank, irms_vpdb, by = c("Batch", "peak")) %>%
  left_join(., select(irms_i45, Batch, peak, i45_intensity))

irms_models <- irms_joined %>%
  group_by(Batch) %>%
  do(model = lm(d13C_vpdb ~ avg_d13C_tank, data = .), d = (.))

fitted_output <- augment(irms_models, model, data = d)

batch_colors <- c("A6_1" = "#7FC97F", 
                  "A7_1" = "#FDD58C", 
                  "A7_2" = "#B2258F", 
                  "B5_1" = "#666666")

fitted_output %>%
  rename(Dd13C = .resid) %>%
  mutate(Batch = as.character(Batch)) %>%
  ggplot(aes(i45_intensity, Dd13C, color = Batch)) +
  geom_point() +
  scale_color_manual(values = batch_colors) +
  xlab("i45 (mV)") +
  ylab(expression(paste("Δδ"^13, "C (‰)"))) +
  theme_bw() +
  ggtitle(expression(paste("Δδ"^13, "C vs. i45 signal intensity"))) 





