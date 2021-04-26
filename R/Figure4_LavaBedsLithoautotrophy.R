pacman::p_load(tidyverse)

ipl_irms <- read_csv("../data/irms_samples_tank.csv") # raw tank values for samples
irms_stds <- read_csv("../data/irms_standards.csv") 
ipl_ug_g_tle <- read_csv("../data/ipl_abun_table_ug_g_tle.csv")
sample_metadata <- read_csv("../data/sample_metadata.csv")

# Convert sample tank d13C_IPL values to VPDB scale

# calculate average IPL d13C tank values and convert to VPDB scale
ipl_tank <- ipl_irms %>%
  pivot_longer(col = 3:ncol(.), names_to = "ipl", values_to = "avg_d13C_tank") %>% #note: not actually mean d13C tank values, just need to call it this for predict() to work correctly.
  group_by(sample_fraction_coded, ipl) %>%
  distinct(sample_fraction_coded, ipl, avg_d13C_tank, .keep_all = T) %>%
  #filter(!sample_fraction %in% c("POS-20180801-F-01-P3", "POS-20180801-F-01-P4")) %>%
  group_by(Batch) %>%
  mutate(ID = as.character(row_number()))

ipl_tank_nested <- ipl_tank %>%
  nest()
ipl_tank_split <- ipl_tank %>%
  split(.$Batch)


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

irms_joined <- left_join(irms_tank, irms_vpdb, by = c("Batch", "peak"))

irms_models <- irms_joined %>%
  group_by(Batch) %>%
  do(model = lm(d13C_vpdb ~ avg_d13C_tank, data = .)) %>%
  left_join(ipl_tank_nested) %>%
  filter(!Batch == "B5_1")

fitted <- irms_models %>%
  mutate(fit = list(as_tibble(predict(model, data)))) %>%
  pull(fit)

a6_1 <- fitted[[1]] %>%
  rownames_to_column("ID") %>%
  left_join(., filter(ipl_tank, Batch == "A6_1")) %>%
  rename(d13C_vpdb_pred = value)
a7_1 <- fitted[[2]] %>%
  rownames_to_column("ID") %>%
  left_join(., filter(ipl_tank, Batch == "A7_1")) %>%
  rename(d13C_vpdb_pred = value)
a7_2 <- fitted[[3]] %>%
  rownames_to_column("ID") %>%
  left_join(., filter(ipl_tank, Batch == "A7_2")) %>%
  rename(d13C_vpdb_pred = value)

ipl_vpdb_pred_data <- full_join(a6_1, a7_1) %>%
  full_join(a7_2)

# determine average d13C (VPDB) and standard error incorporating model and instrument errors

# calculate RMSE for model validation
find_rmse <- function(linear_model){
  rss <- c(crossprod(linear_model$residuals))
  mse <- rss / length(linear_model$residuals)
  sqrt(mse)
}

ipl_vpdb_avg <- ipl_vpdb_pred_data %>% 
  group_by(sample_fraction_coded, ipl) %>%
  mutate(avg_d13C_vpdb_pred = mean(d13C_vpdb_pred),
         sd_d13C_vpdb_pred = sd(d13C_vpdb_pred)) %>%
  distinct(sample_fraction_coded, ipl, avg_d13C_vpdb_pred, .keep_all = T) %>%
  mutate(model_rmse = case_when(Batch == "A6_1" ~ find_rmse(irms_models$model[[1]]),
                                Batch == "A7_1" ~ find_rmse(irms_models$model[[2]]),
                                Batch == "A7_2" ~ find_rmse(irms_models$model[[3]])),
         se_reported = model_rmse + sd_d13C_vpdb_pred) %>%
  filter(!is.na(se_reported)) # remove measurements without duplicates (SE = NA)

# filter out singletons
obs <- ipl_vpdb_avg %>% 
  group_by(ipl) %>% 
  count(ipl) # number of observations for each ipl d13C 
ipl_vpdb_avg <- ipl_vpdb_avg %>%
  left_join(., obs, by = "ipl") %>%
  filter(n > 1)

# merge d13C_IPL data with yield data normalized to TLE

frac_abbrev <- data.frame("glycolipids" = "P2", 
                          "phospholipids" = "P3", 
                          "phosphatidylcholine" = "P4") %>%
  pivot_longer(1:ncol(.), names_to = "fraction", values_to = "frac_abbrev")

tle_long <- ipl_ug_g_tle %>%
  pivot_longer(2:ncol(.), names_to = "ipl_fullname", values_to = "ipl_ug_g_tle") %>%
  mutate(fraction = str_extract(ipl_fullname, "[^.]*$"),
         ipl = str_remove(ipl_fullname, "[^.]*$"),
         ipl = str_remove(ipl, "[.]")) %>%
  left_join(frac_abbrev, by = "fraction") %>%
  mutate(sample_fraction_coded = str_c(sample_name, frac_abbrev, sep = "-"))
  
irms_samples <- ipl_vpdb_avg %>%
  ungroup() %>%
  distinct(sample_fraction_coded) %>%
  pull(sample_fraction_coded)

tle_long_irms <- tle_long %>%
  filter(sample_fraction_coded %in% irms_samples)

# merge d13C and ipl yield data

ipl_ug_g_tle_long <- ipl_ug_g_tle %>%
  pivot_longer(2:ncol(.), names_to = "ipl_fullname", values_to = "ug_g_tle")

yield_plus_d13C <- left_join(ipl_vpdb_avg, tle_long_irms) %>%
  left_join(sample_metadata) %>%
  left_join(ipl_ug_g_tle_long) %>%
  filter(!is.na(sample_name))

# export yield_plus_d13C for use in Figure 5
write_csv(yield_plus_d13C, "../data/ipl_yield_plus_d13C.csv")

##### Visualize stable isotopic composition of C reservoirs in the lava tubes #####

# manually order y-axis by sample type
stype_order <- rev(c("surface soil", 
                     "cave soil", 
                     "ooze", 
                     "polyp/coralloid", 
                     "mineral", 
                     "white biofilm", 
                     "mixed biofilm", 
                     "yellow biofilm", 
                     "tan biofilm"))

# color pallete for plot
point_color <- c("#8F2D56", #dic
                 "#D81159", #doc
                 "#FFBC42", #ipl
                 "#218380") #toc

# plot
yield_plus_d13C %>% 
  filter(!sample_color %in% c("blue-green")) %>%
  ggplot() +
  geom_violin(data = sample_metadata, aes(d13C_toc, factor(sample_class2, levels = stype_order), 
                                          alpha = 0.8, fill = "TOC"), trim = T) +
  geom_violin(aes(avg_d13C_vpdb_pred, factor(sample_class2, levels = stype_order), 
                  alpha = 0.8, fill = "IPL"), trim = T) +
  geom_rect(aes(xmin = -29.79-2.04, xmax = -29.79+2.04, 
                ymin = 0, ymax = Inf, 
                fill = "DOC")) + # measured DOC values +
  geom_rect(aes(xmin = -6.05-1.85, xmax = -6.05+1.85, 
                ymin = 0, ymax = Inf, 
                fill = "DIC")) + # measured DIC values
  geom_vline(aes(xintercept = -29.79), color = "#000000") +
  geom_vline(aes(xintercept = -6.05), color = "#5c1919") +
  geom_violin(aes(avg_d13C_vpdb_pred, factor(sample_class2, levels = stype_order), 
                  alpha = 0.8, fill = "IPL"), trim = T) +
  geom_violin(data = sample_metadata, aes(d13C_toc, factor(sample_class2, levels = stype_order), 
                                          alpha = 0.8, fill = "TOC"), trim = T) +
  scale_fill_manual(values = point_color) +
  guides(alpha = FALSE) +
  labs(fill = "Reservoir") +
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  xlab("d13C (permille)")
