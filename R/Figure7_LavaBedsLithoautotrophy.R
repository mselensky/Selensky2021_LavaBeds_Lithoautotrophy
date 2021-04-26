## ------------------ ##

# --- Script information
# Title: Modeling the influence of in situ lithoautotrophy on IPL fatty acids
# Dataset: Lava Beds Lithoautotrophy
# Date: 26 Apr 2021
# Author: Matt Selensky
# email: mselensky@u.northwestern.edu
# Copyright (c) Matt Selensky, 2021

# --- Summary
# This script ultimately produces a model that estimates the relative fraction of 
# in situ lithoautotrophy on IPL fatty acids extracted from Lava Beds National 
# Monument, CA (Lava Beds). Here, I first visualize the d13C compositions of C 
# reservoirs at Lava Beds. I then model the fraction of lithoautotrophy (fL) that 
# explains observed d13C IPL values, using Equation 1 from the main text. Calculated 
# fL values then exported as a .csv-formatted table and are summarized as both 
# static and interactive figures. 

## ------------------ ##

# load packages and data

pacman::p_load(tidyverse, ggridges)

aqueous <- read_csv("../data/aqueous_data.csv")
yield_plus_d13C <- read_csv("../data/ipl_yield_plus_d13C.csv") %>% # from Figure4_LavaBedsLithoautotrophy.R
  filter(!sample_fraction_coded %in% c("L853-20180801-F-01-P3", "L853-20180801-F-01-P4"))

# define colors for samples 
sample_class2_colorz <- c("tan biofilm" = "#17174a", 
                          "white biofilm" = "#b0b0d4", 
                          "yellow biofilm" = "#5858a6", 
                          "mixed biofilm" = "#191991",
                          "surface soil" = "#F9A304", 
                          "cave soil" = "#c28513", 
                          "cave sludge" = "#c28513",
                          "mineral" = "#f2f2f2", 
                          "polyp/coralloid" = "#bfbfbf", 
                          "ooze" = "#2d7529")

aqueous_global <- aqueous %>%
  pivot_longer(cols = 2:3, names_to = "doc_d13C", values_to = "drip_doc_d13C") %>%
  pivot_longer(cols = 2:3, names_to = "doc_ppm", values_to = "drip_doc_ppm") %>%
  mutate(drip_doc_d13C = replace(drip_doc_d13C, drip_doc_d13C < -45, NA)) %>%
  summarise(mean_drip_doc_d13C = mean(drip_doc_d13C, na.rm = TRUE),
            stdv_drip_doc_d13C = sd(drip_doc_d13C, na.rm = TRUE),
            mean_drip_dic_d13C = mean(drip_dic_d13C, na.rm = TRUE),
            stdv_drip_dic_d13C = sd(drip_dic_d13C, na.rm = TRUE),
            mean_drip_doc_ppm = mean(drip_doc_ppm, na.rm = TRUE), 
            stdv_drip_doc_ppm = sd(drip_doc_ppm, na.rm = TRUE),
            mean_drip_dic_mM = mean(drip_dic_mM, na.rm = TRUE),
            stdv_drip_dic_mM = sd(drip_dic_mM, na.rm = TRUE))

# save aqueous_global summaries as global values

mean_drip_doc_d13C = aqueous_global$mean_drip_doc_d13C
stdv_drip_doc_d13C = aqueous_global$stdv_drip_doc_d13C
mean_drip_dic_d13C = aqueous_global$mean_drip_dic_d13C
stdv_drip_dic_d13C = aqueous_global$stdv_drip_dic_d13C
mean_drip_doc_ppm = aqueous_global$mean_drip_doc_ppm
stdv_drip_doc_ppm = aqueous_global$stdv_drip_doc_ppm
mean_drip_dic_mM = aqueous_global$mean_drip_dic_mM
stdv_drip_dic_mM = aqueous_global$stdv_drip_dic_mM

soil_d13Corg <- yield_plus_d13C %>%
  filter(sample_class2== "surface soil") %>%
  summarise(mean_d13C_toc_soil = mean(d13C_toc, na.rm = TRUE),
            stdv_d13C_toc_soil = sd(d13C_toc, na.rm = TRUE))

model_data <- yield_plus_d13C %>%
  select(sample_fraction_coded, sample_name, sample_class2, ipl, avg_d13C_vpdb_pred, sd_d13C_vpdb_pred, model_rmse, se_reported)

# Assuming Calvin cycle-based fixation,
# (dIPL - dS) / (dL - dS) = fL
# dIPL = avg_d13C_vpdb_pred
# dS = drip_DOC_d13C_avg
# dL = dS - 20
# fL = fraction of isotopic signal attributed to lithoautotrophy

model_output <- model_data %>%
  mutate(fL = if_else(sample_class2 == "surface soil", 
                      (avg_d13C_vpdb_pred - soil_d13Corg$mean_d13C_toc_soil) / ((soil_d13Corg$mean_d13C_toc_soil-20) - soil_d13Corg$mean_d13C_toc_soil),
                      (avg_d13C_vpdb_pred - mean_drip_doc_d13C) / ((mean_drip_doc_d13C-20) - mean_drip_doc_d13C)))

model_output$sample_class2 <- factor(model_output$sample_class2, levels = 
                                       c("surface soil", 
                                         "cave soil", 
                                         "polyp/coralloid", 
                                         "mineral", 
                                         "ooze",
                                         "mixed biofilm", 
                                         "tan biofilm", 
                                         "yellow biofilm", 
                                         "white biofilm"))

fL_binned <- model_output %>%
  filter(grepl("biofilm", sample_class2) | 
           sample_class2 %in% c("surface soil")) %>%
  mutate(sample_class = if_else(str_detect(sample_class2, "biofilm"),
                                str_extract(sample_class2, "biofilm"),
                                str_extract(sample_class2, ".*"))) %>%
  group_by(sample_class, ipl) %>%
  summarize(n_ipls = n(),
            fL_mean_bin = mean(fL),
            fL_sd_bin = sd(fL))

# filter for IPLs discussed in the text
text_ipls <- c("n-C18:0", "n-C22:0", "n-C16:0",
               "i-C16:0", "i-C15:0", "10_14_DiMe_C17", "C16:1 w7 (trans)", "C17:1 w8 (trans)", #"i-C18:0",
               "C16:1 w7 (cis)")

plot <- model_output %>% 
  filter(se_reported < 0.3,
         ipl %in% text_ipls) %>%
  group_by(sample_class2) %>%
  ggplot() +
  geom_density_ridges(aes(fL, 
                          reorder(ipl, fL, mean), 
                          fill = paste(sample_class2))) +
  geom_point(aes(fL,
                 reorder(ipl, fL, mean), 
                 color = sample_class2)) +
  #shape = as.character(ipl_cluster))) +
  scale_fill_manual(values = sample_class2_colorz) +
  scale_color_manual(values = sample_class2_colorz) +
  theme_bw() +
  labs(color = "Sample Type") +
  guides(shape = FALSE,
         fill = FALSE,
         color = guide_legend(override.aes = list(size = 5,
                                                  shape = 16))) +
  ylab("IPL") +
  xlab("fL")

model_table <- model_output %>%
  select(sample_fraction_coded, ipl, fL) %>%
  pivot_wider(names_from = "ipl", values_from = "fL")
#write_csv(model_table, "../data/supplemental_data/TableS4_fL.csv")



