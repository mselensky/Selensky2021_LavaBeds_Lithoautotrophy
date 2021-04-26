## ------------------ ##

# Title: Dimensionality reduction of lipid compositions across soils and lava cave features
# Dataset: Lava Beds Lithoautotrophy
# Date: 26 April 2021
# Author: Matt Selensky
# email: mselensky@u.northwestern.edu
# Copyright (c) Matt Selensky, 2021

# This script produces Figure 3 from the manuscript "Stable Carbon Isotope Depletions in 
# Lipid Biomarkers Suggest Subsurface Carbon Fixation". To do so, 

## ------------------ ##

##### load packages and data #####

# load packages
pacman::p_load(tidyverse, ggpubr, RColorBrewer, vegan, ggrepel)

# load data
ipl_ug_g_tle <- read_csv("../data/ipl_abun_table_ug_g_tle.csv")
ipl_rel_abun <- read_csv("../data/ipl_rel_abun.csv")
sample_metadata <- read_csv("../data/sample_metadata.csv")
ipl_metadata <- read_csv("../data/ipl_metadata.csv")

# color palette data
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
ipl_class_colorz <- c("branched" = "#38514C", # brown/green palette
                      "normal" = "#bf812d",
                      "unsaturated" = "#35978f",
                      "diacid" = "#8c510a",
                      "cyclic" = "#000000",
                      "sterol" = "#000000",
                      "alkanone" = "#000000",
                      "other" = "#000000")

# prepare relative lipid abundance matrix for Bray-Curtis calculation
ipl_rel_abun_mat <- ipl_rel_abun %>%
  column_to_rownames("sample_name") %>%
  as.matrix()

# determine sum of IPLs* extracted for each sample across ipl fractions
# * normalized to TLE (ug/g)
ipl_sum_frac <- ipl_ug_g_tle %>%
  pivot_longer(2:ncol(.), names_to = "ipl", values_to = "ug_g_tle") %>%
  left_join(., select(sample_metadata, sample_name, tle_g), by = "sample_name") %>%
  group_by(sample_name) %>% 
  summarize(ipl_ug_g_tle_sum = sum(ug_g_tle))

ipl_rel_abun_nmds <- metaMDS(log10(ipl_rel_abun_mat+1), "bray", 2, trymax = 1000)
# log10-transformed solution stress: 0.1991953

# extract scores and loadings from NMDS object

frac_scores <- scores(ipl_rel_abun_nmds) %>%
  as.data.frame() %>%
  rownames_to_column("sample_name") %>%
  left_join(., sample_metadata, "sample_name") %>%
  left_join(., ipl_sum_frac) %>%
  distinct(sample_name, .keep_all = TRUE)

# get ipl vector loadings
frac_env <- envfit(ipl_rel_abun_nmds, ipl_rel_abun_mat) 

# function to extract NMDS loadings
get_nmds_vectors <- function(env_output, top, pvalue) {
  vectors <- env_output %>% scores(display = "vectors") %>%
    as.data.frame() %>%
    rownames_to_column(var = "variable")
  pvals <- env_output$vectors$pvals %>%
    as.data.frame() %>%
    rename(., pval = `.`) %>%
    rownames_to_column(var = "variable")
  vectors <- left_join(vectors, pvals, by = "variable") %>%
    rename(NMDS1_vec = "NMDS1", NMDS2_vec = "NMDS2") %>%
    filter(pval < pvalue)
  vectors <- top_n(vectors, -(top))
  vectors
}

frac_vecs <- get_nmds_vectors(frac_env, 100, 0.05) %>%
  rename(ipl = variable) %>%
  left_join(., ipl_metadata) %>%
  mutate(vec_length = abs(NMDS1_vec) + abs(NMDS2_vec)) %>%
  top_n(vec_length, n = 10)

# plot NMDS
nmds_plot <- ggplot() +
  geom_segment(data = frac_vecs, 
               aes(x = 0, xend = NMDS1_vec,
                   y = 0, yend = NMDS2_vec, 
                   color = ipl_class,
                   text = paste0("<br>", 
                                 "Core lipid: ", ipl,
                                 "<br>",
                                 "IPL fraction: ", fraction,
                                 "<br>",
                                 "IPL class: ", ipl_class, 
                                 "<br>",
                                 "p-value: ", pval))) +
  geom_text_repel(data = frac_vecs, aes(x = NMDS1_vec,
                                        y = NMDS2_vec,
                                        color = ipl_class,
                                        label = ipl)) +
  geom_point(data = frac_scores, 
             shape = 21, 
             aes(NMDS1, NMDS2, 
                 size = ipl_ug_g_tle_sum/1000, 
                 fill = sample_class2,
                 text = paste0("<br>",
                               "sample name: ", sample_name,
                               "<br>", 
                               "sample type: ", sample_class2,
                               "<br>", 
                               "TOC d13C (permil): ", d13C_toc))) +
  scale_color_manual(values = ipl_class_colorz) +
  scale_fill_manual(values = sample_class2_colorz) +
  #scale_shape_manual(values = shapes) +
  guides(shape = FALSE, 
         fill = guide_legend(override.aes = list(shape = 21, 
                                                 size = 5))) +
  labs(fill = "Sample Type",
       color = "Fatty Acid Type") +
  xlab("NMDS1") +
  ylab("NMDS2") +
  scale_size(name = "IPL abundance (μg/g TLE)", 
             trans = "log10") +
  theme_bw() 


# visualize distribution of total IPL yields

ipl_sum_frac %>%
  mutate(ipl_mg_g_tle_sum = ipl_ug_g_tle_sum/1000) %>%
  left_join(., select(sample_metadata, sample_name, sample_class2)) %>%
  ggplot() +
  geom_histogram(aes(x = ipl_mg_g_tle_sum,
                     fill = sample_class2, color = "black"), bins = nrow(ipl_sum_frac)) +
  scale_x_log10() +
  scale_fill_manual(values = sample_class2_colorz) +
  scale_color_manual(values = "black") +
  theme_bw()

ipl_sum_frac %>%
  mutate(ipl_mg_g_tle_sum = ipl_ug_g_tle_sum/1000) %>%
  left_join(., select(sample_metadata, sample_name, sample_class2)) %>%
  group_by(sample_class2) %>%
  filter(!sample_name %in% c("V460-20190810-A-32")) %>%
  summarize(ipl_mg_g_tle_bin_mean = mean(ipl_mg_g_tle_sum),
            ipl_mg_g_tle_bin_sd = sd(ipl_mg_g_tle_sum),
            ipl_mg_g_tle_sum = sum(ipl_mg_g_tle_sum)) %>%
  arrange(-ipl_mg_g_tle_bin_mean)

ipl_sum_frac %>%
  mutate(ipl_mg_g_tle_sum = ipl_ug_g_tle_sum/1000) %>%
  left_join(., select(sample_metadata, sample_name, sample_class2)) %>%
  group_by(sample_name) %>%
  filter(sample_class2 %in% c("surface soil", "cave soil")) %>%
  summarize(ipl_mg_g_tle_sum = sum(ipl_mg_g_tle_sum)) %>%
  arrange(-ipl_mg_g_tle_sum)




