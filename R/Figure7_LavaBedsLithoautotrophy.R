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
# Monument, CA (Lava Beds). I model the fraction of lithoautotrophy (fL) that 
# explains observed d13C IPL values, using Equation 1 from the main text. 
# I assume both Calvin cycle (CC)- and Wood Ljungdahl (WL) pathway-based 
# lithoautotrophy. Calculated fL values are then exported as csv-formatted 
# tables and are summarized as static figures (Figure 7 = CC; Figure S5 = WL pathway). 

## ------------------ ##

# load packages and data

pacman::p_load(tidyverse, ggridges)

aqueous <- read_csv("data/aqueous_data.csv")
yield_plus_d13C <- read_csv("data/ipl_yield_plus_d13C.csv") %>% # from Figure4_LavaBedsLithoautotrophy.R
  filter(!sample_fraction_coded %in% c("L853-20180801-F-01-P3", "L853-20180801-F-01-P4"))
ipl_metadata <- read_csv("data/ipl_metadata.csv")

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
  left_join(ipl_metadata, by = c("ipl" = "core_lipid")) %>%
  mutate(d13C_ipl_corr = ((1/(carbon_number+methanol_C_count))*(-38.9)) + (1-(1/(carbon_number+1)))*avg_d13C_vpdb_pred) %>%
  distinct(sample_fraction_coded, sample_name, sample_class2, ipl, avg_d13C_vpdb_pred, d13C_ipl_corr, sd_d13C_vpdb_pred, model_rmse, se_reported)

# Assuming Calvin cycle-based fixation,
# (dIPL - dS) / (dL - dS) = fL_CC
# dIPL = avg_d13C_vpdb_pred
# dS = drip_DOC_d13C_avg
# dL = dS - 20
# fL_CC = fraction of isotopic signal attributed to CC-based lithoautotrophy

# Assuming Wood Ljungdahl-based fixation,
# (dIPL - dS) / (dL - dS) = fL_WL
# dIPL = avg_d13C_vpdb_pred
# dS = drip_DOC_d13C_avg
# dL = dS - 52 (maximum rA-CoA fractionation from Hayes 2001)
# fL_WL = fraction of isotopic signal attributed to WL-based lithoautotrophy

model_output <- model_data %>%
  mutate(fL_CC = if_else(sample_class2 == "surface soil", 
                         (d13C_ipl_corr - soil_d13Corg$mean_d13C_toc_soil) / ((soil_d13Corg$mean_d13C_toc_soil-20) - soil_d13Corg$mean_d13C_toc_soil),
                         (d13C_ipl_corr - mean_drip_doc_d13C) / ((mean_drip_doc_d13C-20) - mean_drip_doc_d13C)),
         fL_WL = if_else(sample_class2 == "surface soil", 
                         (d13C_ipl_corr - soil_d13Corg$mean_d13C_toc_soil) / ((soil_d13Corg$mean_d13C_toc_soil-52) - soil_d13Corg$mean_d13C_toc_soil),
                         (d13C_ipl_corr - mean_drip_doc_d13C) / ((mean_drip_doc_d13C-52) - mean_drip_doc_d13C)))

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
            fL_CC_mean_bin = mean(fL_CC),
            fL_CC_sd_bin = sd(fL_CC),
            fL_WL_mean_bin = mean(fL_WL), 
            fL_WL_sd_bin = sd(fL_WL))

# filter for IPLs specifically discussed in the text for Figure 7
text_ipls <- c("n-C18:0", "n-C22:0", "n-C16:0",
               "i-C16:0", "i-C15:0", "10_14_DiMe_C17", "C16:1 w7 (trans)", "C17:1 w8 (trans)", #"i-C18:0",
               "C16:1 w7 (cis)")

depleted_ipls <- c("i-C16:0", "10_14_DiMe_C17", "C16:1 w7 (trans)", "C17:1 w8 (trans)")

# plot Figure 7
cc_plot <- model_output %>% 
  filter(se_reported < 0.3,
         ipl %in% text_ipls
         ) %>%
  group_by(sample_class2) %>%
  ggplot() +
  geom_density_ridges(aes(fL_CC, 
                          reorder(ipl, fL_CC, mean), 
                          fill = paste(sample_class2)),
                      stat = "density_ridges") +
  geom_point(aes(fL_CC,
                 reorder(ipl, fL_CC, mean), 
                 color = sample_class2, 
                 text = paste("<br>",
                              "IPL: ", ipl,
                              "<br>",
                              "fL: ", round(fL_CC, 3),
                              "<br>",
                              "sample type: ", sample_class2,
                              "<br>",
                              "sample name: ", sample_name))) +
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

# export mass balance calculations assuming CC and WL pathway
model_table_CC <- model_output %>%
  select(sample_fraction_coded, sample_class2, ipl, fL_CC) %>%
  pivot_wider(names_from = "ipl", values_from = "fL_CC")
write_csv(model_table_CC, "data/supplemental_data/TableS5_fL_CC.csv")

model_table_WL <- model_output %>%
  select(sample_fraction_coded, sample_class2, ipl, fL_WL) %>%
  pivot_wider(names_from = "ipl", values_from = "fL_WL")
write_csv(model_table_WL, "data/supplemental_data/TableS6_fL_rACoA.csv")

##### Other code #####

# summarize fL statistics for groups of IPLs discussed in the manuscript
model_output %>%
  #filter(ipl %in% depleted_ipls & se_reported < 0.3) %>%
  filter(ipl %in% c("n-C18:0", "n-C20:0", "n-C22:0", "n-C24:0", "n-C26:0") & se_reported < 0.3) %>%
  left_join(select(yield_plus_d13C, sample_name, sample_class)) %>%
  select(sample_fraction_coded, sample_class, ipl, fL_CC, fL_WL) %>%
  group_by(sample_class) %>%
  #filter(sample_class == "cave biofilm") %>% arrange((fL_WL))
  summarize(fL_CC_mean = mean(fL_CC, na.rm = T), 
            fL_CC_sd = sd(fL_CC, na.rm = T), 
            fL_WL_mean = mean(fL_WL, na.rm = T), 
            fL_WL_sd = sd(fL_WL, na.rm = T))

# Plot again as "wl_plot" (Wood-Ljungdahl/reductive acetyl CoA pathway) for Supplemental Figure S5
wl_plot <- model_output %>% 
  filter(se_reported < 0.3,
         ipl %in% text_ipls) %>%
  group_by(sample_class2) %>%
  ggplot() +
  geom_density_ridges(aes(fL_WL, 
                          reorder(ipl, fL_WL, mean), 
                          fill = paste(sample_class2)),
                      stat = "density_ridges") +
  geom_point(aes(fL_WL,
                 reorder(ipl, fL_WL, mean), 
                 color = sample_class2, 
                 text = paste("<br>",
                              "IPL: ", ipl,
                              "<br>",
                              "fL: ", round(fL_WL, 3),
                              "<br>",
                              "sample type: ", sample_class2,
                              "<br>",
                              "sample name: ", sample_name))) +
  scale_fill_manual(values = sample_class2_colorz) +
  scale_color_manual(values = sample_class2_colorz) +
  theme_bw() +
  labs(color = "Sample Type") +
  guides(shape = FALSE,
         fill = FALSE,
         color = guide_legend(override.aes = list(size = 5,
                                                  shape = 16))) +
  ylab("IPL") +
  xlab("fL_WL")



##### Plot fL model visualization with HTML table of model data #####
my_plotly <- plotly::ggplotly(cc_plot, tooltip = "text") %>%
  hide_legend()
  # note: plotly does not support density ridges; only points show in interactive figure

interactive_model_table <- model_output %>%
  select(sample_name, sample_class2, ipl, fL_CC,  d13C_ipl_corr, se_reported) %>%
  filter(!is.na(fL_CC)) %>%
  mutate(fL_CC = round(fL_CC, 3),
         d13C_ipl_corr = round(d13C_ipl_corr, 1), 
         se_reported = round(se_reported, 1)) %>%
  rename(`Sample Name` = sample_name, 
         `Sample Type` = sample_class2,
         IPL = ipl,
         `d13C (permil)` = d13C_ipl_corr,
         `d13C Standard Error` = se_reported) %>%
  DT::datatable(options = list(lengthMenu = c(15, 20), pageLength = 5))


full_interactive_model <- htmltools::browsable(
  htmltools::tagList(list(
    htmltools::tags$div(
      style = 'width:100%;display:block;float:left;',
      my_plotly
    ),
    htmltools::tags$div(
      style = 'width:100%;display:block;float:left;',
      interactive_model_table
    )
  ))
)

htmltools::save_html(full_interactive_model, "figures/final_figures/full_interactive_fL_model.html")




