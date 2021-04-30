## ------------------ ##

# --- Script information
# Title: Comparing fatty acid d13C values across lava cave samples and surface soils
# Dataset: Lava Beds Lithoautotrophy
# Date: 30 Apr 2021
# Author: Matt Selensky
# email: mselensky@u.northwestern.edu
# Copyright (c) Matt Selensky, 2021

# --- Summary
# This script produces a bubbleplot that compares the d13C values of individual
# fatty acids. Samples are first clustered by log10-transformed relative lipid 
# abundance to inform the x-axis ordering of the final plot. Next, the relative
# lipid abundance matrix is transposed to produce a second dendrogram that will 
# be used to order the IPL-derived fatty acids themselves on the y-axis.
# Points in the bubbleplot are sized by fatty acid abundance relative to the 
# total lipid extract (TLE) and are colored by compound-specific d13C values.

## ------------------ ##

##### load packages and data ######

pacman::p_load(tidyverse, broom, ggdendro, ggpubr, ggtree, treeio)

yield_plus_d13C <- read_csv("../data/ipl_yield_plus_d13C.csv") # from Figure4_LavaBedsLithoautotrophy.R
ipl_rel_abun <- read_csv("../data/ipl_rel_abun.csv")
ipl_ug_g_tle <- read_csv("../data/ipl_abun_table_ug_g_tle.csv")
sample_metadata <- read_csv("../data/sample_metadata.csv")
ipl_metadata <- read_csv("../data/ipl_metadata.csv") %>%
  rename(ipl_fullname = ipl)

# color palette data for plotting
sample_class2_colors <- c("tan biofilm" = "#17174a", 
                          "white biofilm" = "#b0b0d4", 
                          "yellow biofilm" = "#5858a6", 
                          "mixed biofilm" = "#191991",
                          "surface soil" = "#F9A304", 
                          "cave soil" = "#c28513", 
                          "cave sludge" = "#c28513",
                          "mineral" = "#f2f2f2", 
                          "polyp/coralloid" = "#bfbfbf", 
                          "ooze" = "#2d7529")

##### dendrogram construction #####
# We will want to construct two dendrograms to inform the ordering of the 
# columns and rows in the final figure. 
# One dendrogram will represent the hierarchcal clustering of each 
# sample_fraction (the columns), while the other dendrogram will cluster 
# based on the IPLs themselves (the rows). 

# create matrix of relative IPL abundances by sample_fraction. 
# (this will be used for the first dendrogram)

# relative abundance matrix for each individual sample_fraction
frac_abbrev <- data.frame("glycolipids" = "P2", 
                          "phospholipids" = "P3", 
                          "phosphatidylcholine" = "P4") %>%
  pivot_longer(1:ncol(.), names_to = "fraction", values_to = "frac_abbrev")

ipl_rel_long <- ipl_ug_g_tle %>%
  pivot_longer(2:ncol(.), names_to = "ipl_fullname", values_to = "ug_g_tle") %>%
  mutate(fraction = str_extract(ipl_fullname, "[.].*"),
         fraction = str_remove(fraction, "[.]")) %>%
  left_join(frac_abbrev, by = "fraction") %>%
  mutate(sample_fraction_coded = str_c(sample_name, frac_abbrev, sep = "-")) %>%
  group_by(sample_fraction_coded) %>%
  mutate(rel_abun_frac = ug_g_tle/sum(ug_g_tle)) %>%
  filter(sample_fraction_coded %in% pull(distinct(yield_plus_d13C, sample_fraction_coded)),
         rel_abun_frac > 0)

ipl_rel_mat <- ipl_rel_long %>%
  filter(!sample_fraction_coded %in% c("L853-20180801-F-01-P3", "L853-20180801-F-01-P4")) %>%
  left_join(ipl_metadata) %>%
  select(sample_fraction_coded, core_lipid, rel_abun_frac) %>%
  pivot_wider(names_from = "core_lipid", values_from = "rel_abun_frac") %>%
  column_to_rownames("sample_fraction_coded") %>%
  replace(is.na(.), 0) %>% as.matrix() 

# transpose ipl_rel_mat_frac to perform distance calculations on IPLs 
# (this will be used for the second dendrogram)
ipl_rel_mat_ipl <- ipl_rel_mat %>% t()

# create dendrogram based on relative IPL abundance for each sample_fraction (dendrogram 1)
euc_dist1 <- dist(log10(ipl_rel_mat+1), "manhattan") #log10-transformed
dendrogram1 <- hclust(euc_dist1, "complete")

#f lipping dendrogram with ggtree + treeio for aesthetics
ggt1 <- ggtree(dendrogram1) 
ggt1_df <- get.tree(ggt1)$tip.label %>%
  as.data.frame() %>%
  rename(sample_fraction_coded = ".") %>%
  left_join(yield_plus_d13C)
p <- ggt1 %<+% ggt1_df + 
  geom_tippoint(aes(color = sample_class2)) +
  #geom_tiplab(aes(text = label)) +
  scale_color_manual(values = sample_class2_colors) +
  #geom_text(aes(label = node)) +
  theme(legend.position = "none") +
  coord_flip() 

names <- data.frame(sample_fraction_coded = get_taxa_name(p))

# create dendrogram based on relative IPL abundances across sample fractions (dendrogram 2)
euc_dist2 <- dist(log10(ipl_rel_mat_ipl + 1), "manhattan") #log10-transformed
dendrogram2 <- hclust(euc_dist2, "complete")
dendrogram2_data <- dendro_data(as.dendrogram(dendrogram2))

#`labs2` contains the IPL names ordered by hclust() (dendrogram2)
labs2 <- dendrogram2_data[["labels"]]
labs2 <- as.character(labs2$label)

# merge `names` with ggt1$data to join `y` values, which is what we need for ordering subplots in figure
plot_label_data <- ggt1$data %>%
  rename(sample_fraction_coded = label) %>%
  filter(!is.na(sample_fraction_coded))

merged_labs <- names %>% 
  left_join(., plot_label_data) %>%
  rownames_to_column("order") %>%
  select(sample_fraction_coded, order, node) 

# correct for isotopic contribution of methanol C in derivatized FAMEs (d13C = -38.9 permille)
bubble_plot_data <- yield_plus_d13C %>%
  left_join(., merged_labs, "sample_fraction_coded") %>%
  left_join(., ipl_metadata, "ipl_fullname") %>%
  filter(!is.na(avg_d13C_vpdb_pred),
         #!ipl %in% excluded_ipls, 
         !se_reported > 0.3) %>%
  mutate(d13C_ipl_corr = ((1/(carbon_number+methanol_C_count))*(-38.9)) + (1-(1/(carbon_number+1)))*avg_d13C_vpdb_pred)

# export `datatable` as supplemental figure s2
datatable <- bubble_plot_data %>%
  group_by(sample_fraction_coded, ipl) %>%
  mutate(d13C_plus_se = str_c(round(d13C_ipl_corr, 2), 
                              round(se_reported, 2), 
                              sep = " ± ")) %>%
  select(sample_fraction_coded, ipl, d13C_plus_se) %>%
  pivot_wider(names_from = "ipl", values_from = "d13C_plus_se")
write_csv(datatable, "../data/supplemental_data/TableS2_irms_d13C.csv")
##### plotting #####

# plot `sample_fraction` dendrogram

# make bubble plot of IPLs and sample_fraction ordered by hclust() results
p3 <- bubble_plot_data %>%
  
  ggplot(aes(-as.numeric(order), factor(ipl, levels = labs2))) + 
  geom_point(aes(color = d13C_ipl_corr, size = ipl_ug_g_tle/1000,
                 text = paste('<br>', 'IPL: ', ipl,
                              '<br>',
                              'IPL class: ', ipl_class,
                              '<br>',
                              'sample: ', sample_fraction_coded,
                              '<br>',
                              'sample type: ', sample_class2, 
                              '<br>',
                              'd13C (permille): ', format(round(d13C_ipl_corr, 2)),
                              '<br>',
                              'standard error: ', format(round(se_reported, 2))))) +
  scale_color_gradient(low = "blue", high = "#F9A304") +
  theme_bw() + theme(legend.position = "none", 
                     axis.title = element_blank(),
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank()) +
  labs(color = expression(paste("δ"^13, "C (‰)")),
       size = "IPL yield (mg/g TLE)")

# plot dendrogram1 and irms bubble plot together  
ggarrange(p3,
          p,
          nrow = 2, common.legend = T, align = "hv", heights = c(4,1))
