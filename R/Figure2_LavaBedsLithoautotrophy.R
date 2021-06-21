## ------------------ ##

# Title: Lipid distributions across soils and lava cave features
# Dataset: Lava Beds Lithoautotrophy
# Date: 26 April 2021
# Author: Matt Selensky
# email: mselensky@u.northwestern.edu
# Copyright (c) Matt Selensky, 2021

# This script produces Figure 2 from the manuscript "Stable Carbon Isotope Depletions in 
# Lipid Biomarkers Suggest Subsurface Carbon Fixation". To do so, a dendrogram is first 
# constructed using a log10-transformed relative abundance matrix of extracted fatty
# acids from Lava Beds (1). These hierarchical clustering results are used to order 
# samples in the final figure. Next, I create a heatmap that displays relative lipid
# abundance binned by fatty acid type for each sample's IPL fraction (2). I then
# produce an areaplot to display relative abundance of binned fatty acids across 
# IPL fractions (3) and a lineplot to report lipid yields (ug/g TLE) for each
# sample (4). All plots are merged together in a single figure (Figure 2). 

## ------------------ ##

##### load packages and data #####
# load packages
pacman::p_load(tidyverse, ggpubr, RColorBrewer, ggdendro, plotly, ggtree, treeio, vegan, ggrepel)

# load data
ipl_ug_g_tle <- read_csv("../data/ipl_abun_table_ug_g_tle.csv")
ipl_rel_abun <- read_csv("../data/ipl_rel_abun.csv")
sample_metadata <- read_csv("../data/sample_metadata.csv")
ipl_metadata <- read_csv("../data/ipl_metadata.csv")

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

##### 1. Create dendrogram to order samples on y-axis of final figure #####

# create relative abundance matrix for distance calculations
ipl_rel_abun_mat <- ipl_rel_abun %>% 
  column_to_rownames("sample_name") %>%
  as.matrix() 

# dendrogram of ipl_rel_abun_mat
ipl_rel_abun_dist <- dist(log10(ipl_rel_abun_mat+1), "manhattan")
dendrogram <- hclust(ipl_rel_abun_dist, "complete")
ggt1 <- ggtree(dendrogram) 
ggt1_df <- get.tree(ggt1)$tip.label %>%
  as.data.frame() %>%
  rename(sample_name = ".") %>%
  left_join(sample_metadata)

p1 <- ggt1 %<+% ggt1_df + 
  geom_tippoint(aes(color = sample_class2)) +
  #geom_tiplab(aes(text = label)) +
  scale_color_manual(values = sample_class2_colors) +
  #geom_text(aes(label = node)) +
  theme(legend.position = "none")

# flip clusters as needed to be more aesthetically pleasing
dendro_plot <- p1 %>%
  flip(91, 92) %>%
  flip(93, 94) %>%
  flip(104, 105)

# use ggtree to extract dendrogram leaf labels  
names <- data.frame(sample_name = get_taxa_name(dendro_plot))

# merge `names` with ggt1$data to join `y` values, which is what we need for ordering subplots in figure
plot_label_data <- ggt1$data %>%
  rename(sample_name = label) %>%
  filter(!is.na(sample_name))

merged_labs <- names %>% 
  left_join(., plot_label_data) %>%
  rownames_to_column("order") %>%
  select(sample_name, order, node)

# manually reorder IPL_class for plot
ipl_order <- c("branched", "unsaturated", "normal", "diacid", "cyclic", "sterol", 
               "alkanone", "other")

##### 2. Heatmap #####
# make heatmap from relative IPL abundances

# convert relative abundance matrix to longform dataframe and merge with metadata
rel_abun_long <- ipl_rel_abun %>%
  pivot_longer(cols = 2:ncol(.), names_to = "ipl", values_to = "rel_abun") %>%
  left_join(., ipl_metadata, "ipl") %>%
  left_join(., sample_metadata, "sample_name") %>%
  left_join(., merged_labs, "sample_name") %>%
  mutate(fraction = factor(fraction, levels = c("glycolipids", "phospholipids", "residual_polar_lipids")))

# bin by ipl_class and fraction for each sample
heatmap_data <- rel_abun_long %>%
  filter(ipl_class %in% c("branched", "unsaturated", "normal", "diacid")) %>%
  group_by(sample_name, fraction, ipl_class) %>%
  summarize(rel_abun_bin_sum = sum(rel_abun)) %>%
  left_join(., merged_labs, "sample_name") %>%
  left_join(., sample_metadata, "sample_name")

# plot heatmap
heatmap <- heatmap_data %>%
  ggplot() +
  geom_tile(aes(x = factor(ipl_class, levels = ipl_order), 
                y = as.numeric(order),
                fill = rel_abun_bin_sum,
                text = paste("<br>",
                             "sample: ", sample_name,
                             "<br>", 
                             "sample type: ", sample_class,
                             "<br>",
                             "IPL yield (% of total extractable IPLs): ", round(rel_abun_bin_sum, 3)*100,
                             "<br>",
                             "IPL class: ", ipl_class))) +
  scale_fill_gradient(low = "#ffffff", high = "#202020") +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(), 
        strip.text = element_blank(),
        legend.position = "none")

##### 3. Areaplot ##### 

# make areaplot of relative abundance bins ordered by hclust() results (from merged_labs)

areaplot_data <- rel_abun_long %>%
  group_by(sample_name, ipl_class) %>%
  mutate(rel_abun_binned = sum(rel_abun)) %>%
  distinct(sample_name, ipl_class, rel_abun_binned, order) %>%
  pivot_wider(names_from = "ipl_class", values_from = "rel_abun_binned") %>% 
  replace(is.na(.), 0) %>%
  pivot_longer(cols = 3:ncol(.), values_to = "rel_abun_binned", names_to = "ipl_class") %>%
  mutate(ipl_class = factor(ipl_class, levels = ipl_order))

ipl_class_colorz <- c("branched" = "#01665e", # brown/green palette
                      "normal" = "#bf812d",
                      "unsaturated" = "#35978f",
                      "diacid" = "#8c510a",
                      "cyclic" = "#000000",
                      "sterol" = "#000000",
                      "alkanone" = "#000000",
                      "other" = "#000000")

areaplot <- areaplot_data %>%
  ggplot(aes(x = as.numeric(order), 
             y = rel_abun_binned, 
             fill = ipl_class)) +
  geom_area(stat = "identity") +
  scale_fill_manual(values = ipl_class_colorz) +
  scale_x_reverse() +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "blank")

##### 4. Lineplot #####
# make additional lineplot showing ipl yields for each fraction

# convert abundance matrix (ug lipid/g TLE) to long-form dataframe
ipl_ug_g_tle_long <- ipl_ug_g_tle %>%
  pivot_longer(cols = 2:ncol(.), names_to = "ipl", values_to = "ipl_ug_g_tle")

lineplot_data_tle <- rel_abun_long %>%
  left_join(ipl_ug_g_tle_long) %>%
  group_by(sample_name, fraction) %>%
  summarize(total_ipl_ug_g_tle = sum(ipl_ug_g_tle)) %>%
  mutate(total_ipl_ug_g_tle_log10 = log10(total_ipl_ug_g_tle + 1),
         total_ipl_ug_g_tle_check = (10^total_ipl_ug_g_tle_log10)-1,
         total_ipl_ug_g_tle_sqrt = sqrt(total_ipl_ug_g_tle)) %>%
  left_join(., merged_labs) %>%
  mutate(order = as.numeric(order)) %>%
  arrange(desc(order))

lineplot_tle <- lineplot_data_tle %>% 
  ggplot() +
  geom_path(aes(x = total_ipl_ug_g_tle/1000, 
                y = order,
                group = fraction,
                color = fraction)) +
  scale_x_log10() +
  scale_y_reverse() +
  theme_bw() +
  theme(legend.position = "blank",
        axis.text = element_blank(),
        axis.title = element_blank())

##### Plot heatmap, lineplot, dendro_plot, and areaplot together for final figure #####

ggarrange(dendro_plot + 
            theme(panel.grid = element_blank()) +
            scale_x_reverse()+ 
            coord_flip(),
          lineplot_tle + 
            theme(panel.grid = element_blank()) +
            coord_flip(),
          areaplot +
            theme(panel.grid = element_blank()),
          heatmap + 
            theme(panel.grid = element_blank()) +
            facet_wrap(~fraction, ncol = 1) +
            coord_flip() +
            scale_y_reverse(),
          align = "v",
          ncol = 1,
          heights = c(1,1,1,3))





