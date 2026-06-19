# ==============================================================================
# Script Name: 7_Coevolution_Analysis_TreeKnit.R
# Purpose: Visualize genome-wide co-evolutionary landscapes and reassortment  
#          dynamics using lower-triangle heatmaps and fragmentation profiles.
# ==============================================================================

library(ggplot2)
library(tidyr)
library(dplyr)
library(viridis)
library(ggrepel) 
library(patchwork)

# ---- 1. Data Loading & Pre-processing ----
# Load multi-scale evolutionary summary generated from Julia TreeKnit analysis
base_path <- '~'
dat_path  <- '/3_Result/5_TreeKnit/H3N2_Evolutionary_Multi_Scale_Summary.csv'
df_path   <- file.path(base_path, dat_path)

if (!file.exists(df_path)) {
  stop(paste("Target empirical CSV not found at:", df_path))
}

df <- read.csv(df_path)

# Define standard chronological segment order for manuscript consistency
segment_names <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS")
rev_segment_names <- rev(segment_names)

# Tidy raw pair indices into formal segment names and factors
df_tidy <- df %>%
  separate(Pair, into = c("Seg_A_idx", "Seg_B_idx"), sep = "-", remove = FALSE) %>%
  mutate(
    idx_a = as.numeric(gsub("Seg", "", Seg_A_idx)),
    idx_b = as.numeric(gsub("Seg", "", Seg_B_idx)),
    Gene_A = factor(segment_names[idx_a], levels = segment_names),
    Gene_B = factor(segment_names[idx_b], levels = segment_names),
    Gene_Pair = paste0(segment_names[idx_a], "-", segment_names[idx_b])
  )

# ---- 2. Figure 3A: Upper-Triangle Co-evolutionary Interaction Map ----

# [BUILT-IN FIX] Direct data-frame logical filtering instead of matrix conversion
# This guarantees that Seg_X always represents the lower index and Seg_Y the higher index,
# forcing the data into a perfect geometric upper triangle without losing row metadata.

df_triangle <- df_tidy %>%
  # Select only unique pairs where the numeric index of A is strictly less than B
  filter(idx_a < idx_b) %>%
  select(Seg_X = Gene_A, Seg_Y = Gene_B, Value = Size_over_1) %>%
  mutate(
    Seg_X = factor(Seg_X, levels = segment_names),
    # Invert the Y-axis levels so PB2 stays at the top and NS at the bottom
    Seg_Y = factor(Seg_Y, levels = rev_segment_names)
  )

# Verify data presence before plotting (Should print rows to your console)
print("Data check for Upper Triangle plot:")
print(head(df_triangle))

# Generate Panel A Heatmap (Upper Triangle)
p3a_heatmap <- ggplot(df_triangle, aes(x = Seg_X, y = Seg_Y, fill = Value)) +
  geom_tile(color = "white", size = 0.4) +
  # Dynamic text contrast: black labels on light yellow/orange tiles, white on dark purple tiles
  geom_text(aes(label = Value), 
            color = ifelse(df_triangle$Value > 75, "black", "white"), 
            size = 3.5, fontface = "bold") +
  scale_fill_viridis_c(
    option = "magma", 
    name = "Stable MCCs\n(Size > 3)",
    limits = c(min(df_tidy$Size_over_1), max(df_tidy$Size_over_1))
  ) +
  labs(
    title = "A. Co-evolutionary Interaction Map",
    subtitle = "Analysis of stable MCC groups (Size > 3) to filter stochastic noise",
    x = "Genomic Segment A", y = "Genomic Segment B"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
    axis.text.y = element_text(face = "bold", color = "black"),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 13),
    # Anchor the legend in the open, empty bottom-left space
    legend.position = c(0.8, 0.7), 
    legend.background = element_rect(fill = "white", color = "grey90", size = 0.3)
  ) +
  coord_fixed()

# Force print to graphics device
print(p3a_heatmap)


# ---- 3. Figure 3B: Fragmentation vs. Linkage Stability Profile ----
# Re-align factor levels of the scatter plot for flawless multi-panel integration
df_scatter_data <- df_tidy %>%
  mutate(Gene_Pair = factor(Gene_Pair))

# Generate Panel B Scatter Plot
p3b_scatter <- ggplot(df_scatter_data, aes(x = Frag_Rate, y = Size_over_1, label = Gene_Pair)) +
  geom_point(aes(size = Max_Clade_Size, color = Size_over_1), alpha = 0.75) +
  geom_text_repel(
    size = 3.2, 
    fontface = "bold", 
    box.padding = 0.4,
    point.padding = 0.3, 
    max.overlaps = Inf, 
    segment.color = 'grey50', 
    force = 2.5
  ) + 
  scale_color_viridis_c(option = "magma", name = "Stable Group Count") +
  labs(
    title = "B. Genomic Stability vs. Fragmentation Profile",
    subtitle = "Reassortment rates (fragmentation) vs. legacy block linkage stability",
    x = "Fragmentation Rate (Total MCCs / Total Strains)",
    y = "Number of Stable Clades (Size > 3)",
    size = "Max Clade Size"
  ) +
  theme_bw(base_size = 11) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15))) +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.08))) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    strip.background = element_rect(fill = "white"),
    legend.title = element_text(face = "bold", size = 9)
  )
p3b_scatter

# Save high-resolution vectors for manuscript submission

ggsave("Fig3A_TreeKnit_Heatmap.svg", plot = p3a_heatmap, device = "svg",
       path = paste0(base_path, "/4_Figures/3_Figure3"),
       width = 230, height = 250, units = "mm", dpi = 720)

ggsave("Fig3B_TreeKnit_Scatter.svg", plot = p3b_scatter, device = "svg",
       path = paste0(base_path, "/4_Figures/3_Figure3"),
       width = 230, height = 250, units = "mm", dpi = 720)

