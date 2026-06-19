# ==============================================================================
# Script: 4_Dominance_selection_Figure.R
# Purpose: Generate phylogeny, genetic-distance boxplots, and correlation plots
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggtree)
  library(treeio)
  library(ggtreeExtra)
  library(ggnewscale)
  library(ggpubr)
})

# ---- 1. Settings ----
options(ignore.negative.edge = TRUE)

base_path <- "~/"

metadata_path <- file.path(base_path, "1_Data/0_Metadata/H3N2_777strains_Metadata.csv")
ha_tree_path  <- file.path(base_path, "1_Data/2_Tree/MCC/04_HA_treeanno.tre")

fig2_dir <- file.path(base_path, "4_Figures/2_Figure2")
fig3_dir <- file.path(base_path, "4_Figures/3_Figure3")

dir.create(fig2_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig3_dir, recursive = TRUE, showWarnings = FALSE)

segment_order <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS")

segment_colors <- c(
  "PB2" = "#E41A1C",
  "PB1" = "#FF7F00",
  "PA"  = "#FFD92F",
  "HA"  = "#4DAF4A",
  "NP"  = "#377EB8",
  "NA"  = "#984EA3",
  "M"   = "#A65628",
  "NS"  = "#F781BF"
)

status_colors <- c(
  "Persistent" = "#E41A1C",
  "Non-persistent" = "#4DAF4A"
)

# ---- 2. Load data ----
mcc_tree <- treeio::read.beast(ha_tree_path)
metadata <- read.csv(metadata_path, stringsAsFactors = FALSE)

plot_data <- metadata %>%
  mutate(
    StrainName = as.character(StrainName),
    persistence_status = ifelse(dom_category == "Dominant", "Persistent", "Non-persistent")
  )

tree_metadata <- plot_data %>%
  select(
    label = StrainName,
    persistence_status,
    dom_prob,
    HA_Nonsyn
  ) %>%
  distinct(label, .keep_all = TRUE)

distance_long <- plot_data %>%
  select(
    StrainName,
    persistence_status,
    contains("_Nonsyn"),
    contains("_Syn")
  ) %>%
  pivot_longer(
    cols = where(is.numeric),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  separate(Metric, into = c("Segment", "Type"), sep = "_") %>%
  mutate(
    Segment = factor(Segment, levels = segment_order),
    Type = factor(
      Type,
      levels = c("Nonsyn", "Syn"),
      labels = c("Non-synonymous genetic distance", "Synonymous genetic distance")
    ),
    persistence_status = factor(
      persistence_status,
      levels = c("Persistent", "Non-persistent")
    )
  ) %>%
  filter(!is.na(Value), !is.na(persistence_status))

# ---- 3. Figure 2A: HA phylogeny with HA dN tile ----
p_fig2a_base <- ggtree(mcc_tree, mrsd = "2024-01-01", size = 0.4) %<+% tree_metadata +
  geom_tippoint(
    aes(color = persistence_status),
    size = 1.2,
    alpha = 0.7
  ) +
  scale_color_manual(
    values = status_colors,
    name = "Persistence status"
  ) +
  theme_tree2()

p_fig2a <- p_fig2a_base +
  ggnewscale::new_scale_fill() +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = label, fill = HA_Nonsyn),
    width = 3,
    offset = 0.1,
    paxis.draw = FALSE
  ) +
  scale_fill_viridis_c(
    option = "B",
    name = "HA dN",
    na.value = "grey90"
  ) +
  labs(title = "A. HA phylogeny and persistence patterns") +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold", size = 9)
  )

ggsave(
  filename = "Fig2A_PhylogeneticTree.svg",
  plot = p_fig2a,
  path = fig2_dir,
  width = 400,
  height = 1000,
  units = "mm",
  dpi = 720
)

# ---- 4. Figure 2B: Segment-level dN/dS distribution ----
p_fig3a <- distance_long %>%
  ggplot(aes(x = Segment, y = Value, fill = Segment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
  geom_jitter(
    aes(color = Segment),
    width = 0.2,
    size = 0.4,
    alpha = 0.3
  ) +
  facet_wrap(~ Type, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = segment_colors) +
  scale_color_manual(values = segment_colors) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "white"),
    axis.text.x = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "B. Segment-level genetic distances",
    x = "Gene segment",
    y = "Genetic distance"
  )

ggsave(
  filename = "Fig3A_AllSeg_dNdSBox.svg",
  plot = p_fig3a,
  path = fige_dir,
  width = 700,
  height = 400,
  units = "mm",
  dpi = 720
)

# ---- 5. Figure 3B: Persistence-stratified dN/dS distribution ----
p_fig3b <- distance_long %>%
  ggplot(aes(x = persistence_status, y = Value, fill = persistence_status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 0.3, alpha = 0.2) +
  facet_grid(Type ~ Segment, scales = "free_y") +
  scale_fill_manual(values = status_colors) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 10),
    panel.spacing = unit(0.8, "lines")
  ) +
  labs(
    title = "B. Persistence-stratified genetic distances",
    x = NULL,
    y = "Genetic distance"
  )

ggsave(
  filename = "Fig3B_Dom_dNdSBox.svg",
  plot = p_fig3b,
  path = fig3_dir,
  width = 500,
  height = 400,
  units = "mm",
  dpi = 720
)

# ---- 6. Correlation plot helper ----
make_lower_correlation_plot <- function(data, feature_suffix, title, fill_high = "#084594") {
  selected_cols <- paste0(segment_order, feature_suffix)
  
  cor_data <- data %>%
    select(all_of(selected_cols)) %>%
    rename_with(~ gsub(feature_suffix, "", .x))
  
  cor_matrix <- cor(
    cor_data,
    use = "pairwise.complete.obs",
    method = "pearson"
  )
  
  cor_matrix[upper.tri(cor_matrix, diag = TRUE)] <- NA
  
  cor_long <- as.data.frame(as.table(cor_matrix)) %>%
    rename(Var1 = Var1, Var2 = Var2, R = Freq) %>%
    filter(!is.na(R)) %>%
    mutate(
      Var1 = factor(Var1, levels = segment_order),
      Var2 = factor(Var2, levels = segment_order),
      Text_Color = ifelse(R > 0.75, "white", "black")
    )
  
  ggplot(cor_long, aes(x = Var1, y = Var2)) +
    geom_point(
      aes(fill = R, size = R),
      shape = 21,
      color = "grey30",
      stroke = 0.6
    ) +
    geom_text(
      aes(label = sprintf("%.2f", R), color = Text_Color),
      size = 2.8,
      fontface = "bold"
    ) +
    scale_color_identity() +
    scale_size_continuous(range = c(6, 14), guide = "none") +
    scale_fill_gradient(
      low = "#FFFFFF",
      high = fill_high,
      limits = c(0, 1),
      name = "Pearson's R"
    ) +
    scale_y_discrete(limits = rev(segment_order)) +
    coord_fixed() +
    theme_minimal(base_size = 11) +
    theme(
      axis.title = element_blank(),
      axis.text = element_text(face = "bold", size = 10, color = "black"),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      legend.position = c(1.07, 0.8),
      legend.background = element_rect(fill = "white", color = NA)
    ) +
    labs(title = title)
}

# ---- 7. Figure 3C,D: dN dS correlation ----

p_fig3c <- make_lower_correlation_plot(
  data = metadata,
  feature_suffix = "_Nonsyn",
  title = "C. Genome-wide correlation of non-synonymous distances"
)

ggsave(
  filename = "Fig3C_dN_circular_correlation.svg",
  plot = p_fig3c,
  path = fig3_dir,
  width = 450,
  height = 400,
  units = "mm",
  dpi = 720
)


p_fig3d <- make_lower_correlation_plot(
  data = metadata,
  feature_suffix = "_Syn",
  title = "D. Genome-wide correlation of synonymous distances"
)

ggsave(
  filename = "Fig3D_dS_circular_correlation.svg",
  plot = p_fig3D,
  path = fig3_dir,
  width = 450,
  height = 400,
  units = "mm",
  dpi = 720
)
