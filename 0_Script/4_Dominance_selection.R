# ==============================================================================
# Script Name: 4_Dominance_selection.R
# Purpose: Quantify the evolutionary "dominance" of H3N2 strains by analyzing 
#          ancestral relationships with future vaccine strains across 
#          posterior tree distributions.
# ==============================================================================

library(tidyverse)
library(ape)
library(treeio)
library(ggtree)

# ---- 1. Environment Setup & Data Loading ----
# Paths for metadata and pre-processed posterior trees
meta_path <- "~/1_Data/10_Final_MetaData_Submission/H3N2_Testdata_777strains_1_Metadata.csv"
tree_path <- "~/1_Data/2_RNSC/2_Treelog/5_3rd_Verification_777/HA_resample.rds"

full_df <- read.csv(meta_path)
trees   <- read_rds(tree_path)

# Mapping vaccine timeline to numeric indices for successor tracking
v_timeline <- c("Mos99", "Fuj02", "Cal04", "Wis05", "Bris07", "Prth09", "Vic11", 
                "Swtz13", "HK15", "Kan17", "HK19", "Cam20", "Dar21", "TH22", "CR23")
v_map      <- setNames(match(full_df$vaccine_code, v_timeline), full_df$StrainName)

# ---- 2. Core Function: Standardized Dominance Calculation ----
# Calculates dominance scores based on the presence of future successors 
# in sibling lineages and filters out phylogenetic outliers.
calculate_dominance_std <- function(phy, v_map, n_tips) {
  descendants <- prop.part(phy)
  edges <- as.data.frame(phy$edge)
  colnames(edges) <- c("parent", "node")
  
  # Extract branch lengths for outlier detection
  branch_lens <- phy$edge.length[match(1:n_tips, phy$edge[,2])]
  v_indices   <- v_map[phy$tip.label]
  
  scores <- setNames(numeric(n_tips), phy$tip.label)
  
  for (i in 1:n_tips) {
    p_node <- edges$parent[edges$node == i]
    if (length(p_node) == 0) next
    
    # Identify sibling nodes to compare evolutionary trajectories
    siblings <- edges$node[edges$parent == p_node & edges$node != i]
    if (length(siblings) == 0) next
    
    # Get descendant tips of the sibling line
    sib_tips <- unlist(lapply(siblings, function(s) {
      if(s <= n_tips) s else descendants[[s - n_tips]]
    }))
    
    # SCORE: Count how many 'successor' vaccine strains exist in sibling lineages
    succ_count <- sum(v_indices[sib_tips] > v_indices[i], na.rm = TRUE)
    
    # OUTLIER DETECTION: Filter out long branches (potential sequencing errors/dead-ends)
    # Using local average branch length as a baseline for standardization
    sib_avg_len <- mean(branch_lens[sib_tips], na.rm = TRUE)
    
    if(is.na(sib_avg_len) || sib_avg_len == 0) {
      is_outlier <- FALSE
    } else {
      # Standardized threshold (2.5x local average)
      is_outlier <- (branch_lens[i] > sib_avg_len * 2.5)
    }
    
    # Assign score if the strain led to future success and is not an outlier
    if (!is.na(succ_count) && succ_count >= 1 && !is_outlier) {
      # Cap the impact score to prevent bias from large clusters
      scores[i] <- min(succ_count, 5) 
    }
  }
  return(scores)
}

# ---- 3. Multi-tree Analysis Execution ----
n_trees <- length(trees)
target_names <- full_df$StrainName
n_tips <- length(target_names)
impact_matrix <- matrix(0, nrow = n_tips, ncol = n_trees, dimnames = list(target_names, NULL))

cat("Starting multi-tree dominance analysis with standardized logic...\n")
for (t in 1:n_trees) {
  curr_phy <- as.phylo(trees[[t]])
  impact_matrix[, t] <- calculate_dominance_std(curr_phy, v_map, n_tips)[target_names]
  if (t %% 50 == 0) cat(sprintf("\rProgress: %d/%d", t, n_trees))
}
cat("\nAnalysis completed.\n")

# ---- 4. Result Integration & Period-specific Normalization ----
# Apply quantile-based categorization to normalize differences across periods
full_df <- full_df %>%
  mutate(
    dom_prob   = rowMeans(impact_matrix > 0, na.rm = TRUE),
    dom_impact = rowMeans(impact_matrix, na.rm = TRUE)
  ) %>%
  group_by(vaccine_code) %>%
  mutate(
    # Categorize based on relative position within each period's distribution
    dom_category = case_when(
      dom_prob == 0 ~ "Extinct",
      dom_prob >= quantile(dom_prob[dom_prob > 0], 0.7, na.rm = TRUE) & dom_prob >= 0.7 ~ "Dominant",
      dom_prob >= quantile(dom_prob[dom_prob > 0], 0.3, na.rm = TRUE) & dom_prob >= 0.3 ~ "Emerging",
      TRUE ~ "Sporadic"
    ),
    dom_binary = as.factor(ifelse(dom_prob >= 0.5, 1, 0))
  ) %>%
  ungroup()

# Display distribution of dominance categories per period
print(table(full_df$vaccine_code, full_df$dom_category))

# ---- 5. Visualization & Export ----
viz_tree <- as.phylo(trees[[1]])
plot_data <- full_df %>% select(StrainName, dom_prob, dom_category)
row.names(plot_data) <- plot_data$StrainName

p <- ggtree(viz_tree) %<+% plot_data +
  geom_tippoint(aes(color = dom_category, size = dom_prob), alpha = 0.7) +
  scale_color_manual(values = c(
    "Dominant" = "#E41A1C", 
    "Emerging" = "#FF7F00", 
    "Sporadic" = "#4DAF4A", 
    "Extinct"  = "#999999"  
  )) +
  theme_tree2() +
  labs(title = "Phylogenetic Dominance Analysis (H3N2)",
       subtitle = "Standardized by period-specific distributions",
       color = "Category", size = "Dominance Prob")

print(p)

# Save the enriched metadata for Machine Learning steps
write.csv(full_df, '~/1_Data/10_Final_MetaData_Submission/H3N2_Testdata_777strains_1_Metadata_Dominance.csv', row.names = FALSE)