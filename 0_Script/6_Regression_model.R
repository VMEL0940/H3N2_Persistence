# ==============================================================================
# Script Name: 9_Window_Sensitivity_Regression_CV_v2_scheme_level.R
# Purpose:
#   Rolling-window-only comparison of temporal CV designs for H3N2 dominance prediction.
#
# Key question answered by this version:
#   NOT "which window performs best only when the test period is TH22?"
#   BUT:
#     1) For rolling width = 8, average over ALL possible rolling test periods
#        e.g., 1-8 -> test 9, 2-9 -> test 10, ..., 6-13 -> test TH22.
#     2) Repeat for multiple rolling widths, then rank rolling schemes.
#     3) Compare the selected representative model with HA-only, all-dN, and all-dN+all-dS controls.
#
# Important design choices:
#   - dS/synonymous terms are used only in Control 2, not in representative model selection.
#   - Final scheme ranking is based on representative dN-only candidate models only.
#   - Train-fitted scaling is used within each CV split to prevent leakage.
#   - Outputs include both:
#       a) fixed-model-per-scheme summary: fair final scheme/model selection
#       b) oracle-per-test summary: diagnostic upper bound only
# ===============================================================================

suppressPackageStartupMessages({
  library(pROC)
  library(lme4)
  library(tidyverse)
  library(lubridate)
})

# ---- 0. User settings ----
base_path <- "~/"
file_name <- "1_Data/0_Metadata/H3N2_777strains_Metadata.csv"
data_path <- file.path(base_path, file_name)

# Compare these window sizes.
# rolling_width_8 means 1-8, 2-9, 3-10, ...
window_sizes_to_test <- 6:10
last_test_period <- "TH22"   # CR23 is kept for final deployment, not CV model selection.

# Selection score penalty. You can tune this, but keep fixed for all schemes.
# Representative score = AUC - complexity_penalty_per_term * number_of_fixed_terms.
complexity_penalty_per_term <- 0.005

# ---- 1. Load data ----
if (!file.exists(data_path)) stop(paste("Target dataset metadata not found at:", data_path))
full_df <- read.csv(data_path, stringsAsFactors = FALSE)

base_levels  <- c("Mos99", "Fuj02", "Cal04", "Wis05", "Bris07", "Prth09", "Vic11", "Swtz13", "HK15", "Kan17", "HK19", "Cam20", "Dar21", "TH22", "CR23")
unique_codes <- unique(full_df$vaccine_code)
all_periods  <- unique(c(base_levels, unique_codes))

if (!(last_test_period %in% all_periods)) stop(paste("last_test_period not found:", last_test_period))
last_test_idx <- which(all_periods == last_test_period)

raw_df <- full_df %>%
  mutate(
    date = as_date(date),
    decimal_year = decimal_date(date),
    Period = factor(vaccine_code, levels = all_periods)
  )

# ---- 2. Leakage-safe train-fitted scaling ----
dn_raw <- c("PB2_Nonsyn", "PB1_Nonsyn", "PA_Nonsyn", "HA_Nonsyn", "NP_Nonsyn", "NA_Nonsyn", "M_Nonsyn", "NS_Nonsyn")
ds_raw <- c("PB2_Syn",    "PB1_Syn",    "PA_Syn",    "HA_Syn",    "NP_Syn",    "NA_Syn",    "M_Syn",    "NS_Syn")
all_raw_features <- c(dn_raw, ds_raw)

fit_scaler <- function(df, cols) {
  stats <- lapply(cols, function(col) {
    x <- df[[col]]
    mu <- mean(x, na.rm = TRUE)
    sig <- sd(x, na.rm = TRUE)
    if (is.na(mu)) mu <- 0
    if (is.na(sig) || sig == 0) sig <- 1
    list(mean = mu, sd = sig)
  })
  names(stats) <- cols
  stats
}

apply_scaler <- function(df, scaler) {
  for (col in names(scaler)) {
    out_col <- paste0(col, "_std")
    df[[out_col]] <- (df[[col]] - scaler[[col]]$mean) / scaler[[col]]$sd
    df[[out_col]][is.na(df[[out_col]])] <- 0
  }

  df %>%
    transmute(
      dom_binary, dom_prob, dom_category, StrainName, vaccine_code, decimal_year, Period,
      PB2 = PB2_Nonsyn_std, PB1 = PB1_Nonsyn_std, PA = PA_Nonsyn_std,
      HA = HA_Nonsyn_std, NP = NP_Nonsyn_std, NA_seg = NA_Nonsyn_std,
      M = M_Nonsyn_std, NS = NS_Nonsyn_std,
      PB2_syn = PB2_Syn_std, PB1_syn = PB1_Syn_std, PA_syn = PA_Syn_std,
      HA_syn = HA_Syn_std, NP_syn = NP_Syn_std, NA_syn = NA_Syn_std,
      M_syn = M_Syn_std, NS_syn = NS_Syn_std
    )
}

recalc_weights <- function(df) {
  zeros <- sum(df$dom_binary == 0, na.rm = TRUE)
  ones  <- sum(df$dom_binary == 1, na.rm = TRUE)
  if (ones <= 0 || zeros <= 0) return(mutate(df, norm_wts = 1))
  df %>%
    mutate(
      norm_wts = ifelse(dom_binary == 1, zeros / ones, 1),
      norm_wts = norm_wts / mean(norm_wts, na.rm = TRUE)
    )
}
# ---- 3. Formula banks ----
# Model design in this version:
#   1) Representative model candidates: nonsynonymous dN only.
#      - HA is always included as the anchor.
#      - Other dN segments are selected by exhaustive subset search, not pre-fixed.
#      - Optional HA x partner interaction candidates are generated automatically.
#   2) HA dN only model.
#   3) Control 1: all dN segments, PB2~NS.
#   4) Control 2: all dN + all dS segments, PB2~NS.
#
# Important: dS terms are NOT included in representative model selection.

dn_genes_all <- c("PB2", "PB1", "PA", "HA", "NP", "NA_seg", "M", "NS")
ds_genes_all <- c("PB2_syn", "PB1_syn", "PA_syn", "HA_syn", "NP_syn", "NA_syn", "M_syn", "NS_syn")
dn_partners <- setdiff(dn_genes_all, "HA")

# Tune these if the exhaustive model bank is too large.
# max_representative_partners = 7 means all possible dN subset sizes are tested.
max_representative_partners <- length(dn_partners)
representative_include_interactions <- TRUE

clean_gene_name <- function(x) {
  x %>%
    gsub("_seg", "", .) %>%
    gsub("_syn", "", .)
}

make_representative_dn_bank <- function(anchor = "HA",
                                        partners = dn_partners,
                                        prefix = "Representative_dN",
                                        max_partners = length(partners),
                                        include_interactions = TRUE) {
  comb_list <- list()
  for (i in 1:min(length(partners), max_partners)) {
    comb_list <- c(comb_list, combn(partners, i, simplify = FALSE))
  }

  forms <- list()

  for (g in comb_list) {
    g <- unlist(g)

    # Additive dN subset model: HA + selected nonsynonymous partner segments.
    ftxt <- paste0(
      "dom_binary ~ ",
      paste(c(anchor, g), collapse = " + "),
      " + (1 | vaccine_code)"
    )

    model_name <- paste0(
      prefix, "_HA_",
      paste(clean_gene_name(g), collapse = "_")
    )
    forms[[model_name]] <- ftxt

    # Optional interaction model: HA x one selected partner + remaining selected dN terms.
    # This keeps interaction search unbiased across all partner genes rather than fixing NA/PB1/PA/NS.
    if (include_interactions && length(g) >= 1) {
      for (ix in g) {
        additive_rest <- setdiff(g, ix)
        rhs_terms_int <- c(paste0(anchor, " * ", ix), additive_rest)
        ftxt_int <- paste0(
          "dom_binary ~ ",
          paste(rhs_terms_int, collapse = " + "),
          " + (1 | vaccine_code)"
        )

        model_name_int <- paste0(
          prefix, "_HA_x_", clean_gene_name(ix),
          if (length(additive_rest) > 0) {
            paste0("_plus_", paste(clean_gene_name(additive_rest), collapse = "_"))
          } else {
            ""
          }
        )
        forms[[model_name_int]] <- ftxt_int
      }
    }
  }

  forms
}

representative_dn_candidates <- make_representative_dn_bank(
  anchor = "HA",
  partners = dn_partners,
  prefix = "Representative_dN",
  max_partners = max_representative_partners,
  include_interactions = representative_include_interactions
)

ha_dn_only_model <- list(
  "HA_dN_Only" = "dom_binary ~ HA + (1 | vaccine_code)"
)

dn_all_control <- list(
  "Control1_All_dN_PB2_to_NS" = paste0(
    "dom_binary ~ ",
    paste(dn_genes_all, collapse = " + "),
    " + (1 | vaccine_code)"
  )
)

dnds_all_control <- list(
  "Control2_All_dN_All_dS_PB2_to_NS" = paste0(
    "dom_binary ~ ",
    paste(c(dn_genes_all, ds_genes_all), collapse = " + "),
    " + (1 | vaccine_code)"
  )
)

formula_bank <- c(
  representative_dn_candidates,
  ha_dn_only_model,
  dn_all_control,
  dnds_all_control
)

model_class <- tibble(
  Model_Name = names(formula_bank),
  Class = case_when(
    grepl("^Representative_dN", Model_Name) ~ "representative_dN_candidate",
    Model_Name == "HA_dN_Only" ~ "ha_dN_only",
    Model_Name == "Control1_All_dN_PB2_to_NS" ~ "control1_all_dN",
    Model_Name == "Control2_All_dN_All_dS_PB2_to_NS" ~ "control2_all_dN_all_dS",
    TRUE ~ "other"
  )
)

count_fixed_terms <- function(formula_text) {
  rhs <- strsplit(formula_text, "~")[[1]][2]
  rhs <- gsub("\\(1 \\| vaccine_code\\)", "", rhs)
  rhs <- gsub("\\s+", "", rhs)
  if (is.na(rhs) || rhs == "") return(0)
  # Treat A*B as A + B + A:B = 3 fixed terms.
  rhs <- gsub("\\*", "+INTERACTION+", rhs)
  terms <- unlist(strsplit(rhs, "\\+"))
  terms <- terms[terms != ""]
  length(terms)
}

# ---- 4. Rolling-only window grid ----
make_windows <- function(all_periods, last_test_period, window_sizes) {
  max_idx <- which(all_periods == last_test_period)
  windows <- list()

  for (w in window_sizes) {
    if (w >= max_idx) next

    # Rolling fixed width only: 1-w -> test w+1, 2-(w+1) -> test w+2, ...
    for (start_i in 1:(max_idx - w)) {
      end_i <- start_i + w - 1
      windows[[length(windows) + 1]] <- tibble(
        Window_Scheme = paste0("rolling_width_", w),
        Window_Type = "rolling",
        Window_Size = w,
        Train_Start_Index = start_i,
        Train_End_Index = end_i,
        Test_Index = end_i + 1,
        Train_Periods = paste(all_periods[start_i:end_i], collapse = ","),
        Test_Period = all_periods[end_i + 1]
      )
    }
  }

  bind_rows(windows) %>%
    mutate(
      Test_Period = as.character(Test_Period),
      Train_Range_Label = paste0(all_periods[Train_Start_Index], "-", all_periods[Train_End_Index], " -> ", Test_Period)
    ) %>%
    filter(Test_Period != "CR23")
}

window_grid <- make_windows(all_periods, last_test_period, window_sizes_to_test)

# ---- 5. Fit/evaluate one train/test split ----
fit_one_window <- function(win_row) {
  train_periods <- strsplit(win_row$Train_Periods, ",")[[1]]
  test_period <- win_row$Test_Period

  raw_train <- raw_df %>% filter(as.character(Period) %in% train_periods)
  raw_test  <- raw_df %>% filter(as.character(Period) == test_period)

  if (nrow(raw_train) < 10 || nrow(raw_test) < 5) return(NULL)
  if (length(unique(raw_train$dom_binary)) < 2 || length(unique(raw_test$dom_binary)) < 2) return(NULL)

  scaler <- fit_scaler(raw_train, all_raw_features)
  train_df <- apply_scaler(raw_train, scaler) %>% recalc_weights()
  test_df  <- apply_scaler(raw_test, scaler)

  out <- vector("list", length(formula_bank))
  names(out) <- names(formula_bank)

  for (m_name in names(formula_bank)) {
    ftxt <- formula_bank[[m_name]]
    out[[m_name]] <- tryCatch({
      fit <- suppressWarnings(glmer(
        as.formula(ftxt), data = train_df, family = binomial, weights = norm_wts,
        nAGQ = 0,
        control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
      ))

      # New test period should be predicted using fixed effects only.
      probs <- predict(fit, newdata = test_df, type = "response", re.form = ~0, allow.new.levels = TRUE)
      auc_test <- as.numeric(auc(roc(test_df$dom_binary, probs, quiet = TRUE)))

      tibble(
        Window_Scheme = win_row$Window_Scheme,
        Window_Type = win_row$Window_Type,
        Window_Size = win_row$Window_Size,
        Train_Start_Index = win_row$Train_Start_Index,
        Train_End_Index = win_row$Train_End_Index,
        Train_Periods = win_row$Train_Periods,
        Train_Range_Label = win_row$Train_Range_Label,
        Test_Index = win_row$Test_Index,
        Test_Period = test_period,
        Model_Name = m_name,
        Formula = ftxt,
        AUC = auc_test,
        AIC = AIC(fit),
        N_Fixed_Terms = count_fixed_terms(ftxt),
        Error = NA_character_
      )
    }, error = function(e) {
      tibble(
        Window_Scheme = win_row$Window_Scheme,
        Window_Type = win_row$Window_Type,
        Window_Size = win_row$Window_Size,
        Train_Start_Index = win_row$Train_Start_Index,
        Train_End_Index = win_row$Train_End_Index,
        Train_Periods = win_row$Train_Periods,
        Train_Range_Label = win_row$Train_Range_Label,
        Test_Index = win_row$Test_Index,
        Test_Period = test_period,
        Model_Name = m_name,
        Formula = ftxt,
        AUC = NA_real_,
        AIC = NA_real_,
        N_Fixed_Terms = count_fixed_terms(ftxt),
        Error = e$message
      )
    })
  }

  bind_rows(out)
}

# ---- 6. Run all windows/models ----
message("Total train/test splits to evaluate: ", nrow(window_grid))
message("Total models per split: ", length(formula_bank))

all_results <- pmap_dfr(window_grid, function(...) fit_one_window(tibble(...))) %>%
  left_join(model_class, by = "Model_Name")

# ---- 7. Score models ----
# Representative model selection is based on dN-only model performance.
# dS is NOT used as an adjustment term or selection ceiling.
# Control 1 and Control 2 remain in the result tables for comparison only.

control_comparison <- all_results %>%
  filter(Class %in% c("ha_dN_only", "control1_all_dN", "control2_all_dN_all_dS"), !is.na(AUC)) %>%
  select(Window_Scheme, Window_Type, Window_Size, Test_Period, Test_Index, Class, AUC) %>%
  pivot_wider(names_from = Class, values_from = AUC, names_prefix = "AUC_")

scored_results <- all_results %>%
  left_join(
    control_comparison,
    by = c("Window_Scheme", "Window_Type", "Window_Size", "Test_Period", "Test_Index")
  ) %>%
  mutate(
    Delta_vs_HA_dN_Only = AUC - AUC_ha_dN_only,
    Delta_vs_Control1_All_dN = AUC - AUC_control1_all_dN,
    Delta_vs_Control2_All_dN_All_dS = AUC - AUC_control2_all_dN_all_dS,
    Row_Selection_Score = AUC - complexity_penalty_per_term * N_Fixed_Terms,
    Control2_dS_warning = AUC_control2_all_dN_all_dS >= 0.70
  )

# ---- 8A. Fair final model selection: one fixed best biological model per scheme ----
# This answers: "For rolling_width_8 as a whole, what is its average performance across
# all possible test periods, using one selected biological model?"
model_by_scheme_summary <- scored_results %>%
  filter(Class == "representative_dN_candidate", !is.na(AUC)) %>%
  group_by(Window_Scheme, Window_Type, Window_Size, Model_Name, Class, Formula) %>%
  summarise(
    Avg_AUC_All_Tests = mean(AUC, na.rm = TRUE),
    Median_AUC_All_Tests = median(AUC, na.rm = TRUE),
    Min_AUC_All_Tests = min(AUC, na.rm = TRUE),
    SD_AUC_All_Tests = sd(AUC, na.rm = TRUE),
    Avg_Delta_vs_HA_dN_Only_All_Tests = mean(Delta_vs_HA_dN_Only, na.rm = TRUE),
    Median_Delta_vs_HA_dN_Only_All_Tests = median(Delta_vs_HA_dN_Only, na.rm = TRUE),
    Avg_Delta_vs_Control1_All_dN_All_Tests = mean(Delta_vs_Control1_All_dN, na.rm = TRUE),
    Avg_Delta_vs_Control2_All_dN_All_dS_All_Tests = mean(Delta_vs_Control2_All_dN_All_dS, na.rm = TRUE),
    Avg_Row_Selection_Score = mean(Row_Selection_Score, na.rm = TRUE),
    N_Test_Periods = n_distinct(Test_Period),
    Test_Periods = paste(unique(Test_Period[order(Test_Index)]), collapse = ","),
    N_Fixed_Terms = first(N_Fixed_Terms),
    .groups = "drop"
  ) %>%
  mutate(
    Fixed_Model_Selection_Score = Avg_AUC_All_Tests -
      complexity_penalty_per_term * N_Fixed_Terms
  )

best_fixed_model_per_scheme <- model_by_scheme_summary %>%
  group_by(Window_Scheme, Window_Type, Window_Size) %>%
  slice_max(Fixed_Model_Selection_Score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(Fixed_Model_Selection_Score))

fixed_model_per_test <- scored_results %>%
  inner_join(
    best_fixed_model_per_scheme %>% select(Window_Scheme, Best_Fixed_Model = Model_Name),
    by = "Window_Scheme"
  ) %>%
  filter(Model_Name == Best_Fixed_Model) %>%
  arrange(Window_Size, Window_Type, Test_Index)

scheme_ranking_all_available_tests <- best_fixed_model_per_scheme %>%
  transmute(
    Window_Scheme, Window_Type, Window_Size,
    Best_Fixed_Biological_Model = Model_Name,
    Class, Formula,
    Avg_AUC_All_Available_Tests = Avg_AUC_All_Tests,
    Median_AUC_All_Available_Tests = Median_AUC_All_Tests,
    Min_AUC_All_Available_Tests = Min_AUC_All_Tests,
    SD_AUC_All_Available_Tests = SD_AUC_All_Tests,
    Avg_Delta_vs_HA_dN_Only_All_Available_Tests = Avg_Delta_vs_HA_dN_Only_All_Tests,
    Avg_Delta_vs_Control1_All_dN_All_Available_Tests = Avg_Delta_vs_Control1_All_dN_All_Tests,
    Avg_Delta_vs_Control2_All_dN_All_dS_All_Available_Tests = Avg_Delta_vs_Control2_All_dN_All_dS_All_Tests,
    Fixed_Model_Selection_Score,
    N_Test_Periods,
    Test_Periods,
    N_Fixed_Terms
  ) %>%
  arrange(desc(Fixed_Model_Selection_Score))


# ---- 8A-2. Compare the selected representative model against the three fixed controls ----
selected_representative_vs_controls_by_scheme <- scored_results %>%
  left_join(
    best_fixed_model_per_scheme %>%
      select(Window_Scheme, Selected_Representative_Model = Model_Name),
    by = "Window_Scheme"
  ) %>%
  mutate(
    Comparison_Group = case_when(
      Class == "representative_dN_candidate" & Model_Name == Selected_Representative_Model ~ "1_Selected_representative_dN_subset",
      Class == "ha_dN_only" ~ "2_HA_dN_only",
      Class == "control1_all_dN" ~ "3_Control1_all_dN",
      Class == "control2_all_dN_all_dS" ~ "4_Control2_all_dN_all_dS",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Comparison_Group), !is.na(AUC)) %>%
  group_by(Window_Scheme, Window_Type, Window_Size, Comparison_Group, Model_Name, Class, Formula) %>%
  summarise(
    Avg_AUC_All_Tests = mean(AUC, na.rm = TRUE),
    Median_AUC_All_Tests = median(AUC, na.rm = TRUE),
    Min_AUC_All_Tests = min(AUC, na.rm = TRUE),
    SD_AUC_All_Tests = sd(AUC, na.rm = TRUE),
    Avg_Delta_vs_HA_dN_Only = mean(Delta_vs_HA_dN_Only, na.rm = TRUE),
    Avg_Delta_vs_Control1_All_dN = mean(Delta_vs_Control1_All_dN, na.rm = TRUE),
    Avg_Delta_vs_Control2_All_dN_All_dS = mean(Delta_vs_Control2_All_dN_All_dS, na.rm = TRUE),
    Avg_Selection_Score = mean(Row_Selection_Score, na.rm = TRUE),
    N_Test_Periods = n_distinct(Test_Period),
    Test_Periods = paste(unique(Test_Period[order(Test_Index)]), collapse = ","),
    N_Fixed_Terms = first(N_Fixed_Terms),
    .groups = "drop"
  ) %>%
  arrange(Window_Size, Window_Type, Comparison_Group)

# ---- 8A-3. Rolling model ranking: show all representative models better than HA dN only ----
# Table 5A now includes ALL representative dN models whose average AUC exceeds HA dN only
# within the same rolling scheme.

# ---- 8A-3. Table 5A: HA dN only + all representative models better than HA dN only ----
# This table is model-level ranking under each rolling width.
# It includes:
#   1) HA_dN_Only baseline
#   2) all representative dN-subset models whose average AUC is higher than HA_dN_Only

ha_only_by_scheme_all_available <- scored_results %>%
  filter(Class == "ha_dN_only", !is.na(AUC)) %>%
  group_by(Window_Scheme, Window_Type, Window_Size, Model_Name, Class, Formula) %>%
  summarise(
    Avg_AUC_All_Available_Tests = mean(AUC, na.rm = TRUE),
    Median_AUC_All_Available_Tests = median(AUC, na.rm = TRUE),
    Min_AUC_All_Available_Tests = min(AUC, na.rm = TRUE),
    SD_AUC_All_Available_Tests = sd(AUC, na.rm = TRUE),
    Avg_Delta_vs_HA_dN_Only_All_Available_Tests = 0,
    Median_Delta_vs_HA_dN_Only_All_Available_Tests = 0,
    Avg_Delta_vs_Control1_All_dN_All_Available_Tests = mean(Delta_vs_Control1_All_dN, na.rm = TRUE),
    Avg_Delta_vs_Control2_All_dN_All_dS_All_Available_Tests = mean(Delta_vs_Control2_All_dN_All_dS, na.rm = TRUE),
    Fixed_Model_Selection_Score = mean(Row_Selection_Score, na.rm = TRUE),
    N_Test_Periods = n_distinct(Test_Period),
    Test_Periods = paste(unique(Test_Period[order(Test_Index)]), collapse = ","),
    N_Fixed_Terms = first(N_Fixed_Terms),
    .groups = "drop"
  ) %>%
  mutate(
    Ranking_Group = "0_HA_dN_Only_baseline"
  )

representative_better_than_HA_all_available <- model_by_scheme_summary %>%
  filter(
    Class == "representative_dN_candidate",
    !is.na(Avg_AUC_All_Tests),
    !is.na(Avg_Delta_vs_HA_dN_Only_All_Tests),
    Avg_Delta_vs_HA_dN_Only_All_Tests > 0
  ) %>%
  transmute(
    Window_Scheme,
    Window_Type,
    Window_Size,
    Model_Name,
    Class,
    Formula,
    Avg_AUC_All_Available_Tests = Avg_AUC_All_Tests,
    Median_AUC_All_Available_Tests = Median_AUC_All_Tests,
    Min_AUC_All_Available_Tests = Min_AUC_All_Tests,
    SD_AUC_All_Available_Tests = SD_AUC_All_Tests,
    Avg_Delta_vs_HA_dN_Only_All_Available_Tests = Avg_Delta_vs_HA_dN_Only_All_Tests,
    Median_Delta_vs_HA_dN_Only_All_Available_Tests = Median_Delta_vs_HA_dN_Only_All_Tests,
    Avg_Delta_vs_Control1_All_dN_All_Available_Tests = Avg_Delta_vs_Control1_All_dN_All_Tests,
    Avg_Delta_vs_Control2_All_dN_All_dS_All_Available_Tests = Avg_Delta_vs_Control2_All_dN_All_dS_All_Tests,
    Fixed_Model_Selection_Score,
    N_Test_Periods,
    Test_Periods,
    N_Fixed_Terms,
    Ranking_Group = "1_Representative_dN_model_better_than_HA"
  )

scheme_ranking_all_available_tests <- bind_rows(
  ha_only_by_scheme_all_available,
  representative_better_than_HA_all_available
) %>%
  arrange(
    Window_Size,
    Ranking_Group,
    desc(Avg_AUC_All_Available_Tests),
    desc(Avg_Delta_vs_HA_dN_Only_All_Available_Tests),
    N_Fixed_Terms
  )

# ---- 8B. Common-test-period ranking across all schemes ----
# This is stricter when comparing different widths, because width 6 has more test periods
# than width 10. Common tests are the periods that every scheme can evaluate.
common_test_periods <- scored_results %>%
  distinct(Window_Scheme, Test_Period, Test_Index) %>%
  group_by(Test_Period, Test_Index) %>%
  summarise(N_Schemes = n_distinct(Window_Scheme), .groups = "drop") %>%
  filter(N_Schemes == n_distinct(scored_results$Window_Scheme)) %>%
  arrange(Test_Index)

model_by_scheme_common_summary <- scored_results %>%
  semi_join(common_test_periods, by = c("Test_Period", "Test_Index")) %>%
  filter(Class == "representative_dN_candidate", !is.na(AUC)) %>%
  group_by(Window_Scheme, Window_Type, Window_Size, Model_Name, Class, Formula) %>%
  summarise(
    Avg_AUC_Common_Tests = mean(AUC, na.rm = TRUE),
    Median_AUC_Common_Tests = median(AUC, na.rm = TRUE),
    Min_AUC_Common_Tests = min(AUC, na.rm = TRUE),
    SD_AUC_Common_Tests = sd(AUC, na.rm = TRUE),
    Avg_Delta_vs_HA_dN_Only_Common_Tests = mean(Delta_vs_HA_dN_Only, na.rm = TRUE),
    Avg_Delta_vs_Control1_All_dN_Common_Tests = mean(Delta_vs_Control1_All_dN, na.rm = TRUE),
    Avg_Delta_vs_Control2_All_dN_All_dS_Common_Tests = mean(Delta_vs_Control2_All_dN_All_dS, na.rm = TRUE),
    N_Common_Test_Periods = n_distinct(Test_Period),
    Common_Test_Periods = paste(unique(Test_Period[order(Test_Index)]), collapse = ","),
    N_Fixed_Terms = first(N_Fixed_Terms),
    .groups = "drop"
  ) %>%
  mutate(
    Common_Test_Selection_Score = Avg_AUC_Common_Tests -
      complexity_penalty_per_term * N_Fixed_Terms
  )

scheme_ranking_common_tests <- model_by_scheme_common_summary %>%
  group_by(Window_Scheme, Window_Type, Window_Size) %>%
  slice_max(Common_Test_Selection_Score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    Window_Scheme, Window_Type, Window_Size,
    Best_Biological_Model_Common_Tests = Model_Name,
    Class, Formula,
    Avg_AUC_Common_Tests,
    Median_AUC_Common_Tests,
    Min_AUC_Common_Tests,
    SD_AUC_Common_Tests,
    Avg_Delta_vs_HA_dN_Only_Common_Tests,
    Avg_Delta_vs_Control1_All_dN_Common_Tests,
    Avg_Delta_vs_Control2_All_dN_All_dS_Common_Tests,
    Common_Test_Selection_Score,
    N_Common_Test_Periods,
    Common_Test_Periods,
    N_Fixed_Terms
  ) %>%
  arrange(desc(Common_Test_Selection_Score))

# ---- 8B. Common-test-period ranking: show all representative models better than HA dN only ----
# Table 5B now includes ALL representative dN models whose average common-test AUC exceeds
# HA dN only within the same rolling scheme.

# ---- 8B-2. Table 5B: HA dN only + all representative models better than HA dN only under common tests ----

ha_only_by_scheme_common_tests <- scored_results %>%
  semi_join(common_test_periods, by = c("Test_Period", "Test_Index")) %>%
  filter(Class == "ha_dN_only", !is.na(AUC)) %>%
  group_by(Window_Scheme, Window_Type, Window_Size, Model_Name, Class, Formula) %>%
  summarise(
    Avg_AUC_Common_Tests = mean(AUC, na.rm = TRUE),
    Median_AUC_Common_Tests = median(AUC, na.rm = TRUE),
    Min_AUC_Common_Tests = min(AUC, na.rm = TRUE),
    SD_AUC_Common_Tests = sd(AUC, na.rm = TRUE),
    Avg_Delta_vs_HA_dN_Only_Common_Tests = 0,
    Avg_Delta_vs_Control1_All_dN_Common_Tests = mean(Delta_vs_Control1_All_dN, na.rm = TRUE),
    Avg_Delta_vs_Control2_All_dN_All_dS_Common_Tests = mean(Delta_vs_Control2_All_dN_All_dS, na.rm = TRUE),
    Common_Test_Selection_Score = mean(Row_Selection_Score, na.rm = TRUE),
    N_Common_Test_Periods = n_distinct(Test_Period),
    Common_Test_Periods = paste(unique(Test_Period[order(Test_Index)]), collapse = ","),
    N_Fixed_Terms = first(N_Fixed_Terms),
    .groups = "drop"
  ) %>%
  mutate(
    Ranking_Group = "0_HA_dN_Only_baseline"
  )

representative_better_than_HA_common_tests <- model_by_scheme_common_summary %>%
  filter(
    Class == "representative_dN_candidate",
    !is.na(Avg_AUC_Common_Tests),
    !is.na(Avg_Delta_vs_HA_dN_Only_Common_Tests),
    Avg_Delta_vs_HA_dN_Only_Common_Tests > 0
  ) %>%
  transmute(
    Window_Scheme,
    Window_Type,
    Window_Size,
    Model_Name,
    Class,
    Formula,
    Avg_AUC_Common_Tests,
    Median_AUC_Common_Tests,
    Min_AUC_Common_Tests,
    SD_AUC_Common_Tests,
    Avg_Delta_vs_HA_dN_Only_Common_Tests,
    Avg_Delta_vs_Control1_All_dN_Common_Tests,
    Avg_Delta_vs_Control2_All_dN_All_dS_Common_Tests,
    Common_Test_Selection_Score,
    N_Common_Test_Periods,
    Common_Test_Periods,
    N_Fixed_Terms,
    Ranking_Group = "1_Representative_dN_model_better_than_HA"
  )

scheme_ranking_common_tests <- bind_rows(
  ha_only_by_scheme_common_tests,
  representative_better_than_HA_common_tests
) %>%
  arrange(
    Window_Size,
    Ranking_Group,
    desc(Avg_AUC_Common_Tests),
    desc(Avg_Delta_vs_HA_dN_Only_Common_Tests),
    N_Fixed_Terms
  )

# ---- 8C. Oracle per-test upper bound: diagnostic only, not final model selection ----
# This answers: "If the best model could change at every test period, how high could each
# scheme score?" Use this as supplementary diagnostic, not as final selection.
oracle_best_biological_per_scheme_test <- scored_results %>%
  filter(Class == "representative_dN_candidate", !is.na(AUC)) %>%
  group_by(Window_Scheme, Window_Type, Window_Size, Test_Period, Test_Index) %>%
  slice_max(Row_Selection_Score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(Window_Size, Window_Type, Test_Index)

oracle_scheme_summary <- oracle_best_biological_per_scheme_test %>%
  group_by(Window_Scheme, Window_Type, Window_Size) %>%
  summarise(
    Oracle_Avg_AUC = mean(AUC, na.rm = TRUE),
    Oracle_Median_AUC = median(AUC, na.rm = TRUE),
    Oracle_Min_AUC = min(AUC, na.rm = TRUE),
    Oracle_Avg_Delta_vs_HA_dN_Only = mean(Delta_vs_HA_dN_Only, na.rm = TRUE),
    Oracle_Avg_Delta_vs_Control1_All_dN = mean(Delta_vs_Control1_All_dN, na.rm = TRUE),
    Oracle_Avg_Delta_vs_Control2_All_dN_All_dS = mean(Delta_vs_Control2_All_dN_All_dS, na.rm = TRUE),
    Oracle_Avg_Row_Selection_Score = mean(Row_Selection_Score, na.rm = TRUE),
    N_Test_Periods = n_distinct(Test_Period),
    Test_Periods = paste(unique(Test_Period[order(Test_Index)]), collapse = ","),
    .groups = "drop"
  ) %>%
  arrange(desc(Oracle_Avg_Row_Selection_Score))

# ---- Final rolling width selected for downstream comparison ----
final_window_size <- 10
final_window_scheme <- paste0("rolling_width_", final_window_size)

# ---- 8D. Explicit rolling width=8 summary because this is the current setting ----
# ---- 8D. Explicit final rolling width summary ----
# Since rolling_width_10 was selected as the final temporal CV scheme,
# all downstream model-vs-control comparisons are generated for rolling width 10.

current_final_rolling_summary <- selected_representative_vs_controls_by_scheme %>%
  filter(Window_Type == "rolling", Window_Size == final_window_size) %>%
  arrange(Comparison_Group)

current_final_rolling_period_level <- scored_results %>%
  filter(Window_Type == "rolling", Window_Size == final_window_size) %>%
  left_join(
    best_fixed_model_per_scheme %>%
      filter(Window_Type == "rolling", Window_Size == final_window_size) %>%
      select(Window_Scheme, Selected_Representative_Model = Model_Name),
    by = "Window_Scheme"
  ) %>%
  mutate(
    Comparison_Group = case_when(
      Class == "representative_dN_candidate" & Model_Name == Selected_Representative_Model ~ "1_Selected_representative_dN_subset",
      Class == "ha_dN_only" ~ "2_HA_dN_only",
      Class == "control1_all_dN" ~ "3_Control1_all_dN",
      Class == "control2_all_dN_all_dS" ~ "4_Control2_all_dN_all_dS",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Comparison_Group), !is.na(AUC)) %>%
  arrange(Test_Index, Comparison_Group)

# ---- 9. Final table guide ----
final_recommendation <- tibble(
  Recommendation_Type = c(
    "Rolling_model_ranking_all_available_tests",
    "Rolling_model_ranking_common_tests",
    "Final_rolling_width10_model_vs_controls"
  ),
  Recommended_Table = c(
    "Table_5A_Rolling_Models_Better_Than_HA_dN_Only_ALL_AVAILABLE_TESTS.csv",
    "Table_5B_Rolling_Models_Better_Than_HA_dN_Only_COMMON_TESTS.csv",
    "Table_8A_Final_RollingWidth10_Selected_vs_Controls_SUMMARY.csv"
  ),
  Interpretation = c(
    "Includes HA dN only and all representative dN-subset models that outperform HA dN only using all available rolling test periods for each width.",
    "Includes HA dN only and all representative dN-subset models that outperform HA dN only using only test periods shared across all rolling widths.",
    "Compares the selected representative dN-subset model against HA-only, all-dN, and all-dN+all-dS controls under the final rolling width 10 setting."
  )
)

# ---- 10. Export ----
result_dir <- file.path(base_path, "3_Result/7_Rolling_Window_Sensitivity")
if (!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE)

write.csv(
  scheme_ranking_all_available_tests,
  file.path(result_dir, "Table_5A_Rolling_Models_Better_Than_HA_dN_Only_ALL_AVAILABLE_TESTS.csv"),
  row.names = FALSE
)

write.csv(
  scheme_ranking_common_tests,
  file.path(result_dir, "Table_5B_Rolling_Models_Better_Than_HA_dN_Only_COMMON_TESTS.csv"),
  row.names = FALSE
)

write.csv(
  current_final_rolling_summary,
  file.path(result_dir, "Table_8A_Final_RollingWidth10_Selected_vs_Controls_SUMMARY.csv"),
  row.names = FALSE
)

write.csv(
  current_final_rolling_period_level,
  file.path(result_dir, "Table_8B_Final_RollingWidth10_Selected_vs_Controls_PERIOD_LEVEL.csv"),
  row.names = FALSE
)

write.csv(window_grid, file.path(result_dir, "Table_0_Rolling_Window_Grid_ALL_SPLITS.csv"), row.names = FALSE)
write.csv(all_results, file.path(result_dir, "Table_1_Rolling_All_Model_AUCs_Per_Split.csv"), row.names = FALSE)
write.csv(scored_results, file.path(result_dir, "Table_2_Rolling_Scored_Model_AUCs_Per_Split.csv"), row.names = FALSE)
write.csv(model_by_scheme_summary, file.path(result_dir, "Table_3_Rolling_Model_By_Scheme_ALL_AVAILABLE_TESTS.csv"), row.names = FALSE)
write.csv(best_fixed_model_per_scheme, file.path(result_dir, "Table_4_Rolling_Best_Fixed_Model_Per_Scheme.csv"), row.names = FALSE)
write.csv(selected_representative_vs_controls_by_scheme, file.path(result_dir, "Table_4B_Rolling_Selected_Representative_vs_Controls_BY_SCHEME.csv"), row.names = FALSE)
write.csv(common_test_periods, file.path(result_dir, "Table_5C_Rolling_Common_Test_Periods.csv"), row.names = FALSE)
write.csv(fixed_model_per_test, file.path(result_dir, "Table_6_Rolling_Fixed_Model_Per_Test_Performance.csv"), row.names = FALSE)
write.csv(oracle_best_biological_per_scheme_test, file.path(result_dir, "Table_9A_Rolling_Oracle_Best_Model_Per_Scheme_Test_DIAGNOSTIC.csv"), row.names = FALSE)
write.csv(oracle_scheme_summary, file.path(result_dir, "Table_9B_Rolling_Oracle_Scheme_Summary_DIAGNOSTIC.csv"), row.names = FALSE)

