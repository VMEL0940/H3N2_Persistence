# ==============================================================================
# Script Name: 10_Final_Rolling_Fixed_Models_TrainTest_AUROC.R
# Purpose:
#   After selecting the final rolling-window model, refit ONLY the final model and
#   predefined controls on the same rolling CV splits, and report train/test AUROC.
#
# Models compared:
#   1) Final representative dN model
#      dom_binary ~ HA * PB1 + PB2 + NA_seg + (1 | vaccine_code)
#   2) HA dN only
#      dom_binary ~ HA + (1 | vaccine_code)
#   3) Control 1: all dN segments
#      dom_binary ~ PB2 + PB1 + PA + HA + NP + NA_seg + M + NS + (1 | vaccine_code)
#   4) Control 2: all dN + all dS segments
#      dom_binary ~ PB2 + PB1 + PA + HA + NP + NA_seg + M + NS +
#                   PB2_syn + PB1_syn + PA_syn + HA_syn + NP_syn + NA_syn + M_syn + NS_syn +
#                   (1 | vaccine_code)
#
# Notes:
#   - Scaling is fit on the training split only, then applied to both train/test.
#   - Test AUROC uses fixed effects only: re.form = ~0
#     because the next vaccine_code level should not borrow a fitted random intercept.
#   - Train AUROC is reported in two ways:
#       a) Train_AUC_FixedOnly: same prediction rule as test; preferred for fair comparison.
#       b) Train_AUC_ConditionalRE: includes fitted random intercepts; useful to inspect overfit.
# ==============================================================================

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

# Final validation setting
final_window_type <- "rolling"
final_window_width <- 10
last_test_period <- "TH22"   # CR23 is kept for final deployment, not CV model selection.

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

safe_auc <- function(y, p) {
  if (length(unique(y[!is.na(y)])) < 2) return(NA_real_)
  as.numeric(auc(roc(y, p, quiet = TRUE)))
}

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

# ---- 3. Final fixed models ----
final_model_bank <- c(
  "1_Final_Rolling_Best_dN_HA_x_PB1_plus_PB2_NA" =
    "dom_binary ~ HA * PB1 + PB2 + NA_seg + (1 | vaccine_code)",

  "2_HA_dN_Only" =
    "dom_binary ~ HA + (1 | vaccine_code)",

  "3_Control1_All_dN_PB2_to_NS" =
    "dom_binary ~ PB2 + PB1 + PA + HA + NP + NA_seg + M + NS + (1 | vaccine_code)",

  "4_Control2_All_dN_All_dS_PB2_to_NS" =
    "dom_binary ~ PB2 + PB1 + PA + HA + NP + NA_seg + M + NS + PB2_syn + PB1_syn + PA_syn + HA_syn + NP_syn + NA_syn + M_syn + NS_syn + (1 | vaccine_code)"
)

model_class <- tibble(
  Model_Name = names(final_model_bank),
  Model_Role = c(
    "Selected representative dN model",
    "HA dN only baseline",
    "Control 1: all dN segments",
    "Control 2: all dN + all dS segments"
  ),
  Formula = unname(final_model_bank),
  N_Fixed_Terms = sapply(unname(final_model_bank), count_fixed_terms)
)

# ---- 4. Rolling window grid ----
make_rolling_windows <- function(all_periods, last_test_period, rolling_width) {
  max_idx <- which(all_periods == last_test_period)
  if (rolling_width >= max_idx) stop("rolling_width must be smaller than last_test_period index.")

  windows <- list()
  for (start_i in 1:(max_idx - rolling_width)) {
    end_i <- start_i + rolling_width - 1
    windows[[length(windows) + 1]] <- tibble(
      Window_Scheme = paste0("rolling_width_", rolling_width),
      Window_Type = "rolling",
      Window_Size = rolling_width,
      Train_Start_Index = start_i,
      Train_End_Index = end_i,
      Test_Index = end_i + 1,
      Train_Periods = paste(all_periods[start_i:end_i], collapse = ","),
      Test_Period = all_periods[end_i + 1],
      Train_Range_Label = paste0(all_periods[start_i], "-", all_periods[end_i], " -> ", all_periods[end_i + 1])
    )
  }

  bind_rows(windows) %>%
    mutate(Test_Period = as.character(Test_Period)) %>%
    filter(Test_Period != "CR23")
}

window_grid <- make_rolling_windows(all_periods, last_test_period, final_window_width)

# ---- 5. Fit/evaluate one rolling split ----
fit_one_final_window <- function(win_row) {
  train_periods <- strsplit(win_row$Train_Periods, ",")[[1]]
  test_period <- win_row$Test_Period

  raw_train <- raw_df %>% filter(as.character(Period) %in% train_periods)
  raw_test  <- raw_df %>% filter(as.character(Period) == test_period)

  if (nrow(raw_train) < 10 || nrow(raw_test) < 5) return(NULL)
  if (length(unique(raw_train$dom_binary)) < 2 || length(unique(raw_test$dom_binary)) < 2) return(NULL)

  scaler <- fit_scaler(raw_train, all_raw_features)
  train_df <- apply_scaler(raw_train, scaler) %>% recalc_weights()
  test_df  <- apply_scaler(raw_test, scaler)

  out <- vector("list", length(final_model_bank))
  names(out) <- names(final_model_bank)

  for (m_name in names(final_model_bank)) {
    ftxt <- final_model_bank[[m_name]]

    out[[m_name]] <- tryCatch({
      fit <- suppressWarnings(glmer(
        as.formula(ftxt), data = train_df, family = binomial, weights = norm_wts,
        nAGQ = 0,
        control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
      ))

      # Fair primary comparison: fixed-effect-only predictions for both train and test.
      train_probs_fixed <- predict(fit, newdata = train_df, type = "response", re.form = ~0, allow.new.levels = TRUE)
      test_probs_fixed  <- predict(fit, newdata = test_df,  type = "response", re.form = ~0, allow.new.levels = TRUE)

      # Diagnostic in-sample fit: includes fitted random intercepts for training data.
      train_probs_conditional <- predict(fit, newdata = train_df, type = "response", re.form = NULL, allow.new.levels = TRUE)

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
        Train_N = nrow(train_df),
        Train_Pos = sum(train_df$dom_binary == 1, na.rm = TRUE),
        Train_Neg = sum(train_df$dom_binary == 0, na.rm = TRUE),
        Test_N = nrow(test_df),
        Test_Pos = sum(test_df$dom_binary == 1, na.rm = TRUE),
        Test_Neg = sum(test_df$dom_binary == 0, na.rm = TRUE),
        Train_AUC_FixedOnly = safe_auc(train_df$dom_binary, train_probs_fixed),
        Train_AUC_ConditionalRE = safe_auc(train_df$dom_binary, train_probs_conditional),
        Test_AUC_FixedOnly = safe_auc(test_df$dom_binary, test_probs_fixed),
        Train_minus_Test_AUC_FixedOnly = Train_AUC_FixedOnly - Test_AUC_FixedOnly,
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
        Train_N = nrow(train_df),
        Train_Pos = sum(train_df$dom_binary == 1, na.rm = TRUE),
        Train_Neg = sum(train_df$dom_binary == 0, na.rm = TRUE),
        Test_N = nrow(test_df),
        Test_Pos = sum(test_df$dom_binary == 1, na.rm = TRUE),
        Test_Neg = sum(test_df$dom_binary == 0, na.rm = TRUE),
        Train_AUC_FixedOnly = NA_real_,
        Train_AUC_ConditionalRE = NA_real_,
        Test_AUC_FixedOnly = NA_real_,
        Train_minus_Test_AUC_FixedOnly = NA_real_,
        AIC = NA_real_,
        N_Fixed_Terms = count_fixed_terms(ftxt),
        Error = e$message
      )
    })
  }

  bind_rows(out)
}

# ---- 6. Run final model/control comparison ----
message("Rolling splits to evaluate: ", nrow(window_grid))
message("Fixed models per split: ", length(final_model_bank))

final_train_test_auc_by_split <- pmap_dfr(window_grid, function(...) fit_one_final_window(tibble(...))) %>%
  left_join(model_class %>% select(Model_Name, Model_Role), by = "Model_Name") %>%
  relocate(Model_Role, .after = Model_Name)

final_train_test_auc_summary <- final_train_test_auc_by_split %>%
  filter(is.na(Error)) %>%
  group_by(Model_Name, Model_Role, Formula, N_Fixed_Terms) %>%
  summarise(
    N_Splits = n(),
    Test_Periods = paste(Test_Period[order(Test_Index)], collapse = ","),
    Mean_Train_AUC_FixedOnly = mean(Train_AUC_FixedOnly, na.rm = TRUE),
    Median_Train_AUC_FixedOnly = median(Train_AUC_FixedOnly, na.rm = TRUE),
    Min_Train_AUC_FixedOnly = min(Train_AUC_FixedOnly, na.rm = TRUE),
    SD_Train_AUC_FixedOnly = sd(Train_AUC_FixedOnly, na.rm = TRUE),
    Mean_Train_AUC_ConditionalRE = mean(Train_AUC_ConditionalRE, na.rm = TRUE),
    Median_Train_AUC_ConditionalRE = median(Train_AUC_ConditionalRE, na.rm = TRUE),
    Mean_Test_AUC_FixedOnly = mean(Test_AUC_FixedOnly, na.rm = TRUE),
    Median_Test_AUC_FixedOnly = median(Test_AUC_FixedOnly, na.rm = TRUE),
    Min_Test_AUC_FixedOnly = min(Test_AUC_FixedOnly, na.rm = TRUE),
    SD_Test_AUC_FixedOnly = sd(Test_AUC_FixedOnly, na.rm = TRUE),
    Mean_Train_minus_Test_AUC_FixedOnly = mean(Train_minus_Test_AUC_FixedOnly, na.rm = TRUE),
    Median_Train_minus_Test_AUC_FixedOnly = median(Train_minus_Test_AUC_FixedOnly, na.rm = TRUE),
    Mean_AIC = mean(AIC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Mean_Test_AUC_FixedOnly))

# Pairwise comparison versus the selected final representative model.
selected_model_name <- "1_Final_Rolling_Best_dN_HA_x_PB1_plus_PB2_NA"
final_pairwise_vs_selected_by_split <- final_train_test_auc_by_split %>%
  filter(is.na(Error)) %>%
  select(Test_Period, Test_Index, Train_Range_Label, Model_Name, Model_Role,
         Train_AUC_FixedOnly, Train_AUC_ConditionalRE, Test_AUC_FixedOnly) %>%
  pivot_wider(
    id_cols = c(Test_Period, Test_Index, Train_Range_Label),
    names_from = Model_Name,
    values_from = c(Train_AUC_FixedOnly, Train_AUC_ConditionalRE, Test_AUC_FixedOnly)
  )

# Long-format deltas are safer than very wide names for downstream plotting.
selected_auc_by_split <- final_train_test_auc_by_split %>%
  filter(Model_Name == selected_model_name, is.na(Error)) %>%
  select(Test_Period, Test_Index,
         Selected_Train_AUC_FixedOnly = Train_AUC_FixedOnly,
         Selected_Train_AUC_ConditionalRE = Train_AUC_ConditionalRE,
         Selected_Test_AUC_FixedOnly = Test_AUC_FixedOnly)

final_delta_vs_selected_by_split <- final_train_test_auc_by_split %>%
  filter(Model_Name != selected_model_name, is.na(Error)) %>%
  left_join(selected_auc_by_split, by = c("Test_Period", "Test_Index")) %>%
  mutate(
    Delta_Test_AUC_vs_Selected = Test_AUC_FixedOnly - Selected_Test_AUC_FixedOnly,
    Delta_Train_AUC_FixedOnly_vs_Selected = Train_AUC_FixedOnly - Selected_Train_AUC_FixedOnly,
    Delta_Train_AUC_ConditionalRE_vs_Selected = Train_AUC_ConditionalRE - Selected_Train_AUC_ConditionalRE
  ) %>%
  select(Window_Scheme, Test_Period, Test_Index, Train_Range_Label,
         Model_Name, Model_Role, Formula,
         Train_AUC_FixedOnly, Selected_Train_AUC_FixedOnly, Delta_Train_AUC_FixedOnly_vs_Selected,
         Train_AUC_ConditionalRE, Selected_Train_AUC_ConditionalRE, Delta_Train_AUC_ConditionalRE_vs_Selected,
         Test_AUC_FixedOnly, Selected_Test_AUC_FixedOnly, Delta_Test_AUC_vs_Selected)

final_delta_vs_selected_summary <- final_delta_vs_selected_by_split %>%
  group_by(Model_Name, Model_Role, Formula) %>%
  summarise(
    N_Splits = n(),
    Mean_Delta_Test_AUC_vs_Selected = mean(Delta_Test_AUC_vs_Selected, na.rm = TRUE),
    Median_Delta_Test_AUC_vs_Selected = median(Delta_Test_AUC_vs_Selected, na.rm = TRUE),
    N_Test_AUC_Wins_vs_Selected = sum(Delta_Test_AUC_vs_Selected > 0, na.rm = TRUE),
    N_Test_AUC_Losses_vs_Selected = sum(Delta_Test_AUC_vs_Selected < 0, na.rm = TRUE),
    Wilcoxon_p_Test_AUC_vs_Selected = tryCatch(
      wilcox.test(Test_AUC_FixedOnly, Selected_Test_AUC_FixedOnly, paired = TRUE, exact = FALSE)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  arrange(desc(Mean_Delta_Test_AUC_vs_Selected))

# ---- 7. Export ----
result_dir <- file.path(base_path, "3_Result/8_Final_Rolling_TrainTest_AUROC")
if (!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE)

write.csv(window_grid, file.path(result_dir, "Table_F0_Final_Rolling_Window_Grid.csv"), row.names = FALSE)
write.csv(model_class, file.path(result_dir, "Table_F1_Final_Model_Definitions.csv"), row.names = FALSE)
write.csv(final_train_test_auc_by_split, file.path(result_dir, "Table_F2_Train_Test_AUROC_By_Rolling_Split.csv"), row.names = FALSE)
write.csv(final_train_test_auc_summary, file.path(result_dir, "Table_F3_Train_Test_AUROC_Summary.csv"), row.names = FALSE)
write.csv(final_delta_vs_selected_by_split, file.path(result_dir, "Table_F4_Delta_vs_Selected_By_Split.csv"), row.names = FALSE)
write.csv(final_delta_vs_selected_summary, file.path(result_dir, "Table_F5_Delta_vs_Selected_Summary.csv"), row.names = FALSE)

cat("\n### Done: final rolling fixed-model train/test AUROC comparison complete. ###\n")
cat("\nPrimary summary table: Table_F3_Train_Test_AUROC_Summary.csv\n")
print(final_train_test_auc_summary)

cat("\nDelta vs selected final model: Table_F5_Delta_vs_Selected_Summary.csv\n")
print(final_delta_vs_selected_summary)


# ============================================================
# 6. Publication-style Figures: Final rolling test AUROC
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
})

# ------------------------------------------------------------
# Directory
# ------------------------------------------------------------

out_dir <- "~/3_Result/8_Final_Rolling_TrainTest_AUROC"

fig_dir <- file.path(out_dir, "Figures3")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# Load results
# ------------------------------------------------------------

split_auc <- read_csv(
  file.path(out_dir, "Table_F2_Train_Test_AUROC_By_Rolling_Split.csv"),
  show_col_types = FALSE
)

auc_summary <- read_csv(
  file.path(out_dir, "Table_F3_Train_Test_AUROC_Summary.csv"),
  show_col_types = FALSE
)

# ------------------------------------------------------------
# Model labels and colors
# ------------------------------------------------------------

model_labels <- c(
  "Selected_HA_x_PB1_plus_PB2_NA" = "Selected\nHA × PB1 + PB2 + NA",
  "HA_dN_Only" = "HA dN\nonly",
  "Control1_All_dN_PB2_to_NS" = "Control 1\nall dN",
  "Control2_All_dN_All_dS_PB2_to_NS" = "Control 2\nall dN + all dS"
)

model_map <- tibble::tribble(
  ~Model_Name, ~Model_Label,
  "1_Final_Rolling_Best_dN_HA_x_PB1_plus_PB2_NA", "Selected\nHA × PB1 + PB2 + NA",
  "2_HA_dN_Only", "HA dN\nonly",
  "3_Control1_All_dN_PB2_to_NS", "Control 1\nall dN",
  "4_Control2_All_dN_All_dS_PB2_to_NS", "Control 2\nall dN + all dS"
)

model_levels <- model_map$Model_Label

model_colors <- c(
  "Selected\nHA × PB1 + PB2 + NA" = "#D55E00",
  "HA dN\nonly" = "#0072B2",
  "Control 1\nall dN" = "#009E73",
  "Control 2\nall dN + all dS" = "#CC79A7"
)

split_auc <- split_auc %>%
  select(-any_of("Model_Label")) %>%
  left_join(model_map, by = "Model_Name") %>%
  mutate(
    Model_Label = factor(Model_Label, levels = model_levels),
    Test_Period = factor(Test_Period, levels = unique(Test_Period))
  )

auc_summary <- auc_summary %>%
  select(-any_of("Model_Label")) %>%
  left_join(model_map, by = "Model_Name") %>%
  mutate(
    Model_Label = factor(Model_Label, levels = model_levels)
  )

# ------------------------------------------------------------
# Helper: summary from split-level test AUROC
# ------------------------------------------------------------

test_summary <- split_auc %>%
  group_by(Model_Name, Model_Label) %>%
  summarise(
    N_Splits = sum(!is.na(Test_AUC_FixedOnly)),
    Mean_Test_AUC = mean(Test_AUC_FixedOnly, na.rm = TRUE),
    SD_Test_AUC = sd(Test_AUC_FixedOnly, na.rm = TRUE),
    SE_Test_AUC = SD_Test_AUC / sqrt(N_Splits),
    CI95_Low = Mean_Test_AUC - qt(0.975, df = N_Splits - 1) * SE_Test_AUC,
    CI95_High = Mean_Test_AUC + qt(0.975, df = N_Splits - 1) * SE_Test_AUC,
    .groups = "drop"
  )

# ============================================================
# Figure 1A. Main figure: mean test AUROC with split-level dots
# ============================================================

# ============================================================
# Figure 1A. Main figure: mean test AUROC ± SD with split-level dots
# ============================================================

test_summary <- split_auc %>%
  group_by(Model_Name, Model_Label) %>%
  summarise(
    N_Splits = sum(!is.na(Test_AUC_FixedOnly)),
    Mean_Test_AUC = mean(Test_AUC_FixedOnly, na.rm = TRUE),
    SD_Test_AUC = sd(Test_AUC_FixedOnly, na.rm = TRUE),
    Low = Mean_Test_AUC - SD_Test_AUC,
    High = Mean_Test_AUC + SD_Test_AUC,
    .groups = "drop"
  )

p_main_A <- ggplot() +
  geom_hline(
    yintercept = 0.5,
    linetype = "dashed",
    linewidth = 0.45,
    color = "gray55"
  ) +
  geom_jitter(
    data = split_auc,
    aes(x = Model_Label, y = Test_AUC_FixedOnly, color = Model_Label),
    width = 0.10,
    height = 0,
    size = 2.5,
    alpha = 0.60
  ) +
  geom_errorbar(
    data = test_summary,
    aes(x = Model_Label, ymin = Low, ymax = High, color = Model_Label),
    width = 0.18,
    linewidth = 0.9
  ) +
  geom_point(
    data = test_summary,
    aes(x = Model_Label, y = Mean_Test_AUC, fill = Model_Label),
    shape = 21,
    size = 5.2,
    color = "black",
    stroke = 0.7
  ) +
  geom_text(
    data = test_summary,
    aes(
      x = Model_Label,
      y = pmin(Mean_Test_AUC + 0.08, 0.98),
      label = sprintf("%.3f", Mean_Test_AUC)
    ),
    size = 4.0,
    fontface = "bold"
  ) +
  scale_color_manual(values = model_colors, guide = "none") +
  scale_fill_manual(values = model_colors, guide = "none") +
  coord_cartesian(ylim = c(0.35, 1.02)) +
  labs(
    title = "Final model performance under rolling-window validation",
    subtitle = "Dots indicate individual rolling test periods; points and bars indicate mean test AUROC ± SD",
    x = NULL,
    y = "Test AUROC"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_text(face = "bold"),
    axis.line = element_line(linewidth = 0.5)
    # plot.margin = margin(10, 20, 10, 20)
  )

ggsave(
  file.path(fig_dir, "Main_Figure_1A_Test_AUROC_Mean_SD_Dots.svg"),
  p_main_A,
  width = 8.8,
  height = 5.2,
  dpi = 400
)

ggsave(
  file.path(fig_dir, "Main_Figure_1A_Test_AUROC_Mean_SD_Dots.pdf"),
  p_main_A,
  width = 8.8,
  height = 5.2
)

# ============================================================
# Figure 1B. Rolling test AUROC trajectory
# ============================================================

p_main_B <- ggplot(
  split_auc,
  aes(
    x = Test_Period,
    y = Test_AUC_FixedOnly,
    group = Model_Label,
    color = Model_Label
  )
) +
  geom_hline(
    yintercept = 0.5,
    linetype = "dashed",
    linewidth = 0.45,
    color = "gray50"
  ) +
  geom_line(linewidth = 1.0, alpha = 0.9) +
  geom_point(size = 2.5, alpha = 0.95) +
  scale_color_manual(values = model_colors) +
  coord_cartesian(ylim = c(0.4, 1.0)) +
  labs(
    title = "Rolling test AUROC across future periods",
    subtitle = "Performance was evaluated only on held-out future test periods",
    x = "Future test period",
    y = "Test AUROC",
    color = "Model"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank()
  )

ggsave(
  file.path(fig_dir, "Main_Figure_1B_Rolling_Test_AUROC_Trajectory.svg"),
  p_main_B,
  width = 9.5,
  height = 5.4,
  dpi = 400
)

ggsave(
  file.path(fig_dir, "Main_Figure_1B_Rolling_Test_AUROC_Trajectory.pdf"),
  p_main_B,
  width = 9.5,
  height = 5.4
)

# ============================================================
# Figure 1C. Delta test AUROC versus selected model
# ============================================================

selected_model_name <- "1_Final_Rolling_Best_dN_HA_x_PB1_plus_PB2_NA"

selected_by_split <- split_auc %>%
  filter(Model_Name == selected_model_name) %>%
  select(
    Test_Period,
    Selected_Test_AUC = Test_AUC_FixedOnly
  )

delta_df <- split_auc %>%
  left_join(selected_by_split, by = "Test_Period") %>%
  mutate(
    Delta_vs_Selected = Test_AUC_FixedOnly - Selected_Test_AUC
  ) %>%
  filter(Model_Name != selected_model_name)

delta_df %>%
  select(Model_Name, Test_Period, Test_AUC_FixedOnly, Selected_Test_AUC, Delta_vs_Selected)

delta_summary <- delta_df %>%
  group_by(Model_Name, Model_Label) %>%
  summarise(
    N_Splits = sum(!is.na(Delta_vs_Selected)),
    Mean_Delta = mean(Delta_vs_Selected, na.rm = TRUE),
    SD_Delta = sd(Delta_vs_Selected, na.rm = TRUE),
    SE_Delta = SD_Delta / sqrt(N_Splits),
    CI95_Low = Mean_Delta - qt(0.975, df = N_Splits - 1) * SE_Delta,
    CI95_High = Mean_Delta + qt(0.975, df = N_Splits - 1) * SE_Delta,
    .groups = "drop"
  )

p_main_C <- ggplot() +
  geom_hline(
    yintercept = 0,
    linewidth = 0.55,
    color = "gray35"
  ) +
  geom_jitter(
    data = delta_df,
    aes(x = Model_Label, y = Delta_vs_Selected, color = Model_Label),
    width = 0.10,
    height = 0,
    size = 2.2,
    alpha = 0.55
  ) +
  geom_errorbar(
    data = delta_summary,
    aes(x = Model_Label, ymin = CI95_Low, ymax = CI95_High, color = Model_Label),
    width = 0.18,
    linewidth = 0.9
  ) +
  geom_point(
    data = delta_summary,
    aes(x = Model_Label, y = Mean_Delta, fill = Model_Label),
    shape = 21,
    size = 5.0,
    color = "black",
    stroke = 0.6
  ) +
  geom_text(
    data = delta_summary,
    aes(
      x = Model_Label,
      y = Mean_Delta,
      label = sprintf("%.3f", Mean_Delta)
    ),
    nudge_y = -0.045,
    size = 3.7,
    fontface = "bold"
  ) +
  scale_color_manual(values = model_colors, guide = "none") +
  scale_fill_manual(values = model_colors, guide = "none") +
  labs(
    title = "Difference in test AUROC relative to the selected model",
    subtitle = "Negative values indicate lower performance than the selected dN-subset model",
    x = NULL,
    y = expression(Delta~"test AUROC vs selected model")
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11),
    axis.text.x = element_text(size = 11),
    axis.title.y = element_text(face = "bold")
  )

ggsave(
  file.path(fig_dir, "Main_Figure_1C_Delta_AUROC_vs_Selected.svg"),
  p_main_C,
  width = 7.6,
  height = 5.2,
  dpi = 400
)

ggsave(
  file.path(fig_dir, "Main_Figure_1C_Delta_AUROC_vs_Selected.pdf"),
  p_main_C,
  width = 7.6,
  height = 5.2
)

# ============================================================
# Optional Supplementary Figure. Train-test AUROC gap
# ============================================================

gap_df <- split_auc %>%
  mutate(
    FixedOnly_Gap = Train_AUC_FixedOnly - Test_AUC_FixedOnly,
    ConditionalRE_Gap = Train_AUC_ConditionalRE - Test_AUC_FixedOnly
  ) %>%
  select(Model_Label, Test_Period, FixedOnly_Gap, ConditionalRE_Gap) %>%
  pivot_longer(
    cols = c(FixedOnly_Gap, ConditionalRE_Gap),
    names_to = "Gap_Type",
    values_to = "Gap"
  ) %>%
  mutate(
    Gap_Type = dplyr::recode(
      Gap_Type,
      "FixedOnly_Gap" = "Train fixed-only - Test",
      "ConditionalRE_Gap" = "Train conditional RE - Test"
    ),
    Gap_Type = factor(
      Gap_Type,
      levels = c("Train fixed-only - Test", "Train conditional RE - Test")
    )
  )

gap_colors <- c(
  "Train fixed-only - Test" = "#0072B2",
  "Train conditional RE - Test" = "#D55E00"
)

p_supp_gap <- ggplot(
  gap_df,
  aes(x = Test_Period, y = Gap, group = Gap_Type, color = Gap_Type)
) +
  geom_hline(yintercept = 0, linewidth = 0.45, color = "gray35") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.0) +
  facet_wrap(~ Model_Label, ncol = 2) +
  scale_color_manual(values = gap_colors) +
  labs(
    title = "Train-test AUROC gap across rolling splits",
    subtitle = "Conditional random-effect train AUROC is shown only as an in-sample diagnostic",
    x = "Future test period",
    y = "Train AUROC - test AUROC",
    color = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

ggsave(
  file.path(fig_dir, "Supplementary_Train_Test_AUROC_Gap.svg"),
  p_supp_gap,
  width = 10,
  height = 7,
  dpi = 400
)

ggsave(
  file.path(fig_dir, "Supplementary_Train_Test_AUROC_Gap.pdf"),
  p_supp_gap,
  width = 10,
  height = 7
)

message("Saved publication-style figures to: ", fig_dir)



# ============================================================
# Supplementary Figure: Train vs Test AUROC by test period
# ============================================================

train_test_long <- split_auc %>%
  select(
    Model_Name,
    Model_Label,
    Test_Period,
    Train_AUC_FixedOnly,
    Test_AUC_FixedOnly
  ) %>%
  pivot_longer(
    cols = c(Train_AUC_FixedOnly, Test_AUC_FixedOnly),
    names_to = "Dataset",
    values_to = "AUROC"
  ) %>%
  mutate(
    Dataset = recode(
      Dataset,
      "Train_AUC_FixedOnly" = "Train",
      "Test_AUC_FixedOnly" = "Test"
    ),
    Dataset = factor(Dataset, levels = c("Train", "Test"))
  )

dataset_colors <- c(
  "Train" = "#4D4D4D",
  "Test" = "#D55E00"
)

p_train_test_bar <- ggplot(
  train_test_long,
  aes(x = Test_Period, y = AUROC, fill = Dataset)
) +
  geom_hline(
    yintercept = 0.5,
    linetype = "dashed",
    linewidth = 0.4,
    color = "gray55"
  ) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65,
    color = "black",
    linewidth = 0.25,
    alpha = 0.9
  ) +
  facet_wrap(~ Model_Label, ncol = 2) +
  scale_fill_manual(values = dataset_colors) +
  coord_cartesian(ylim = c(0.35, 1.0)) +
  labs(
    title = "Train versus test AUROC across rolling validation periods",
    subtitle = "Train AUROC was calculated using fixed effects only, matching the test prediction setting",
    x = "Future test period",
    y = "AUROC",
    fill = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "bottom"
  )

ggsave(
  file.path(fig_dir, "Supplementary_Train_vs_Test_AUROC_Barplot_By_Period.svg"),
  p_train_test_bar,
  width = 10,
  height = 7,
  dpi = 400
)

ggsave(
  file.path(fig_dir, "Supplementary_Train_vs_Test_AUROC_Barplot_By_Period.pdf"),
  p_train_test_bar,
  width = 10,
  height = 7
)


# ============================================================
# Figure: Top rolling width 10 dN-subset models better than HA dN only
# ============================================================

rolling_result_dir <- "~/3_Result/7_Rolling_Window_Sensitivity"

table5a_path <- file.path(
  rolling_result_dir,
  "Table_5A_Rolling_Models_Better_Than_HA_dN_Only_ALL_AVAILABLE_TESTS.csv"
)

# 만약 예전 파일명으로 저장했다면 fallback
if (!file.exists(table5a_path)) {
  table5a_path <- file.path(
    rolling_result_dir,
    "Table_5A_Rolling_Scheme_Ranking_ALL_AVAILABLE_TESTS.csv"
  )
}

rolling_rank <- read_csv(table5a_path, show_col_types = FALSE)

final_window_size <- 10
top_n_models <- 12

pretty_model_name <- function(x) {
  x %>%
    str_replace("^Representative_dN_", "") %>%
    str_replace("^HA_dN_Only$", "HA dN only") %>%
    str_replace_all("_plus_", " + ") %>%
    str_replace_all("_x_", " × ") %>%
    str_replace_all("_", " ") %>%
    str_replace("^HA ", "HA + ") %>%
    str_replace("HA ×", "HA ×") %>%
    str_wrap(width = 42)
}

width10_rank <- rolling_rank %>%
  filter(Window_Size == final_window_size) %>%
  mutate(
    Plot_Group = case_when(
      Class == "ha_dN_only" ~ "HA dN only baseline",
      Class == "representative_dN_candidate" ~ "dN-subset model",
      TRUE ~ Class
    )
  )

ha_baseline_auc <- width10_rank %>%
  filter(Class == "ha_dN_only") %>%
  summarise(x = mean(Avg_AUC_All_Available_Tests, na.rm = TRUE)) %>%
  pull(x)

top_width10_models <- width10_rank %>%
  filter(Class == "representative_dN_candidate") %>%
  arrange(desc(Avg_AUC_All_Available_Tests)) %>%
  slice_head(n = top_n_models) %>%
  bind_rows(
    width10_rank %>% filter(Class == "ha_dN_only")
  ) %>%
  mutate(
    Model_Label_Long = pretty_model_name(Model_Name),
    Model_Label_Long = fct_reorder(Model_Label_Long, Avg_AUC_All_Available_Tests),
    Plot_Group = factor(
      Plot_Group,
      levels = c("HA dN only baseline", "dN-subset model")
    )
  )

model_bar_colors <- c(
  "HA dN only baseline" = "#0072B2",
  "dN-subset model" = "#D55E00"
)

p_top_width10_auc <- ggplot(
  top_width10_models,
  aes(
    x = Model_Label_Long,
    y = Avg_AUC_All_Available_Tests,
    fill = Plot_Group
  )
) +
  geom_hline(
    yintercept = 0.5,
    linetype = "dashed",
    linewidth = 0.4,
    color = "gray55"
  ) +
  geom_hline(
    yintercept = ha_baseline_auc,
    linetype = "dotted",
    linewidth = 0.7,
    color = "#0072B2"
  ) +
  geom_col(
    width = 0.72,
    color = "black",
    linewidth = 0.25,
    alpha = 0.92
  ) +
  geom_text(
    aes(
      label = sprintf("%.3f", Avg_AUC_All_Available_Tests)
    ),
    hjust = -0.12,
    size = 3.4
  ) +
  coord_flip(ylim = c(0.45, min(1.02, max(top_width10_models$Avg_AUC_All_Available_Tests, na.rm = TRUE) + 0.08))) +
  scale_fill_manual(values = model_bar_colors) +
  labs(
    title = "Top dN-subset models under rolling width 10",
    subtitle = "Models shown are representative dN-subset models outperforming HA dN only; dotted line indicates HA dN only AUROC",
    x = NULL,
    y = "Mean test AUROC across rolling splits",
    fill = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11),
    axis.title.x = element_text(face = "bold"),
    axis.text.y = element_text(size = 9.5),
    legend.position = "bottom"
  )

ggsave(
  file.path(fig_dir, "Supplementary_Top_RollingWidth10_dN_Subset_Models_AUROC.svg"),
  p_top_width10_auc,
  width = 10,
  height = 7.5,
  dpi = 400
)

ggsave(
  file.path(fig_dir, "Supplementary_Top_RollingWidth10_dN_Subset_Models_AUROC.pdf"),
  p_top_width10_auc,
  width = 10,
  height = 7.5
)

# ============================================================
# 6B. Extract fixed-effect odds ratios for the selected final model
# ============================================================

selected_model_name <- "1_Final_Rolling_Best_dN_HA_x_PB1_plus_PB2_NA"
selected_formula <- final_model_bank[[selected_model_name]]

extract_selected_model_or_one_window <- function(win_row) {
  train_periods <- strsplit(win_row$Train_Periods, ",")[[1]]
  test_period <- win_row$Test_Period
  
  raw_train <- raw_df %>% filter(as.character(Period) %in% train_periods)
  raw_test  <- raw_df %>% filter(as.character(Period) == test_period)
  
  if (nrow(raw_train) < 10 || nrow(raw_test) < 5) return(NULL)
  if (length(unique(raw_train$dom_binary)) < 2 || length(unique(raw_test$dom_binary)) < 2) return(NULL)
  
  scaler <- fit_scaler(raw_train, all_raw_features)
  train_df <- apply_scaler(raw_train, scaler) %>% recalc_weights()
  
  tryCatch({
    fit <- suppressWarnings(glmer(
      as.formula(selected_formula),
      data = train_df,
      family = binomial,
      weights = norm_wts,
      nAGQ = 0,
      control = glmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5)
      )
    ))
    
    beta <- fixef(fit)
    se <- sqrt(diag(vcov(fit)))
    
    tibble(
      Window_Scheme = win_row$Window_Scheme,
      Window_Type = win_row$Window_Type,
      Window_Size = win_row$Window_Size,
      Test_Index = win_row$Test_Index,
      Test_Period = test_period,
      Train_Range_Label = win_row$Train_Range_Label,
      Term = names(beta),
      Beta = as.numeric(beta),
      SE = as.numeric(se),
      CI95_Low_Beta = Beta - 1.96 * SE,
      CI95_High_Beta = Beta + 1.96 * SE,
      OR = exp(Beta),
      OR_CI95_Low = exp(CI95_Low_Beta),
      OR_CI95_High = exp(CI95_High_Beta),
      Error = NA_character_
    )
  }, error = function(e) {
    tibble(
      Window_Scheme = win_row$Window_Scheme,
      Window_Type = win_row$Window_Type,
      Window_Size = win_row$Window_Size,
      Test_Index = win_row$Test_Index,
      Test_Period = test_period,
      Train_Range_Label = win_row$Train_Range_Label,
      Term = NA_character_,
      Beta = NA_real_,
      SE = NA_real_,
      CI95_Low_Beta = NA_real_,
      CI95_High_Beta = NA_real_,
      OR = NA_real_,
      OR_CI95_Low = NA_real_,
      OR_CI95_High = NA_real_,
      Error = e$message
    )
  })
}

selected_model_or_by_split <- pmap_dfr(
  window_grid,
  function(...) extract_selected_model_or_one_window(tibble(...))
)

selected_model_or_summary <- selected_model_or_by_split %>%
  filter(is.na(Error), !is.na(Term), Term != "(Intercept)") %>%
  group_by(Term) %>%
  summarise(
    N_Splits = n(),
    Mean_Beta = mean(Beta, na.rm = TRUE),
    SD_Beta = sd(Beta, na.rm = TRUE),
    Median_Beta = median(Beta, na.rm = TRUE),
    Mean_OR = exp(Mean_Beta),
    Median_OR = exp(Median_Beta),
    OR_Low_Approx = exp(Mean_Beta - SD_Beta),
    OR_High_Approx = exp(Mean_Beta + SD_Beta),
    Mean_SE = mean(SE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Term_Label = case_when(
      Term == "HA" ~ "HA dN",
      Term == "PB1" ~ "PB1 dN",
      Term == "PB2" ~ "PB2 dN",
      Term == "NA_seg" ~ "NA dN",
      Term == "HA:PB1" ~ "HA dN × PB1 dN",
      TRUE ~ Term
    ),
    Term_Label = factor(
      Term_Label,
      levels = rev(c("HA dN", "PB1 dN", "PB2 dN", "NA dN", "HA dN × PB1 dN"))
    )
  )

write.csv(
  selected_model_or_by_split,
  file.path(result_dir, "Table_F6_Selected_Model_OR_By_Rolling_Split.csv"),
  row.names = FALSE
)

write.csv(
  selected_model_or_summary,
  file.path(result_dir, "Table_F7_Selected_Model_OR_Summary.csv"),
  row.names = FALSE
)

# ============================================================
# Figure: Selected model effect sizes as odds ratios
# ============================================================

or_path <- file.path(out_dir, "Table_F6_Selected_Model_OR_By_Rolling_Split.csv")

selected_or <- read_csv(or_path, show_col_types = FALSE) %>%
  filter(is.na(Error), !is.na(Term), Term != "(Intercept)") %>%
  mutate(
    Term_Label = case_when(
      Term == "HA" ~ "HA dN",
      Term == "PB1" ~ "PB1 dN",
      Term == "PB2" ~ "PB2 dN",
      Term == "NA_seg" ~ "NA dN",
      Term == "HA:PB1" ~ "HA dN × PB1 dN",
      TRUE ~ Term
    ),
    Term_Label = factor(
      Term_Label,
      levels = rev(c("HA dN", "PB1 dN", "PB2 dN", "NA dN", "HA dN × PB1 dN"))
    )
  )

or_summary <- selected_or %>%
  group_by(Term, Term_Label) %>%
  summarise(
    N_Splits = n(),
    Mean_Beta = mean(Beta, na.rm = TRUE),
    SD_Beta = sd(Beta, na.rm = TRUE),
    Mean_OR = exp(Mean_Beta),
    OR_Low = exp(Mean_Beta - SD_Beta),
    OR_High = exp(Mean_Beta + SD_Beta),
    .groups = "drop"
  )

p_or <- ggplot() +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.5,
    color = "gray45"
  ) +
  geom_jitter(
    data = selected_or,
    aes(x = OR, y = Term_Label),
    height = 0.10,
    width = 0,
    size = 2.2,
    alpha = 0.45,
    color = "#4D4D4D"
  ) +
  geom_errorbarh(
    data = or_summary,
    aes(
      y = Term_Label,
      xmin = OR_Low,
      xmax = OR_High
    ),
    height = 0.18,
    linewidth = 0.9,
    color = "#D55E00"
  ) +
  geom_point(
    data = or_summary,
    aes(x = Mean_OR, y = Term_Label),
    shape = 21,
    size = 5.0,
    fill = "#D55E00",
    color = "black",
    stroke = 0.7
  ) +
  geom_text(
    data = or_summary,
    aes(
      x = Mean_OR,
      y = Term_Label,
      label = sprintf("%.2f", Mean_OR)
    ),
    nudge_y = 0.24,
    size = 3.7,
    fontface = "bold"
  ) +
  scale_x_log10(
    breaks = c(0.25, 0.5, 1, 2, 4, 8, 16),
    labels = c("0.25", "0.5", "1", "2", "4", "8", "16")
  ) +
  labs(
    title = "Effect sizes of the selected rolling width 10 model",
    subtitle = "Points show split-specific odds ratios; orange points and bars show mean OR ± 1 SD on the log-odds scale",
    x = "Odds ratio per 1-SD increase in predictor",
    y = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11),
    axis.title.x = element_text(face = "bold"),
    axis.text.y = element_text(size = 11)
  )

ggsave(
  file.path(fig_dir, "Supplementary_Selected_Model_OR_Effect_Size.svg"),
  p_or,
  width = 8.5,
  height = 5.4,
  dpi = 400
)

ggsave(
  file.path(fig_dir, "Supplementary_Selected_Model_OR_Effect_Size.pdf"),
  p_or,
  width = 8.5,
  height = 5.4
)


# ============================================================
# Add p-values to existing OR summary without refitting models
# Uses the original Table_F6 generated by the old script
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

base_path <- "~/"

result_dir <- file.path(
  base_path,
  "3_Result/8_Final_Rolling_TrainTest_AUROC"
)

f6_path <- file.path(
  result_dir,
  "Table_F6_Selected_Model_OR_By_Rolling_Split.csv"
)

if (!file.exists(f6_path)) {
  stop("Existing Table_F6_Selected_Model_OR_By_Rolling_Split.csv not found: ", f6_path)
}

or_by_split <- readr::read_csv(
  f6_path,
  show_col_types = FALSE
)

# ------------------------------------------------------------
# Standardize column names if needed
# ------------------------------------------------------------

or_by_split <- or_by_split %>%
  mutate(
    Term_Label = case_when(
      Term == "HA" ~ "HA dN",
      Term == "PB1" ~ "PB1 dN",
      Term == "PB2" ~ "PB2 dN",
      Term == "NA_seg" ~ "NA dN",
      Term == "HA:PB1" ~ "HA dN × PB1 dN",
      TRUE ~ Term
    )
  )

# If Wald_Z / Wald_p are absent, compute them from Beta and SE
or_by_split <- or_by_split %>%
  mutate(
    Wald_Z = if ("Wald_Z" %in% names(.)) {
      Wald_Z
    } else {
      Beta / SE
    },
    Wald_p = if ("Wald_p" %in% names(.)) {
      Wald_p
    } else {
      2 * pnorm(abs(Wald_Z), lower.tail = FALSE)
    },
    OR = if ("OR" %in% names(.)) {
      OR
    } else {
      exp(Beta)
    },
    OR_CI95_Low = if ("OR_CI95_Low" %in% names(.)) {
      OR_CI95_Low
    } else {
      exp(Beta - 1.96 * SE)
    },
    OR_CI95_High = if ("OR_CI95_High" %in% names(.)) {
      OR_CI95_High
    } else {
      exp(Beta + 1.96 * SE)
    }
  )

# ------------------------------------------------------------
# Summary with p-values
# Important:
# Summary_OR is exp(mean Beta), matching the old figure logic
# ------------------------------------------------------------

or_summary_with_p <- or_by_split %>%
  group_by(Term, Term_Label) %>%
  summarise(
    N_Splits = n(),
    Test_Periods = paste(unique(Test_Period[order(Test_Index)]), collapse = ","),
    
    Mean_Beta = mean(Beta, na.rm = TRUE),
    Median_Beta = median(Beta, na.rm = TRUE),
    SD_Beta = sd(Beta, na.rm = TRUE),
    SE_Beta_Across_Splits = SD_Beta / sqrt(N_Splits),
    
    Summary_OR = exp(Mean_Beta),
    Median_OR = median(OR, na.rm = TRUE),
    Mean_OR = mean(OR, na.rm = TRUE),
    
    Summary_OR_SDlog_Low = exp(Mean_Beta - SD_Beta),
    Summary_OR_SDlog_High = exp(Mean_Beta + SD_Beta),
    
    Mean_SE_Model = mean(SE, na.rm = TRUE),
    Mean_Wald_Z = mean(Wald_Z, na.rm = TRUE),
    Median_Wald_p = median(Wald_p, na.rm = TRUE),
    
    # 1) Across-split t-test of beta against 0
    TTest_p_Beta_NE_0 = tryCatch(
      t.test(Beta, mu = 0)$p.value,
      error = function(e) NA_real_
    ),
    
    # 2) Inverse-variance weighted meta-analytic Wald p-value
    Meta_Beta_IVW = {
      ww <- 1 / (SE^2)
      sum(ww * Beta, na.rm = TRUE) / sum(ww, na.rm = TRUE)
    },
    Meta_SE_IVW = sqrt(1 / sum(1 / (SE^2), na.rm = TRUE)),
    Meta_Z_IVW = Meta_Beta_IVW / Meta_SE_IVW,
    Meta_p_Wald_IVW = 2 * pnorm(abs(Meta_Z_IVW), lower.tail = FALSE),
    
    Meta_OR_IVW = exp(Meta_Beta_IVW),
    Meta_OR_CI95_Low = exp(Meta_Beta_IVW - 1.96 * Meta_SE_IVW),
    Meta_OR_CI95_High = exp(Meta_Beta_IVW + 1.96 * Meta_SE_IVW),
    
    # 3) Fisher combined p-value across split-level Wald p-values
    Fisher_ChiSq = -2 * sum(log(Wald_p), na.rm = TRUE),
    Fisher_df = 2 * sum(!is.na(Wald_p)),
    Fisher_Combined_p_Wald = pchisq(
      Fisher_ChiSq,
      df = Fisher_df,
      lower.tail = FALSE
    ),
    
    N_Positive_Beta = sum(Beta > 0, na.rm = TRUE),
    N_Negative_Beta = sum(Beta < 0, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  mutate(
    Term_Label = factor(
      Term_Label,
      levels = c(
        "HA dN",
        "PB1 dN",
        "PB2 dN",
        "NA dN",
        "HA dN × PB1 dN"
      )
    )
  ) %>%
  arrange(Term_Label)

# ------------------------------------------------------------
# Save without overwriting the old F7 unless you want to
# ------------------------------------------------------------

readr::write_csv(
  or_by_split,
  file.path(result_dir, "Table_F6_Selected_Model_OR_By_Rolling_Split_with_Wald_p.csv")
)

readr::write_csv(
  or_summary_with_p,
  file.path(result_dir, "Table_F7_Selected_Model_OR_Summary_with_pvalues_FROM_EXISTING_F6.csv")
)

print(or_summary_with_p, n = Inf)
