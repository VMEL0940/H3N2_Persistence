# ==============================================================================
# Script Name: 3_Renaissance_Calculator.R
# Purpose: Calculate the phylogenetic distance (Nonsynonymous and Synonymous 
#          substitutions) between reference vaccine strains and circulating 
#          strains across the posterior distribution of trees.
# ==============================================================================

library(treeio)
library(ape)
library(phytools)
library(future.apply)

# 1. Path Configuration
baseDir       <- "/mnt/shared610/H3N2_DomPred/1_Data/1_Phylogeny/2_MCC/5_3rd_CrossValidation_777"
outputBaseDir <- "/mnt/shared610/H3N2_DomPred/1_Data/2_RNSC/3_DistanceLog/4_3rd_CrossValidation_777/v2"
metaPath      <- "/mnt/shared610/H3N2_DomPred/1_Data/10_Final_MetaData_Submission/H3N2_Testdata_777strains_1_Metadata.csv"

# 2. Parallel Processing Setup
# Utilize future.apply for efficient multi-core computation
plan(multisession, workers = parallel::detectCores() - 1)

# 3. Load Metadata and Parameters
metadata <- read.csv(metaPath, stringsAsFactors = FALSE)
segments <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS")
target_periods <- 1:15  # Corresponds to vaccine groups Mos99 to CR23

# Reference Vaccine Strain List (Cleaning IDs for matching)
vaccineStrains_raw <- c(
  "EPI103320_A_Moscow_10_1999", "EPI358781_A_Fujian_411_2002", "EPI367109_A_California_7_2004",
  "EPI502253_A_Wisconsin_67_2005", "EPI577980_A_Brisbane_10_2007", "EPI577969_A_Perth_16_2009",
  "EPI417234_A_Victoria_361_2011", "EPI614441_A_Switzerland_9715293_2013", "EPI686117_A_Hong_Kong_15611_2015",
  "MG974447_A_Kansas_14_2017", "1592032_A_Hong_Kong_2671_2019", "1841681_A_Cambodia_e0826360_2020",
  "EPI2415906_A_Darwin_9_2021", "19194107_A_Thailand_8_2022", "19296516_A_Croatia_10136RV_2023"
)
vaccineStrains <- gsub("_NA", "", vaccineStrains_raw)

# FASTA Template Names
baseInFileNames <- c(
  "1_SEG_Mos99_76seq.fasta", "2_SEG_Fuj02_54seq.fasta", "3_SEG_Cal04_7seq.fasta", 
  "4_SEG_Wis05_17seq.fasta", "5_SEG_Bris07_51seq.fasta", "6_SEG_Prth09_54seq.fasta", 
  "7_SEG_Vic11_135seq.fasta", "8_SEG_Swtz13_25seq.fasta", "9_SEG_HK15_193seq.fasta", 
  "10_SEG_Kan17_24seq.fasta", "11_SEG_HK19_24seq.fasta", "12_SEG_Cam20_26seq.fasta", 
  "13_SEG_Dar21_58seq.fasta", "14_SEG_TH22_12seq.fasta", "15_SEG_CR23_21seq.fasta"
)

# 4. Main Calculation Loop
for (currentSegment in segments) {
  cat(paste0("\n>>> Analyzing Segment: ", currentSegment, " <<<\n"))
  
  segOutputDir <- file.path(outputBaseDir, currentSegment)
  if(!dir.exists(segOutputDir)) dir.create(segOutputDir, recursive = TRUE)
  
  # Load the RDS object containing posterior trees (pre-processed in Step 2)
  rds_path <- paste0("/mnt/shared610/H3N2_DomPred/1_Data/2_RNSC/2_Treelog/5_3rd_Verification_777/", currentSegment, "_resample.rds")
  if(!file.exists(rds_path)) {
    cat(paste0("RDS file not found: ", rds_path, "\n"))
    next
  }
  
  trees <- readRDS(rds_path)
  
  # Clean tip labels to match reference IDs
  for(t in 1:length(trees)) {
    trees[[t]]@phylo$tip.label <- gsub("_NA", "", trees[[t]]@phylo$tip.label)
  }
  
  inFileNames <- gsub("SEG", currentSegment, baseInFileNames)
  
  # Loop through each vaccine group/period
  for (v_period in target_periods) {
    vaccineStrain <- vaccineStrains[v_period]
    fastaName     <- inFileNames[v_period]
    fastaDir      <- paste0("/mnt/shared610/H3N2_DomPred/1_Data/1_Phylogeny/1_Alignment/5_3rd_Verification_777/", currentSegment, "_separate")
    fastaPath     <- file.path(fastaDir, fastaName)
    
    if(!file.exists(fastaPath)) next
    
    # Load sequence IDs to identify tips in the group
    fasta <- scan(fastaPath, what="", sep="\n", quiet=TRUE)
    seq_IDs <- gsub('>', '', fasta[grepl(">", fasta)])
    full_tip_names <- seq_IDs[trimws(seq_IDs) != trimws(vaccineStrain)]
    
    cat(paste0("  - Period ", v_period, " [", vaccineStrain, "] | Tips: ", length(full_tip_names), "\n"))
    
    # Calculate substitution distances for each circulating strain
    results_list <- lapply(full_tip_names, function(originalTip) {
      # Retrieve sampling date from metadata
      match_idx <- which(metadata$StrainName == originalTip)
      tipdate <- if(length(match_idx) > 0) metadata$date[match_idx[1]] else "Unknown"
      
      # Vectorized distance calculation across posterior trees
      dist_stats <- sapply(trees, function(tr) {
        phy <- tr@phylo
        idx_v <- which(phy$tip.label == vaccineStrain)
        idx_t <- which(phy$tip.label == originalTip)
        
        if(length(idx_v) == 0 || length(idx_t) == 0) return(c(N=NA, S=NA))
        
        # Identify nodes along the shortest path between vaccine and isolate
        path_nodes <- nodepath(phy, from = idx_v, to = idx_t)
        d_idx <- match(path_nodes, tr@data$node)
        d_idx <- d_idx[!is.na(d_idx)]
        
        # Accumulate Nonsynonymous (N) and Synonymous (S) counts
        return(c(N = sum(as.numeric(tr@data$N[d_idx]), na.rm = TRUE),
                 S = sum(as.numeric(tr@data$S[d_idx]), na.rm = TRUE)))
      })
      
      all_N <- dist_stats["N", ]
      all_S <- dist_stats["S", ]
      
      # Return summarized statistics for the strain
      return(data.frame(
        vaccineStrain = vaccineStrain,
        compareStrain = originalTip,
        tipDate = tipdate,
        mean_N = mean(all_N, na.rm=TRUE), median_N = median(all_N, na.rm=TRUE), 
        min_N = min(all_N, na.rm=TRUE), max_N = max(all_N, na.rm=TRUE),
        mean_S = mean(all_S, na.rm=TRUE), median_S = median(all_S, na.rm=TRUE), 
        min_S = min(all_S, na.rm=TRUE), max_S = max(all_S, na.rm=TRUE)
      ))
    })
    
    # Save the summarized results to a text file
    if (length(results_list) > 0) {
      finalDF <- do.call(rbind, results_list)
      outFileName <- paste0(currentSegment, "_G", v_period, "_summary_stats.txt")
      write.table(finalDF, file = file.path(segOutputDir, outFileName), sep = "\t", row.names = FALSE, quote = FALSE)
    }
  }
}

cat("\n### Renaissance Counting statistics calculated for all segments. ###\n")

# ==============================================================================
# 5. Post-Processing: Multi-Segment Data Consolidation
# Purpose: 1) Aggregate group-specific summaries into segment-wide logs.
#          2) Consolidate all segments into a single master dataset for 
#             downstream evolutionary analysis and machine learning.
# ==============================================================================

cat("\n>>> Initiating Multi-Segment Data Consolidation <<<\n")

# Define the formal segment order for the final publication-ready table
ordered_segments <- c("HA", "M", "NA", "NP", "NS", "PA", "PB1", "PB2")
segment_data_list <- list()

for (seg in segments) {
  segDir <- file.path(outputBaseDir, seg)
  
  # A. Aggregate individual period files into a single segment log ([Seg]_Gs.txt)
  # ---------------------------------------------------------------------------
  pattern <- paste0("^", seg, "_G[0-9]+_summary_stats\\.txt$")
  file_list <- list.files(segDir, pattern = pattern, full.names = TRUE)
  
  if (length(file_list) == 0) {
    warning(paste("No statistical logs found for segment:", seg))
    next
  }
  
  # Import and bind all group-level statistics
  combined_seg_df <- do.call(rbind, lapply(file_list, function(f) {
    read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  }))
  
  # Export the consolidated segment-wide file
  segMergeFileName <- paste0(seg, "_Gs.txt")
  write.table(combined_seg_df, file = file.path(segDir, segMergeFileName), 
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # B. Prepare for Master Merging (Extraction of substitution metrics)
  # ---------------------------------------------------------------------------
  # Retain unique strain identifiers and mean substitution counts
  temp_df <- combined_seg_df[, c("vaccineStrain", "compareStrain", "tipDate", "mean_N", "mean_S")]
  
  # Rename columns using the [Segment]_[Metric] convention
  colnames(temp_df)[4:5] <- c(paste0(seg, "_Nonsyn"), paste0(seg, "_Syn"))
  
  segment_data_list[[seg]] <- temp_df
}

# C. Construct the Final Consolidated Master Table
# ---------------------------------------------------------------------------
# Perform a horizontal merge (Inner Join) across all segments using key identifiers
master_df <- Reduce(function(x, y) {
  merge(x, y, by = c("vaccineStrain", "compareStrain", "tipDate"), all = TRUE)
}, segment_data_list[intersect(ordered_segments, names(segment_data_list))])

# D. Final Column Reordering and Formatting
# ---------------------------------------------------------------------------
# Requirement: [Identifiers] + [All Nonsynonymous Columns] + [All Synonymous Columns]
id_cols      <- c("vaccineStrain", "compareStrain", "tipDate")
nonsyn_cols  <- paste0(ordered_segments, "_Nonsyn")
syn_cols     <- paste0(ordered_segments, "_Syn")

# Verify column existence in the master dataframe before reindexing
final_cols <- c(id_cols, 
                nonsyn_cols[nonsyn_cols %in% colnames(master_df)], 
                syn_cols[syn_cols %in% colnames(master_df)])

master_df <- master_df[, final_cols]

# Export the final master dataset as a CSV for manuscript submission
finalOutPath <- file.path(outputBaseDir, "H3N2_Renaissance_Consolidated_Master.csv")
write.csv(master_df, file = finalOutPath, row.names = FALSE)

cat(paste0("\n### Post-Processing Complete ###\n"))
cat(paste0("Master consolidated dataset generated at: ", finalOutPath, "\n"))