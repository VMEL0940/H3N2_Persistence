# ==============================================================================
# Script Name: 1_Dataset_Separation.R
# Purpose: This script partitions global H3N2 influenza sequence data into 
#          subsets based on specific vaccine-referenced time periods.
#          It ensures the reference vaccine strain is positioned at the top
#          of each FASTA file for subsequent alignment and phylogenetic analysis.
# ==============================================================================

# Load required libraries
library(Biostrings)
library(dplyr)

# 1. Path Configuration
# Note: Update these paths to match your local or server environment
base_path <- "~/"
meta_path <- "1_Data/0_Metadata/H3N2_777strains_Metadata.csv"

# Set Working Directory
setwd(base_path)

# 2. Data Loading
# Load metadata containing strain names and corresponding vaccine codes
metadata <- read.csv(meta_path, stringsAsFactors = FALSE)

# 3. Define Analysis Parameters
# Define the 8 influenza segments and the chronological vaccine reference codes
segments <- c("PB2", "PB1", "PA", "NP", "HA", "NA", "M", "NS")
v_codes <- c("Mos99", "Fuj02", "Cal04", "Wis05", "Bris07", "Prth09", 
             "Vic11", "Swtz13", "HK15", "Kan17", "HK19", "Cam20", "Dar21", "TH22", "CR23")

# 4. Processing Loop by Segment
for (seg in segments) {
  cat("\n--- Processing Segment:", seg, "---\n")
  
  # Identify the target master FASTA file for the current segment
  pattern <- paste0("_", seg, "_work3")
  all_files <- list.files(base_path, full.names = TRUE)
  target_file <- all_files[grepl(pattern, all_files)][1]
  
  if (is.na(target_file)) {
    cat("Warning: Could not find master file for segment", seg, "- skipping...\n")
    next
  }
  
  # Load Master FASTA file
  fasta_all <- readDNAStringSet(target_file)
  
  # Create output directory for the current segment
  out_dir <- file.path(base_path, paste0(seg, "_separate"))
  if (!dir.exists(out_dir)) dir.create(out_dir)
  
  # 5. Subset Generation by Vaccine Code
  for (i in seq_along(v_codes)) {
    current_code <- v_codes[i]
    
    # Filter metadata for the current vaccine group
    current_meta <- metadata %>% filter(vaccine_code == current_code)
    
    # Extract unique strain names associated with this group
    target_strains <- unique(current_meta$StrainName)
    target_strains <- target_strains[target_strains != "" & !is.na(target_strains)]
    
    # Step A: Filter sequences from the master FASTA
    matched_fasta <- fasta_all[names(fasta_all) %in% target_strains]
    
    if (length(matched_fasta) > 0) {
      # Step B: Reorder - Place the Reference Vaccine Strain at the first position
      # This is crucial for bioinformatics tools that require a reference at top
      v_strain_name <- unique(current_meta$vaccineStrain)[1] 
      
      if (v_strain_name %in% names(matched_fasta)) {
        # Separate the reference strain
        v_seq <- matched_fasta[v_strain_name]
        # Separate the rest
        other_seqs <- matched_fasta[names(matched_fasta) != v_strain_name]
        # Re-combine with the reference strain at index 1
        final_fasta <- c(v_seq, other_seqs)
      } else {
        # Keep original subset if reference strain name is not found in FASTA headers
        final_fasta <- matched_fasta
        cat("Notice: Reference strain", v_strain_name, "not found in FASTA for", current_code, "\n")
      }
      
      # Step C: Save the partitioned FASTA file
      num_seq <- length(final_fasta)
      # Naming convention: [Index]_[Segment]_[VaccineCode]_[SequenceCount].fasta
      file_name <- paste0(i, "_", seg, "_", current_code, "_", num_seq, "seq.fasta")
      full_out_path <- file.path(out_dir, file_name)
      
      writeXStringSet(final_fasta, full_out_path)
      cat("Saved:", file_name, "(Reference strain prioritized)\n")
    }
  }
}

cat("\n### Dataset partitioning completed for all segments. ###\n")