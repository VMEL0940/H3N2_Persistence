library(foreach)
library(doParallel)
library(treeio)

home_path <- "~/"
meta_path <- "1_Data/2_Tree/Resample"

base_path <- paste0(home_path, meta_path)

setwd(base_path)

seg_list <- c("HA", "NA", "PB2", "PB1", "PA", "NP", "M", "NS")
log_file <- file.path(base_path, "conversion_log.txt")

# Remove previous log file
if (file.exists(log_file)) file.remove(log_file)

# 1. Set up parallel processing
num_cores <- min(length(seg_list), parallel::detectCores() - 2)
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Store the main library paths
main_lib_paths <- .libPaths()

cat(sprintf("\n[Start] Using %d cores. Log file: %s\n", num_cores, log_file))

# 2. Run conversion in parallel
foreach(seg = seg_list) %dopar% {
  .libPaths(main_lib_paths)
  library(treeio)
  
  # Log-writing function
  write_log <- function(msg) {
    cat(
      sprintf("[%s] %s: %s\n", Sys.time(), seg, msg),
      file = log_file,
      append = TRUE
    )
  }
  
  tree_path <- file.path(base_path, paste0(seg, "_resample.trees"))
  rds_path <- file.path(base_path, paste0(seg, "_resample.rds"))
  
  if (file.exists(tree_path)) {
    write_log("Loading data started (.trees parsing in progress).")
    
    # Main conversion step
    trees_obj <- treeio::read.beast(tree_path)
    write_log(sprintf("Loading completed (%d trees). Saving RDS file.", length(trees_obj)))
    
    saveRDS(trees_obj, file = rds_path)
    write_log("RDS file saved. Cleaning memory.")
    
    rm(trees_obj)
    gc()
  } else {
    write_log("Input file not found. Skipping this segment.")
  }
}

stopCluster(cl)

cat("\n[Done] RDS conversion completed for all segments.\n")