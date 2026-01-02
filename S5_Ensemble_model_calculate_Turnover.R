# Load required packages
library(terra)

### 1. Basic Settings
# --------------------------------------------------------------------------
# Set base path
ipath_base <- "D:/Biomod2_current_model/"

# Load species list (used for the all-species group)
sp.n0 <- read.csv(file.path(ipath_base, "target_spp.csv"))
sp.n0 <- sp.n0[, 1]

# Create top-level output directory
output_dir <- file.path(ipath_base, "Community_Analysis")
dir.create(output_dir, showWarnings = FALSE)

# Define all scenarios to process
time_periods <- c("T1") 
future_periods <- c("T22", "T28", "T32", "T38")
gcm_suffixes <- c("a", "b", "c", "d")
future_scenarios <- as.vector(outer(future_periods, gcm_suffixes, paste0))
all_scenarios <- c(time_periods, future_scenarios)


### 2. Load Subgroup Info and Generate Lists
# --------------------------------------------------------------------------
cat("--- Loading and Preparing Species Groups ---\n")

# Load group information CSV file
group_file <- file.path(ipath_base, "BIOMOD_species_list_group_2025Jul14.csv")
group_data <- read.csv(group_file, header = TRUE) # Assuming the file has a header

# Explicitly name columns (adjust to actual column names if necessary)
# E.g., Col 1 = Species Name, Col 2 = Life Form, Col 3 = Endemism
colnames(group_data)[1:3] <- c("Species", "Life_form", "Endemism")

# Create subgroup names (e.g., "tree-endemic")
group_data$SubGroup <- paste(group_data$Life_form, group_data$Endemism, sep = "-")

# Create list of groups to analyze
group_list <- list()

# 1. Add all-species group
group_list[["All_Species"]] <- sp.n0

# 2. Add specific subgroups (6 groups)
unique_subgroups <- unique(group_data$SubGroup)
for (subgroup in unique_subgroups) {
  group_list[[subgroup]] <- group_data$Species[group_data$SubGroup == subgroup]
}

### 3. Start Main Analysis (Iterate by Group)
# ==========================================================================
for (group_name in names(group_list)) {
  
  cat(paste("\n\n========== Processing Group:", group_name, "==========\n"))
  
  # Species list for the current loop
  current_species_list <- group_list[[group_name]]
  
  # Create output directory for the specific group
  group_output_dir <- file.path(output_dir, group_name)
  dir.create(group_output_dir, showWarnings = FALSE)
  
  
  ### 3A. Calculate Species Richness
  # --------------------------------------------------------------------------
  cat(paste("--- Part 1: Calculating Species Richness for group", group_name, "---\n"))
  
  for (scen in all_scenarios) {
    cat("  Processing richness for scenario:", scen, "\n")
    
    # Use the species list for the current group
    richness_files <- sapply(current_species_list, function(sp) {
      file_path <- file.path(ipath_base, sp, "individual_GCM_terra", 
                             paste0(sp, "_proj", scen, "_binary_1d.tif"))
      if (file.exists(file_path)) return(file_path)
      return(NA)
    })
    
    richness_files <- richness_files[!is.na(richness_files)]
    
    if (length(richness_files) > 0) {
      community_stack <- rast(richness_files)
      richness_map <- sum(community_stack, na.rm = TRUE)
      # Save to the group folder
      output_filename <- file.path(group_output_dir, paste0("Richness_", scen, ".tif"))
      writeRaster(richness_map, output_filename, overwrite = TRUE)
    } else {
      cat("  No files found for scenario:", scen, "\n")
    }
  }
  
  
  ### 3B. Calculate Species Turnover and Change Rates
  # --------------------------------------------------------------------------
  cat(paste("\n--- Part 2: Calculating Turnover for group", group_name, "---\n"))
  
  # Use the species list for the current group
  current_files <- sapply(current_species_list, function(sp) {
    file_path <- file.path(ipath_base, sp, "individual_GCM_terra", 
                           paste0(sp, "_projT1_binary_1d.tif"))
    if (file.exists(file_path)) return(file_path)
    return(NA)
  })
  current_files <- current_files[!is.na(current_files)]
  
  if (length(current_files) > 0) {
    current_community_stack <- rast(current_files)
    current_richness_map <- rast(file.path(group_output_dir, "Richness_T1.tif"))
    
    for (scen in future_scenarios) {
      cat("  Processing turnover for scenario:", scen, "vs. T1\n")
      
      # Use the species list for the current group
      future_files <- sapply(current_species_list, function(sp) {
        file_path <- file.path(ipath_base, sp, "individual_GCM_terra", 
                               paste0(sp, "_proj", scen, "_binary_1d.tif"))
        if (file.exists(file_path)) return(file_path)
        return(NA)
      })
      future_files <- future_files[!is.na(future_files)]
      
      if (length(future_files) > 0) {
        future_community_stack <- rast(future_files)
        
        losses_map <- sum(current_community_stack == 1 & future_community_stack == 0, na.rm = TRUE)
        gains_map <- sum(current_community_stack == 0 & future_community_stack == 1, na.rm = TRUE)
        
        perc_loss_map <- 100 * losses_map / current_richness_map
        perc_loss_map[current_richness_map == 0] <- 0
        
        perc_gain_map <- 100 * gains_map / current_richness_map
        perc_gain_map[current_richness_map == 0] <- 0
        
        turnover_denominator <- current_richness_map + gains_map
        turnover_map <- 100 * (losses_map + gains_map) / turnover_denominator
        turnover_map[turnover_denominator == 0] <- 0
        
        # Save to the group folder
        writeRaster(perc_loss_map, file.path(group_output_dir, paste0("PercLoss_", scen, "_vs_T1.tif")), overwrite = TRUE)
        writeRaster(perc_gain_map, file.path(group_output_dir, paste0("PercGain_", scen, "_vs_T1.tif")), overwrite = TRUE)
        writeRaster(turnover_map, file.path(group_output_dir, paste0("Turnover_", scen, "_vs_T1.tif")), overwrite = TRUE)
        
      } else {
        cat("  No future files found for scenario:", scen, "\n")
      }
    }
  } else {
    cat("  Could not find baseline T1 files for this group. Skipping analysis.\n")
  }
}

cat("\n\n===== ALL GROUPS PROCESSED! =====\n")