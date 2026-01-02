# Load terra package
library(terra)

# Set base path (same as provided code)
ipath_base <- "D:/Biomod2_current_model/"
setwd(ipath_base)

# Load species list
sp.n0 <- read.csv("target_spp.csv")

sp.n0 <- sp.n0[,1]

# Initialize list to store results
all_results <- list()

# Analysis start message
cat("Starting distribution change analysis...\n")

# Iterate through each species
for (sp.n in sp.n0) {
  
  cat("Processing species:", sp.n, "\n")
  
  # Folder path where result files are stored
  output_dir <- file.path(ipath_base, sp.n, "individual_GCM_terra")
  
  # Load current (T1) distribution raster (1d scenario)
  current_raster_path <- file.path(output_dir, paste0(sp.n, "_projT1_binary_1d.tif"))
  
  # Check if file exists
  if (!file.exists(current_raster_path)) {
    cat("  Warning: Current (T1) raster file for", sp.n, "not found. Skipping to next species.\n")
    next
  }
  
  current_raster <- rast(current_raster_path)
  
  # Calculate only non-NA pixels (exclude background)
  # Create a mask to target only pixels where current_raster values are 0 or 1
  analysis_mask <- sum(current_raster, rast(current_raster)) # Force all pixels with values to be 0, others NA
  
  # Calculate current distribution range size (number of pixels)
  current_range_size_df <- global(current_raster, "sum", na.rm = TRUE)
  current_range_size <- current_range_size_df[1,1]
  
  # If there is no current distribution, skip to the next species
  if (current_range_size == 0) {
    cat("  Info:", sp.n, "has no currently modeled distribution area. Excluded from analysis.\n")
    next
  }
  
  # List of future periods and GCMs to analyze
  future_periods <- c("projT22", "projT28", "projT32", "projT38")
  gcm_suffixes <- c("a", "b", "c", "d")
  
  # Iterate through each future scenario
  for (period in future_periods) {
    for (gcm in gcm_suffixes) {
      
      # Generate future projection raster file path
      future_raster_path <- file.path(output_dir, paste0(sp.n, "_", period, gcm, "_binary_1d.tif"))
      
      # Check if file exists
      if (!file.exists(future_raster_path)) {
        cat("  Warning: File", future_raster_path, "not found.\n")
        next
      }
      
      future_raster <- rast(future_raster_path)
      
      # Calculate pixel-wise changes
      # Current: 1, Future: 0 => Disa (Loss)
      disa_raster <- (current_raster == 1 & future_raster == 0)
      disa_count <- global(disa_raster, "sum", na.rm = TRUE)[1,1]
      
      # Current: 1, Future: 1 => Stable1 (Maintained)
      stable1_raster <- (current_raster == 1 & future_raster == 1)
      stable1_count <- global(stable1_raster, "sum", na.rm = TRUE)[1,1]
      
      # Current: 0, Future: 1 => Gain
      gain_raster <- (current_raster == 0 & future_raster == 1)
      # Gain area is calculated within the total analysis area (mask)
      gain_raster_masked <- mask(gain_raster, analysis_mask)
      gain_count <- global(gain_raster_masked, "sum", na.rm = TRUE)[1,1]
      
      # Current: 0, Future: 0 => Stable0 (Non-habitat maintained)
      stable0_raster <- (current_raster == 0 & future_raster == 0)
      # Stable0 area is calculated within the total analysis area (mask)
      stable0_raster_masked <- mask(stable0_raster, analysis_mask)
      stable0_count <- global(stable0_raster_masked, "sum", na.rm = TRUE)[1,1]
      
      # Calculate derived variables
      # Prevent division by zero
      if (current_range_size > 0) {
        perc_loss <- (disa_count / current_range_size) * 100
        perc_gain <- (gain_count / current_range_size) * 100
      } else {
        perc_loss <- 0
        perc_gain <- 0
      }
      
      species_range_change <- perc_gain - perc_loss
      
      # Create result data frame
      result_row <- data.frame(
        Species = sp.n,
        Period = period,
        GCM = gcm,
        CurrentRangeSize = current_range_size,
        Disa = disa_count,
        Stable1 = stable1_count,
        Gain = gain_count,
        Stable0 = stable0_count,
        PercLoss = perc_loss,
        PercGain = perc_gain,
        SpeciesRangeChange = species_range_change
      )
      
      # Add to the all-results list
      all_results[[length(all_results) + 1]] <- result_row
    }
  }
}

# Combine all results into a single data frame
final_summary <- do.call(rbind, all_results)

# Save results to a CSV file
output_csv_path <- file.path(ipath_base, "distribution_change_summary_1d_v2.csv")
write.csv(final_summary, output_csv_path, row.names = FALSE, fileEncoding = "UTF-8")

cat("\nAnalysis complete. Results saved to:\n", output_csv_path, "\n")

