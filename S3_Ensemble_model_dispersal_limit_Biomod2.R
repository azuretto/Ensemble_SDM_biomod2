library(sp)
library(sf)
library(raster)
library(rgeos) 
library(adehabitatHR)
library(snowfall)

# -------------------------------------------------------------------------
# Global Settings
# -------------------------------------------------------------------------

### Set working directory
# Define base path globally or ensure it is correct within the function
base_path = "D:/Biomod2_current_model/"
setwd(base_path)

## Load species occupancy data
# Verify that the file name and path are correct
dataSpecies <- read.csv("occ_table_L50_NA_scname_working_nameDec24_normalonly_part1.csv")

# Make species name list
# sp.n0 <- colnames(dataSpecies)
# sp.n0 <- sp.n0 [-1:-3]
sp.n0 <- read.csv("target_spp.csv")
sp.n0 <- sp.n0[,1]

# -------------------------------------------------------------------------
# Define Function for Parallel Processing (T1 ~ T38)
# -------------------------------------------------------------------------
MyBiomod_projection <- function(sp.n){ 
  
  # Load libraries within the parallel node
  require(raster)
  require(sp)
  require(rgeos)
  require(adehabitatHR)
  
  # Reset path
  ipath = "D:/Biomod2_current_model/"
  myRespName = sp.n  
  
  # -------------------------------------------------------
  # 1. Species Data & CRS Setup
  # -------------------------------------------------------
  # Exception handling: Check if species exists in data
  if(!(myRespName %in% colnames(dataSpecies))) return(NULL) 
  
  myRespXY <- dataSpecies[,c("X_WGS84","Y_WGS84", myRespName)]
  colnames(myRespXY) <- c("X_WGS84", "Y_WGS84", "occ")
  myRespXY <- myRespXY[myRespXY$occ == 1, ]
  myRespXY <- myRespXY[,-3]
  
  # Define CRS (WGS84) 
  sp::coordinates(myRespXY) <- ~ X_WGS84 + Y_WGS84
  sp::proj4string(myRespXY) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
  
  # -------------------------------------------------------
  # 2. Make MCP & Buffer
  # -------------------------------------------------------
  bgExt0 <- mcp(myRespXY, percent = 100)
  bgExt0 <- gBuffer(bgExt0, width = 0.03) # Approx. 3km (varies by latitude)
  bgExt1 <- gBuffer(bgExt0, width = 1)    # Approx. 111km (varies by latitude)
  
  # -------------------------------------------------------
  # 3. Load & Process Probability Rasters (T1 ~ T38)
  # -------------------------------------------------------
  ipath_sp = paste0(ipath, sp.n) 
  files.probability <- list.files(ipath_sp, pattern= 'ensemble[.]grd$', full.names = TRUE, recursive = TRUE)
  
  # [IMPORTANT] Check file order: Remove the first file (algorithm) then index
  if(length(files.probability) > 0) files.probability <- files.probability[-1] 
  
  # Ensure indices [c(2)], [c(3,4,5,6)] match the actual file order in the folder.
  HS.raster.T1  <- stack(files.probability[c(2)])
  HS.raster.T22 <- stack(files.probability[c(3,4,5,6)])
  HS.raster.T28 <- stack(files.probability[c(7,8,9,10)])
  HS.raster.T32 <- stack(files.probability[c(11,12,13,14)])
  HS.raster.T38 <- stack(files.probability[c(15,16,17,18)])
  
  # Subset (Remove CV raster logic: Keep TRUE, Remove FALSE)
  HS.raster.T1  <- HS.raster.T1[[which(c(TRUE, FALSE))]]
  
  # GCM Mean Calculation
  gcm_filter <- c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE)
  HS.raster.T22 <- calc(HS.raster.T22[[which(gcm_filter)]], fun = mean, na.rm = T)
  HS.raster.T28 <- calc(HS.raster.T28[[which(gcm_filter)]], fun = mean, na.rm = T)
  HS.raster.T32 <- calc(HS.raster.T32[[which(gcm_filter)]], fun = mean, na.rm = T)
  HS.raster.T38 <- calc(HS.raster.T38[[which(gcm_filter)]], fun = mean, na.rm = T)
  
  # -------------------------------------------------------
  # 4. Load & Process Binary Rasters (Calculate Threshold)
  # -------------------------------------------------------
  files.binary <- list.files(ipath_sp, pattern= 'TSSbin[.]grd$', full.names = TRUE, recursive = TRUE)
  if(length(files.binary) > 0) files.binary <- files.binary[-1]
  
  binary.T1 <- stack(files.binary[c(2)])
  binary.T1 <- binary.T1[[which(c(TRUE, FALSE))]]
  
  # Calculate Threshold based on T1
  binary.temp <- binary.T1
  binary.temp[binary.temp == 0] <- NA
  binary.probability <- HS.raster.T1 * binary.temp
  
  # Extract min value using cellStats 
  TSS.threshold <- cellStats(binary.probability, stat = 'min', na.rm = TRUE)
  
  # Reclassify Future Scenarios based on T1 Threshold
  rcl_mat <- c(-Inf, TSS.threshold, 0, TSS.threshold, +Inf, 1)
  
  binary.T22 <- reclassify(HS.raster.T22, rcl_mat)
  binary.T28 <- reclassify(HS.raster.T28, rcl_mat)
  binary.T32 <- reclassify(HS.raster.T32, rcl_mat)
  binary.T38 <- reclassify(HS.raster.T38, rcl_mat)
  
  # -------------------------------------------------------
  # 5. Stack & Clip Loop
  # -------------------------------------------------------
  # Stack 5 scenarios 
  HS.stack <- stack(HS.raster.T1, HS.raster.T22, HS.raster.T28, HS.raster.T32, HS.raster.T38)
  binary.stack <- stack(binary.T1, binary.T22, binary.T28, binary.T32, binary.T38)
  
  # Rasterize Buffers (Base: HS.raster.T1)
  bgExt0_r <- rasterize(bgExt0, HS.raster.T1, field = 1, background = NA) # No dispersal
  bgExt1_r <- rasterize(bgExt1, HS.raster.T1, field = 1, background = NA) # 1 degree dispersal
  
  # Name vector (5 scenarios)
  binaryTnames <- c("projT1", "projT22", "projT28", "projT32", "projT38")
  save_dir <- paste0(ipath_sp, "/GeoLimit_v2/")
  dir.create(save_dir, showWarnings = F)
  
  # Loop: 1 to 5 (Matches number of scenarios)
  for (n in 1:5){
    current_name <- binaryTnames[n]
    
    # (1) Original
    cliped <- binary.stack[[n]]
    writeRaster(cliped, paste0(save_dir, sp.n, "_", current_name, "_binary_all.tif"), overwrite = TRUE)
    
    cliped.HS <- HS.stack[[n]]
    writeRaster(cliped.HS, paste0(save_dir, sp.n, "_", current_name, "_prob_all.tif"), overwrite = TRUE)
    
    # (2) No dispersal (clipped by bgExt0_r)
    cliped <- cliped * bgExt0_r
    writeRaster(cliped, paste0(save_dir, sp.n, "_", current_name, "_binary_nd.tif"), overwrite = TRUE)
    
    cliped.HS <- cliped.HS * bgExt0_r
    writeRaster(cliped.HS, paste0(save_dir, sp.n, "_", current_name, "_prob_nd.tif"), overwrite = TRUE)
    
    # (3) 1 degree dispersal (clipped by bgExt1_r)
    # Caution: Reload variable from stack to avoid using previously clipped data
    cliped <- binary.stack[[n]] 
    cliped <- cliped * bgExt1_r
    writeRaster(cliped, paste0(save_dir, sp.n, "_", current_name, "_binary_1d.tif"), overwrite = TRUE)
    
    cliped.HS <- HS.stack[[n]]
    cliped.HS <- cliped.HS * bgExt1_r
    writeRaster(cliped.HS, paste0(save_dir, sp.n, "_", current_name, "_prob_1d.tif"), overwrite = TRUE)
  }
  
  return(paste(sp.n, "Finished"))
}

# -------------------------------------------------------------------------
# Parallel Execution (Snowfall)
# -------------------------------------------------------------------------

# Initialize snowfall (Adjust 'cpus' based on your hardware)
sfInit(parallel=TRUE, cpus=48) 

# Export libraries to slave nodes
sfLibrary('sp', character.only=TRUE)
sfLibrary('sf', character.only=TRUE)
sfLibrary('raster', character.only=TRUE)
sfLibrary('rgeos', character.only=TRUE)
sfLibrary('adehabitatHR', character.only=TRUE)

# Export variables to slave nodes
sfExport('dataSpecies') 
sfExport('sp.n0') 

# Run parallel processing
mySFModelsOut <- sfLapply(sp.n0, MyBiomod_projection)

# Stop snowfall cluster
sfStop()