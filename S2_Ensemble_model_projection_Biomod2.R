##### #Ensemble modeling with biomod2 (part 2. projection) #########

###############################################################################
## Install versions of packages for biomod2 v3.5.1 compatibility
###############################################################################

# 1. Install 'remotes' package
if (!require("remotes")) install.packages("remotes")

# 2. Install archived spatial packages required by biomod2 3.5.x
if (!require("rgdal")) remotes::install_version("rgdal", version = "1.6-7", repos = "https://cloud.r-project.org")
if (!require("rgeos")) remotes::install_version("rgeos", version = "0.6-4", repos = "https://cloud.r-project.org")
if (!require("maptools")) remotes::install_version("maptools", version = "1.1-8", repos = "https://cloud.r-project.org")

# 3. Install other core dependencies available on CRAN
req_packages <- c("raster", "sp", "itertools", "foreach", "ggplot2", 
                  "reshape2", "nnet", "doParallel", "snowfall",
                  "ENMeval", "ade4", "Hmisc", "gam", "ecospat", 
                  "maxnet", "randomForest", "caret", "dplyr", "sf") 

new_packages <- req_packages[!(req_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# 4. Install biomod2 version 3.5.1 specifically
if (packageVersion("biomod2") != "3.5.1") {
  remotes::install_version("biomod2", version = "3.5.1", repos = "https://cloud.r-project.org")
}

# Check the installed version
print(packageVersion("biomod2"))


###############################################
## load packages
library(biomod2)
library(randomForest)
library(caret)
library(raster)
library(dplyr)
library(sf) 

sessionInfo()


####### Create biomod2 objects to use later #########
# preparation of species, environment data

## need to change according to the number of scenarios
## set directories
path_folder = "D:/Biomod2_current_model"

setwd(path_folder) # Set working directory

# Load species occupancy data
dataSpecies <- read.csv("occ_table_L50_NA_scname_working_nameDec24_normalonly_part1.csv")

dataSpecies_xy <- dataSpecies[,c(1:3)]
dataSpecies <- dataSpecies[,c(4:1336)]

colnames(dataSpecies)
dataSpecies <- dataSpecies[,order(colnames(dataSpecies))]

dataSpecies <- cbind(dataSpecies_xy, dataSpecies)

#check order
colnames(dataSpecies)[378]

# make species name list
sp.n <- colnames(dataSpecies) 
sp.n <- sp.n [-1:-3]
head(sp.n) 

# check name (test)
match(c("Rubus.palmatus.coptophyllus"), sp.n)

# Look at structure of dataSpecies dataframe
str(dataSpecies)

################################
## Load in current & future environmental covariates 

# fixed environmental covariates
env_fix <- list.files(path=paste0(path_folder,"/","env","/Fix"), pattern = "tif$", full.names = TRUE)
myExpl_fix <- raster::stack(env_fix)

# T1 1981-2000 (current)
env1 <- list.files(path=paste0(path_folder, "/", "env", "/T1"), pattern = "tif$", full.names = TRUE)
myExpl1 <- raster::stack(append(env1, env_fix))
myExpl1$Soil_type <- as.factor(myExpl1$Soil_type)

#######################################
# T2 2031-2050 RCP2.6
env22a <- list.files(path=paste0(path_folder, "/", "env", "/T2R25/GFDL"), pattern = "tif$", full.names = TRUE)
myExpl22a <- raster::stack(append(env22a, env_fix))
myExpl22a$Soil_type <- as.factor(myExpl22a$Soil_type)

env22b <- list.files(path=paste0(path_folder, "/", "env", "/T2R25/HadGEM"), pattern = "tif$", full.names = TRUE)
myExpl22b <- raster::stack(append(env22b, env_fix))
myExpl22b$Soil_type <- as.factor(myExpl22b$Soil_type)

env22c <- list.files(path=paste0(path_folder, "/", "env", "/T2R25/MIROC"), pattern = "tif$", full.names = TRUE)
myExpl22c <- raster::stack(append(env22c, env_fix))
myExpl22c$Soil_type <- as.factor(myExpl22c$Soil_type)

env22d <- list.files(path=paste0(path_folder, "/", "env", "/T2R25/MRI"), pattern = "tif$", full.names = TRUE)
myExpl22d <- raster::stack(append(env22d, env_fix))
myExpl22d$Soil_type <- as.factor(myExpl22d$Soil_type)

#######################################
# T2 2031-2050 RCP8.5
env28a <- list.files(path=paste0(path_folder, "/", "env", "/T2R85/GFDL"), pattern = "tif$", full.names = TRUE)
myExpl28a <- raster::stack(append(env28a, env_fix))
myExpl28a$Soil_type <- as.factor(myExpl28a$Soil_type)

env28b <- list.files(path=paste0(path_folder, "/", "env", "/T2R85/HadGEM"), pattern = "tif$", full.names = TRUE)
myExpl28b <- raster::stack(append(env28b, env_fix))
myExpl28b$Soil_type <- as.factor(myExpl28b$Soil_type)

env28c <- list.files(path=paste0(path_folder, "/", "env", "/T2R85/MIROC"), pattern = "tif$", full.names = TRUE)
myExpl28c <- raster::stack(append(env28c, env_fix))
myExpl28c$Soil_type <- as.factor(myExpl28c$Soil_type)

env28d <- list.files(path=paste0(path_folder, "/", "env", "/T2R85/MRI"), pattern = "tif$", full.names = TRUE)
myExpl28d <- raster::stack(append(env28d, env_fix))
myExpl28d$Soil_type <- as.factor(myExpl28d$Soil_type)

#######################################
# T3 2081-2100 RCP2.6
env32a <- list.files(path=paste0(path_folder, "/", "env", "/T3R25/GFDL"), pattern = "tif$", full.names = TRUE)
myExpl32a <- raster::stack(append(env32a, env_fix))
myExpl32a$Soil_type <- as.factor(myExpl32a$Soil_type)

env32b <- list.files(path=paste0(path_folder, "/", "env", "/T3R25/HadGEM"), pattern = "tif$", full.names = TRUE)
myExpl32b <- raster::stack(append(env32b, env_fix))
myExpl32b$Soil_type <- as.factor(myExpl32b$Soil_type)

env32c <- list.files(path=paste0(path_folder, "/", "env", "/T3R25/MIROC"), pattern = "tif$", full.names = TRUE)
myExpl32c <- raster::stack(append(env32c, env_fix))
myExpl32c$Soil_type <- as.factor(myExpl32c$Soil_type)

env32d <- list.files(path=paste0(path_folder, "/", "env", "/T3R25/MRI"), pattern = "tif$", full.names = TRUE)
myExpl32d <- raster::stack(append(env32d, env_fix))
myExpl32d$Soil_type <- as.factor(myExpl32d$Soil_type)

#######################################
# T3 2081-2100 RCP8.5
env38a <- list.files(path=paste0(path_folder, "/", "env", "/T3R85/GFDL"), pattern = "tif$", full.names = TRUE)
myExpl38a <- raster::stack(append(env38a, env_fix))
myExpl38a$Soil_type <- as.factor(myExpl38a$Soil_type)

env38b <- list.files(path=paste0(path_folder, "/", "env", "/T3R85/HadGEM"), pattern = "tif$", full.names = TRUE)
myExpl38b <- raster::stack(append(env38b, env_fix))
myExpl38b$Soil_type <- as.factor(myExpl38b$Soil_type)

env38c <- list.files(path=paste0(path_folder, "/", "env", "/T3R85/MIROC"), pattern = "tif$", full.names = TRUE)
myExpl38c <- raster::stack(append(env38c, env_fix))
myExpl38c$Soil_type <- as.factor(myExpl38c$Soil_type)
names(myExpl38c) = names(myExpl38b)

env38d <- list.files(path=paste0(path_folder, "/", "env", "/T3R85/MRI"), pattern = "tif$", full.names = TRUE)
myExpl38d <- raster::stack(append(env38d, env_fix))
myExpl38d$Soil_type <- as.factor(myExpl38d$Soil_type)
names(myExpl38d) = names(myExpl38c)


envName <- c("T1", "T22a", "T22b", "T22c", "T22d", "T28a", "T28b", "T28c", "T28d",
             "T32a", "T32b", "T32c", "T32d", "T38a", "T38b", "T38c", "T38d")

envTime <- c(myExpl1, myExpl22a, myExpl22b, myExpl22c, myExpl22d, myExpl28a, myExpl28b, myExpl28c, myExpl28d,
             myExpl32a, myExpl32b, myExpl32c, myExpl32d, myExpl38a, myExpl38b, myExpl38c, myExpl38d)


###############################################
## biomod2 function ##

MyBiomod_projection <- function(sp.n) {
  
  # parallel working check
  library(dplyr)
  library(raster)
  
  myRespName = sp.n
  
  # create Temp folder 
  dir.create(paste('D:/Temp_files_', myRespName, sep=""), showWarnings = F, recursive = T)
  rasterOptions(tmpdir = paste('D:/Temp_files_', myRespName, sep=""))
  
  # delete user app local temp folder
  temp_del_folder <- tempdir()
  
  cat('\n',myRespName,'future projection CCCA...')
  
  # Load model output
  path_to_model <- paste0(myRespName, "/", myRespName, ".m0.models.out")
  
  if (!file.exists(path_to_model)) {
    return(paste("Error: Model file not found for", myRespName))
  }
  
  myBiomodModelOut_name <- load(path_to_model)  
  myBiomodModelOut <- get(myBiomodModelOut_name)
  
  # Get evaluation statistics
  myBiomodModelEval <- get_evaluations(myBiomodModelOut)
  
  # select TSS threshold for ensemble 
  TSS.threshold <- myBiomodModelEval["TSS", "Testing.data" ,,,]
  TSS.threshold[is.na(TSS.threshold)] <- 0
  TSS.threshold <- as.data.frame(TSS.threshold)
  
  # TSS threshold check
  TSS.temp <- data.frame(TSS = c(t(TSS.threshold)))
  TSS.top24 <- top_n(TSS.temp, 24, TSS) 
  TSS.val <- min(TSS.top24$TSS) # 상위 24개 중 최소값
  TSS.threshold <- max(TSS.val, 0.2) ## 최소 0.2 보장
  
  if (TSS.threshold > 0.6) {
    TSS.threshold <- 0.6
  }
  
  # Ensemble Modeling
  myBiomodEM <- NULL 
  
  try({
    myBiomodEM <- BIOMOD_EnsembleModeling(myBiomodModelOut,
                                          chosen.models = 'all',
                                          em.by = 'all',
                                          eval.metric = c('TSS'),
                                          eval.metric.quality.threshold = TSS.threshold,
                                          prob.mean = TRUE,
                                          prob.cv = TRUE,
                                          prob.ci = FALSE,
                                          prob.median = FALSE,
                                          committee.averaging = FALSE,
                                          prob.mean.weight = FALSE,
                                          VarImport = 1)
  })
  
  if (!is.null(myBiomodEM)) {
    for (i in 1:17) {
      try({
        bm.ef.all <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM
                                                , new.env = envTime[[i]]
                                                , output.format = ".grd"
                                                , proj.name = (paste0(envName[[i]]))
                                                , selected.models = "all"
                                                , binary.meth = 'TSS')
      })
    }
  } else {
    cat(paste("\n Ensemble Model Failed for:", myRespName))
  }
  
  ## delete temp folder
  unlink(paste("./",myRespName,"/Temp", sep=""), recursive = TRUE)
  
  #remove all temporary rasters
  unlink(paste('D:/Temp_files_',myRespName, sep=""), recursive = T, force = FALSE)
  unlink(paste("./",myRespName,"/proj_", sep=""), recursive = T, force = FALSE)
  
  return(paste("Done:", myRespName))
  
} # biomod2 iteration

###############################################
###############################################

## parallel calculation part
## Init snowfall
library(snowfall)

# Core number
n_cores <- 44 # physical cores
sfInit(parallel=TRUE, cpus=n_cores) 

## Export packages
sfLibrary('biomod2', character.only=TRUE)
sfLibrary('dplyr', character.only=TRUE)
sfLibrary('tidyr', character.only=TRUE)
sfLibrary('raster', character.only=TRUE)
sfLibrary('randomForest', character.only=TRUE)
sfLibrary('caret', character.only=TRUE)
sfLibrary('sf', character.only=TRUE) 

## Export variables
sfExport('dataSpecies')
sfExport('myExpl1')
sfExport('myExpl22a', 'myExpl22b', 'myExpl22c', 'myExpl22d')
sfExport('myExpl28a', 'myExpl28b', 'myExpl28c', 'myExpl28d')
sfExport('myExpl32a', 'myExpl32b', 'myExpl32c', 'myExpl32d')
sfExport('myExpl38a', 'myExpl38b', 'myExpl38c', 'myExpl38d')

sfExport('sp.n')
sfExport('envName')
sfExport('envTime')
sfExport('MyBiomod_projection')

# Run parallel
mySFModelsOut <- sfLapply(sp.n, MyBiomod_projection)

## stop snowfall
sfStop()

###############################################
###############################################
