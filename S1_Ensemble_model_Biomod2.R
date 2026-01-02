##### #Ensemble modeling with biomod2 #########

###############################################################################
## Install versions of packages for biomod2 v3.5.1 compatibility
###############################################################################

# 1. Install 'remotes' package to allow installation of specific versions
if (!require("remotes")) install.packages("remotes")

# 2. Install archived spatial packages required by biomod2 3.5.x
if (!require("rgdal")) remotes::install_version("rgdal", version = "1.6-7", repos = "https://cloud.r-project.org")
if (!require("rgeos")) remotes::install_version("rgeos", version = "0.6-4", repos = "https://cloud.r-project.org")
if (!require("maptools")) remotes::install_version("maptools", version = "1.1-8", repos = "https://cloud.r-project.org")

# 3. Install other core dependencies available on CRAN
# Added essential packages from your original list that are commonly used
req_packages <- c("raster", "sp", "itertools", "foreach", "ggplot2", 
                  "reshape2", "nnet", "doParallel", "snowfall",
                  "ENMeval", "ade4", "Hmisc", "gam", "ecospat", 
                  "maxnet", "randomForest", "caret", "dplyr")

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
library(blockCV) 
library(sf)
library(dplyr)
library(tidyr)
#library(doParallel)

sessionInfo()

###############################################
###############################################
## parallel part
## Init snowfall
library(snowfall)
detectCores(logical = TRUE)
sfInit(parallel=TRUE, cpus=48)  # use only physical cores
#sfInit(parallel=FALSE)  ## sequential mode

## Export packages
sfLibrary('biomod2', character.only=TRUE)
sfLibrary('dplyr', character.only=TRUE)
sfLibrary('tidyr', character.only=TRUE)
sfLibrary('raster', character.only=TRUE)
sfLibrary('randomForest', character.only=TRUE)
sfLibrary('caret', character.only=TRUE)
sfLibrary('sf', character.only=TRUE)

## Export variables  (export cluster for all of variables)
# Export variable list
vars_to_export <- c("dataSpecies", "sp.n", "env0", "env1",
                    "env22a", "env22b", "env22c", "env22d",
                    "env28a", "env28b", "env28c", "env28d",
                    "env32a", "env32b", "env32c", "env32d",
                    "env38a", "env38b", "env38c", "env38d")

# Export variables 
#Use sfExportAll() to export all your workspace variables

sfExport(list = vars_to_export)
#sfExport('myBiomodOption')

# Do the run
mySFModelsOut <- sfLapply(sp.n, MyBiomod_pl)

## stop snowfall
sfStop( nostop=FALSE )
###############################################


###############################################
## serial run lapply

# myLapply_SFModelsOut <- lapply(sp.n, MyBiomod_pl)

###############################################

###################
##### BIOMOD2 #####
###################

####### Create biomod2 objects to use later #########
# preparation of species, environment data


## set directories
## need to change according to the number of scenarios and drvie setting

path_folder = "G:/Biomod2_work"  

setwd(path_folder) # Set working directory

# Save graphics defaults
#par.defaults <- par(no.readonly=TRUE)
#save(par.defaults, file="R.default.par.RData")

# Load species occupancy data
dataSpecies <- read.csv("occ_table_L50_NA_scname_working_nameDec24_normalonly_part1.csv")

# make species name list
sp.n <- colnames(dataSpecies)
sp.n <- sp.n [-1:-3]
sp.n

# Look at structure of dataSpecies dataframe
str(dataSpecies)

# Load in environmental covariates 
# Environment rasters from each institution (see Method section) should be saved in /current folder.
# Create stack of all environmental covariate layers
env1 <- list.files(path=paste0(path_folder,"/","current"), pattern = "tif$", full.names = TRUE)
myExpl <- raster::stack(env1)
myExpl$Soil_type <- as.factor(myExpl$Soil_type)
is.factor(myExpl$Soil_type) # check category variables as factor


# Set the options that you have chosen for the different model algorithms

#GLM : Generalized Linear Model (glm)
#GAM : Generalized Additive Model (gam, gam or bam, see BIOMOD_ModelingOptions for details on algorithm selection)
#GBM : Generalized Boosting Model or usually called Boosted Regression Trees (gbm)
#CTA: Classification Tree Analysis (rpart)
#ANN: Artificial Neural Network (nnet)
#SRE: Surface Range Envelop or usually called BIOCLIM
#FDA: Flexible Discriminant Analysis (fda)
#MARS: Multiple Adaptive Regression Splines (earth)
#RF: Random Forest (randomForest)

myBiomodOption <- BIOMOD_ModelingOptions(
  GLM = list( type = 'quadratic',
              interaction.level = 0,
              myFormula = NULL,
              test = 'AIC',
              family = binomial(link = 'logit'),
              mustart = 0.5,
              control = glm.control(epsilon = 1e-08, maxit = 50, trace = FALSE
              ) ),
  
  GBM = list( distribution = 'bernoulli',
              n.trees = 2500,
              interaction.depth = 7,
              n.minobsinnode = 5,
              shrinkage = 0.001,
              bag.fraction = 0.5,
              train.fraction = 1,
              cv.folds = 3,
              keep.data = FALSE,
              verbose = FALSE,
              perf.method = 'cv'),
  
  GAM = list( algo = 'GAM_mgcv',
              type = 's_smoother',
              k = -1,
              interaction.level = 0,
              myFormula = NULL,
              family = binomial(link = 'logit'),
              method = 'GCV.Cp',
              optimizer = c('outer','newton'),
              select = FALSE,
              knots = NULL,
              paraPen = NULL,
              control = list(nthreads = 1, irls.reg = 0, epsilon = 1e-07
                             , maxit = 200, trace = FALSE, mgcv.tol = 1e-07, mgcv.half = 15
                             , rank.tol = 1.49011611938477e-08
                             , nlm = list(ndigit=7, gradtol=1e-06, stepmax=2, steptol=1e-04, iterlim=200, check.analyticals=0)
                             , optim = list(factr=1e+07)
                             , newton = list(conv.tol=1e-06, maxNstep=5, maxSstep=2, maxHalf=30, use.svd=0)
                             , outerPIsteps = 0, idLinksBases = TRUE, scalePenalty = TRUE, keepData = FALSE
                             , edge.correct = FALSE) ),
  
  CTA = list( method = 'class',
              parms = 'default',
              cost = NULL,
              control = list(xval = 5, minbucket = 5, minsplit = 5, cp = 0.001
                             , maxdepth = 25) ),
  
  ANN = list( NbCV = 5,
              size = NULL,
              decay = NULL,
              rang = 0.1,
              maxit = 200),
  
  # FDA = list( method = 'mars',
  #              add_args = NULL),
  
  # MARS = list( type = 'simple',
  #              interaction.level = 0,
  #              myFormula = NULL,
  #              nk = NULL,
  #              penalty = 2,
  #              thresh = 0.001,
  #              nprune = NULL,
  #              pmethod = 'backward'),
  
  RF = list( do.classif = TRUE,
             ntree = 500,
             mtry = 'default',
             nodesize = 5,
             maxnodes = NULL))

###############################################
## biomod2 function ##

MyBiomod_pl <- function(sp.n) {
  
  myRespName = sp.n
  
  cat('\n',myRespName,'modeling...')
  
  
  init.proc.time = proc.time()
  
  ###############################################
  #test dummy 
  #myRespName <- "Abelia.spathulata.stenophylla"
  ###############################################
  
  myResp <- as.numeric(dataSpecies[,myRespName])
  myRespXY <- dataSpecies[,c("X_WGS84","Y_WGS84")]
  
  
  ###############################################
  ########## cross-validation using blockCV #############
  
  ##?? different number of records in resulting CV table (need to check issue)
  #check optimal range of block
  #sac <- spatialAutoRange(rasterLayer = myExpl,
  #                         sampleNumber = 5000,
  #                         border = NULL,
  #                         showPlots = TRUE,
  #                         plotVariograms = FALSE,
  #                         doParallel = FALSE)
  
  #PA <- dataSpecies
  #PA_data <- st_as_sf(PA, coords = c("X_WGS84", "Y_WGS84"), crs = crs(myExpl))
  
  # spatial blocking by specified range with random assignment
  #sblock <- spatialBlock(speciesData = PA_data,
  #                       species = 'Sp_FagusDtest',
  #                       rasterLayer = myExpl,
  #                       theRange = 114000, # size of the blocks
  #                       k = 3, # change based on 
  #                       selection = "systematic",
  #                       iteration = 100, # find evenly dispersed folds
  #                       biomod2Format = TRUE,
  #                       xOffset = 0, # shift the blocks horizontally
  #                       yOffset = 0)
  # 2. Defining the folds for DataSplitTable 
  # note that biomodTable should be used here not folds
  # DataSplitTable <- sblock$biomodTable 
  ###############################################################
  
  
  ############# Choosing biomod options and running models ##############
  # Create object (myBiomodData) to contain all the previous objects within it, 
  # formatted correctly
  myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                       expl.var = myExpl,
                                       resp.xy = myRespXY,
                                       resp.name = myRespName,
                                       PA.nb.rep = 0)
  
  ###############################################
  ## DataSplit_CV
  
  ###############################################
  #test dummy
  # DataSplitTable1 <- BIOMOD_cv(myBiomodData, k = 2, do.full.models = FALSE,
  #           stratified.cv = TRUE, stratify = "x", balance = "presences")
  
  # DataSplitTable.y <- BIOMOD_cv(myBiomodData, k = 2, do.full.models = FALSE,
  #                               stratified.cv = TRUE, stratify = "y", balance = "presences")
  # DataSplitTable.y[1, 1] <- c("NA")
  # head(DataSplitTable.y)
  ###############################################
  
  myResp.sum <- sum(myResp)
  
  if (myResp.sum >= 200) {
    DataSplitTable.y <- BIOMOD_cv(myBiomodData, k = 5, do.full.models = FALSE,
                                  stratified.cv = TRUE, stratify = "y", balance = "presences")
    DataSplitTable.x <- BIOMOD_cv(myBiomodData, k = 5, do.full.models = FALSE,
                                  stratified.cv = TRUE, stratify = "x", balance = "presences")
    DataSplitTable <- cbind(DataSplitTable.y, DataSplitTable.x)
    colnames(DataSplitTable) <- c("RUN1", "RUN2", "RUN3", "RUN4", "RUN5", 
                                  "RUN6", "RUN7", "RUN8", "RUN9", "RUN10")
    head(DataSplitTable)
    
  }
  else{
    DataSplitTable.y <- BIOMOD_cv(myBiomodData, k = 2, do.full.models = FALSE,
                                  stratified.cv = TRUE, stratify = "y", balance = "presences")
    
    DataSplitTable.x <- BIOMOD_cv(myBiomodData, k = 2, do.full.models = FALSE,
                                  stratified.cv = TRUE, stratify = "x", balance = "presences")
    
    DataSplit.all <- list(as.data.frame(DataSplitTable.x[, 1]),
                          as.data.frame(DataSplitTable.x[, 1]),
                          as.data.frame(DataSplitTable.x[, 2]),
                          as.data.frame(DataSplitTable.x[, 2]),
                          as.data.frame(DataSplitTable.y[, 1]), 
                          as.data.frame(DataSplitTable.y[, 1]),
                          as.data.frame(DataSplitTable.y[, 1]),
                          as.data.frame(DataSplitTable.y[, 2]),
                          as.data.frame(DataSplitTable.y[, 2]),
                          as.data.frame(DataSplitTable.y[, 2]))
    
    cv.resample <- function (x) {
      cv.group.total <- sum(x[ , 1])
      cv.p = 0.7
      random.sub <- round(cv.group.total*cv.p) ## desired total proportion of group
      nr <- nrow(x); nc <- ncol(x)
      ina <- unlist(x) %in% TRUE  ## logical vector, TRUE corresponds to TRUE positions
      ina[sample(which(ina), random.sub)] <- FALSE
      x[matrix(ina, nr=nr, nc=nc)]<- FALSE ## using matrix indexing
      x <- as.data.frame(x)
    }
    
    DataSplit.Ran <- lapply(DataSplit.all,  cv.resample)
    summary(DataSplit.Ran[[10]])
    DataSplitTable <- cbind(DataSplit.Ran[[1]], DataSplit.Ran[[2]],
                            DataSplit.Ran[[3]], DataSplit.Ran[[4]],
                            DataSplit.Ran[[5]], DataSplit.Ran[[6]],
                            DataSplit.Ran[[7]], DataSplit.Ran[[8]],
                            DataSplit.Ran[[9]], DataSplit.Ran[[10]])
    colnames(DataSplitTable) <- c("RUN1", "RUN2", "RUN3", "RUN4", "RUN5", 
                                  "RUN6", "RUN7", "RUN8", "RUN9", "RUN10")
    DataSplitTable <- data.matrix(DataSplitTable)
    DataSplitTable <- DataSplitTable>0
  }
  
  # head(DataSplitTable1)
  # summary(DataSplitTable1)
  # 
  # head(DataSplitTable)
  # summary(DataSplitTable)
  
  
  # Run the biomod2 models that you have chosen on the data provided.
  beforeModeling.proc.time = proc.time()
  
  # 'GLM', 'GAM', 'CTA', 'ANN', 'FDA', 'GBM', 'MARS', 'RF'
  # 'GLM', 'GAM', 'MARS', 'FDA' regression models
  # 'GBM', 'CTA', 'ANN', 'RF' ML
  # 'GLM', 'GAM', 'ANN', 'GBM', 'RF', 'CTA'
  
  ###############################################
  ###############################################
  ##Tuning Option
  ###############################################
  
  # cl<-makeCluster(4);registerDoParallel(cl)
  # Biomod.tuning <- BIOMOD_tuning(myBiomodData, models = c('GLM', 'GAM', 'ANN', 'GBM', 'RF', 'CTA'),
  #                               metric='TSS', type.GLM = 'quadratic')
  #stopCluster(cl)
  
  ##############################################
  ##############################################
  ##############################################
  # with tuning option
  # cl<-makeCluster(4);registerDoParallel(cl)
  # 
  # myBiomodModelOut <- BIOMOD_Modeling(
  #   myBiomodData,
  #   models = c('GLM', 'GAM', 'ANN', 'GBM', 'RF', 'CTA'),
  #   models.options = Biomod.tuning$models.options,
  #   #  NbRunEval = 3,
  #   DataSplitTable=DataSplitTable, #block cross-validation
  #   #  DataSplit = 70,
  #   Prevalence = 0.5,
  #   VarImport = 3,
  #   models.eval.meth = c('ROC', 'TSS', 'KAPPA'),
  #   SaveObj = TRUE,
  #   #  rescal.all.models = FALSE,
  #   #  do.full.models = FALSE,
  #   modeling.id = paste("sdm0", sep=""))
  # 
  # stopCluster(cl)
  
  ###############################################
  
  myBiomodModelOut <- BIOMOD_Modeling(
    myBiomodData,
    models = c('GLM', 'GAM', 'ANN', 'GBM', 'RF', 'CTA'),
    models.options = myBiomodOption,
    #  NbRunEval = 3,
    DataSplitTable=DataSplitTable, #block cross-validation
    #  DataSplit = 70,
    Prevalence = 0.5,
    VarImport = 3,
    models.eval.meth = c('ROC', 'TSS', 'KAPPA'),
    SaveObj = TRUE,
    #  rescal.all.models = FALSE, # turn on when GBM gives error
    #  do.full.models = FALSE,
    modeling.id = paste("sdm0", sep=""))
  
  #The available evaluations methods are :
  #??ROC?? : Relative Operating Characteristic
  #??KAPPA?? : Cohen's Kappa (Heidke skill score)
  #??TSS?? : True skill statistic (Hanssen and Kuipers discriminant, Peirce's skill score)
  #??FAR?? : False alarm ratio
  #??SR?? : Success ratio
  #??ACCURANCY?? : Accuracy (fraction correct)
  #??BIAS?? : Bias score (frequency bias)
  #??POD?? : Probability of detection (hit rate)
  #??CSI?? : Critical success index (threat score)
  #??ETS?? : Equitable threat score (Gilbert skill score)
  
  ################# Evaluating models ################
  
  # Get evaluation statistics for each model
  myBiomodModelEval <- get_evaluations(myBiomodModelOut)
  
  # Display the model evaluation statistics
  #myBiomodModelEval
  #myBiomodModelEval["ROC", "Testing.data",,,]
  #myBiomodModelEval["TSS", "Testing.data" ,,,]
  #myBiomodModelEval["KAPPA", "Testing.data",,,]
  
  # select TSS threshold for ensemble 
  TSS.threshold <- myBiomodModelEval["TSS", "Testing.data" ,,,]
  TSS.threshold[is.na(TSS.threshold)] <- 0
  TSS.threshold <- as.data.frame(TSS.threshold)
  TSS.max <- max(TSS.threshold) # chekcer for ensemble
  TSS.threshold <- data.frame(TSS = c(t(TSS.threshold)))
  TSS.threshold <- max(min(top_n(TSS.threshold, 20, TSS)), 0.3) ## select larget threshold 0.3
  TSS.threshold
  
  #save csv table of evaluation statistics
  write.table(myBiomodModelEval, file = sprintf("EVAL_%s.csv",myRespName), row.names=FALSE)
  
  # evaluation statistics summary by algorithm 
  eval.mean <- apply(myBiomodModelEval, c(1,2,3), mean)
  eval.sd <- apply(myBiomodModelEval, c(1,2,3), sd)
  write.table(eval.mean, file = sprintf("EVALmean_%s.csv",myRespName), row.names=FALSE)
  write.table(eval.sd, file = sprintf("EVALsd_%s.csv",myRespName), row.names=FALSE)
  
  # Display relative importance of each environmental covariate to each model 
  VarImp <- get_variables_importance(myBiomodModelOut)
  
  #save csv table of variable importance
  write.table(VarImp, file = sprintf("VarImp_%s.csv",myRespName), row.names=TRUE)
  
  ###############################################

  beforeEnsemble.proc.time = proc.time()
  
  ###############################################
  ## check for minimum TSS threshod before ensemble / prevent ensemble error
  ## minimum TSS threshod check, prevent ensemble error
  ## if TSS_max < 0.3  no ensemble  
  ##         TSS_max >= 0.3 do ensemble 
  
  if (TSS.max >= 0.3){
    
    ###############################################
    ######### Project and ensemble models #########
    
    # Build ensemble model from selected models.
    try({
      myBiomodEM <- BIOMOD_EnsembleModeling(myBiomodModelOut,
                                            chosen.models = 'all',
                                            em.by = 'all',
                                            eval.metric = c('TSS'),  # TSS 0.3
                                            eval.metric.quality.threshold = TSS.threshold, ## testing threshold, 0.2-0.4 moderate,  >0.4 good / 0.7
                                            prob.mean = FALSE, #all model mean
                                            prob.cv = FALSE,
                                            prob.ci = TRUE,
                                            prob.ci.alpha = 0.05,
                                            prob.median = TRUE,
                                            committee.averaging = FALSE,
                                            prob.mean.weight = TRUE, #weighted model mean
                                            prob.mean.weight.decay = 'proportional',
                                            VarImport = 1)
    })
    # Uncomment to plot selected models, if you want
    # plot(myBiomodProj)
    
    ######################### Project ensemble model###############################
    beforeProjection.proc.time = proc.time()
    
    #Project and plot predictions from the ensemble model 
    #EMplot_current <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEMprojection, output = myBiomodProj)
    
    
    # Project selected models for later building ensemble model
    try({
      
      myBiomodProj <- BIOMOD_Projection(
        modeling.output = myBiomodModelOut,
        new.env = myExpl,
        proj.name = 'current',
        selected.models = "all",
        binary.meth = 'TSS',
        compress = 'xz',
        build.clamping.mask = F,
        output.format = '.grd')
      
      EMplot_current <- BIOMOD_EnsembleForecasting(myBiomodEM, projection.output = myBiomodProj, compress = 'xz')
      
     
      ########################### Evaluating ensemble models #############################
      
      # Get evaluation statistics for each model
      myBiomodModelEval_EM <- get_evaluations(myBiomodEM)
      #myBiomodModelEval_EM
      
      #save csv table of evaluation statistics
      write.csv(myBiomodModelEval_EM, file = sprintf("%s_emEVAL.csv",myRespName), row.names=TRUE)
      
      # Display relative importance of each environmental covariate to each model 
      emVarImp <- get_variables_importance(myBiomodEM)
      
      #save csv table of variable importance
      write.csv(emVarImp, file = sprintf("%s_emVarImp.csv",myRespName), row.names=TRUE)
    })
    

    ###############################################
    # time checker
    #path_data = paste0(path_folder,"/", myRespName)
    path_save = paste0(path_folder, "/", myRespName)
    
    end.proc.time = proc.time()
    PROC.TIME = list(init.proc.time
                     , beforeModeling.proc.time
                     , beforeEnsemble.proc.time
                     , beforeProjection.proc.time
                     , end.proc.time)
    names(PROC.TIME) = c("init", "beforeModeling", "beforeEnsemble", "beforeProjection", "end")
    
    total.time <- end.proc.time[3] - init.proc.time[3]
    total.time
    write.csv(total.time, file = sprintf("%s_totaltime.csv",myRespName), row.names=TRUE)
    
  }
  else {
    
    ###############################################
    # time checker
    #path_data = paste0(path_folder,"/", myRespName)
    path_save = paste0(path_folder, "/", myRespName)
    
    end.proc.time = proc.time()
    PROC.TIME = list(init.proc.time
                     , beforeModeling.proc.time
                     , beforeEnsemble.proc.time
                     , end.proc.time)
    
    names(PROC.TIME) = c("init", "beforeModeling", "beforeEnsemble", "end")
    
    total.time <- end.proc.time[3] - init.proc.time[3]
    total.time
    write.csv(total.time, file = sprintf("%s_totaltime.csv",myRespName), row.names=TRUE)
    
    ###############################################  
  } ## else part (ensemble check)
  
} ## end of biomod2 function



#Create an unique temp_file folder so we can delete it at the end of the Tile loop and avoid crashing because of HD full of temp files with little raster stacks
dir.create(paste('Temp_files_', SpecieName, sep=""), showWarnings = F, recursive = T)
rasterOptions(tmpdir = paste('Temp_files_', SpecieName, sep=""))


#remove all temporary rasters (100% necessary when huge predictor rasters)
unlink(paste('Temp_files_',SpecieName, sep=""), recursive = T, force = FALSE)
unlink(paste("./",SpecieName,"/proj_", sep=""), recursive = T, force = FALSE)

