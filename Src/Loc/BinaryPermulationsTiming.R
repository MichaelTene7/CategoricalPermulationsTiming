#This script is being used to time-trial multiple methods of binary permulations on similar datasets. 
clusterRun = F
clusterRun = T
if(clusterRun){.libPaths("/share/ceph/wym219group/shared/libraries/R4")} #add path to custom libraries to searched locations


# r = filePrefix                      This is a prefix used to organize and separate files by analysis run. Always required. 
# m = mainTreeFilename.txt or .rds    This sets the location of the maintrees file
# n = numberOfPermulations            This is the number of permulations to run in the script 
# i = runInstanceValue                This is used to generate unique filenames for each instance of the script. Used in parrallelization. 
# l = relaxationValue <int, 0-1>      This is a value between 0 and 1, the percentage off of exact a permulation is allowed to be. Increases runspeed, but reduces permulation accuracy. Recommend 0.1 if more than 3 categories
# t=rootSpeciesName                    This is the name of the root species, if not using REFERENCE(human)


library(RERconverge)
library(tools)
library(data.table)
source("Src/Reu/cmdArgImport.R")
source("Src/Reu/convertForegroundVector.R")
source("Src/Reu/fast_bin_perm.R")
source("Src/Reu/permPValCorReport.R")
source("Src/Reu/GetPermsBinaryFudged.R")
source("Src/Reu/RelaxedRejectionPermFuncs2.R")
source("Src/Reu/CategoricalPermulationsParallelFunctions.R")
source("Src/Reu/scriptsForTestingPermulationFunctions/simBinPhenoSSM_fromMasterTree.r")
source("Src/Reu/convertPathsToPhenVector.R")
source("Src/Reu/convertBinaryToCategoricalCorrelationsFormat.R")
source("Src/Reu/scriptsForTestingPermulationFunctions/neededPermulationFunctions.r")
source("Src/Reu/RERConvergeFunctions.R")

args = c('r=setupTest', 'm=../RunRER/data/zoonomiaAllMammalsTrees.rds', 'v=F', 't=vs_HLornAna3', 'n=5', 'l=0.05')
args = c('r=emilyPhen', 'm=Data/emilyMultiphylo.rds', 'v=F', 't=Lophuromys_woosnami_LSUMZ37793', 'n=2', 'l=0.05', 'c=T')





#-------------------------
# --- Standard start-up code ---
#-------------------------
if(clusterRun)args = commandArgs(trailingOnly = TRUE)
{  # Bracket used for collapsing purposes
  #File Prefix
  if(!is.na(cmdArgImport('r'))){                                                #This cmdArgImport script is a way to import arguments from the command line. 
    filePrefix = cmdArgImport('r')
  }else{
    stop("THIS IS AN ISSUE MESSAGE; SPECIFY FILE PREFIX")
  }
  
  #  Output Directory 
  if(!dir.exists("Output")){                                      #Make output directory if it does not exist
    dir.create("Output")
  }
  outputFolderNameNoSlash = paste("Output/",filePrefix, sep = "") #Set the prefix sub directory
  if(!dir.exists(outputFolderNameNoSlash)){                       #create that directory if it does not exist
    dir.create(outputFolderNameNoSlash)
  }
  outputFolderName = paste("Output/",filePrefix,"/", sep = "")
  
  #  Force update argument
  forceUpdate = FALSE
  if(!is.na(cmdArgImport('v'))){                                 #Import if update being forced with argument 
    forceUpdate = cmdArgImport('v')
    forceUpdate = as.logical(forceUpdate)
  }else{
    message("Force update not specified, not forcing update")
  }
}


#-------------------------
#-- Arugment Imports ---
#-------------------------
#defaults: 
mainTreesLocation = "../RunRERBinaryMT/Data/zoonomiaAllMammalsTrees.rds"  #this is the standard location on Sol 
permulationAmount = 100
runInstanceValue = NULL
useRelaxation = FALSE
relaxationValue = NULL
PathsObject = NA
rootSpeciesValue = "REFERENCE"
speciesFilter = NULL
runDaniel = T
runFudged = T
runCategorical = T
runEmily = T


#Settings without an arugment but centralized here 
min.sp = 10
min.pos = 2
modelType = "ER"
rootProbability = "auto"

#MainTrees Location
if(!is.na(cmdArgImport('m'))){
  mainTreesLocation = cmdArgImport('m')
}else{
  message("No maintrees arg, using default")
}


#Number of permulations
if(!is.na(cmdArgImport('n'))){
  permulationAmount = cmdArgImport('n')
  permulationAmount = as.numeric(permulationAmount)
}else{
  paste("Number of permulations not specified, using 100")
}
#Run instance value
if(!is.na(cmdArgImport('i'))){
  runInstanceValue = cmdArgImport('i')
}else{
  paste("This script does not have a run instance value")
}
#Relaxation
if(!all(is.na(cmdArgImport('l')))){
  useRelaxation = TRUE
  relaxationValue = cmdArgImport('l')
  relaxationValue = as.numeric(relaxationValue)
}else{
  paste("Not relaxing permulations.")
}
#Root species
if(!is.na(cmdArgImport('t'))){
  rootSpeciesValue = cmdArgImport('t')
}else{
  message("No root species specified, using 'REFERENCE'")
}

#Run Daniel
if(!is.na(cmdArgImport('d'))){
  runDaniel = cmdArgImport('d')
}else{
  message("Default, running Daniel")
}

#Run fudged
if(!is.na(cmdArgImport('f'))){
  runFudged = cmdArgImport('f')
}else{
  message("default, running fudged")
}

#Run categorical
if(!is.na(cmdArgImport('c'))){
  runCategorical = cmdArgImport('c')
}else{
  message("default, running categorical")
}

#Run emily
if(!is.na(cmdArgImport('e'))){
  runEmily = cmdArgImport('e')
}else{
  message("default, running emily")
}

#-------------------------
#-- Import data  ---
#-------------------------

#maintrees
if(file_ext(mainTreesLocation) == "rds"){
  if(!exists("mainTrees")){mainTrees = readRDS(mainTreesLocation)}
}else{
  if(!exists("mainTrees")){mainTrees = readTrees(mainTreesLocation)} 
}

#speciesFilter 
speciesFilterFileName = paste(outputFolderName, filePrefix, "SpeciesFilter.rds",sep="") #Make the name of the location a pre-made filter would have to test for it
if (file.exists(paste(speciesFilterFileName))){                  
  speciesFilter = readRDS(speciesFilterFileName)                       #if so, use it 
  message("Pre-made filter found, using pre-made filter.")
}else{                                                    
  message("No speciesFilter arg, using all species")                           #if not, use no filter
}
if(all(is.null(speciesFilter))){ # If the species filter is meant to be empty, that is, all of the species should be used
  speciesFilter = mainTrees$masterTree$tip.label #include all of the species on the tree in the filter 
}


#PathsFile
PathsFileFilename = paste(outputFolderName, filePrefix, "PathsFile.rds", sep="")
if(file.exists(PathsFileFilename)){
  PathsObject = readRDS(PathsFileFilename)
}else{
  if(runEmily){
    stop("THIS IS AN ISSUE MESSAGE, GENERATE A PATHSFILE")  
  }else{
    message("No Paths found, not running emily, not required.")
  }
}

#phenotypeVector
phenotypeVectorFilename = paste(outputFolderName, filePrefix, "PhenotypeVector.rds", sep="")
if(file.exists(phenotypeVectorFilename)){
  phenotypeVector = readRDS(phenotypeVectorFilename)
}else{
  if(!all(is.na(PathsObject))){
    message("Generating phenotype vector from paths")
    phenotypeVector = convertPathsToPhenVector(PathsObject)
    phenotypeVector= phenotypeVector[which(names(phenotypeVector) %in% mainTrees$masterTree$tip.label)]
  }else{
    stop("THIS IS AN ISSUE MESSAGE, GENERATE A PHENTOTYPEVECTOR")
  }
}
phenotypeVector = convertForegroundVector(phenotypeVector)

#RERFile
RERFileName = paste(outputFolderName, filePrefix, "RERFile.rds", sep= "")
if(!file.exists(paste(RERFileName)) | forceUpdate){
  stop("THIS IS AN ISSUE MESSAGE,  GENERATE AN RER FILE")
}else{
  RERObject = readRDS(RERFileName)
}

#CorrelationFile 
correlationFileName = paste(outputFolderName, filePrefix, "CorrelationFile.rds", sep= "")
if(!file.exists(paste(correlationFileName)) | forceUpdate){
  stop("THIS IS AN ISSUE MESSAGE,  GENERATE AN Basic Correlation FILE")
}else{
  CorrelationObject = readRDS(correlationFileName)
}

#CombinedCorrelationFile 
combinedCorrelationFileName = paste(outputFolderName, filePrefix, "CombinedCategoricalCorrelationFile.rds", sep= "")
if(!file.exists(paste(combinedCorrelationFileName)) | forceUpdate){
  if(runCategorical){
    message("Generating categorical correlaiton format from binary correlation.")
    combinedCorrelationObject = convertBinaryToCategoricalCorrelationsFormat(CorrelationObject)
  }
}else{
  combinedCorrelationObject = readRDS(combinedCorrelationFileName)
}



masterTree = mainTrees$masterTree
rootNode = which(masterTree$tip.label %in% rootSpeciesValue)
rootedMasterTree = multi2di(masterTree)

timeListNameSet = c("SimulationTimes", "PathTimes", "CorrelationTimes", "PermulationAmount", "PvalueCalculationTime")

#-------------------------
#-- Functions  ---
#-------------------------
{
  # ---- Timer ----
  
  functionTimer = function(timedFunction, repeatNumber, count = T, report = F, ...){
    functionToTime = get(timedFunction, mode = "function", inherits = T)
    TimingStart = Sys.time()
    durationVector = NULL
    
    for(i in 1:repeatNumber){
      instanceStart = Sys.time()
      out = functionToTime(...)
      
      if(count){message(i)}
      if(report){message(out)}
      instanceEnd = Sys.time()
      instanceDuration = instanceEnd - instanceStart
      durationVector = append(durationVector, instanceDuration)
      
    }
    if(report){cat(durationVector)}
    return(durationVector)
  }
  
  
  # ---- Record species in foreground ----
  
  makeForegroundSpeciesStorage = function(namePrefix){
    permulatedForegroundEdgeList = list()
    tipInForegroundBinary= c()
    storage = list(permulatedForegroundEdgeList, tipInForegroundBinary)
    storageNames = c("PermulatedForegroundEdgeList", "TipInForegroundBinary")
    names(storage) = paste0(namePrefix, storageNames)
    storage
  }
  
  recordForegroundSpecies = function(permTree, storage, categorical = F){
    permFG_list = storage[[1]]
    inFG = storage[[2]]
    if(!categorical){
      fgEdges<-permTree$edge[which(permTree$edge.length==1),2]
      permFgs<-permTree$tip.label[fgEdges]
    }else{
      permFgs = names(permTree$tips[which(permTree$tips ==2)])
    }
    
    permFG_list[[length(permFG_list)+1]]<-permFgs
    inFG<-cbind(inFG, unlist(sapply(masterTree$tip.label, function(x) if(x %in% permFgs){1} else{0})))
    
    storage[[1]] = permFG_list
    storage[[2]] = inFG
    storage
  }
  
  
  
  # ----- Categorical --------
  
  categoricalPermulationPipline = function(message = F, permulationAmount){
    if(message){cat("Simulating Phenotype \n")}
    simulationTime = suppressWarnings(
      system.time({
        permulationsData = categoricalPermulations(mainTrees, phenotypeVector, rm = modelType, rp = rootProbability, ntrees = permulationAmount, percent_relax = currentRelaxation)
      })
    )
    if(message){cat("Total Simulation time: ", simulationTime["elapsed"], "\n")}
    categoricalSimulationTimes <<- append(categoricalSimulationTimes, simulationTime["elapsed"])
    
    #record the foregrounds
    categoricalTrees = permulationsData$trees
    categoricalForegroundStoragePrefix = paste0("categorical", currentRelaxation)
    categoricalTempForegroundStorage = makeForegroundSpeciesStorage(categoricalForegroundStoragePrefix)
    for(i in 1:length(categoricalTrees)){
      categoricalTempForegroundStorage = recordForegroundSpecies(categoricalTrees[[i]], categoricalTempForegroundStorage, categorical = T)
    }
    categoricalTempForegroundStorage <<- categoricalTempForegroundStorage
    
    
    if(!is.null(permulationsData)){
      
      if(message){cat("Calculating Paths \n")}
      pathsTime = system.time({}) #nothing here, no operation required 
      if(message){cat("Paths time: ", pathsTime["elapsed"], "\n")}
      categoricalPathTimes <<- append(categoricalPathTimes, pathsTime["elapsed"])
      
      
      if(message){cat("Performing Correlations \n")}
      correlationTime = system.time({permulatedCorrelations = CategoricalPermulationGetCor(combinedCorrelationObject, permulationsData$trees, phenotypeVector, mainTrees, RERObject, report=T)})  
      if(message){cat("\n Total Correlation time: ", correlationTime["elapsed"], "\n")}
      categoricalCorrelationTimes <<- append(categoricalCorrelationTimes, correlationTime["elapsed"])
      
      
      return(permulatedCorrelations)
    }else{
      message("Perm failed, skipping")
      return(NULL)
    }
  }
  
  # ---- Daniel ---- 
  #Convert format function
  convertPermulationFormat = function(permulationCorList, RERObj = RERObject, permulationNum = permulationAmount){
    permulationCorList
    permPvals = data.frame(matrix(ncol = permulationNum, nrow = nrow(RERObj)))
    rownames(permPvals) = rownames(RERObj)
    permRhovals = data.frame(matrix(ncol = permulationNum, nrow = nrow(RERObj)))
    rownames(permRhovals) = rownames(RERObj)
    permStatvals = data.frame(matrix(ncol = permulationNum, nrow = nrow(RERObj)))
    rownames(permStatvals) = rownames(RERObj)
    for (i in 1:length(permulationCorList)) {
      permPvals[, i] = permulationCorList[[i]]$P
      permRhovals[, i] = permulationCorList[[i]]$Rho
      permStatvals[, i] = sign(permulationCorList[[i]]$Rho) * -log10(permulationCorList[[i]]$P)
    }
    output = vector("list", 3)
    output[[1]] = permPvals
    output[[2]] = permRhovals
    output[[3]] = permStatvals
    names(output) = c("corP", "corRho", "corStat")
    output
  }
  
  danielSinglePermulation = function(message = F){
    if(message){cat("Simulating Phenotype \n")}
    simulationTime = suppressWarnings(
      system.time({
        permulatedForeground = fastSimBinPhenoVec(tree=rootedMasterTree, phenvec=phenotypeVector, internal=internalNumberValue)
        })
      )
    if(message){cat("Simulation time: ", simulationTime["elapsed"], "\n")}
    danielSimulationTimes <<- append(danielSimulationTimes, simulationTime["elapsed"])
    
    permulatedTree = NULL
    if(message){cat("Creating Tree \n")}
    treeTime = system.time({tryCatch({permulatedTree = foreground2Tree(permulatedForeground, mainTrees, plotTree=F, clade="all", transition="bidirectional", useSpecies=speciesFilter)}, error = function(cond){permulatedTree = NULL; message("Original error:"); message(cond)})}) #generate a tree using that foregound
    if(message){cat("Tree time: ", treeTime["elapsed"], "\n")}
    
    if(!is.null(permulatedTree)){
      
      danielForegroundStorage <<- recordForegroundSpecies(permulatedTree, danielForegroundStorage)
      
      if(message){cat("Calculating Paths \n")}
      pathsTime = system.time({permulatedPaths = tree2Paths(permulatedTree, mainTrees, binarize=T, useSpecies=speciesFilter)})                                                    #generate a path from that tree
      if(message){cat("Paths time: ", pathsTime["elapsed"], "\n")}
      danielPathTimes <<- append(danielPathTimes, pathsTime["elapsed"])
      
      if(message){cat("Performing Correlations \n")}
      correlationTime = system.time({permulatedCorrelations = correlateWithBinaryPhenotype(RERObject, permulatedPaths, min.sp=min.sp, min.pos = min.pos)})                                                 #Use that path to get a coreelation of the null phenotype to genes (this is the outbut of a Get PermsBinary run)
      if(message){cat("\n Correlation time: ", correlationTime["elapsed"], "\n")}
      danielCorrelationTimes <<- append(danielCorrelationTimes, correlationTime["elapsed"])
      
      
      return(permulatedCorrelations)
    }else{
      pathsTime["elapsed"] = NA
      danielPathTimes <<- append(danielPathTimes, pathsTime["elapsed"])
      correlationTime["elapsed"] = NA
      danielCorrelationTimes <<- append(danielCorrelationTimes, correlationTime["elapsed"])
      message("Perm failed, skipping")
      return(NULL)
    }
  }
  
  danielPermulationPipeline= function(permulationAmount, message = F){
    
    danielSimulationTimes = NULL
    danielPathTimes = NULL
    danielCorrelationTimes = NULL
    danielPValueTime = NULL
    
    {cat("Generating Correlation Data \n")}
    dataTime = system.time({
      correlationList = list()
      for(i in 1:permulationAmount){                                                  #Repeat for the number of permulations
        singlePermCorrelation = danielSinglePermulation(message = message)
        correlationList = append(correlationList, list(singlePermCorrelation))        #add it to a growing list of the dataframes outputted from CorrelateWithBinaryPhenotype
        if(message){message("Completed permulation: ", i)}                                         #report completed the permulation
      }
      convertedPermulations = convertPermulationFormat(correlationList)
    })
    cat("\n -------------- \n Permulation Data time: ", dataTime["elapsed"], "\n ----------- \n")
    
    {cat("Calculating P values \n")}
    danielPValueTime = system.time({  permulationPValues = permPValCorReport(CorrelationObject, convertedPermulations)})
    cat("Pvalue time: ", danielPValueTime["elapsed"], "\n")
    danielPValueTime <<- danielPValueTime["elapsed"]
    
    return(permulationPValues)
    
  }
  
  
  
  
  # --- Emily ----
  
  
  emilySinglePermulation = function(message = F, useMidpoint = T){
    if(message){cat("Simulating Phenotype \n")}
    simulationTime = suppressWarnings(
      system.time({
        if(useMidpoint){
          permulatedTree = simBinPhenoSSM_fromMasterTree(tree=masterTree, trees=mainTrees, fg_vec=foregroundSpecies, pathvec=PathsObject)
        }else{
          permulatedTree = simBinPhenoSSM_fromMasterTree_noMidpointRoot(tree=masterTree, trees=mainTrees, fg_vec=foregroundSpecies, pathvec=PathsObject)
        }
        })
      )
    if(message){cat("Simulation time: ", simulationTime["elapsed"], "\n")}
    emilySimulationTimes <<- append(emilySimulationTimes, simulationTime["elapsed"])
    
    if(message){cat("Calculating Paths \n")}
    permulatedPaths = NULL
    pathsTime = system.time({ 
      tryCatch({
        permulatedPaths = tree2Paths(permulatedTree, mainTrees, binarize=T, useSpecies=speciesFilter)
      }, error = function(cond){permulatedPaths = NULL; message("Original error:"); message(cond)})
    })                                                    #generate a path from that tree
    if(message){cat("Paths time: ", pathsTime["elapsed"], "\n")}
    emilyPathTimes <<- append(emilyPathTimes, pathsTime["elapsed"])
    
    if(!is.null(permulatedPaths)){

      if(useMidpoint){
        emilyMidpointForegroundStorage <<- recordForegroundSpecies(permulatedTree, emilyMidpointForegroundStorage)
      }else{
        emilyNoMidpointForegroundStorage <<- recordForegroundSpecies(permulatedTree, emilyNoMidpointForegroundStorage)
      }
      
      if(message){cat("Performing Correlations \n")}
      correlationTime = system.time({permulatedCorrelations = correlateWithBinaryPhenotype(RERObject, permulatedPaths, min.sp=min.sp, min.pos = min.pos)})                                                 #Use that path to get a coreelation of the null phenotype to genes (this is the outbut of a Get PermsBinary run)
      if(message){cat("\n Correlation time: ", correlationTime["elapsed"], "\n")}
      emilyCorrelationTimes <<- append(emilyCorrelationTimes, correlationTime["elapsed"])
      
      
      return(permulatedCorrelations)
    }else{
      correlationTime = system.time({1+1})
      correlationTime["elapsed"] = NA
      emilyCorrelationTimes <<- append(emilyCorrelationTimes, correlationTime["elapsed"])
      message("Perm failed, skipping")
      return(NULL)
    }
  }
  
  emilyPermulationPipeline= function(permulationAmount, message = F, midpoint = T){
    
    emilySimulationTimes = NULL
    emilyPathTimes = NULL
    emilyCorrelationTimes = NULL
    emilyPValueTime = NULL
    
    {cat("Generating Correlation Data \n")}
    dataTime = system.time({
      correlationList = list()
      for(i in 1:permulationAmount){                                                  #Repeat for the number of permulations
        singlePermCorrelation = emilySinglePermulation(message = message, useMidpoint = midpoint)
        correlationList = append(correlationList, list(singlePermCorrelation))        #add it to a growing list of the dataframes outputted from CorrelateWithBinaryPhenotype
        if(message){message("Completed permulation: ", i)}                                         #report completed the permulation
      }
      convertedPermulations = convertPermulationFormat(correlationList)
      #convertedPermulations = correlationList
    })
    cat("\n -------------- \n Permulation Data time: ", dataTime["elapsed"], "\n ----------- \n")
    
    {cat("Calculating P values \n")}
    emilyPValueTime = system.time({  permulationPValues = permPValCorReport(CorrelationObject, convertedPermulations)})
    cat("Pvalue time: ", emilyPValueTime["elapsed"], "\n")
    emilyPValueTime <<- emilyPValueTime["elapsed"]
    
    return(permulationPValues)
    
  }
}
#-------------------------
#-- Basic Permulations ---
#-------------------------

# Basic GetPermsBinary is not being used due to the need for a sisterList 
# If desired, this can be re-added.
# getPermsBinary()


#-------------------------
#-- Categorical Permulations ---
#-------------------------

#create Timing storage objects
categoricalSimulationTimes = NULL
categoricalPathTimes = NULL
categoricalCorrelationTimes = NULL
categoricalPValueTime = NULL
categoricalTempForegroundStorage = NULL

if(runCategorical){
  for(i in 1:length(relaxationValue)){
    currentRelaxation = relaxationValue[i]
    message(paste("Running categorical at relaxtion:", currentRelaxation))
    
    categoricalResult = categoricalPermulationPipline(T, permulationAmount)
    
    categoricalTimes = list(categoricalSimulationTimes, categoricalPathTimes, categoricalCorrelationTimes, permulationAmount, categoricalPValueTime )
    names(categoricalTimes) = paste0("categorical", currentRelaxation, timeListNameSet)
    categoricalTimesFilename = paste(outputFolderName, filePrefix, "Categorical", currentRelaxation, "-", "TimesFile",  runInstanceValue, ".rds", sep= "")
    saveRDS(categoricalTimes, categoricalTimesFilename)
    
    categoricalForegroundsFilename = paste(outputFolderName, filePrefix, "Categorical", currentRelaxation, "-", "ForegroundsFile",  runInstanceValue, ".rds", sep= "")
    saveRDS(categoricalTempForegroundStorage, categoricalForegroundsFilename)
    
    categoricalResultsFilename = paste(outputFolderName, filePrefix, "Categorical", currentRelaxation, "-", "ResultsFile",  runInstanceValue, ".rds", sep= "")
    saveRDS(categoricalResult, categoricalResultsFilename)
  }
  
}

#-------------------------
#-- Emily Permulations ---
#-------------------------

#create Timing storage objects
emilySimulationTimes = NULL
emilyPathTimes = NULL
emilyCorrelationTimes = NULL
emilyPValueTime = NULL

emilyMidpointForegroundStorage = makeForegroundSpeciesStorage("emilyMidpoint")
emilyNoMidpointForegroundStorage = makeForegroundSpeciesStorage("emilyNoMidpoint")

foregroundSpecies = names(phenotypeVector)[phenotypeVector == 1]

if(runEmily){
  emilySimulationTimes = NULL
  emilyPathTimes = NULL
  emilyCorrelationTimes = NULL
  emilyPValueTime = NULL
  cat("----------------- \n Midpoint Analysis \n ----------- \n")
  
  emilyResults = emilyPermulationPipeline(permulationAmount,T)
  
  emilyTimes = list(emilySimulationTimes, emilyPathTimes, emilyCorrelationTimes, permulationAmount, emilyPValueTime)
  names(emilyTimes) = paste0("emilyMidpoint", timeListNameSet)
  emilyTimesFilename = paste(outputFolderName, filePrefix, "EmilyMidpointTimesFile", runInstanceValue, ".rds", sep= "")
  saveRDS(emilyTimes, emilyTimesFilename)
  
  emilyMidForegroundsFilename = paste(outputFolderName, filePrefix, "EmilyMidpointForegroundsFile", runInstanceValue, ".rds", sep= "")
  saveRDS(emilyMidpointForegroundStorage, emilyMidForegroundsFilename)
  
  emilyResultsFilename = paste(outputFolderName, filePrefix, "EmilyMidpointResultsFile", runInstanceValue, ".rds", sep= "")
  saveRDS(emilyResults, emilyResultsFilename)
  
  
  
  
  # --- no midpoint version --- 
  message("Completed Emily Midpoint Analysis")
  cat("----------------- \n No Midpoint Analysis \n ----------- \n")
  
  #reset the containers 
  emilySimulationTimes = NULL
  emilyPathTimes = NULL
  emilyCorrelationTimes = NULL
  emilyPValueTime = NULL
  
  emilyNoMidResults = emilyPermulationPipeline(permulationAmount,T, F)
  
  emilyNoMidTimes = list(emilySimulationTimes, emilyPathTimes, emilyCorrelationTimes, permulationAmount, emilyPValueTime)
  names(emilyTimes) = paste0("emilyNoMidpoint", timeListNameSet)
  emilyNoMidTimesFilename = paste(outputFolderName, filePrefix, "EmilyNoMidpointTimesFile", runInstanceValue, ".rds", sep= "")
  saveRDS(emilyNoMidTimes, emilyNoMidTimesFilename)
  
  emilyNoMidForegroundsFilename = paste(outputFolderName, filePrefix, "EmilyNoMidpointForegroundsFile", runInstanceValue, ".rds", sep= "")
  saveRDS(emilyNoMidpointForegroundStorage, emilyNoMidForegroundsFilename)
  
  emilyNoMidResultsFilename = paste(outputFolderName, filePrefix, "EmilyNoMidpointResultsFile", runInstanceValue, ".rds", sep= "")
  saveRDS(emilyNoMidResults, emilyNoMidResultsFilename)
}



#-------------------------
#-- BinaryPermulationsFudged ---
#-------------------------


#Convert percentage fudge to absolute fudge 
foregroundSpecies = names(phenotypeVector)[which(phenotypeVector == 1)]
originalTree = foreground2Tree(foregroundSpecies, mainTrees, plotTree = F, clade = "all", useSpecies = speciesFilter)
totalOriginalForeground = sum(originalTree$edge.length)


#Create timing storage objects
fudgedSimulationTimes = NULL
fudgedPathTimes = NULL
fudgedCorrelationTimes = NULL
fudgedPValueTime = NULL

if(runFudged){
  for(i in 1:length(relaxationValue)){
    currentRelaxation = relaxationValue[i]
    message(paste("Running fudge at relaxtion:", currentRelaxation))
    fudgeNumber = as.integer(totalOriginalForeground * currentRelaxation)
    
    fudgedForegroundPrefixname = paste0("fudged", currentRelaxation)
    fudgedTempForegroundStorage = makeForegroundSpeciesStorage(fudgedForegroundPrefixname)
    
    fudgedResult = getPermsBinaryFudgedReport(foregroundSpecies, RERObject, mainTrees, speciesFilter, permulationAmount, rootNode, fudge = fudgeNumber, CorrelationObject, phenotypeVector, foregroundStorage = fudgedTempForegroundStorage)
    
    fudgedForegroundsFilename = paste(outputFolderName, filePrefix, "Fudged", currentRelaxation, "-", "ForegroundsFile",  runInstanceValue, ".rds", sep= "")
    saveRDS(fudgedTempForegroundStorage, fudgedForegroundsFilename)
    
    
    fudgedTimes = list(fudgedSimulationTimes, fudgedPathTimes, fudgedCorrelationTimes,permulationAmount, fudgedPValueTime )
    names(fudgedTimes) = paste0("fudged",currentRelaxation, timeListNameSet)
    fudgedTimesFilename = paste(outputFolderName, filePrefix, "Fudged", currentRelaxation, "-", "TimesFile",  runInstanceValue, ".rds", sep= "")
    saveRDS(fudgedTimes, fudgedTimesFilename)
    
    fudgedResultsFilename = paste(outputFolderName, filePrefix, "Fudged", currentRelaxation, "-", "ResultsFile",  runInstanceValue, ".rds", sep= "")
    saveRDS(fudgedResult, fudgedResultsFilename)
  }
  
}



#-------------------------
#-- Daniel Fast Binary Permulations ---
#-------------------------



#Single-time setup
bitMap = makeLeafMap(rootedMasterTree)
internalNumberValue = countInternal(rootedMasterTree, bitMap,fg=names(phenotypeVector)[which(phenotypeVector==1)])


#create Timing storage objects
danielSimulationTimes = NULL
danielPathTimes = NULL
danielCorrelationTimes = NULL
danielPValueTime = NULL
danielForegroundStorage = makeForegroundSpeciesStorage("daniel")

if(runDaniel){
  danielResults = danielPermulationPipeline(permulationAmount,T)
  
  danielTimes = list(danielSimulationTimes, danielPathTimes, danielCorrelationTimes, permulationAmount, danielPValueTime)
  danielTimesFilename = paste(outputFolderName, filePrefix, "DanielTimesFile", runInstanceValue, ".rds", sep= "")
  saveRDS(danielTimes, danielTimesFilename)
  
  danielForegroundsFilename = paste(outputFolderName, filePrefix, "DanielForegroundsFile", runInstanceValue, ".rds", sep= "")
  saveRDS(danielForegroundStorage, danielForegroundsFilename)
  
  danielResultsFilename = paste(outputFolderName, filePrefix, "DanielResultsFile", runInstanceValue, ".rds", sep= "")
  saveRDS(danielResults, danielResultsFilename)
}

# This is an old line which gives better diagnositics on jus the permulations aspects. 
# if(runDaniel){danielTimes = functionTimer("danielPermulation", permulationAmount, message = T)}









