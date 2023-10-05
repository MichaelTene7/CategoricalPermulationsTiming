.libPaths("/share/ceph/wym219group/shared/libraries/R4")
library(RERconverge)
phenotypeVectors = readRDS("Data/CategoricalPermuationsTimingPhenotypes.rds")
mainTrees = readRDS("Data/CategoricalPermulationsTimingTrees.rds")
premadeTimes = readRDS("Data/CategoricalPermulationsTimes.rds")
source("Src/Reu/cmdArgImport.R")
source("Src/Reu/RelaxedRejectionPermFuncs.R")
args = c("m=cstm")

# -- Argument IMports ---
args = commandArgs(trailingOnly = TRUE)
relaxation = 0 
customRateMatrix = matrix(c(1,2,3,2,4,5,6,7,8),3)

#Rate Model
if(!is.na(cmdArgImport('m'))){
  rateModel = cmdArgImport('m')
  if(rateModel == "cstm"){rateModel = customRateMatrix}
}else{
  stop("Specify rate model")
}

#Relaxation Level
if(!is.na(cmdArgImport('r'))){
  relaxation = cmdArgImport('r')
  relaxation = as.numeric(relaxation)
}else{
  message("Relaxation level not specified, using 0")
}

fileNameIdentifier = paste(rateModel[1], "Relax", relaxation, sep="")
"Output/CategoricalPermuationsTimingPhenotypesER.rds"
phenotypeOutFilename = paste("Output/CategoricalPermulationsTimingPhenotypes", fileNameIdentifier, ".rds", sep="")
timesOutFilename = paste("Output/CategoricalPermulationsTimes", fileNameIdentifier, ".rds", sep="")

# --- Example code used to generate permualtions times ---
phenotypeVector = phenotypeVectors$phenotypeCH_OG

#Code used to get the times
permulationTimeTest = function(phenotypeVector){
  permsStartTime = Sys.time()                                                     #get the time before start of permulations
  permulationData = categoricalPermulations(mainTrees, phenotypeVector, rm = rateModel, rp = "auto", ntrees = 5, percent_relax = relaxation)
  permsEndTime = Sys.time()                                                       #get time at end of permulations
  timetaken = permsEndTime-permsStartTime
  timetaken
}
timePhenotype = permulationTimeTest(phenotypeVector)

#
#Function to run additional time tests
phenotypeVectorMain = phenotypeVectors$phenotypeVectorAllCategories
phenotypeVectorsOut = phenotypeVectors
premadeTimesOut = premadeTimes

#Specify the desired phenotype names; the abbreviation to use, and any substitutions in the form of list(c("replace1","with1"), c("replace2","With2"))
timeTestSave= function(phenNames, abbreviation, substitutions=NULL, phenotypesMain = phenotypeVectors$phenotypeVectorAllCategories){
  phenotypeList = phenNames
  phenotypeName = paste("phenotype", abbreviation, sep = "")
  phenotypeVec = phenotypesMain[which(phenotypesMain %in% phenNames)]
  if(!is.null(substitutions)){
    for( i in 1:length(substitutions)){
      substitutePhenotypes = substitutions[[i]]
      message(paste("replacing", substitutePhenotypes[1], "with", substitutePhenotypes[2]))
      phenotypeVec = gsub(substitutePhenotypes[1], substitutePhenotypes[2], phenotypeVec)
    }
  }
  assign(phenotypeName, phenotypeVec)
  assign(phenotypeName, phenotypeVec, envir = .GlobalEnv)
  
  position = (length(phenotypeVectorsOut)+1)
  phenotypeVectorsOut[position] <<- list(phenotypeVec)
  names(phenotypeVectorsOut)[position] <<- phenotypeName
  
  timeName = paste("time", abbreviation, sep = "")
  timeVal = permulationTimeTest(phenotypeVec)
  assign(timeName, timeVal)
  assign(timeName, timeVal, envir = .GlobalEnv)
  
  position2 = length(premadeTimesOut)+1
  premadeTimesOut[position2] <<- list(timeVal)
  names(premadeTimesOut)[position2] <<- timeName
}

{ #SYM
  phenotypeVectorsOut = phenotypeVectors
  premadeTimesOut = premadeTimes
  rateModel = "SYM"
  
  #HOGC
  timeTestSave(c("Herbivore", "Omnivore", "Generalist", "Carnivore"), "HOGC")
  
  #HOGP
  timeTestSave(c("Herbivore", "Omnivore", "Piscivore", "Generalist"), "HOGP")
  timeHOGP
  
  #HOGI
  timeTestSave(c("Herbivore", "Omnivore", "Generalist", "Insectivore"), "HOGI")
  
  #HOCP
  timeTestSave(c("Herbivore", "Omnivore", "Carnivore", "Piscivore"), "HOCP")
  
  #HOCI 
  timeTestSave(c("Herbivore", "Omnivore", "Carnivore", "Insectivore"), "HOCI")
  
  #HOPI
  timeTestSave(c("Herbivore", "Omnivore", "Piscivore", "Insectivore"), "HOPI")
  
  #HGCP
  timeTestSave(c("Herbivore", "Generalist", "Carnivore", "Piscivore"), "HGCP")
  
  #HGCI
  timeTestSave(c("Herbivore", "Generalist", "Carnivore", "Insectivore"), "HGCI")
  
  #HGPI
  timeTestSave(c("Herbivore", "Generalist", "Piscivore", "Insectivore"), "HGPI")
  
  #HCPI
  timeTestSave(c("Herbivore", "Carnivore", "Piscivore", "Insectivore"), "HCPI")
  
  #OGCP
  timeTestSave(c("Omnivore", "Generalist", "Carnivore", "Piscivore"), "OGCP")
  
  #OGCI
  timeTestSave(c("Omnivore", "Generalist", "Carnivore", "Insectivore"), "OGCI")
  
  #OGPI
  timeTestSave(c("Omnivore", "Generalist", "Piscivore", "Insectivore"), "OGPI")
  
  #OCPI
  timeTestSave(c("Omnivore", "Carnivore", "Piscivore", "Insectivore"), "OCPI")
  
  #GCPI
  timeTestSave(c("Generalist", "Carnivore", "Piscivore", "Insectivore"), "GCPI")
  
  saveRDS(phenotypeVectorsOut, phenotypeOutFilename)
  saveRDS(premadeTimesOut, timesOutFilename)
}

