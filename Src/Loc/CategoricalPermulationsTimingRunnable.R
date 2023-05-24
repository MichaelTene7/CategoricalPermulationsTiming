library(RERconverge)
phenotypeVectors = readRDS("Data/CategoricalPermuationsTimingPhenotypes.rds")
mainTrees = readRDS("Data/CategoricalPermulationsTimingTrees.rds")
premadeTimes = readRDS("Data/CategoricalPermulationsTimes.rds")
source("Src/Reu/cmdArgImport.R")

# --- Example code used to generate permualtions times ---
phenotypeVector = phenotypeVectors$phenotypeCH_OG
rateModel = "SYM"

#Code used to get the times
permulationTimeTest = function(phenotypeVector){
  permsStartTime = Sys.time()                                                     #get the time before start of permulations
  permulationData = categoricalPermulations(mainTrees, phenotypeVector, rm = rateModel, rp = "auto", ntrees = 5)
  permsEndTime = Sys.time()                                                       #get time at end of permulations
  timetaken = permsEndTime-permsStartTime
  timetaken
}
timePhenotype = permulationTimeTest(phenotypeVector)


# ---- Results display code --- 
{
# - Length effect - 
length(phenotypeVectors$phenotypeGPC)
premadeTimes$timeGPC
length(phenotypeVectors$phenotypeCHO)
premadeTimes$timeCHO
length(phenotypeVectors$phenotypeCH_OG)
premadeTimes$timeCH_OG


# - Clustering Effect -
percentCategories = function(input){
  categories = unique(input)
  output = NULL
  for(i in 1:length(categories)){
    percent = length(which(input ==categories[i]))/length(input)
    result = paste(categories[i], paste(substr(as.character(percent*100), 1, 2), "%", sep=''))
    output = append(output, result)
  }
  output
}

#All of these are fairly similar
premadeTimes$timeGPH
premadeTimes$timeGPO
premadeTimes$timeGPC
#Note the increase in this case
percentCategories(phenotypeVectors$phenotypeGPI) 
premadeTimes$timeGPI
#
char2TreeCategorical(phenotypeVectors$phenotypeGPC, mainTrees, plot=T)
premadeTimes$timeGPC

char2TreeCategorical(phenotypeVectors$phenotypeGPI, mainTrees, plot=T)
premadeTimes$timeGPI
#
#Note the decrease in time when percentage is reduced
#Even though the overall tree is larger
length(phenotypeVectors$phenotypeGPI)
percentCategories(phenotypeVectors$phenotypeGPI) 
premadeTimes$timeGPI

length(phenotypeVectors$phenotypeHOI)
percentCategories(phenotypeVectors$phenotypeHOI) 
char2TreeCategorical(phenotypeVectors$phenotypeHOI, mainTrees, plot=T)
premadeTimes$timeHOI


# - Phenotype number effect - 
premadeTimes$timeGPC
premadeTimes$timeGPO
percentCategories(phenotypeVectors$phenotypeGPCO) 
char2TreeCategorical(phenotypeVectors$phenotypeGPCO, mainTrees, plot=T)
premadeTimes$timeGPCO

premadeTimes$timeCHO
premadeTimes$timeCH_OG
percentCategories(phenotypeVectors$phenotypeCHOG) 
char2TreeCategorical(phenotypeVectors$phenotypeCHOG, mainTrees, plot=T)
premadeTimes$timeCHOG

percentCategories(phenotypeVectors$phenotypeCHOI)
char2TreeCategorical(phenotypeVectors$phenotypeCHOI, mainTrees, plot=T)
premadeTimes$timeCHOI


#However, this one appears to be an exception, it does not experience the slowdown
premadeTimes$timeHPI
premadeTimes$timeHOI
premadeTimes$timeHOP
percentCategories(phentotypeHOPI) 
char2TreeCategorical(phenotypeVectors$phenotypeHOPI, mainTrees, plot=T)
premadeTimes$timeHOPI
}

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
  
  saveRDS(phenotypeVectorsOut, "Output/CategoricalPermuationsTimingPhenotypesSYM.rds")
  saveRDS(premadeTimesOut, "Output/CategoricalPermulationsTimesSYM.rds")
}

{ #ER
  phenotypeVectorsOut = phenotypeVectors
  premadeTimesOut = premadeTimes
  rateModel = "ER"
  
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
  
  saveRDS(phenotypeVectorsOut, "Output/CategoricalPermuationsTimingPhenotypesER.rds")
  saveRDS(premadeTimesOut, "Output/CategoricalPermulationsTimesER.rds")
}

{ #SYM
  phenotypeVectorsOut = phenotypeVectors
  premadeTimesOut = premadeTimes
  rateModel = "ARD"
  
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
  
  saveRDS(phenotypeVectorsOut, "Output/CategoricalPermuationsTimingPhenotypesARD.rds")
  saveRDS(premadeTimesOut, "Output/CategoricalPermulationsTimesARD.rds")
}