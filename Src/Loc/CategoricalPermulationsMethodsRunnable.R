.libPaths("/share/ceph/wym219group/shared/libraries/R4")
library(RERconverge)
phenotypeVectorMain = readRDS("Data/CategoricalPermulationsTimingHillerPhenotypes.rds")
mainTrees = readRDS("Data/CategoricalPermulationsTimingHillerTrees.rds")
#premadeTimes = readRDS("Data/CategoricalPermulationsTimes.rds")
source("Src/Reu/cmdArgImport.R")
source("Src/Reu/RelaxedRejectionPermFuncs.R")
args = c("m=SYM", "r=.20")

# -- Argument Imports ---
args = commandArgs(trailingOnly = TRUE)
relaxation = 0 
customRateMatrix = matrix(c(1,2,3,2,4,5,6,7,8),3)

#Rate Model
if(!is.na(cmdArgImport('m'))){
  rateModel = cmdArgImport('m')
  names(rateModel) = rateModel[1]
  if(rateModel == "cstm"){rateModel = customRateMatrix; names(rateModel) = "custom"}
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

fileNameIdentifier = paste(names(rateModel)[1], "Relax", relaxation, sep="")
"Output/CategoricalPermuationsTimingPhenotypesER.rds"
phenotypeOutFilename = paste("Output/Hiller/CategoricalPermulationsTimingHillerPhenotypes", fileNameIdentifier, ".rds", sep="")
timesOutFilename = paste("Output/Hiller/CategoricalPermulationsHillerTimes", fileNameIdentifier, ".rds", sep="")



# --- Example code used to generate permualtions times ---
#phenotypeVector = phenotypeVectors$phenotypeCH_OG

#Code used to get the times
permulationTimeTest = function(phenotypeVector){
  permsStartTime = Sys.time()                                                     #get the time before start of permulations
  permulationData = categoricalPermulations(mainTrees, phenotypeVector, rm = rateModel, rp = "auto", ntrees = 5, percent_relax = relaxation)
  permsEndTime = Sys.time()                                                       #get time at end of permulations
  timetaken = permsEndTime-permsStartTime
  timetaken
}
#timePhenotype = permulationTimeTest(phenotypeVector)

# -----  Function to run additional time tests  ---------
phenotypeVectorsOut = list()
timesOut = list()

makePhenotypeVector = function(phenNames, abbreviation, substitutions=NULL, phenotypesMain = phenotypeVectorMain){
  #Create a phenotypeVector for the Output  list environment
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
  #assign(phenotypeName, phenotypeVec, envir = .GlobalEnv)
  
  position = (length(phenotypeVectorsOut)+1)
  phenotypeVectorsOut[position] <<- list(phenotypeVec)
  names(phenotypeVectorsOut)[position] <<- phenotypeName
  phenotypeVec
}

testPermTime = function(phenotypeVec, phenCode){
  timeName = paste(phenCode, "-time",  sep = "")
  timeVal = permulationTimeTest(phenotypeVec)
  units(timeVal) <- "secs"
  assign(timeName, timeVal)
  assign(timeName, timeVal, envir = .GlobalEnv)
  
  position2 = length(timesOut)+1
  timesOut[position2] <<- list(timeVal)
  names(timesOut)[position2] <<- timeName
}


timeTrials = function(phenSet, outMat, numTrials = 5, subs = NULL, subOnly = F, phenOnly = F){
  skipRun = F
  for(i in 1:ncol(phenSet)){
    phenCode = colnames(phenSet)[i]
    message(phenCode)
    if(!is.null(subs)){
      for( j in 1:length(subs)){
        substitutePhenotypes = subs[[j]]
        substitutePhenotypesCode = substring(substitutePhenotypes, 1, 1)
        if(subOnly){
          if(!length(grep(substitutePhenotypesCode[1], phenCode) == 0))
            skipRun = T
        }
        #message(paste("replacing", substitutePhenotypesCode[1], "with", substitutePhenotypesCode[2]))
        phenCode = gsub(substitutePhenotypesCode[1], "", phenCode)
        if(!length(grep(substitutePhenotypesCode[2], phenCode) == 0)){
          phenCode = paste(phenCode, substitutePhenotypesCode[2], sep="")
        }
      }
      if(skipRun){
        skipRun = F
        message("Only running on subs, no subs found, skipping")
        next()
      }
    }
    phenotypesUsed = phenSet[,1]
    phenVec = makePhenotypeVector(phenotypesUsed, phenCode, substitutions = subs)
    if(!phenOnly){
      for(k in 1:numTrials){
        testPermTime(phenVec, paste(phenCode, "_", k, sep=""))
      }
    }
    message(paste("completed", phenCode, ";", i, "of", ncol(phenSet)))
  }
}

addBreakToOutputs = function(breakName){
  position = (length(phenotypeVectorsOut)+1)
  phenotypeVectorsOut[position] <<- NA
  names(phenotypeVectorsOut)[position] <<- breakName
  
  position2 = length(timesOut)+1
  timesOut[position2] <<- NA
  names(timesOut)[position2] <<- breakName
}

# ---- Generate the phenotype sets ----
phenotypeOptions = unique(phenotypeVectorMain)
phenotypeLetterOptions = substring(phenotypeOptions, 1, 1)

combinations2Phen = combn(phenotypeOptions, 2)
comb2PhenCodes = substring(combinations2Phen[,], 1, 1)
comb2PhenCodes = paste(comb2PhenCodes[1,], comb2PhenCodes[2,], sep = "")

combinations3Phen = combn(phenotypeOptions, 3)
comb3PhenCodes = substring(combinations3Phen[,], 1, 1)
comb3PhenCodes = paste(comb3PhenCodes[1,], comb3PhenCodes[2,], comb3PhenCodes[3,], sep = "")

combinations4Phen = combn(phenotypeOptions, 4)
comb4PhenCodes = substring(combinations4Phen[,], 1, 1)
comb4PhenCodes = paste(comb4PhenCodes[1,], comb4PhenCodes[2,], comb4PhenCodes[3,],comb4PhenCodes[4,], sep = "")

combinations5Phen = combn(phenotypeOptions, 5)
comb5PhenCodes = substring(combinations5Phen[,], 1, 1)
comb5PhenCodes = paste(comb5PhenCodes[1,], comb5PhenCodes[2,], comb5PhenCodes[3,], comb5PhenCodes[4,], comb5PhenCodes[5,],sep = "")

combinations6Phen = combn(phenotypeOptions, 6)
comb6PhenCodes = substring(combinations6Phen[,], 1, 1)
comb6PhenCodes = paste(comb6PhenCodes,collapse = "")


colnames(combinations2Phen) = comb2PhenCodes
colnames(combinations3Phen) = comb3PhenCodes
colnames(combinations4Phen) = comb4PhenCodes
colnames(combinations5Phen) = comb5PhenCodes
colnames(combinations6Phen) = comb6PhenCodes

# ----- Create the combined phenotypes ------
#M = Meativore, Carnivore+Piscivore
#G = Generalivore, Omnivore+Anthropivore 

# - Meativore - 
meativoreSubs = list(c("Piscivore", "Meativore"), c("Carnivore", "Meativore"))

phen3meatCols = which(grep("C", colnames(combinations3Phen)) %in% grep("P", colnames(combinations3Phen)))
phen4meatCols = which(grep("C", colnames(combinations4Phen)) %in% grep("P", colnames(combinations4Phen)))
phen5meatCols = which(grep("C", colnames(combinations5Phen)) %in% grep("P", colnames(combinations5Phen)))
phen6meatCols = which(grep("C", colnames(combinations6Phen)) %in% grep("P", colnames(combinations6Phen)))


# - Generalivore - 
generalivoreSubs = list(c("Omnivore", "Generalivore"), c("Anthropivore", "Generalivore"))

phen3GeneralCols = which(grep("O", colnames(combinations3Phen)) %in% grep("A", colnames(combinations3Phen)))
phen4GeneralCols = which(grep("O", colnames(combinations4Phen)) %in% grep("A", colnames(combinations4Phen)))
phen5GeneralCols = which(grep("O", colnames(combinations5Phen)) %in% grep("A", colnames(combinations5Phen)))
phen6GeneralCols = which(grep("O", colnames(combinations6Phen)) %in% grep("A", colnames(combinations6Phen)))

# - Both - 
bothSubs = append(meativoreSubs, generalivoreSubs)




# ----- Run the permulations to time ----
#outputs are saved after each number of phenotypes, to allow for data collection even if the higher phenotype numbers cause the script to time out

#clear the outputs
phenotypeVectorsOut = list()
timesOut = list()

# - Run timing Tests - 
addBreakToOutputs("2Phenotypes")
phen2Times = data.frame()
timeTrials(combinations2Phen, phen2Times)
timeTrials(combinations3Phen, phen2Times, subs = meativoreSubs, subOnly = T)
timeTrials(combinations3Phen, phen2Times, subs = generalivoreSubs, subOnly = T)
timeTrials(combinations4Phen, phen2Times, subs = bothSubs, subOnly = T)

saveRDS(phenotypeVectorsOut, phenotypeOutFilename)
saveRDS(timesOut, timesOutFilename)


addBreakToOutputs("3Phenotypes")
phen3Times = data.frame()
timeTrials(combinations3Phen, phen3Times)
timeTrials(combinations4Phen, phen3Times, subs = meativoreSubs, subOnly = T)
timeTrials(combinations4Phen, phen3Times, subs = generalivoreSubs, subOnly = T)
timeTrials(combinations5Phen, phen3Times, subs = bothSubs, subOnly = T)

saveRDS(phenotypeVectorsOut, phenotypeOutFilename)
saveRDS(timesOut, timesOutFilename)

addBreakToOutputs("4Phenotypes")
phen4Times = data.frame()
timeTrials(combinations4Phen, phen4Times)
timeTrials(combinations5Phen, phen4Times, subs = meativoreSubs, subOnly = T)
timeTrials(combinations5Phen, phen4Times, subs = generalivoreSubs, subOnly = T)
timeTrials(combinations6Phen, phen4Times, subs = bothSubs, subOnly = T)

saveRDS(phenotypeVectorsOut, phenotypeOutFilename)
saveRDS(timesOut, timesOutFilename)

addBreakToOutputs("5Phenotypes")
phen5Times = data.frame()
timeTrials(combinations5Phen, phen5Times)
timeTrials(combinations6Phen, phen5Times, subs = meativoreSubs, subOnly = T)
timeTrials(combinations6Phen, phen5Times, subs = generalivoreSubs, subOnly = T)

saveRDS(phenotypeVectorsOut, phenotypeOutFilename)
saveRDS(timesOut, timesOutFilename)

addBreakToOutputs("6Phenotypes")
phen6Times = data.frame()
timeTrials(combinations6Phen, phen6Times)

saveRDS(phenotypeVectorsOut, phenotypeOutFilename)
saveRDS(timesOut, timesOutFilename)


# -- Save the outputs -- 
saveRDS(phenotypeVectorsOut, phenotypeOutFilename)
saveRDS(timesOut, timesOutFilename)






































