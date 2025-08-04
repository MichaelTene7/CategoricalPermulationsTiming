#-------------------------------------------------------
# --- Check on paths ----
#-------------------------------------------------------
library(RERconverge)
emilyPaths = readRDS("Output/emilyPhen/emilyPhenPathsFile.rds")
setupTestPaths = readRDS("Output/setupTest/setupTestPathsFile.rds")

categoricalTree = readRDS("Output/setupTest/setupTestCategoricalTree.rds")
mainTrees = readRDS("../RunRERBinaryMT/data/zoonomiaAllMammalsTrees.rds")
speciesFilter = readRDS("Output/setupTest/setupTestSpeciesFilter.rds")

setupTestPathsBinary = tree2Paths(categoricalTree, mainTrees, useSpecies = speciesFilter)
saveRDS(setupTestPathsBinary, "Output/setupTest/setupTestPathsFile.rds")
?tree2Paths
#-------------------------------------------------------
# --- analyze combined data ----
#-------------------------------------------------------

fileprefix = "emilyPhen"
outputFolder = paste0("Output/", fileprefix, "/")
relaxableTypes = c("Fudged", "Categorical")
unrelaxableTypes = c("EmilyMidpoint", "EmilyNoMidpoint")
relaxLevels = c("0-", "0.05-", "0.1-", "0.2-")
relaxCombos = expand.grid(a = relaxLevels, b = relaxableTypes)
runtypes = append(paste0(relaxCombos$b, relaxCombos$a), unrelaxableTypes)
k=1
processTimes = function(timeObject, runType = NA, totaled = F){
  
  simulationTime = mean(timeObject[[1]], na.rm = T)
  simulationDeviation = sd(timeObject[[1]], na.rm = T)
  
  pathTime = mean(timeObject[[2]], na.rm = T)
  pathDeviation = sd(timeObject[[2]], na.rm = T)
  
  correlationTime = mean(timeObject[[3]], na.rm = T)
  correlationDeviation = sd(timeObject[[3]], na.rm = T)
  
  agregatePvalTime = mean(timeObject[[5]], na.rm = T)
  agregatePvalDeviation = sd(timeObject[[5]], na.rm = T)
  permsPerAgregate = mean(timeObject[[4]], na.rm = T)
  
  pvalTime = agregatePvalTime/permsPerAgregate
  pvalDeviation = agregatePvalDeviation/permsPerAgregate
  
  fulltimeMatirx = rbind(timeObject[[1]], timeObject[[2]], timeObject[[3]])
  fulltimeValues = colSums(fulltimeMatirx, na.rm = T)
  
  subfullTime = mean(fulltimeValues, na.rm = T)
  subfullDeviation = sd(fulltimeValues, na.rm = T)
    
  if(totaled){
    simulationTime = simulationTime / permsPerAgregate
    simulationDeviation = simulationDeviation / permsPerAgregate
    
    pathTime = pathTime / permsPerAgregate
    pathDeviation = pathDeviation / permsPerAgregate
    
    correlationTime = correlationTime / permsPerAgregate
    correlationDeviation = correlationDeviation / permsPerAgregate
    
    subfullTime = subfullTime / permsPerAgregate
    subfullDeviation = subfullDeviation / permsPerAgregate
    
  }
  
  if(!is.na(pvalTime)){
    fullTime = subfullTime + pvalTime    
  }else{fullTime = subfullTime}
  if(!is.na(pvalDeviation)){
    fullDeviation = subfullDeviation + pvalDeviation
  }else{fullDeviation = subfullDeviation}
  
  
  outMessage = paste("\n-------- \n Run type:", runType, "| Totaled Type:", totaled, 
                     "\n Time per permulation:", fullTime, "±", fullDeviation, 
                     "\n Based on: simulation time", simulationTime, "±", simulationDeviation, "; Path time:", pathTime, "±", pathDeviation, "; \n Correlation time:", correlationTime, "±", correlationDeviation, "; pValue time:", pvalTime, "±", pvalDeviation)
  cat(outMessage)
  
  output = list(outMessage, fullTime, fullDeviation, subfullTime, subfullDeviation, pvalTime, pvalDeviation, correlationTime, correlationDeviation, pathTime, pathDeviation, simulationTime, simulationDeviation)
  names(output) = c("outMessage", "fullTime", "fullDeviation", "noPvalTime", "noPvalDeviation", "pvalTime", "pvalDeviation", "correlationTime", "correlationDeviation", "pathTime", "pathDeviation", "simulationTime", "simulationDeviation")
  output
}

summaries = list()
foregrounds = list(edgeList = list(), tipBinaries = list())

for(k in 1:length(runtypes)){

  currentRunType = runtypes[k]
  
  if(grepl("Categorical", currentRunType)){
    categoricalRun = T
  }else{categoricalRun = F}
  
  currentCombinedFileprefix = paste0(outputFolder, currentRunType, "Combined")
  currentCombinedTimesFile = paste0(currentCombinedFileprefix, "Times.rds")
  currentCombinedTimes = readRDS(currentCombinedTimesFile)
  
  typeSummary = processTimes(currentCombinedTimes, currentRunType, categoricalRun)

  summaries[[length(summaries)+1]] = typeSummary
  names(summaries)[length(summaries)] = currentRunType
  
  currentCombinedForegroundsFile = paste0(currentCombinedFileprefix, "Foreground.rds")
  currentCombinedForegrounds = readRDS(currentCombinedForegroundsFile)
  
  foregrounds[1][k] = currentCombinedForegrounds[1]
  foregrounds[2][k] = currentCombinedForegrounds[2]
  
  
  
}

summariesDataframe = do.call(rbind, lapply(summaries, as.data.frame))

View(currentCombinedForegrounds[[2]])

timeSumarry = readRDS(summariesFilename)
foregroundsSummary = readRDS(foregroundSummariesFilename)

foregroundsSummary[2]
View(timeSumarry)

cat(timeSumarry$outMessage)

summariesFilenameCsv = paste0(outputFolder, fileprefix, "CombinedTimingSummaries.csv")
write.csv(timeSumarry, summariesFilenameCsv)

#-------------------------------------------------------
# --- Combine outputs from cluster ----
#-------------------------------------------------------

fileprefix = "emilyPhen"
outputFolder = paste0("Output/", fileprefix, "/")
relaxableTypes = c("Fudged", "Categorical")
unrelaxableTypes = c("EmilyMidpoint", "EmilyNoMidpoint")
relaxLevels = c("0-", "0.05-", "0.1-", "0.2-")
relaxCombos = expand.grid(a = relaxLevels, b = relaxableTypes)
runtypes = append(paste0(relaxCombos$b, relaxCombos$a), unrelaxableTypes)

currentRuntype = "Categorical0.1-"

instances = 1:100

for(k in 1:length(runtypes)){
  
  currentRuntype = runtypes[k]

  fileStarter = paste0(outputFolder, fileprefix, currentRuntype)
  combinedTimes = list()
  combinedForeground = list()
  

  for(i in instances){
    fileSuffix = paste0(i, ".rds")
    
    # -- Time File -- 
    
    currentTimeFile = paste0(fileStarter, "TimesFile", fileSuffix)
    currentTime = readRDS(currentTimeFile)
    if(i == 1){
      for(j in 1:length(currentTime)){
        if(!is.null(currentTime[[j]])){
          combinedTimes[[j]] = currentTime[[j]]
        }else{
          combinedTimes[[j]] = NA
        }
        names(combinedTimes) = names(currentTime)
      }
    }else{
      for(j in 1:length(currentTime)){
        combinedTimes[[j]] = append(combinedTimes[[j]], currentTime[[j]])
        
      }
      
    }
    
    # -- foreground File -- 
    currentForegroundFile = paste0(fileStarter, "ForegroundsFile", fileSuffix)
    currentForeground = readRDS(currentForegroundFile)
    
    if(i == 1){
      for(j in 1:length(currentForeground)){
        combinedForeground[[j]] = currentForeground[[j]]
        
      }
    }else{
      for(j in 1:length(currentForeground)){
        combinedForeground[[j]] = append(combinedForeground[[j]], currentForeground[[j]])
        
      }
    }
  }

  
  combineTimesName = paste0(currentRuntype, "CombinedTimes")
  combineForegroundName = paste0(currentRuntype, "CombinedForeground")
  
  assign(combineTimesName, combinedTimes)
  saveRDS(combinedTimes, paste0(outputFolder, combineTimesName, ".rds"))
  
  assign(combineForegroundName, combinedForeground)
  saveRDS(combinedForeground, paste0(outputFolder, combineForegroundName, ".rds"))
  
}


currentRuntype = "Fudged0-"

test = currentForeground[[2]]

test2 = test

test3 = cbind(test, test2)
#-------------------------------------------------------
# --- Old code ----
#-------------------------------------------------------

load("Src/Reu/scriptsForTestingPermulationFunctions/RERconverge_output.logRTM_binary.OUmodel.RTMspeciesOnly.rds")

myTrees$trees[[4]]

saveRDS(myTrees, "Data/emilyMultiphylo.rds")

plotTree(mainTrees$masterTree)

length(mainTrees$masterTree$tip.label)

saveRDS(pathvec, PathsFileFilename)

saveRDS(myRER, RERFileName)
saveRDS(res, correlationFileName)

fakeCombinedCorrelation = list(NA, NA)

fakePairwiseCorrelation = list(res, res)
fakeCombinedCorrelation = list(res, fakePairwiseCorrelation)

saveRDS(fakeCombinedCorrelation, combinedCorrelationFileName)

fgspec<-c("Pseudomys_novaehollandiae_ABTC08140", "Pseudomys_delicatulus_U1509", "Notomys_alexis_U1308", "Notomys_fuscus_M22830", "Notomys_mitchellii_M21518", "Zyzomys_pedunculatus_Z34925", "Bandicota_indica_ABTC119185", "Nesokia_indica_ABTC117074", "Hyomys_goliath_ABTC42697", "Pseudomys_shortridgei_Z25113", "Paruromys_dominator_JAE4870", "Eropeplus_canus_NMVZ21733")

unique(names(phenotypeVector[which(phenotypeVector == 1)]))


test = readRDS("Output/setupTest/setupTestPhenotypeVector.rds")

paths = phenotypeVector
test = readRDS("Output/setupTest/setupTestCombinedCategoricalCorrelationFile.rds")
test2 = readRDS("Output/setupTest/setupTestRERFile.rds")

tree = masterTree
trees = mainTrees
fg_vec = foregroundSpecies
pathvec = PathsObject


fgspec<-c("Pseudomys_novaehollandiae_ABTC08140", "Pseudomys_delicatulus_U1509", "Notomys_alexis_U1308", "Notomys_fuscus_M22830", "Notomys_mitchellii_M21518", "Zyzomys_pedunculatus_Z34925", "Bandicota_indica_ABTC119185", "Nesokia_indica_ABTC117074", "Hyomys_goliath_ABTC42697", "Pseudomys_shortridgei_Z25113", "Paruromys_dominator_JAE4870", "Eropeplus_canus_NMVZ21733")

(names(foregroundSpecies) %in% fgspec)
(fgspec  %in% names(foregroundSpecies))

testTree = mainTrees$trees[1]
testTree$`ENSMUSP00000000001-mafft-cds.filter.AAreplace`$tip.label %in% fg_vec

unique(names(PathsObject)) %in% myTrees$masterTree$tip.label
unique(names(PathsObject)) %in% mainTrees$masterTree$tip.label
names(fg_vec) %in% unique(names(PathsObject))
