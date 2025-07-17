#-------------------------------------------------------
# --- Combine outputs from cluster ----
#-------------------------------------------------------

fileprefix = "emilyPhen"
outputFolder = paste0("D:/Output/", fileprefix, "/")
relaxableTypes = c("Fudged", "Categorical")
unrelaxableTypes = c("EmilyMidpoint", "EmilyNoMidpoint")
relaxLevels = c("0-", "0.05-", "0.1-", "0.2-")
instances = 1:100

relaxCombos = expand.grid(a = relaxLevels, b = relaxableTypes)
runtypes = append(paste0(relaxCombos$b, relaxCombos$a), unrelaxableTypes)


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


currentRuntype = "Categorical0-"



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
