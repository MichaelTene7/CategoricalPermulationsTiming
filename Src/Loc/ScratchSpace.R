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
