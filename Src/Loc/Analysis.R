library(ggplot2)
library(cowplot)
library(dplyr)

# -- read data times ---
SYMTimes = readRDS("Output/Hiller/CategoricalPermulationsHillerTimesSYMRelax0.rds")
SYMPhenotypes = readRDS("Output/Hiller/CategoricalPermulationsTimingHillerPhenotypesSYMRelax0.rds")

ERDTimes = readRDS("Output/Hiller/CategoricalPermulationsHillerTimesERRelax0.rds")
ERPhenotypes = readRDS("Output/Hiller/CategoricalPermulationsTimingHillerPhenotypesERRelax0.rds")


SYM20Times = readRDS("Output/Hiller/CategoricalPermulationsHillerTimesSYMRelax0.2.rds")
SYM20Phenotypes = readRDS("Output/Hiller/CategoricalPermulationsTimingHillerPhenotypesSYMRelax0.2.rds")

# -- set run parameters -- 
inUsePhenotypes = SYM20Phenotypes
inUseTimes = SYM20Times
comboletters = c("M", "G")
replaceLetters = c("PC", "AO")
savingPrefix = "SYM20"




# -- organize the data into a dataframe ---
{
phenotypes = data.frame(names(inUsePhenotypes))
names(phenotypes) = "phenotypeName"
for(i in 1:nrow(phenotypes)){
  phenotypes$length[i] = length(inUsePhenotypes[[i]])
}

compareData = data.frame(names(inUseTimes))
names(compareData) = "runName"
compareData$time = as.numeric(inUseTimes)
compareData$phen = substring(compareData$runName, 1, nchar(compareData$runName) - 7)

compareData = compareData[-which(is.na(compareData$time)),]



uniquePhens = unique(compareData$phen)
compareData$speciesNum = rep(NA)
for(i in 1:length(uniquePhens)){
  currentPhenotype = uniquePhens[i]
  currentPhenotypeNames2 = paste("phenotype", currentPhenotype, sep='')
  phenotypeNumber = which(phenotypes$phenotypeName == currentPhenotypeNames2)
  if(length(phenotypeNumber) != 0){
    compareData$speciesNum[which(compareData$phen == currentPhenotype)] = phenotypes$length[phenotypeNumber]
  }
}
compareData$speciesNum = as.numeric(compareData$speciesNum)

compareData$categoryNumber = NA
for(i in 1:nrow(compareData)){
  compareData$categoryNumber[i] = nchar(compareData$phen[i])
}
compareData$categoryNumber = nchar(compareData$phen)
#compareData$categoryNumber = as.character(compareData$categoryNumber)

compareData$combo = F
for(i in 1:length(comboletters)){
  compareData$combo[grep(comboletters[i], compareData$phen)] = T
}

compareData$categoryChar = as.character(compareData$categoryNumber)
}
# -- Make subsets of the datafor each category number --- 
{
twoCategory = compareData[compareData$categoryNumber ==2, ]
threeCategory = compareData[compareData$categoryNumber ==3, ]
fourCategory = compareData[compareData$categoryNumber ==4, ]
fiveCategory = compareData[compareData$categoryNumber ==5, ]
sixCategory = compareData[compareData$categoryNumber ==6, ]

zscores = function(dataVector){
  zscores = (dataVector - mean(dataVector)) / sd (dataVector)
  zscores
}

twoCategory$zScore = zscores(twoCategory$time)
threeCategory$zScore = zscores(threeCategory$time)
fourCategory$zScore = zscores(fourCategory$time)
fiveCategory$zScore = zscores(fiveCategory$time)
sixCategory$zScore = zscores(sixCategory$time)
}

# -- remove outliers from each set -- 
allCategory = compareData

getOutlierPhenotypes = function(data, outlierCutoff){
  outlierPhens = character()
  uniquePhens = unique(data$phen)
  for(i in 1:length(uniquePhens)){
    currentPhen = uniquePhens[i]
    currentRows = data[which(data$phen == currentPhen),]
    averageZ = mean(currentRows$zScore)
    if(averageZ > outlierCutoff){
      message(paste( "Phenotype", currentPhen, "is an outlier phenotype, average zScore of", averageZ))
      outlierPhens = append(outlierPhens, currentPhen)
    }
  }
  if(length(outlierPhens) == 0){
    message("No outlier phenotypes")
  }else{
    message("Outliers:", paste(outlierPhens))
  }
  outlierPhens
}

outlierCutoff = 2

twoCategoryOutlier = getOutlierPhenotypes(twoCategory, outlierCutoff)
threeCategoryOutlier = getOutlierPhenotypes(threeCategory, outlierCutoff)
fourCategoryOutlier = getOutlierPhenotypes(fourCategory, outlierCutoff)
fiveCategoryOutlier = getOutlierPhenotypes(fiveCategory, outlierCutoff)
sixCategoryOutlier = getOutlierPhenotypes(sixCategory, outlierCutoff)

#manually add any remaining outliers
#fiveCategoryOutlier = append(fiveCategoryOutlier, "AOHIM")

twoCategoryNoOutliers = subset(twoCategory, !phen %in% twoCategoryOutlier)
threeCategoryNoOutliers = subset(threeCategory, !phen %in% threeCategoryOutlier)
fourCategoryNoOutliers = subset(fourCategory, !phen %in% fourCategoryOutlier)
fiveCategoryNoOutliers = subset(fiveCategory, !phen %in% fiveCategoryOutlier)
sixCategoryNoOutliers = subset(sixCategory, !phen %in% sixCategoryOutlier)

allCategoryNoOutliers = rbind(twoCategoryNoOutliers, threeCategoryNoOutliers, fourCategoryNoOutliers, fiveCategoryNoOutliers, sixCategoryNoOutliers)

outlieredDataList = list(allCategory, twoCategory, threeCategory, fourCategory, fiveCategory, sixCategory)
noOutlierList = list(allCategoryNoOutliers, twoCategoryNoOutliers, threeCategoryNoOutliers, fourCategoryNoOutliers, fiveCategoryNoOutliers, sixCategoryNoOutliers)
allDataList = list(outlieredDataList, noOutlierList)

dataFilename = paste("Output/Analysis/RawData", savingPrefix, ".rds", sep="")
save.rds(allDataList, dataFilename)


# --- end organization of data ---



# -- plot the data -- 
plotDataCompare = function(dataSet, xAxis = "speciesNum", yAxis = "time", colorVar = "phen", modelType = "lm"){
  dataName = deparse(substitute(dataSet))
  message(dataName)
  linearModel = lm(paste(yAxis, "~", xAxis), data = dataSet)
  plot1 = ggplot(dataSet, aes_string(x = xAxis, y = yAxis)) +
    geom_point(aes_string(color = colorVar))+
    labs(title = dataName)+
    geom_smooth(method = "lm", show.legend = F)+
    annotate("text",
      x = min(dataSet[,xAxis]+ ((max(dataSet[,xAxis]) - min(dataSet[,xAxis]))/4)), y = max(dataSet[, yAxis]), 
      label = paste(
        "y =", 
        round(coef(linearModel)[2], 2), 
        "x +", 
        round(coef(linearModel)[1], 2), 
        "   R^2 = ", 
        format(summary(linearModel)$r.squared, digits = 3)
        )
      )
  print(plot1)

  dataNoOutliersName = paste(dataName, "NoOutliers", sep='')
  dataNoOutliers = get(dataNoOutliersName)
  if(!nrow(dataSet) == nrow(dataNoOutliers)){
    plot2 = ggplot(dataNoOutliers, aes_string(x = xAxis, y = yAxis)) +
      geom_point(aes_string(color = colorVar))+
      labs(title = dataNoOutliersName)+
      geom_smooth(method = "lm", show.legend = F)+
      annotate("text",
               x = min(dataSet[,xAxis]+ ((max(dataSet[,xAxis]) - min(dataSet[,xAxis]))/4)), y = max(dataNoOutliers[, yAxis]), 
               label = paste(
                 "y =", 
                 round(coef(linearModel)[2], 2), 
                 "x +", 
                 round(coef(linearModel)[1], 2), 
                 "   R^2 = ", 
                 format(summary(linearModel)$r.squared, digits = 3)
               )
      )
    
    
    print(plot2)
  }else{
    message("Dataset has no outliers.")
    return(plot1)
  }
}
# -- Species number effect --
cat2SpecNumPlot = plotDataCompare(twoCategory)
cat3SpecNumPlot = plotDataCompare(threeCategory)
cat4SpecNumPlot = plotDataCompare(fourCategory)
cat5SpecNumPlot = plotDataCompare(fiveCategory)
cat6SpecNumPlot = plotDataCompare(sixCategory)

plot_grid(cat2SpecNumPlot, cat3SpecNumPlot, cat4SpecNumPlot, cat5SpecNumPlot, cat6SpecNumPlot, nrow =2, ncol = 3)

pdfName = paste("Output/Analysis/SpeciesNumberPlots", savingPrefix, ".pdf", sep="")
pdf(file = pdfName, height = 10, width = 16)
plot_grid(cat2SpecNumPlot, cat3SpecNumPlot, cat4SpecNumPlot, cat5SpecNumPlot, cat6SpecNumPlot, nrow =2, ncol = 3)
dev.off()


# -- Category number effect table --
CompareCategoryMeans= c(mean(twoCategory$time), mean(threeCategory$time), mean(fourCategory$time), mean(fiveCategory$time), mean(sixCategory$time))
CompareCatgeorySds = c(sd(twoCategory$time), sd(threeCategory$time), sd(fourCategory$time), sd(fiveCategory$time), sd(sixCategory$time))
compareCategories = data.frame(CompareCategoryMeans, CompareCatgeorySds)
colnames(compareCategories) = c("Mean", "SD")
rownames(compareCategories) = c("twoCategory", "threeCategory", "fourCategory", "fiveCategory", "sixCategory")

compareCategories$noOutlierMean = c(mean(twoCategoryNoOutliers$time), mean(threeCategoryNoOutliers$time), mean(fourCategoryNoOutliers$time), mean(fiveCategoryNoOutliers$time), mean(sixCategoryNoOutliers$time))
compareCategories$noOutlierSD = c(sd(twoCategoryNoOutliers$time), sd(threeCategoryNoOutliers$time), sd(fourCategoryNoOutliers$time), sd(fiveCategoryNoOutliers$time), sd(sixCategoryNoOutliers$time))

compareCategories
write.csv(compareCategories, paste("Output/Analysis/CompareCategoryMeans", savingPrefix, ".csv", sep=""))


# -- category number effect plot --


plotDataCompare(allCategory, x = "categoryNumber", color = "categoryChar")

dataSet = allCategoryNoOutliers
xAxis = "categoryNumber"
colorVar = "categoryChar"
yAxis = "time"
dataName = deparse(substitute(allCategoryNoOutliers))
message(dataName)
linearModel = lm(time ~ exp(categoryNumber), data = dataSet)
plot1 = ggplot(dataSet, aes(x = categoryNumber, y = time)) +
  geom_point(aes_string(color = colorVar))+
  labs(title = dataName)+
  geom_smooth(method = "gam", formula = (y ~ exp(x)), show.legend = F)+
  annotate("text",
           x = min(dataSet[,xAxis]+ ((max(dataSet[,xAxis]) - min(dataSet[,xAxis]))/4)), y = max(dataSet[, yAxis]), 
           label = paste(
             "y = x^", 
             round(coef(linearModel)[2], 2), 
             "+", 
             round(coef(linearModel)[1], 2), 
             "   R^2 = ", 
             format(summary(linearModel)$r.squared, digits = 3)
           )
  )
print(plot1)



# -- direct category comparisons -- 
orderAnamgrams = function(inString){
  paste(sort(strsplit(inString, NULL)[[1]]), collapse = "")
}

phenoTypes = unique(allCategory$phen)


unmergedPhenotypes = NULL
mergedAndUnmergedPhenotypes = NULL
comboComparisions = data.frame()

for(i in 1:length(comboletters)){

  letterInUse = comboletters[i]
  replacementInUse = replaceLetters[i]
  currentComboPhens = phenoTypes[grep(letterInUse, phenoTypes)]
  decomboedPhens = sub(letterInUse, replacementInUse, currentComboPhens)
  for(j in 1:length(currentComboPhens)){
    comboPhen = currentComboPhens[j]
    comboData = allCategory[allCategory$phen == comboPhen,]
    comboAverage = mean(comboData$time)
    comboSD = sd(comboData$time)
    
    decomboedPhenInput = decomboedPhens[j]
    categoryNumber = nchar(decomboedPhenInput)
    possibleDecomboPhens = phenoTypes[grep(orderAnamgrams(decomboedPhenInput), sapply(phenoTypes, orderAnamgrams ))]
    decomboPhenOutput = possibleDecomboPhens[which(nchar(possibleDecomboPhens) == categoryNumber)]
    
    decomboedData = allCategory[allCategory$phen == decomboPhenOutput,]
    deComboAverage = mean(decomboedData$time)
    deComboSD = sd(decomboedData$time)
    
    timeDifference = deComboAverage - comboAverage 
    diffPercent = ((deComboAverage - comboAverage) / comboAverage) * 100
    
    ComboPhen = letterInUse
    
    output = data.frame(categoryNumber, ComboPhen, comboSD, deComboSD, comboAverage, deComboAverage, timeDifference, diffPercent)
    rownames(output) = paste(comboPhen, "/", decomboPhenOutput, sep="")
    comboComparisions = rbind(comboComparisions, output)
    
    names(decomboPhenOutput) = comboPhen
    names(comboPhen) = comboPhen
    unmergedPhenotypes = append(unmergedPhenotypes, decomboPhenOutput)
    mergedAndUnmergedPhenotypes = append(mergedAndUnmergedPhenotypes, comboPhen)
    mergedAndUnmergedPhenotypes = append(mergedAndUnmergedPhenotypes, decomboPhenOutput)
  }
}
options(scipen = 5)
comboComparisions
unmergedPhenotypes
mergedAndUnmergedPhenotypes

mergeUnmergeData = allCategory[which(allCategory$phen %in% mergedAndUnmergedPhenotypes),]
phenNames = mergeUnmergeData$phen
mergeNames = NULL
for(i in 1:length(phenNames)){
  mergeNames[i] = names(which(mergedAndUnmergedPhenotypes == phenNames[i]))
}
mergeUnmergeData$mergeNames = mergeNames

mergeUnmergeDataNoOutliers = allCategoryNoOutliers[which(allCategoryNoOutliers$phen %in% mergedAndUnmergedPhenotypes),]
phenNames = mergeUnmergeDataNoOutliers$phen
mergeNames = NULL
for(i in 1:length(phenNames)){
  mergeNames[i] = names(which(mergedAndUnmergedPhenotypes == phenNames[i]))
}
mergeUnmergeDataNoOutliers$mergeNames = mergeNames




plotDataSimple = function(dataSet, xAxis = "speciesNum", yAxis = "time", colorVar = "phen", modelType = "lm"){
  dataName = deparse(substitute(dataSet))
  message(dataName)
  #linearModel = lm(paste(yAxis, "~", xAxis), data = dataSet)
  plot1 = ggplot(dataSet, aes_string(x = xAxis, y = yAxis)) +
    geom_point(aes_string(color = colorVar))
  print(plot1)
  
  #dataNoOutliersName = paste(dataName, "NoOutliers", sep='')
  #dataNoOutliers = get(dataNoOutliersName)
  #if(!nrow(dataSet) == nrow(dataNoOutliers)){
  #  plot2 = ggplot(dataNoOutliers, aes_string(x = xAxis, y = yAxis)) +
  #    geom_point(aes_string(color = colorVar))+
  #    labs(title = dataNoOutliersName)+
  #  print(plot2)
  #}
}
mergePlot2 = plotDataSimple(mergeUnmergeData[nchar(mergeUnmergeData$mergeNames)==2,], xAxis = "mergeNames", colorVar = "categoryChar")
mergePlot3 = plotDataSimple(mergeUnmergeData[nchar(mergeUnmergeData$mergeNames)==3,], xAxis = "mergeNames", colorVar = "categoryChar")
mergePlot4 = plotDataSimple(mergeUnmergeData[nchar(mergeUnmergeData$mergeNames)==4,], xAxis = "mergeNames", colorVar = "categoryChar")
mergePlot5 = plotDataSimple(mergeUnmergeData[nchar(mergeUnmergeData$mergeNames)==5,], xAxis = "mergeNames", colorVar = "categoryChar")

plot_grid(mergePlot2, mergePlot3, mergePlot4, mergePlot5, nrow = 2)

mergePlotNoOL2 = plotDataSimple(mergeUnmergeDataNoOutliers[nchar(mergeUnmergeDataNoOutliers$mergeNames)==2,], xAxis = "mergeNames", colorVar = "categoryChar")
mergePlotNoOL3 = plotDataSimple(mergeUnmergeDataNoOutliers[nchar(mergeUnmergeDataNoOutliers$mergeNames)==3,], xAxis = "mergeNames", colorVar = "categoryChar")
mergePlotNoOL4 = plotDataSimple(mergeUnmergeDataNoOutliers[nchar(mergeUnmergeDataNoOutliers$mergeNames)==4,], xAxis = "mergeNames", colorVar = "categoryChar")
mergePlotNoOL5 = plotDataSimple(mergeUnmergeDataNoOutliers[nchar(mergeUnmergeDataNoOutliers$mergeNames)==5,], xAxis = "mergeNames", colorVar = "categoryChar")


pdf2Name = paste("Output/Analysis/mergeUnmergePlots", savingPrefix, ".pdf", sep="")
pdf(file = pdf2Name, height = 10, width = 16)
plot_grid(mergePlot2, mergePlot3, mergePlot4, mergePlot5, nrow = 2)
plot_grid(mergePlotNoOL2, mergePlotNoOL3, mergePlotNoOL4, mergePlotNoOL5, nrow = 2)
dev.off()



combo3 = comboComparisions[which(comboComparisions$categoryNumber ==3),]
mean(combo3$diffPercent)
sd(combo3$diffPercent)
mean(combo3$diffPercent[combo3$ComboPhen == "M"])
mean(combo3$diffPercent[combo3$ComboPhen == "G"])

combo4 = comboComparisions[which(comboComparisions$categoryNumber ==4),]
mean(combo4$diffPercent)
sd(combo4$diffPercent)
mean(combo4$diffPercent[combo4$ComboPhen == "M"])
mean(combo4$diffPercent[combo4$ComboPhen == "G"])

mean(combo4$diffPercent[-12])
sd(combo4$diffPercent[-12])
mean(combo4$diffPercent[combo4$ComboPhen == "M"])
mean(combo4$diffPercent[combo4$ComboPhen == "G"][-4])


combo5 = comboComparisions[which(comboComparisions$categoryNumber ==5),]
mean(combo5$diffPercent)
sd(combo5$diffPercent)
mean(combo5$diffPercent[combo5$ComboPhen == "M"])
mean(combo5$diffPercent[combo5$ComboPhen == "G"])

mean(c(combo5$diffPercent[combo5$ComboPhen == "M"][-3], combo5$diffPercent[combo5$ComboPhen == "G"][-5]))
sd(c(combo5$diffPercent[combo5$ComboPhen == "M"][-3], combo5$diffPercent[combo5$ComboPhen == "G"][-5]))
mean(combo5$diffPercent[combo5$ComboPhen == "M"][-3])
mean(combo5$diffPercent[combo5$ComboPhen == "G"][-5])



combo6 = comboComparisions[which(comboComparisions$categoryNumber ==6),]
mean(combo6$diffPercent)
mean(combo6$diffPercent[combo6$ComboPhen == "M"])
mean(combo6$diffPercent[combo6$ComboPhen == "G"])




# -- double combo comparissions 
letter1phens = phenoTypes[grep(comboletters[1], phenoTypes)]
letter2phens = phenoTypes[grep(comboletters[2], phenoTypes)]
doubleComboPhens =  letter1phens[letter1phens %in% letter2phens]

doubleComboComparisions = data.frame()
doubleUnmergedPhenotypes = NULL
doubleMergedAndUnmergedPhenotypes = NULL


decomboedDoubles = doubleComboPhens

for(j in 1:length(comboletters)){
  letterInUse = comboletters[j]
  replacementInUse = replaceLetters[j]
  decomboedDoubles = sub(letterInUse, replacementInUse, decomboedDoubles)
}

for(j in 1:length(decomboedDoubles)){
  comboPhen = doubleComboPhens[j]
  comboData = allCategory[allCategory$phen == comboPhen,]
  comboAverage = mean(comboData$time)
  comboSD = sd(comboData$time)
  
  decomboedPhenInput = decomboedDoubles[j]
  categoryNumber = nchar(decomboedPhenInput)
  possibleDecomboPhens = phenoTypes[grep(orderAnamgrams(decomboedPhenInput), sapply(phenoTypes, orderAnamgrams ))]
  decomboPhenOutput = possibleDecomboPhens[which(nchar(possibleDecomboPhens) == categoryNumber)]
  
  decomboedData = allCategory[allCategory$phen == decomboPhenOutput,]
  deComboAverage = mean(decomboedData$time)
  deComboSD = sd(decomboedData$time)
  
  timeDifference = deComboAverage - comboAverage 
  diffPercent = ((deComboAverage - comboAverage) / comboAverage) * 100
  
  ComboPhen = paste(comboletters, collapse = "")
  
  output = data.frame(categoryNumber, ComboPhen, comboSD, deComboSD, comboAverage, deComboAverage, timeDifference, diffPercent)
  rownames(output) = paste(comboPhen, "/", decomboPhenOutput, sep="")
  doubleComboComparisions = rbind(doubleComboComparisions, output)
  
  names(decomboPhenOutput) = comboPhen
  names(comboPhen) = comboPhen
  doubleUnmergedPhenotypes = append(doubleUnmergedPhenotypes, decomboPhenOutput)
  doubleMergedAndUnmergedPhenotypes = append(doubleMergedAndUnmergedPhenotypes, comboPhen)
  doubleMergedAndUnmergedPhenotypes = append(doubleMergedAndUnmergedPhenotypes, decomboPhenOutput)
}
doubleComboComparisions  


doubleMergeUnmergeData = allCategory[which(allCategory$phen %in% doubleMergedAndUnmergedPhenotypes),]
phenNames = doubleMergeUnmergeData$phen
mergeNames = NULL
for(i in 1:length(phenNames)){
  mergeNames[i] = names(which(doubleMergedAndUnmergedPhenotypes == phenNames[i]))
}
doubleMergeUnmergeData$mergeNames = mergeNames

doubleMergePlot = plotDataSimple(doubleMergeUnmergeData, xAxis = "mergeNames", colorVar = "categoryChar")




spacer = data.frame(NA, NA, NA, NA, NA, NA, NA, NA)
colnames(spacer) = colnames(comboComparisions)
rownames(spacer) = "-"
allComboComparisions = rbind(comboComparisions, spacer)
allComboComparisions = rbind(allComboComparisions, doubleComboComparisions)

comboComparisionsFilename = paste("Output/Analysis/DirectComboComparisions", savingPrefix, ".csv")
write.csv(allComboComparisions, comboComparisionsFilename)



