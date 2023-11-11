library(ggplot2)
library(cowplot)
library(gridExtra)
library(grid)
library(pdftools)
library(stringr)
library(english)
# -- Setup arguments --
comboletters = c("M", "G")
replaceLetters = c("PC", "AO")
outlierCutoff = 2

# -- Read the Data -- 
# - Phenotypes - 
PhenotypesData = readRDS("Output/Hiller/CategoricalPermulationsTimingHillerPhenotypesSYMRelax0.2.rds")
phenotypes = data.frame(names(PhenotypesData))

names(phenotypes) = "phenotypeName"
for(i in 1:nrow(phenotypes)){
  phenotypes$length[i] = length(PhenotypesData[[i]])
}
phenotypes$phenotypeCodes = substring(phenotypes$phenotypeName, 10, nchar(phenotypes$phenotypeName))
phenotypes$phenotypeCodes[phenotypes$length == 1] = NA

uniquePhens = unique(phenotypes$phenotypeCodes)
uniquePhens = uniquePhens[-1]

# - Times - 
SYM20Times = readRDS("Output/Hiller/CategoricalPermulationsHillerTimesSYMRelax0.2.rds")
SYM10Times = readRDS("Output/Hiller/CategoricalPermulationsHillerTimesSYMRelax0.1.rds")


SYM5oneToFourTimes = readRDS("Output/Hiller/CategoricalPermulationsHillerTimesSYMRelax0.05OneToFour.rds")
SYM5Times = SYM5oneToFourTimes
SYM5FiveTimes = readRDS("Output/Hiller/CategoricalPermulationsHillerTimesSYMRelax0.05FiveOnly.rds")
SYM5SixTimes = readRDS("Output/Hiller/CategoricalPermulationsHillerTimesSYMRelax0.05SixOnly.rds")
SYM5Times = append(SYM5oneToFourTimes, SYM5FiveTimes, SYM5SixTimes)

SYM0oneToFourTimes = readRDS("Output/Hiller/CategoricalPermulationsHillerTimesSYMRelax0OneToFour.rds")
SYM0Times = SYM0oneToFourTimes
SYM0FiveTimes = readRDS("Output/Hiller/CategoricalPermulationsHillerTimesSYMRelax0FiveOnly.rds")
SYM0SixTimes = readRDS("Output/Hiller/CategoricalPermulationsHillerTimesSYMRelax0SixOnly.rds")
SYM0Times = append(SYM0oneToFourTimes, SYM0FiveTimes, SYM0SixTimes)



# --- organize the data into a dataframe ---
{
# - outlier removal functions - 
getOutlierPhen = function(data){
  uniquePhens = unique(data$phen)
  zScores = NULL
  for(i in 1:length(uniquePhens)){
    currentPhen = uniquePhens[i]
    currentRows = data[which(data$phen == currentPhen),]
    averageZ = mean(currentRows$zScore)
    names(averageZ) = currentPhen
    zScores = append(zScores, averageZ)
  }
  highestZ = zScores[which.max(zScores)]
  outlierPhen = highestZ
  outlierPhen
}

removeOutliers = function(yesOutliers, outlierFactor =1, outlierCutoff = 3){
  noOutliers = yesOutliers
  outlier = (getOutlierPhen(noOutliers))
  if(outlier>outlierCutoff){
    message(paste( "Phenotype", names(outlier), "is an extreme outlier phenotype, average zScore of", outlier, ". Removing regardless of mean vs SD."))
    noOutliers = subset(noOutliers, !phen %in% names(outlier))
  }else if(!sd(noOutliers$time) > (mean(noOutliers$time)/outlierFactor)){
    message("No Outlier Phenotypes")
  }
  while(sd(noOutliers$time) > (mean(noOutliers$time)/outlierFactor)){
    outlier = NULL
    outlier = (getOutlierPhen(noOutliers))
    message(paste( "Phenotype", names(outlier), "is an outlier phenotype, average zScore of", outlier))
    #if(outlier > outlierCutoff){
    noOutliers = subset(noOutliers, !phen %in% names(outlier))
    #}else{
    #  message("Z-score too low, not removing phenotype")
    #  break()
    #}
  }
  noOutliers
}

# - MAIN FUNCTION - 
dataframeTimes = function(inUseTimes, prefix){
  # -- Make the main dataframe --
  compareData = data.frame(names(inUseTimes))
  names(compareData) = "runName"
  compareData$time = as.numeric(inUseTimes)
  compareData$phen = substring(compareData$runName, 1, nchar(compareData$runName) - 7)
  compareData = compareData[-which(is.na(compareData$time)),]
  

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
  compareData$categoryNumber = nchar(compareData$phen)
  compareData$categoryChar = as.character(compareData$categoryNumber)

  
  compareData$combo = F
  for(i in 1:length(comboletters)){
    compareData$combo[grep(comboletters[i], compareData$phen)] = T
  }
  #return(compareData)
  
  # -- Make the category data frames --
  allCategory = compareData
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
  allCategory = rbind(twoCategory, threeCategory, fourCategory, fiveCategory, sixCategory)
  # -- Remove outliers --
  message("TwoCategory")
  twoCategoryNoOutliers = removeOutliers(twoCategory)
  message("ThreeCategory")
  threeCategoryNoOutliers = removeOutliers(threeCategory)
  message("FourCategory")
  fourCategoryNoOutliers = removeOutliers(fourCategory)
  if(!nrow(fiveCategory)==0){
    message("FiveCategory")
    fiveCategoryNoOutliers = removeOutliers(fiveCategory)
  }else{
    fiveCategoryNoOutliers = fiveCategory 
  }
  if(!nrow(sixCategory)==0){
    message("SixCategory")
    sixCategoryNoOutliers = removeOutliers(sixCategory)
  }else{
    sixCategoryNoOutliers = sixCategory
  }
  
  allCategoryNoOutliers = rbind(twoCategoryNoOutliers, threeCategoryNoOutliers, fourCategoryNoOutliers, fiveCategoryNoOutliers, sixCategoryNoOutliers)
  
  outlieredDataList = list(allCategory, twoCategory, threeCategory, fourCategory, fiveCategory, sixCategory)
  names(outlieredDataList) = c("allCategory", "twoCategory", "threeCategory", "fourCategory", "fiveCategory", "sixCategory")
  noOutlierList = list(allCategoryNoOutliers, twoCategoryNoOutliers, threeCategoryNoOutliers, fourCategoryNoOutliers, fiveCategoryNoOutliers, sixCategoryNoOutliers)
  names(noOutlierList) = c("allCategory", "twoCategory", "threeCategory", "fourCategory", "fiveCategory", "sixCategory")
  allDataList = list(outlieredDataList, noOutlierList)
  names(allDataList) = c("yesOutlier", "noOutlier")
  
  allCatName = paste(prefix, "AllCategory", sep="")
  assign(allCatName, allCategory, envir = .GlobalEnv)
  
  AllCatNoOutlierName = paste(prefix, "NoOutliersAllCategory", sep="")
  assign(AllCatNoOutlierName, allCategoryNoOutliers, envir = .GlobalEnv)
  
  dataListName = paste(prefix, "DataList", sep="")
  assign(dataListName, allDataList, envir = .GlobalEnv)
  
  return()
}
}
#---

# --- plotting functions ---
{
# -- core plotting function -- 

plotDataSimple = function(dataSet, xAxis = "speciesNum", yAxis = "time", colorVar = "phen"){
  dataName = deparse(substitute(dataSet))
  message(dataName)
  plot1 = ggplot(dataSet, aes_string(x = xAxis, y = yAxis)) +
    geom_point(aes_string(color = colorVar))
  print(plot1)
}  
  
  
  
plotLM = function(dataSet, dataName, xAxis = "speciesNum", yAxis = "time", colorVar = "phen", modelType = paste(yAxis, "~", xAxis)){
  linearModel = lm(modelType, data = dataSet)
  plot = ggplot(dataSet, aes_string(x = xAxis, y = yAxis)) +
    geom_point(aes_string(color = colorVar))+
    labs(title = paste(dataName))+
    geom_smooth(method = "lm", show.legend = F)+
    if(!is.na(coef(linearModel)[2])){
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
    }
  return(plot)
}

plotDataLinear = function(dataSet, xAxis = "speciesNum", yAxis = "time", colorVar = "phen", modelType = "lm", print = F){
  dataName = deparse(substitute(dataSet))
  dataName = sub(".*\\$", "", dataName)
  message(dataName)
  plot1 = plotLM(dataSet, dataName)
  if(print){print(plot1)}
  return(plot1)
}


plotDataExponential = function(dataSet, graphTitle, print = F){
  xAxis = "categoryNumber"
  colorVar = "categoryChar"
  yAxis = "time"
  dataName = deparse(substitute(dataSet))
  message(dataName)
  linearModel = lm(time ~ exp(categoryNumber), data = dataSet)
  plot1 = ggplot(dataSet, aes(x = categoryNumber, y = time)) +
    geom_point(aes_string(color = colorVar))+
    labs(title = graphTitle)+
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
  if(print){print(plot1)}
  return(plot1)
}
#--


# - species time plot - 
speciesTimePlot = function(dataList, outlierBoolean, printSubplots = F){
  dataListName = deparse(substitute(dataList))
  dataListName = sub("DataList", "", dataListName)
  yesOutlierData = dataList[[1]]
  names(yesOutlierData) = "List"
  noOutlierData = dataList[[2]]
  names(yesOutlierData) = "List"
  
  if(outlierBoolean){
    dataUsed = dataList[[1]]
    if(identical(dataList[[1]],dataList[[2]])){
      outlierString = "No Outliers In Dataset"
    }else{
      outlierString = "Outliers Included"
    }
  }else{
    if(identical(dataList[[1]],dataList[[2]])){
      message("No Outliers in Dataset")
      return()
    }else{
    dataUsed = dataList[[2]]
    outlierString = "Outliers Excluded"
    }
  }

cat2SpecNumPlot = plotDataLinear(dataUsed$twoCategory, print = printSubplots)
cat3SpecNumPlot = plotDataLinear(dataUsed$threeCategory, print = printSubplots)
cat4SpecNumPlot = plotDataLinear(dataUsed$fourCategory, print = printSubplots)
if(!nrow(dataUsed$fiveCategory)==0){cat5SpecNumPlot = plotDataLinear(dataUsed$fiveCategory, print = printSubplots)}else{cat5SpecNumPlot = NULL}
if(!nrow(dataUsed$sixCategory)==0){cat6SpecNumPlot = plotDataLinear(dataUsed$sixCategory, print = printSubplots)}else{cat6SpecNumPlot = NULL}

gridplot = arrangeGrob(cat2SpecNumPlot, cat3SpecNumPlot, cat4SpecNumPlot, cat5SpecNumPlot, cat6SpecNumPlot, nrow =2, ncol = 3)
titleGrob = textGrob(paste(dataListName, outlierString), gp = gpar(fontface = "bold", fontsize = 16))
finalPlot = grid.arrange(titleGrob, gridplot, ncol = 1, heights = c(1, 10))
print(finalPlot)
}

#- Category effect plot - 

categoryEffectPlot = function(dataList, printSubplots = F){
  
  dataListName = deparse(substitute(dataList))
  dataListName = sub("DataList", "", dataListName)
  
  
  catEffectPlotOutlier = plotDataExponential(dataList[[1]]$allCategory, "Outliers Included", printSubplots)
  catEffectPlotNoOutlier = plotDataExponential(dataList[[2]]$allCategory, "Outliers Excluded", printSubplots)
  gridplot = arrangeGrob(catEffectPlotOutlier, catEffectPlotNoOutlier, nrow =1)
  titleGrob = textGrob(paste(dataListName), gp = gpar(fontface = "bold", fontsize = 16))
  finalPlot = grid.arrange(titleGrob, gridplot, ncol = 1, heights = c(1, 10))
  print(finalPlot)
}


# - relaxation effect plots - 

plotDataRelax = function(dataSet, xAxis = "phen", yAxis = "Mean", colorVar = "relaxLevel", title = "", print = F){
  
  dataName = deparse(substitute(dataSet))
  message(dataName)
  plot1 = ggplot(dataSet, aes_string(x = xAxis, y = yAxis)) +
    geom_point(aes_string(color = colorVar))+
    labs(title = title)+
    scale_color_manual(values = c("0" = "red", "5" = "orange", "10" = "cyan", "20"= "blue"), limits = c("0", "5", "10", "20"))
    #scale_color_gradient2(high = "blue", mid = "yellow", low = "red")
  if(print){print(plot1)}
}  

relaxationEffectPlot = function(meansSet, printSubplots = F){
  for(i in min(meansSet$categoryNumber):max(meansSet$categoryNumber)){
    currentMeans = meansSet[meansSet$categoryNumber == i,]
    categoryText = paste(english(i), "Category", sep ="")
    
    if(any(currentMeans$isOutlier)){
      categoryplot = plotDataRelax(currentMeans[currentMeans$isOutlier == F,], title = paste(categoryText, "Outliers Excluded"), print = printSubplots)
      
      plotVar = paste("plot", i, sep="" )
      assign(plotVar, categoryplot)
    }else{
      categoryplot = plotDataRelax(currentMeans, title = paste(categoryText, "No Outliers In Dataset"), print = printSubplots)
      plotVar = paste("plot", i, sep="" )
      assign(plotVar, categoryplot)
    }
  }
  outlierPlot = plotDataRelax(meansSet[meansSet$isOutlier == T,], title = paste("All category Outliers"), print = printSubplots)
  gridplot = plot_grid(plot2,plot3,plot4,plot5,plot6,outlierPlot, nrow =3)
  print(gridplot)
}




}
# ---

# --- Collect means for comparissions ---
collecPhenMeans = function(dataListSet){
  output = NULL
  for(i in 1:length(dataListSet)){
    currentDataset = dataListSet[[i]][[1]][[1]]
    currentNoOutlierDataset = dataListSet[[i]][[2]][[1]]
    outliers = unique(currentDataset$phen)[which(!unique(currentDataset$phen) %in% unique(currentNoOutlierDataset$phen) )]
    
    
    datasetName = names(dataListSet)[i]
    datasetName = sub("DataList", "", datasetName)
    modelType = sub("\\d+.*", "", datasetName)
    relaxLevel = str_extract_all(datasetName, "\\d+") %>% unlist()
    relaxLevelNumeric = as.numeric(relaxLevel)
    
    for(j in 1:length(uniquePhens)){
      phen = uniquePhens[j]
      categoryNumber = nchar(phen)
      categoryChar = as.character(categoryNumber)
      if(phen %in% outliers){
        isOutlier = T
      }else{
        isOutlier = F
      }
      relevantRows = currentDataset[currentDataset$phen == phen,]
      
      Mean = mean(relevantRows$time)
      Sd = sd(relevantRows$time)
      
      PhenID = paste(phen, relaxLevel, modelType, sep="_")
      
      currentPhenResult = data.frame(PhenID, phen, Mean, Sd, categoryNumber, isOutlier, modelType, relaxLevel, datasetName, categoryChar, relaxLevelNumeric)
      rownames(currentPhenResult) = PhenID
      output = rbind(output, currentPhenResult)
    }
  }
  output
}  
  
getCatgeoryMeans = function(meansSet, outlierInclusion = F){
  output = NULL
  for(i in 2:6){
    if(outlierInclusion){
      currentRows = meansSet[meansSet$categoryNumber == i,]
    }else{
      currentRows = meansSet[meansSet$categoryNumber == i & meansSet$isOutlier ==F,]
    }
    Mean_0 = mean(currentRows[currentRows$relaxLevelNumeric == 0,]$Mean)
    Sd_0 = sd(currentRows[currentRows$relaxLevelNumeric == 0,]$Mean)
    Mean_5 = mean(currentRows[currentRows$relaxLevelNumeric == 5,]$Mean)
    Sd_5 = sd(currentRows[currentRows$relaxLevelNumeric == 5,]$Mean)
    PercentDif_5 = ((Mean_0 - Mean_5) / Mean_5) * 100
    Mean_10 = mean(currentRows[currentRows$relaxLevelNumeric == 10,]$Mean)
    Sd_10 = sd(currentRows[currentRows$relaxLevelNumeric == 10,]$Mean)
    PercentDif_10 = ((Mean_0 - Mean_10) / Mean_10) * 100
    Mean_20 = mean(currentRows[currentRows$relaxLevelNumeric == 20,]$Mean)
    Sd_20 = sd(currentRows[currentRows$relaxLevelNumeric == 20,]$Mean)
    PercentDif_20 = ((Mean_0 - Mean_20) / Mean_20) * 100
    categoryOuput = data.frame(Mean_0, Sd_0, Mean_5, Sd_5, PercentDif_5, Mean_10, Sd_10, PercentDif_10,Mean_20, Sd_20, PercentDif_20)
    if(outlierInclusion){rownames(categoryOuput) = paste(english(i), "CategoryAll", sep="")}else{rownames(categoryOuput) = paste(english(i), "CategoryNoOutliers", sep="")}
    output = rbind(output, categoryOuput)
  }
  output
}
  
 
  
  









dataframeTimes(SYM20Times, "SYM20")
dataframeTimes(SYM10Times, "SYM10")
dataframeTimes(SYM5Times, "SYM5")
dataframeTimes(SYM0Times, "SYM0")





SYM20SpPlots = speciesTimePlot(SYM20DataList, T)
speciesTimePlot(SYM20DataList, F)
SYM10SpPlots = speciesTimePlot(SYM10DataList, T)
speciesTimePlot(SYM10DataList, F)
SYM5SpPlots = speciesTimePlot(SYM5DataList, T)
SYM0SpPlots = speciesTimePlot(SYM0DataList, T)


categoryEffectPlot(SYM20DataList)
categoryEffectPlot(SYM10DataList)
categoryEffectPlot(SYM5DataList)
categoryEffectPlot(SYM0DataList)


dataListSet = list(SYM20DataList, SYM10DataList, SYM5DataList, SYM0DataList)
names(dataListSet) = c("SYM20DataList", "SYM10DataList", "SYM5DataList", "SYM0DataList")

phenotypeTimeMeans = collecPhenMeans(dataListSet)
meansSet = phenotypeTimeMeans
relaxationEffectPlot(meansSet, printSubplots = T)

getCatgeoryMeans(phenotypeTimeMeans, outlierInclusion = F) 
getCatgeoryMeans(phenotypeTimeMeans, outlierInclusion = T) 




# -- file output code -- 
a = b #prevent accidental running 
savingPrefix = "Multirelax"

# - species number plots - 
pdfName = paste("Output/Analysis/SpeciesNumberPlots", savingPrefix, ".pdf", sep="")
pdf(file = pdfName, height = 10, width = 20)
speciesTimePlot(SYM20DataList, T)
speciesTimePlot(SYM20DataList, F)
speciesTimePlot(SYM10DataList, T)
speciesTimePlot(SYM10DataList, F)
speciesTimePlot(SYM5DataList, T)
speciesTimePlot(SYM0DataList, T)
dev.off()

# - category effect plots - 
pdfName = paste("Output/Analysis/CategoryEffectPlots", savingPrefix, ".pdf", sep="")
pdf(file = pdfName, height = 10, width = 20)
categoryEffectPlot(SYM20DataList)
categoryEffectPlot(SYM10DataList)
categoryEffectPlot(SYM5DataList)
categoryEffectPlot(SYM0DataList)
dev.off()


# - Relaxation effect plots - 
pdfName = paste("Output/Analysis/RelaxationEffectPlots", savingPrefix, ".pdf", sep="")
pdf(file = pdfName, height = 10, width = 20)
relaxationEffectPlot(meansSet, printSubplots = F)
dev.off()

# - Relaxation Means comparison table - 

catMeansNoOutliers = getCatgeoryMeans(phenotypeTimeMeans, outlierInclusion = F) 
catMeansYesOutliers = getCatgeoryMeans(phenotypeTimeMeans, outlierInclusion = T) 
allCategorywideMeans = rbind(catMeansNoOutliers, catMeansYesOutliers)

catMeansFilename = paste("Output/Analysis/RelaxationEffectTable.csv")
write.csv(allCategorywideMeans, catMeansFilename)


# -- debug code ---
dataList = SYM10DataList
dataListName = deparse(substitute(SYM10DataList))
dataListName = sub("DataList", "", dataListName)


SYM10DataList$yesOutlier$fiveCategory
SYM5DataList$yesOutlier$fourCategory
SYM0DataList$yesOutlier$fourCategory
