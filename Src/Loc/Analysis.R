library(ggplot2)
library(cowplot)

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


# -- Category number effect --
CompareCategoryMeans= c(mean(twoCategory$time), mean(threeCategory$time), mean(fourCategory$time), mean(fiveCategory$time), mean(sixCategory$time))
CompareCatgeorySds = c(sd(twoCategory$time), sd(threeCategory$time), sd(fourCategory$time), sd(fiveCategory$time), sd(sixCategory$time))
compareCategories = data.frame(CompareCategoryMeans, CompareCatgeorySds)
colnames(compareCategories) = c("Mean", "SD")
rownames(compareCategories) = c("twoCategory", "threeCategory", "fourCategory", "fiveCategory", "sixCategory")

compareCategories$noOutlierMean = c(mean(twoCategoryNoOutliers$time), mean(threeCategoryNoOutliers$time), mean(fourCategoryNoOutliers$time), mean(fiveCategoryNoOutliers$time), mean(sixCategoryNoOutliers$time))
compareCategories$noOutlierSD = c(sd(twoCategoryNoOutliers$time), sd(threeCategoryNoOutliers$time), sd(fourCategoryNoOutliers$time), sd(fiveCategoryNoOutliers$time), sd(sixCategoryNoOutliers$time))

compareCategories
write.csv(compareCategories, paste("Output/Analysis/CompareCategoryMeans", savingPrefix, ".csv", sep=""))


# -- category number effect --




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







plotDataCompare(allCategory, color = "categoryChar")



timeSpeciesCorrelation = function(dataSet){
  model = lm(time ~ speciesNum, data = dataSet)
  summary(model)
}

timeCateogryCorrelation = function(dataSet){
  model = lm(time ~ categoryNumber, data = dataSet)
  summary(model)
}

# -- plot with outliers -- 

timeSpeciesCorrelation(compareData) # no correlation between speceis number and time when mutliple categories 
timeCateogryCorrelation(compareData)

timeSpeciesCorrelation(twoCategory)





compareDataOutliers = which(compareData$time >2300)
compareDataNoOutliers = compareData[-compareDataOutliers,]
ggplot(compareDataNoOutliers, aes(x = speciesNum, y = time, color = categoryNumber)) +
  geom_point()




cat4 = compareData[compareData$categoryNumber ==4, ]
cat4Outliers = which(cat4$time > 2500)
cat4NoOutliers = cat4[-cat4Outliers,]


ggplot(cat4NoOutliers, aes(x = speciesNum, y = time, color = phen)) +
  geom_jitter()



plot(compareData$categoryNumber, compareData$time, color = compareData$combo)




