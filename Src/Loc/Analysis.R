library(ggplot2)

# -- read data times ---
SYMTimes = readRDS("Output/Hiller/CategoricalPermulationsHillerTimesSYMRelax0.rds")
SYMPhenotypes = readRDS("Output/Hiller/CategoricalPermulationsTimingHillerPhenotypesSYMRelax0.rds")

ERDTimes = readRDS("Output/Hiller/CategoricalPermulationsHillerTimesERRelax0.rds")
ERPhenotypes = readRDS("Output/Hiller/CategoricalPermulationsTimingHillerPhenotypesERRelax0.rds")


# -- set run parameters -- 
inUsePhenotypes = SYMPhenotypes
inUseTimes = SYMTimes
comboletters = c("M", "G")


# -- organize the data into a dataframe ---
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

# --- end organization of data ---
dataSet = compareData

timeSpeciesCorrelation = function(dataSet){
  model = lm(time ~ speciesNum, data = dataSet)
  summary(model)
}

timeCateogryCorrelation = function(dataSet){
  model = lm(time ~ categoryNumber, data = dataSet)
  summary(model)
}


# -- plot with outliers -- 
ggplot(compareData, aes(x = speciesNum, y = time, color = categoryNumber)) +
  geom_point()

timeSpeciesCorrelation(compareData) # no correlation between speceis number and time when mutliple categories 
timeCateogryCorrelation(compareData)







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




