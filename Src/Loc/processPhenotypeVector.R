originalPhenotypeVector = readRDS("Data/HillerPhenotypesVectorAllPhenotypes.rds")
merges =list(
  c("Carnivore", "Piscivore") #The carnivore/Piscivores are primarily piscivores
  ,c("Anthropivore", "Omnivore") #Anthropivore is a subset of Omnivores, when in doubt, use the more general term 
  ,c("Planktivore", "Carnivore") #This planktivore is hybrid Carnivore, and was considered carnivores in the main analysis
  )
substitutions = list(
  c("Planktivore", "Insectivore") #Planktivores consume primarily chitinized species, similar to insectivores. 
  ,c("Hematophagy", "Carnivore") #Hematophagy is a specialized form of carnivory. 
  ,c("Ambiguous", "Insectivore") # The two ambiguous species, while slightly unclear if characterized would be done as an insectivore. 
  ,c("Insectivore/Carnivore", "Insectivore") #These two species are primarily insectivorous. For inclusion in the main dataset which combines the two categories, these were assigned insectivore
  ,c("Insectivore/herbivore", "Herbivore") #These species are primarily Herbivorous, and were considered herbivorous in the main analysis
  ,c("Insectivore/Herbivore", "Herbivore") #These species are primarily Herbivorous, and were considered herbivorous in the main analysis
  
)

mergedPhenotypVector = originalPhenotypeVector

# - Merge any phenotypes which are hybrids of mergable phenotypes - 
if(!is.null(merges)){                                                    #Consider species with multiple combined categories as the merged category
  for( i in 1:length(merges)){                                           #Eg if [X] is replaced with [Y], [X/Y] becomes [Y]
    substitutePhenotypes = merges[[i]]
    message(paste("Merging Hybrids of", substitutePhenotypes[1], "/", substitutePhenotypes[2], "to", substitutePhenotypes[2]))
    entriesWithPhen1 = grep(substitutePhenotypes[1],  mergedPhenotypVector)
    entriesWithPhen2 = grep(substitutePhenotypes[2],  mergedPhenotypVector)
    combineEntries = which(entriesWithPhen1 %in% entriesWithPhen2)
    combineIndexes = entriesWithPhen1[combineEntries]
    mergedPhenotypVector[combineIndexes] = substitutePhenotypes[2]
  }
}

if(!is.null(substitutions)){
  for( i in 1:length(substitutions)){
    substitutePhenotypes = substitutions[[i]]
    message(paste("replacing", substitutePhenotypes[1], "with", substitutePhenotypes[2]))
    mergedPhenotypVector = gsub(substitutePhenotypes[1], substitutePhenotypes[2], mergedPhenotypVector)
  }
}


# - Reduce the vector to only species in the 6 phenotypes being examined - 
phenotypes = c("Herbivore", "Omnivore", "Insectivore", "Carnivore", "Piscivore", "Anthropivore")

cleanedPhenotypVector = mergedPhenotypVector[which(mergedPhenotypVector %in% phenotypes)] #This only remove two species which are not in the main analysis, and have a phenotype of NA

mainAnalysisSpecies = readRDS("data/Old/MainRubyPhenotypes.rds")
speciesNotInMainAnalysis = names(cleanedPhenotypVector)[which(!names(cleanedPhenotypVector) %in% names(mainAnalysisSpecies))]

cleanedPhenotypVector = cleanedPhenotypVector[-which(!names(cleanedPhenotypVector) %in% names(mainAnalysisSpecies))]

saveRDS(cleanedPhenotypVector, "Data/CategoricalPermulationsTimingHillerPhenotypes.rds")


table(cleanedPhenotypVector)



# - Breakdown of species removed - 
removedSpecies = mergedPhenotypVector[-which(mergedPhenotypVector %in% phenotypes)]

mainAnalysisSpecies = readRDS("data/Old/MainRubyPhenotypes.rds")
speciesNotInMainAnalysis = which(!names(removedSpecies) %in% names(mainAnalysisSpecies))
removedSpecies = removedSpecies[-speciesNotInMainAnalysis]


nameConversionTable = readRDS("Data/HillerZoonomPhenotypeTable.rds")
removedSpeciesTable = nameConversionTable[which(nameConversionTable$Hiller %in% names(removedSpecies)),c(3,1,4,5)]
removedSpeciesTable[,3]





