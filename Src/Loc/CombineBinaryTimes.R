#-------------------------------------------------------
# --- Combine outputs from cluster ----
#-------------------------------------------------------

fileprefix = "emilyPhen"
outputFolder = paste0("Output/", fileprefix, "/")
relaxableTypes = c("Fudged", "Categorical")
unrelaxableTypes = c("EmilyMidpoint", "EmilyNoMidpoint")
relaxLevels = c("0-", "0.05-", "0.1-", "0.2-")
instances = 1:100

relaxCombos = expand.grid(a = relaxLevels, b = relaxableTypes)
runtypes = append(paste0(relaxCombos$b, relaxCombos$a), unrelaxableTypes)

summaries = list()
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


for(k in 1:length(runtypes)){
  
  currentRuntype = runtypes[k]
  
  message(k)
  message(currentRuntype)
  
  if(grepl("Categorical", currentRuntype)){
    categoricalRun = T
  }else{categoricalRun = F}
  
  fileStarter = paste0(outputFolder, fileprefix, currentRuntype)
  combinedTimes = list()
  combinedForeground = list()
  
  
  for(i in instances){
    message(i)
    
    fileSuffix = paste0(i, ".rds")
    
    # -- Time File -- 
    message("Times")
    currentTimeFile = paste0(fileStarter, "TimesFile", fileSuffix)
    currentTime = readRDS(currentTimeFile)
    if(i == 1 | is.null(combinedTimes)){
      for(j in 1:length(currentTime)){
        if(!is.null(currentTime[[j]])){
          combinedTimes[[j]] = currentTime[[j]]
        }else{
          combinedTimes[[j]] = NA
        }
      }
      names(combinedTimes) = names(currentTime)
    }else{
      for(j in 1:length(currentTime)){
        combinedTimes[[j]] = append(combinedTimes[[j]], currentTime[[j]])
        
      }
      
    }
    
    # -- foreground File -- 
    currentForegroundFile = paste0(fileStarter, "ForegroundsFile", fileSuffix)
    currentForeground = readRDS(currentForegroundFile)
    message("Foreground")
    if(i == 1 | is.null(combinedForeground)){
      for(j in 1:length(currentForeground)){
        if(!is.null(currentForeground[[j]])){
          combinedForeground[[j]] = currentForeground[[j]]
        }else{
          currentForeground[[j]] = NA
        }
        
      }
      names(combinedForeground) = names(currentForeground)
    }else{
      for(j in 1:length(currentForeground)){
        if(is.double(currentForeground[[j]]) | j=2){
          combinedForeground[[j]] = rbind(combinedForeground[[j]], currentForeground[[j]])
        }
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
 
  
  typeSummary = processTimes(combinedTimes, currentRuntype, categoricalRun)
  summaries[[length(summaries)+1]] = typeSummary
  names(summaries)[length(summaries)] = currentRunType 
}

summariesDataframe = do.call(rbind, lapply(summaries, as.data.frame))
summariesFilename = paste0(outputFolder, fileprefix, "CombinedTimingSummaries.rds")
saveRDS(summariesDataframe, summariesFilename)