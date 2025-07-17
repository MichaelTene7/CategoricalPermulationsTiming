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


for(k in 1:length(runtypes)){
  
  currentRuntype = runtypes[k]
  
  message(k)
  message(currentRuntype)
  
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
