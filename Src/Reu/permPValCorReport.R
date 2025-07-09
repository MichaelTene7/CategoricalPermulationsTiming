#Usage: This function calculates p values from inputed permulations 
# Arguments: 
  # Start number: the gene number to start at. Used for parrallelization
  # Gene Number: The number of genes to iterate across. If <NA>, will calculate for all genes
  # Report: Determines if the script should report times for all genes. If FALSE, only reports for the first and last 10 genes. 
  # plusOne: Determines if the denominator should be equal to the number of valid permulations, or number of valid permulations +1
  # signedDenominator: Determines if the denominator should include all permulations, or only permulations with values in the same direction (sign) as the real value. 
    # This prevents, in effect, half of the permulations always being counted 



permPValCorReport = function (realcor, permvals, startNumber=1, geneNumber = NA, report = FALSE, plusOne = FALSE, signedDenominator = TRUE) {
  timeStart = Sys.time()
  permcor = permvals$corRho                                                     #Rho values from the permulations
  timePermExtractEnd = Sys.time()
  message("Time to extract permulation values: ", timePermExtractEnd - timeStart, attr(timePermExtractEnd - timeStart, "units"))
  
  realstat = realcor$Rho                                                        #Rho values from the real data 
  timeStatExtractEnd = Sys.time()
  message("Time to extract real values: ", timeStatExtractEnd - timePermExtractEnd, attr(timeStatExtractEnd - timePermExtractEnd, "units"))
  names(realstat) = rownames(realcor)                                           #Name the values after the row names 
  
  permcor = permcor[match(names(realstat), rownames(permcor)),]                 #trim the permulations' Rho values to only ones in the real data 
  timeMatchEnd = Sys.time()
  message("Time to match names: ", timeMatchEnd - timeStatExtractEnd, attr(timeMatchEnd - timeStatExtractEnd, "units"))
  
  permpvals = vector(length = length(realstat))                                 #Make a vector to store the pValues of a length equal to the number of realstat genes
  names(permpvals) = names(realstat)                                            #Copy the names from the real Rho values onto that vector
  count = startNumber                                                           #Set a count equal to the gene to start at  
  
  timeBeforeLoopStart = Sys.time()
  message("Total time before loop start: ", timeBeforeLoopStart - timeStart, attr(timeBeforeLoopStart - timeStart, "units"))
  
  if(is.na(geneNumber)){endpoint = length(realstat)}else{endpoint = startNumber + geneNumber -1}  # set the endpoint to be all if gene number is not specified; or start value + number of genes to do if both specified. 
  while (count <= endpoint) {                                                   #While the count hasn't done all of the entries                                         
    timeLoopBegin = Sys.time()
    if(report){message("Gene number: ", count)}
    if (is.na(realstat[count])) {                                               #if the real correlation is NA
      permpvals[count] = NA                                                     #Make the permulation value NA
      if(report){message("is NA")}
    }
    else if((count <= 10 | count >= (endpoint-10)) & report) {                                    #If not, send detailed time messages for the last 10  
      permRow = permcor[count, ]
      permCol = t(permRow)
      permColMakeTime = Sys.time()
      message("Time to make permulation Col: ", permColMakeTime - timeLoopBegin, attr(permColMakeTime - timeLoopBegin, "units"))
      
      realRow = realstat[count]
      realRowMakeTime = Sys.time()
      message("Time to make real row: ", realRowMakeTime - permColMakeTime, attr(realRowMakeTime - permColMakeTime, "units"))
      
      timeNumeratorCalc = Sys.time()
      message("Numerator Calculate time: ", timeNumeratorCalc - realRowMakeTime, attr(timeNumeratorCalc - realRowMakeTime, "units"))
      if(signedDenominator){
        correctSign = sign(realstat[count])
        signFilteredPermCol = permCol[sign(permCol) == correctSign]
        num = sum(abs(signFilteredPermCol) >= abs(realRow),                   #Make a numerator which is the sum of the permulated values greater than the actual values; after removing NA entires. 
                  na.rm = T)
        denom = sum(!is.na(signFilteredPermCol)) #Make a denominator which is the sum of the non-NA permulation values
        if(plusOne){denom = denom+1; num = num+1} #If using complex denominator, set denominator = number of permulations +1
      }else{
        num = sum(abs(permCol) >= abs(realRow),                   #Make a numerator which is the sum of the permulated values greater than the actual values; after removing NA entires. 
                  na.rm = T)
        denom = sum(!is.na(permCol)) #Make a denominator which is the sum of the non-NA permulation values
        if(plusOne){denom = denom+1; num = num+1} #If using simple denominator, set denominator = number of permulations +1
      }
      timeDenomSum = Sys.time()
      message("Denominator sum time: ", timeDenomSum - timeNumeratorCalc, attr(timeDenomSum - timeNumeratorCalc, "units"))
      permpvals[count] = num/denom                                              #p-value = numerator/denominator
    }
    else{                                                                       #Stop sending detailed time messages
      permRow = permcor[count, ]
      permCol = t(permRow)
     
      realRow = realstat[count]

      if(signedDenominator){
        correctSign = sign(realstat[count])
        signFilteredPermCol = permCol[sign(permCol) == correctSign]
        num = sum(abs(signFilteredPermCol) >= abs(realRow),                   #Make a numerator which is the sum of the permulated values greater than the actual values; after removing NA entires. 
                  na.rm = T)
        denom = sum(!is.na(signFilteredPermCol)) #Make a denominator which is the sum of the non-NA permulation values
        if(plusOne){denom = denom+1; num = num+1} #If using complex denominator, set denominator = number of permulations +1
      }else{
        num = sum(abs(permCol) >= abs(realRow),                   #Make a numerator which is the sum of the permulated values greater than the actual values; after removing NA entires. 
                  na.rm = T)
        denom = sum(!is.na(permCol)) #Make a denominator which is the sum of the non-NA permulation values
        if(plusOne){denom = denom+1; num = num+1} #If using simple denominator, set denominator = number of permulations +1
      }
      permpvals[count] = num/denom 
    }
    timeLoopEnd = Sys.time()
    if(report){message("Total loop time: ", timeLoopEnd - timeLoopBegin, attr(timeLoopEnd - timeLoopBegin, "units"))}
    count = count + 1                                                           #increase count by one
    }
    message("calculation compelte.")
    timeEnd = Sys.time()
    message("Total calculation time: ", timeEnd - timeStart, attr(timeEnd - timeStart, "units"))
  permpvals = permpvals[startNumber:endpoint]                          #remove any entires from the vector that are below the start value or above the end value 
  permpvals                                                                     #Return the vector of pValues 
}
