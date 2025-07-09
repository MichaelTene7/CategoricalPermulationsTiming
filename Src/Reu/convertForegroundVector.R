convertForegroundVector = function(phenotypeVector){
  foregroundValue = NULL
  backgroundValue = NULL
  #Determine how to assess foreground vector 
  if(!is.null(foregroundValue)){
    #use the one supplied; no code here
  }else if(!is.null(backgroundValue)){
    foregroundValue = convertBackgroundValueToForegroundValue(phenotypeVector, backgroundValue)
  }else if(isZeroOneVector(phenotypeVector)){
    foregroundValue = 1
  }else if("foreground" %in% phenotypeVector){
    foregroundValue = "foreground"
  }else if("Foreground" %in% phenotypeVector){
    foregroundValue = "Foreground"
  }else if("background" %in% phenotypeVector){
    backgroundValue = "background"
    foregroundValue = convertBackgroundValueToForegroundValue(phenotypeVector, backgroundValue)
  }else if("Background" %in% phenotypeVector){
    backgroundValue = "Background"
    foregroundValue = convertBackgroundValueToForegroundValue(phenotypeVector, backgroundValue)
  }
  
  convertedVector = phenotypeVector
  toForeground = which(convertedVector == foregroundValue)
  toBackground = which(!convertedVector == foregroundValue)
  
  convertedVector[toForeground] = 1
  convertedVector[toBackground] = 0
  namesSaving = names(convertedVector)
  convertedVector = as.numeric(convertedVector)
  names(convertedVector) = namesSaving
  return(convertedVector)
  
}

isZeroOneVector <- function(x) {
  # Try to coerce to numeric safely
  suppressWarnings({
    x_num = as.numeric(x)
  })
  # Check for NAs (which include failed coercion or actual NA)
  if (any(is.na(x_num))) return(FALSE)
  # Check if only values are 0 and 1
  all(sort(unique(x_num)) == c(0, 1))
}

convertBackgroundValueToForegroundValue = function(foregroundVector, backgroundValue){
  if(length(unique(foregroundVector))==2){
    foregroundValue = foregroundVector[which(!unique(foregroundVector) == backgroundValue)]
  }else{
    stop("Background value only valid for 2-state foreground vectors. Use foregroundValue instead.")
  }
  return(foregroundValue)
}
