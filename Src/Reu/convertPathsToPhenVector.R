convertPathsToPhenVector = function(paths){
  tips = unique(names(paths))
  tips = tips[-which(is.na(tips))]
  phenotypeVector = rep(NA, length(tips))
  names(phenotypeVector) = tips
  phenotypeVector= paths[match(tips, names(paths))]
}