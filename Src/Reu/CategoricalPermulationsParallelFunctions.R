CategoricalPermulationGetCor =  function (realCors, nullPhens, phenvals, treesObj, RERmat, method = "kw", 
                                          min.sp = 10, min.pos = 2, winsorizeRER = NULL, winsorizetrait = NULL, 
                                          weighted = F, extantOnly = FALSE, report=F) 
{
  tree = treesObj$masterTree
  keep = intersect(names(phenvals), tree$tip.label)
  tree = pruneTree(tree, keep)
  if (is.rooted(tree)) {
    tree = unroot(tree)
  }
  if(report){pathStartTime = Sys.time()}
  message("Generating null paths")
  nullPaths = lapply(nullPhens, function(x) {
    if(report){message("One path complete")}
    tr = tree
    tr$edge.length = c(x$tips, x$nodes)[tr$edge[,2]]
    tree2Paths(tr, treesObj, categorical = TRUE, useSpecies = names(phenvals))
  })
  if(report){pathsEndTime = Sys.time(); pathsDuration = pathsEndTime - pathStartTime; message(paste("Completed paths;","Duration", pathsDuration, attr(pathsDuration, "units")))}
  
  message("Calculating correlation statistics")
  corsMatPvals = matrix(nrow = nrow(RERmat), ncol = length(nullPhens), dimnames = list(rownames(RERmat), NULL))
  corsMatEffSize = matrix(nrow = nrow(RERmat), ncol = length(nullPhens), dimnames = list(rownames(RERmat), NULL))
  if(report){message("Matrixes")}
  Ppvals = lapply(1:length(realCors[[2]]), matrix, data = NA, nrow = nrow(RERmat), ncol = length(nullPhens), dimnames = list(rownames(RERmat), NULL))
  names(Ppvals) = names(realCors[[2]])
  Peffsize = lapply(1:length(realCors[[2]]), matrix, data = NA, nrow = nrow(RERmat), ncol = length(nullPhens), dimnames = list(rownames(RERmat), NULL))
  names(Peffsize) = names(realCors[[2]])
  if(report){message("pVals")}
  for (i in 1:length(nullPaths)) {
    if(report){corStartTime = Sys.time()}
    cors = getAllCor(RERmat, nullPaths[[i]], method = method, 
                     min.sp = min.sp, min.pos = min.pos, winsorizeRER = winsorizeRER, 
                     winsorizetrait = winsorizetrait, weighted = weighted)
    if(report){corEndTime = Sys.time(); corDuration = corEndTime - corStartTime; message(paste("Completed Correlation", i, "Duration", corDuration, attr(corDuration, "units")))}
    corsMatPvals[, i] = cors[[1]]$P
    corsMatEffSize[, i] = cors[[1]]$Rho
    for (j in 1:length(cors[[2]])) {
      Ppvals[[names(cors[[2]])[j]]][, i] = cors[[2]][[j]]$P
      Peffsize[[names(cors[[2]])[j]]][, i] = cors[[2]][[j]]$Rho
    }
    #if(report){message(paste("compelted", i))}
    gc()
  }
  output = list(corsMatEffSize, Peffsize, corsMatPvals, Ppvals)
  names(output) = c("corsMatEffSize", "Peffsize", "corsMatPvals", "Ppvals")
  return(output)
}


CategoricalCalculatePermulationPValues = function(realCors, intermediateList, start=1, end=NULL, report=F){
  {totalStart = Sys.time()}
  corsMatEffSize = intermediateList[[1]]
  Peffsize = intermediateList[[2]]
  corsMatPvals = intermediateList[[3]]
  Ppvals = intermediateList[[4]]
  message("Obtaining permulations p-values")
  N = nrow(realCors[[1]]) #
  if(start == 1){ #Only do this if start = 1, because otherwise it's already made and you'll overwrite the old script's results 
    realCors[[1]]$permP = rep(NA, N) #Make a column for permP values in all of the dataframes 
    for (j in 1:length(realCors[[2]])) {
      realCors[[2]][[j]]$permP = rep(NA, N) #Make a column for permP values in all of the dataframes 
    }
  }
  
  #Start updating the correlations
  if(is.null(end)){ #if no end specified
    stop = N
  }else{
    stop = end
  }
  for (gene in start:stop) {
    if(report){geneStart = Sys.time()}
    if (is.na(realCors[[1]]$Rho[gene])) {
      p = NA
    }
    else {
      signVal = sign(realCors[[1]]$Rho[gene])
      MatEffSizes = corsMatEffSize[gene, ]
      signedMatEffSizes = MatEffSizes[which(sign(MatEffSizes) == signVal)]
      p = sum(abs(signedMatEffSizes) > abs(realCors[[1]]$Rho[gene]), na.rm = TRUE)/(sum(!is.na(signedMatEffSizes))+1)
    }
    realCors[[1]]$permP[gene] = p
    for (j in 1:length(realCors[[2]])) {
      if (is.na(realCors[[2]][[j]]$Rho[gene])) {
        p = NA
      }
      else {
        realValue = realCors[[2]][[j]]$Rho[gene]
        signValue = sign(realValue)
        peffValues = Peffsize[[names(realCors[[2]][j])]][gene, ]
        signedPeffValues = peffValues[which( sign(peffValues) == signValue)]
        p = (sum(abs(signedPeffValues) > abs(realValue), na.rm = TRUE)+1)/ (sum(!is.na(signedPeffValues))+1)
      }
      realCors[[2]][[j]]$permP[gene] = p
    }
    if(report){geneEnd = Sys.time(); geneDuration = geneEnd - geneStart;message(paste("Completed Gene", gene, "Duration", geneDuration, attr(geneDuration, "units")))}
  }
  message("Done")
  {totalEnd = Sys.time(); totalDuration = totalEnd - totalStart;message(paste("Completed p-Values; Duration", totalDuration, attr(totalDuration, "units")))}
  return(list(res = realCors, pvals = list(corsMatPvals, Ppvals), effsize = list(corsMatEffSize, Peffsize)))
  
}

#