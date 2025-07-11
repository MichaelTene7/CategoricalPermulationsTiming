#'Calculates permuted correlation and enrichment statistics for binary phenotype (does not enforce matching structure of sister species, but maintains internal & tips foreground counts within a fudge of the original tree)
#'@param fgdspecs A vector containing the tip foreground species
#'@param RERs An RER matrix calculated using \code{\link{getAllResiduals}}.
#'@param trees treesObj from \code{\link{readTrees}}
#'@param useSpecies A vector of species to include
#'@param ntrees An integer number of permulations
#'@param root The species to root the tree on
#'@param fudge The number up to which the permulated tree can differ in total foreground species
#'@param cors The gene correlation results from \code{\link{correlateWithBinaryPhenotype}}
#'@param phenvec A named vector of the phenotype
#'@return A list object with enrichment statistics, correlation p-val, rho, and correlation effect size
#'@export
getPermsBinaryFudged <- function(fgdspecs, RERs, trees, useSpecies, ntrees, root, fudge = 5, cors,
                                 phenvec) {
  # get counts on original tree
  t = foreground2Tree(fgdspecs, trees, plotTree = TRUE, clade = "all",
                      useSpecies = useSpecies)
  fgnum = sum(t$edge.length)
  tips = length(fgdspecs)
  internal = fgnum - tips
  
  # print summary
  print(paste("fgnum:", fgnum))
  print(paste("tips:", tips))
  print(paste("internal:", internal))
  
  # drop species in the tree that we don't want to use
  # this is the tree passed to the simulation function
  drop = trees$masterTree$tip.label[!(trees$masterTree$tip.label %in% useSpecies)]
  t=drop.tip(trees$masterTree, drop)
  t=root.phylo(t, root, resolve.root = T)
  
  # get the ratematrix
  rm=ratematrix(t, phenvec)
  
  # make the data frames 
  statdf=data.frame(matrix(data=NA, nrow=nrow(RERs),ncol=ntrees),row.names=rownames(RERs))
  pvaldf=data.frame(matrix(data=NA, nrow=nrow(RERs),ncol=ntrees),row.names=rownames(RERs))
  
  # generate the trees
  count=1
  while(count<=ntrees){
    
    #get phenotype:
    blsum=0
    
    while(blsum>(fgnum+fudge) | blsum<(fgnum-fudge)){
      ###########################################
      sims=sim.char(t, rm, nsim = 1)[,,1] #sim.char returns a weird array data structure, [,,1] is the named vector we want
      top=names(sort(sims, decreasing = TRUE))[1:tips]
      ###########################################
      tf=foreground2Tree(top, trees, clade="all", plotTree = F)
      blsum=sum(tf$edge.length)
    }
    
    #get path:
    p=tree2Paths(tf, trees, useSpecies = useSpecies)
    
    #run correlation:
    c=correlateWithBinaryPhenotype(RERs, p)
    
    ###########################################
    # this assumes rownames will always match
    statdf[,count]=c$Rho
    pvaldf[,count]=c$P
    ###########################################
    
    print(paste0("finished perm: ", count))
    count=count+1
  }
  
  # get perm p-val:
  corswithpermp=cors
  rows=nrow(corswithpermp)
  corswithpermp$permP=rep(NA,rows)
  # this assumes rownames will always match
  for(g in 1:rows){
    if(is.na(corswithpermp$Rho[g])){
      p=NA
    }
    else{
      p=sum(abs(statdf[g,])>abs(cors$Rho[g]),na.rm=TRUE)/sum(!is.na(statdf[g,]))
    }
    corswithpermp$permP[g]=p
  }
  corswithpermp$permP.adj=p.adjust(corswithpermp$permP, method="BH")
  
  # return results
  return(list(res=corswithpermp, stat=statdf, pval=pvaldf))
}



#'Calculates permuted correlation and enrichment statistics for binary phenotype (does not enforce matching structure of sister species, but maintains internal & tips foreground counts within a fudge of the original tree)
#'@param fgdspecs A vector containing the tip foreground species
#'@param RERs An RER matrix calculated using \code{\link{getAllResiduals}}.
#'@param trees treesObj from \code{\link{readTrees}}
#'@param useSpecies A vector of species to include
#'@param ntrees An integer number of permulations
#'@param root The species to root the tree on
#'@param fudge The number up to which the permulated tree can differ in total foreground species
#'@param cors The gene correlation results from \code{\link{correlateWithBinaryPhenotype}}
#'@param phenvec A named vector of the phenotype
#'@return A list object with enrichment statistics, correlation p-val, rho, and correlation effect size
#'@export
getPermsBinaryFudgedReport <- function(fgdspecs, RERs, trees, useSpecies, ntrees, root, fudge = 5, cors,
                                 phenvec, runCorrelation = T, foregroundStorage = NULL) {
  
  
  # get counts on original tree
  t = foreground2Tree(fgdspecs, trees, plotTree = TRUE, clade = "all",
                      useSpecies = useSpecies)
  fgnum = sum(t$edge.length)
  tips = length(fgdspecs)
  internal = fgnum - tips
  
  # print summary
  print(paste("fgnum:", fgnum))
  print(paste("tips:", tips))
  print(paste("internal:", internal))
  
  # drop species in the tree that we don't want to use
  # this is the tree passed to the simulation function
  drop = trees$masterTree$tip.label[!(trees$masterTree$tip.label %in% useSpecies)]
  t=drop.tip(trees$masterTree, drop)
  t=root.phylo(t, root, resolve.root = T)
  
  # get the ratematrix
  rm=ratematrix(t, phenvec)
  
  # make the data frames 
  statdf=data.frame(matrix(data=NA, nrow=nrow(RERs),ncol=ntrees),row.names=rownames(RERs))
  pvaldf=data.frame(matrix(data=NA, nrow=nrow(RERs),ncol=ntrees),row.names=rownames(RERs))
  
  # generate the trees
  count=1
  while(count<=ntrees){
    
    fudgedSimulationTime = system.time({
      #get phenotype:
      blsum=0
      
      while(blsum>(fgnum+fudge) | blsum<(fgnum-fudge)){
        ###########################################
        sims=sim.char(t, rm, nsim = 1)[,,1] #sim.char returns a weird array data structure, [,,1] is the named vector we want
        top=names(sort(sims, decreasing = TRUE))[1:tips]
        ###########################################
        tf=foreground2Tree(top, trees, clade="all", plotTree = F)
        blsum=sum(tf$edge.length)
      }
    })
    fudgedSimulationTimes <<- append(fudgedSimulationTimes, fudgedSimulationTime["elapsed"])
    
    if(!all(is.null(foregroundStorage))){
      foregroundStorage = recordForegroundSpecies(tf, foregroundStorage)
    }
    fudgedTempForegroundStorage <<- foregroundStorage
    
    
    #get path:
    fudgedPathTime = system.time({p=tree2Paths(tf, trees, useSpecies = useSpecies)})
    fudgedPathTimes <<- append(fudgedPathTimes, fudgedPathTime["elapsed"])
    
    
    
    #run correlation:
    fugedCorrelationTime = system.time({c=correlateWithBinaryPhenotype(RERs, p)})
    fudgedCorrelationTimes <<- append(fudgedCorrelationTimes, fugedCorrelationTime["elapsed"])
    
    
    ###########################################
    # this assumes rownames will always match
    statdf[,count]=c$Rho
    pvaldf[,count]=c$P
    ###########################################
    
    print(paste0("finished perm: ", count))
    count=count+1
  }
  
  # get perm p-val:
  fudgedPValueTime = system.time({
    corswithpermp=cors
    rows=nrow(corswithpermp)
    corswithpermp$permP=rep(NA,rows)
    # this assumes rownames will always match
    for(g in 1:rows){
      if(is.na(corswithpermp$Rho[g])){
        p=NA
      }
      else{
        p=sum(abs(statdf[g,])>abs(cors$Rho[g]),na.rm=TRUE)/sum(!is.na(statdf[g,]))
      }
      corswithpermp$permP[g]=p
    }
    corswithpermp$permP.adj=p.adjust(corswithpermp$permP, method="BH")
  })
  fudgedPValueTime <<- fudgedPValueTime["elapsed"]
  # return results
  return(list(res=corswithpermp, stat=statdf, pval=pvaldf))
}
