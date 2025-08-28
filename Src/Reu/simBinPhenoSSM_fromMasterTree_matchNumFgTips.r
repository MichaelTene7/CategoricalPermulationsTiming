#A modification of the original simBinPhenoSSM function that runs faster and reduces bias in which species end up in the simulated foregrounds
#Changes from original simBinPhenoSSM:
  #fixed variable naming bug; output of foreground2Tree is now "t_iter" to be consistent with simBinPhenoCC
  #hard-coded variance as 0.02 for simulations, instead of trying to calculate variance for a simulated continuous trait based on binary data
  #subset taxa from the master tree based on the gene tree, instead of giving it the real gene tree, to account for the fact that many branches have length zero in the gene trees
  #midpoint root the master tree for running simulations over the tree; leads to a more even distribution of branches that end up in the simulated foregrounds
  #relaxed foreground structure requirements for simulated trees (match number of fg tips ONLY, ignoring number and structure of internal fg nodes)
  #added a counter to the while loops such that it only tries 50 times to find a permulated tree the matches the conditions (same number of foreground branches); if it cannot find one after 50 tries it returns a NULL tree
#'Produces one SSM binary permulation for a gene
#' @param tree Tree of the gene of interest
#' @param trees treesObj from \code{\link{readTrees}}
#' @param fg_vec A vector containing the foreground species
#' @param pathvec A path vector generated from the real set of foreground animals
#' @param plotTreeBool Boolean indicator for plotting the output tree (default=FALSE)
#' @return A SSM binary permulated tree
#' @export
simBinPhenoSSM_fromMasterTree_matchNumFgTips=function(tree, trees, fg_vec, pathvec, plotTreeBool=F){
  require(phytools)
  tip.labels = tree$tip.label # the set of species that exist in the gene tree
  ind_fg = which(tip.labels %in% fg_vec) # indices of the observed foreground animals that exist in the gene tree

  if (length(ind_fg) == 0){ #If no foregrounds, return NULL tree
    t_iter = tree
    t_iter$edge = NULL
    t_iter$edge.length = NULL
    t_iter$Nnode = NULL
    t_iter$tip.label = NULL
  } else {
    #Get the number of foreground tips present in this gene tree
    fg_k = tip.labels[ind_fg] # the list of the observed foreground animals that exist in the gene tree
    #print(paste("foregrounds:", fg_k))
    tips=length(fg_k) #Number of foreground tips in real data
    print(paste("Number of foreground tips in real data:", tips))
    
    #Generate the tree on which simulations will be run
    #Emily Kopania's modifications from original simBinPhenoSSM function, which reduce bias in which species end up in simulated foregrounds:
      #keeps the topology and branch lengths of the master tree but ONLY keeps the tips present in this gene tree
      #midpoint roots the tree
    t = midpoint.root(keep.tip(trees$masterTree, tip.labels))
    #print(t)
    #write.tree(t, "midpoint_root.speciesSubset_masterTree.tre", append=TRUE)
    rm = ratematrix(t, pathvec)
    #print(rm)

    #Simulates a tree with the same number of foreground tips as the real data; continues if it can't simulate a tree matching that condition in 50 tries
    #Note: It should get it on the first try, since it's just taking the top n tips based on simulated phenotype values, where n is the number of fg tips in the real data
    num_fg_tips = 0
    try_count = 0
    while( (num_fg_tips != tips) && (try_count < 50) ){
        #Simulate continuous values using a Brownian motion model
        sims = sim.char(t, rm, nsim = 1, model="BM")
        #print(sims)
        #Get top n species based on simulated data, where n is the number of foreground tips in te real data
        nam = rownames(sims)
        s = as.data.frame(sims)
        simulatedvec = s[,1]
        names(simulatedvec) = nam
        top.all = names(sort(simulatedvec, decreasing = TRUE))
        top.tree_k = top.all[top.all %in% tip.labels]
        top = top.tree_k[1:tips]
        #print("Top:")
        #print(top)
        #Generate a simulated foreground tree
        t_iter = foreground2Tree(top, trees, clade = "all", plotTree = F, useSpecies=tip.labels)
        print("Simulated foreground tree:")
        write.tree(t_iter, file = "simFGtrees.tre", append=TRUE)
        #Get all foreground edges, regardless of whether they are internal or tips
        fgEdges = t_iter$edge[which(t_iter$edge.length==1),2]
        #Get all foreground tips
        permFgs = t_iter$tip.label[fgEdges]
        #Number of foreground tips for checking match; replaces "blsum" in original simBinPhenoSSM function
        num_fg_tips = length(permFgs)
        print(paste("number of simulated foreground tips:", num_fg_tips))
        try_count = try_count+1
    }
    if(try_count==50){ #This shouldn't be necessary, but leaving it in so the function has something to return just in case
        print("Assigning null tree")
        t_iter = tree
        t_iter$edge = NULL
        t_iter$edge.length = NULL
        t_iter$Nnode = NULL
        t_iter$tip.label = NULL
    }
  }
  if (plotTreeBool){
    if(!(is.null(t_iter$tip.label))){
        print(t_iter)
        plot(t_iter)
        write.tree(t_iter, "temp.tre", append=T)
    } else{
        write("NULL", "temp.tre", append=T)
    }
  }
  return(t_iter)
}
