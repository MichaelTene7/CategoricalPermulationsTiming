#EK - fixed variable naming bug; "t" is still the output of foreground2Tree and "sim.t" is now the midpoint rooted tree for simulating trait values across
#EK - hard-coded variance as 0.02 for simulations, instead of trying to calculate variance for a simulated continuous trait based on binary data
#EK - subset taxa from the master tree based on the gene tree, instead of giving it the real gene tree, to account for the fact that many branches have length zero in the gene trees
#Ek - midpoint root the master tree for running simulations over the tree; leads to a more even distribution of branches that end up in the simulated foregrounds
#EK - removed the requirement that permulated tree structure matches true structure for the gene tree (commented out testcondition in bottom while loop)
#EK - added a counter to the while loops such that it only tries 50 times to find a permulated tree the matches the conditions (same number of foreground branches); if it cannot find one after 50 tries it returns a NULL tree
#'Produces one SSM binary permulation for a gene
#' @param tree Tree of the gene of interest
#' @param trees treesObj from \code{\link{readTrees}}
#' @param fg_vec A vector containing the foreground species
#' @param sisters_list A list containing pairs of "sister species" in the foreground set (put NULL if empty)
#' @param pathvec A path vector generated from the real set of foreground animals
#' @param plotTreeBool Boolean indicator for plotting the output tree (default=FALSE)
#' @return A SSM binary permulated tree
#' @export
simBinPhenoSSM_fromMasterTree=function(tree, trees, fg_vec, sisters_list=NULL, pathvec, plotTreeBool=F){
  require(phytools)
  tip.labels = tree$tip.label # the set of species that exist in the gene tree
  ind_fg = which(tip.labels %in% fg_vec) # indices of the observed foreground animals that exist in the gene tree

  if (length(ind_fg) == 0){
    t = tree
    t$edge = NULL
    t$edge.length = NULL
    t$Nnode = NULL
    t$tip.label = NULL
  } else {
    fg_k = tip.labels[ind_fg] # the list of the observed foreground animals that exist in the gene tree
    res = getForegroundInfoClades(fg_k,sisters_list,trees,plotTree=F,useSpecies=tip.labels)
    fg_tree = res$tree
    fg.table = res$fg.sisters.table

    #EK - keeps the topology and branch lengths of the master tree but ONLY keeps the tips present in this gene tree
    #EK - also midpoint roots the tree
    #EK - variable is named "sim.t" to distinguish it from "t" (output of foreground2Trees below)
    sim.t=midpoint.root(keep.tip(trees$masterTree, tip.labels))
    #print(sim.t)
    #write.tree(sim.t, "midpoint_root.subset_masterTree.10loci.tre", append=TRUE)
    #EK - hard code a constant rate matrix variance value
    rm=0.02
    #print(rm)

    if (!is.null(sisters_list)){
      fg_tree_info = getBinaryPermulationInputsFromTree(fg_tree)
      num_tip_sisters_true = unlist(fg_tree_info$sisters_list)
      num_tip_sisters_true = num_tip_sisters_true[which(num_tip_sisters_true %in% tip.labels)]
      num_tip_sisters_true = length(num_tip_sisters_true)
      #fg_tree_depth_order = getDepthOrder(fg_tree)
    }

    #EK - I think this needs to be out of the if(!is.null) because later fg_tree_depth_order needs to exist regardless of if sisters_list is null of not
    fg_tree_depth_order = getDepthOrder(fg_tree)
    fgnum = length(which(fg_tree$edge.length == 1))
    if (!is.null(sisters_list)){
      internal = nrow(fg.table)
    } else {
      internal = 0
    }
    tips=fgnum-internal # the number of tips
    #print(paste("Number of foregrounds:", fgnum))
    testcondition=FALSE
    while(!testcondition){
      blsum=0
      #EK - Modified so this tries 50 times to generate a tree with the same number of foregrounds; if it doesn't in that time, move on
      #EK - Based on a previous run, <1% of gene trees took 50 times or more, so this shouldn't get rid of that many
      try_count=0
      #while(blsum!=fgnum){
      while( (blsum!=fgnum) && (try_count < 50) ){
        sims=sim.char(sim.t, rm, nsim = 1, model="BM")
        #print(sims)
        nam=rownames(sims)
        s=as.data.frame(sims)
        simulatedvec=s[,1]
        names(simulatedvec)=nam
        top.all=names(sort(simulatedvec, decreasing = TRUE))
        top.tree_k = top.all[top.all %in% tip.labels]
        top = top.tree_k[1:tips]
        #print("Top:")
        #print(top)
        t=foreground2Tree(top, trees, clade="all", plotTree = F, useSpecies=tip.labels)
        blsum=sum(t$edge.length)
        try_count=try_count+1
        #print(paste("blsum and try_count:", blsum, try_count))
      }
      if(try_count==50){
        #print("Assigning null tree")
        t = tree
        t$edge = NULL
        t$edge.length = NULL
        t$Nnode = NULL
        t$tip.label = NULL
        testcondition=TRUE
        #print("NULL TREE")
      } else{
        t_info = getBinaryPermulationInputsFromTree(t)
        if (!is.null(sisters_list)){
          num_tip_sisters_fake = unlist(t_info$sisters_list)
num_tip_sisters_fake = num_tip_sisters_fake[which(num_tip_sisters_fake %in% tip.labels)]
          num_tip_sisters_fake = length(num_tip_sisters_fake)
          t_depth_order = getDepthOrder(t)
          #testcondition = setequal(sort(t_depth_order), sort(fg_tree_depth_order)) &&
          #  (num_tip_sisters_fake == num_tip_sisters_true)
          testcondition = TRUE
        } else {
          t_depth_order = getDepthOrder(t)
          #testcondition = setequal(sort(t_depth_order), sort(fg_tree_depth_order))
          testcondition = TRUE
        }
      }
    }
  }
  if (plotTreeBool){
    if(!(is.null(t$tip.label))){
        print(t)
        plot(t)
        write.tree(t, "temp.tre", append=T)
    } else{
        write("NULL", "temp.tre", append=T)
    }
  }
  return(t)
}


###SAME AS ABOVE BUT WITHOUT MIDPOINT ROOTING###
#Use to test effect of midpoint root
simBinPhenoSSM_fromMasterTree_noMidpointRoot=function(tree, trees, fg_vec, sisters_list=NULL, pathvec, plotTreeBool=F){
  require(phytools)
  tip.labels = tree$tip.label # the set of species that exist in the gene tree
  ind_fg = which(tip.labels %in% fg_vec) # indices of the observed foreground animals that exist in the gene tree

  if (length(ind_fg) == 0){
    t = tree
    t$edge = NULL
    t$edge.length = NULL
    t$Nnode = NULL
    t$tip.label = NULL
  } else {
    fg_k = tip.labels[ind_fg] # the list of the observed foreground animals that exist in the gene tree
    res = getForegroundInfoClades(fg_k,sisters_list,trees,plotTree=F,useSpecies=tip.labels)
    fg_tree = res$tree
    fg.table = res$fg.sisters.table

    #EK - keeps the topology and branch lengths of the master tree but ONLY keeps the tips present in this gene tree
    #EK - arbitrarily roots the tree on the first species listed in tip.label vector
    #EK - variable is named "sim.t" to distinguish it from "t" (output of foreground2Trees below)
    sim.t=root.phylo(keep.tip(trees$masterTree, tip.labels), tree$tip.label[1], resolve.root=T)
    #print(sim.t)
    #EK - hard code a constant rate matrix variance value
    rm=0.02
    #print(rm)

    if (!is.null(sisters_list)){
      fg_tree_info = getBinaryPermulationInputsFromTree(fg_tree)
      num_tip_sisters_true = unlist(fg_tree_info$sisters_list)
      num_tip_sisters_true = num_tip_sisters_true[which(num_tip_sisters_true %in% tip.labels)]
      num_tip_sisters_true = length(num_tip_sisters_true)
      #fg_tree_depth_order = getDepthOrder(fg_tree)
    }

    #EK - I think this needs to be out of the if(!is.null) because later fg_tree_depth_order needs to exist regardless of if sisters_list is null of not
    fg_tree_depth_order = getDepthOrder(fg_tree)
    fgnum = length(which(fg_tree$edge.length == 1))
    if (!is.null(sisters_list)){
      internal = nrow(fg.table)
    } else {
      internal = 0
    }
    tips=fgnum-internal # the number of tips
    #print(paste("Number of foregrounds:", fgnum))
    testcondition=FALSE
    while(!testcondition){
      blsum=0
      #EK - Modified so this tries 50 times to generate a tree with the same number of foregrounds; if it doesn't in that time, move on
      #EK - Based on a previous run, <1% of gene trees took 50 times or more, so this shouldn't get rid of that many
      try_count=0
      #while(blsum!=fgnum){
      while( (blsum!=fgnum) && (try_count < 50) ){
        sims=sim.char(sim.t, rm, nsim = 1, model="BM")
        #print(sims)
        nam=rownames(sims)
        s=as.data.frame(sims)
        simulatedvec=s[,1]
        names(simulatedvec)=nam
        top.all=names(sort(simulatedvec, decreasing = TRUE))
        top.tree_k = top.all[top.all %in% tip.labels]
        top = top.tree_k[1:tips]
        #print("Top:")
        #print(top)
        t=foreground2Tree(top, trees, clade="all", plotTree = F, useSpecies=tip.labels)
        blsum=sum(t$edge.length)
        try_count=try_count+1
        #print(paste("blsum and try_count:", blsum, try_count))
      }
      if(try_count==50){
        #print("Assigning null tree")
        t = tree
        t$edge = NULL
        t$edge.length = NULL
        t$Nnode = NULL
        t$tip.label = NULL
        testcondition=TRUE
        #print("NULL TREE")
      } else{
        t_info = getBinaryPermulationInputsFromTree(t)
        if (!is.null(sisters_list)){
          num_tip_sisters_fake = unlist(t_info$sisters_list)
num_tip_sisters_fake = num_tip_sisters_fake[which(num_tip_sisters_fake %in% tip.labels)]
          num_tip_sisters_fake = length(num_tip_sisters_fake)
          t_depth_order = getDepthOrder(t)
          #testcondition = setequal(sort(t_depth_order), sort(fg_tree_depth_order)) &&
          #  (num_tip_sisters_fake == num_tip_sisters_true)
          testcondition = TRUE
        } else {
          t_depth_order = getDepthOrder(t)
          #testcondition = setequal(sort(t_depth_order), sort(fg_tree_depth_order))
          testcondition = TRUE
        }
      }
    }
  }
  if (plotTreeBool){
    if(!(is.null(t$tip.label))){
        print(t)
        plot(t)
        write.tree(t, "temp.tre", append=T)
    } else{
        write("NULL", "temp.tre", append=T)
    }
  }
  return(t)
}


###SAME AS ABOVE BUT WITHOUT USING THE MASTER TREE###
#Use to test effect of using master tree and midpoint root together
simBinPhenoSSM_noMidpointRoot=function(tree, trees, fg_vec, sisters_list=NULL, pathvec, plotTreeBool=F){
  require(phytools)
  tip.labels = tree$tip.label # the set of species that exist in the gene tree
  ind_fg = which(tip.labels %in% fg_vec) # indices of the observed foreground animals that exist in the gene tree

  if (length(ind_fg) == 0){
    t = tree
    t$edge = NULL
    t$edge.length = NULL
    t$Nnode = NULL
    t$tip.label = NULL
  } else {
    fg_k = tip.labels[ind_fg] # the list of the observed foreground animals that exist in the gene tree
    res = getForegroundInfoClades(fg_k,sisters_list,trees,plotTree=F,useSpecies=tip.labels)
    fg_tree = res$tree
    fg.table = res$fg.sisters.table

    #EK - arbitrarily roots the tree on the first species listed in tip.label vector
    #EK - variable is named "sim.t" to distinguish it from "t" (output of foreground2Trees below)
    sim.t=root.phylo(tree, tree$tip.label[1], resolve.root=T)
    #print(sim.t)
    #EK - hard code a constant rate matrix variance value
    rm=0.02
    #print(rm)

    if (!is.null(sisters_list)){
      fg_tree_info = getBinaryPermulationInputsFromTree(fg_tree)
      num_tip_sisters_true = unlist(fg_tree_info$sisters_list)
      num_tip_sisters_true = num_tip_sisters_true[which(num_tip_sisters_true %in% tip.labels)]
      num_tip_sisters_true = length(num_tip_sisters_true)
      #fg_tree_depth_order = getDepthOrder(fg_tree)
    }

    #EK - I think this needs to be out of the if(!is.null) because later fg_tree_depth_order needs to exist regardless of if sisters_list is null of not
    fg_tree_depth_order = getDepthOrder(fg_tree)
    fgnum = length(which(fg_tree$edge.length == 1))
    if (!is.null(sisters_list)){
      internal = nrow(fg.table)
    } else {
      internal = 0
    }
    tips=fgnum-internal # the number of tips
    #print(paste("Number of foregrounds:", fgnum))
    testcondition=FALSE
    while(!testcondition){
      blsum=0
      #EK - Modified so this tries 50 times to generate a tree with the same number of foregrounds; if it doesn't in that time, move on
      #EK - Based on a previous run, <1% of gene trees took 50 times or more, so this shouldn't get rid of that many
      try_count=0
      #while(blsum!=fgnum){
      while( (blsum!=fgnum) && (try_count < 50) ){
        sims=sim.char(sim.t, rm, nsim = 1, model="BM")
        #print(sims)
        nam=rownames(sims)
        s=as.data.frame(sims)
        simulatedvec=s[,1]
        names(simulatedvec)=nam
        top.all=names(sort(simulatedvec, decreasing = TRUE))
        top.tree_k = top.all[top.all %in% tip.labels]
        top = top.tree_k[1:tips]
        #print("Top:")
        #print(top)
        t=foreground2Tree(top, trees, clade="all", plotTree = F, useSpecies=tip.labels)
        blsum=sum(t$edge.length)
        try_count=try_count+1
        #print(paste("blsum and try_count:", blsum, try_count))
      }
      if(try_count==50){
        #print("Assigning null tree")
        t = tree
        t$edge = NULL
        t$edge.length = NULL
        t$Nnode = NULL
        t$tip.label = NULL
        testcondition=TRUE
        #print("NULL TREE")
      } else{
        t_info = getBinaryPermulationInputsFromTree(t)
        if (!is.null(sisters_list)){
          num_tip_sisters_fake = unlist(t_info$sisters_list)
num_tip_sisters_fake = num_tip_sisters_fake[which(num_tip_sisters_fake %in% tip.labels)]
          num_tip_sisters_fake = length(num_tip_sisters_fake)
          t_depth_order = getDepthOrder(t)
          #testcondition = setequal(sort(t_depth_order), sort(fg_tree_depth_order)) &&
          #  (num_tip_sisters_fake == num_tip_sisters_true)
          testcondition = TRUE
        } else {
          t_depth_order = getDepthOrder(t)
          #testcondition = setequal(sort(t_depth_order), sort(fg_tree_depth_order))
          testcondition = TRUE
        }
      }
    }
  }
  if (plotTreeBool){
    if(!(is.null(t$tip.label))){
        print(t)
        plot(t)
        write.tree(t, "temp.tre", append=T)
    } else{
        write("NULL", "temp.tre", append=T)
    }
  }
  return(t)
}

