
#' @keywords internal
getForegroundsFromBinaryTree=function(tree){
  nameEdgesPerms.tree = nameEdgesPerms(tree)
  edge.length = as.logical(tree$edge.length)
  foregrounds = nameEdgesPerms.tree[edge.length]
  ind.tip = which(foregrounds != "")
  foregrounds = foregrounds[ind.tip]
  return(foregrounds)
}

#' @keywords internal
nameEdgesPerms=function(tree){
  if (is.null(tree$tip.label)) {
    nn = NULL
  } else {
    nn=character(nrow(tree$edge))
    iim=match(1:length(tree$tip.label), tree$edge[,2])
    nn[iim]=tree$tip.label
  }
  nn
}

#NOTE: I didn't change anything in getForegoundInfoClades; I just copied it here to add print statements to figure out where permulations were stalling when I was troubleshooting
#'Generates a binary phenotype tree and foreground clades information using the list of tip foreground animals, the presence of foreground common ancestors, and their phylogenetic relationships
#' @param fg_vec A vector containing the tip foreground species
#' @param sisters_list A list containing pairs of "sister species" in the foreground set (put NULL if empty)
#' @param trees treesObj from \code{\link{readTrees}}
#' @param plotTree A boolean indicator for plotting the output tree (default=FALSE)
#' @param useSpecies An array containing the tip labels in the output tree
#' @return output.list A list containing 1) "tree" = a binary phenotype tree corresponding to the input information, 2) "fg.sisters.table" = a table containing all sister species in the foreground set
#' @export
getForegroundInfoClades=function(fg_vec,sisters_list=NULL,trees,plotTree=T,useSpecies=NULL){
  if (length(useSpecies)==0){
    useSpecies = trees$masterTree$tip.label
  }

  if (is.null(sisters_list)){
    fg.sisters.table=NULL
    fg_tree = foreground2Tree(fg_vec, trees, plotTree=F, clade="terminal", useSpecies=useSpecies)
  } else {
    # Start with a temp phenotype tree assuming that all internal nodes are foregrounds
    #print("Starting foreground2Tree")
    fg_tree = foreground2Tree(fg_vec,trees,plotTree=F,clade="all",useSpecies=useSpecies)
    #print("Done with foreground2Tree")
    #write.tree(fg_tree, "temp.tre", append=T)
    edge = fg_tree$edge
    edge.length=fg_tree$edge.length

    ind.fg.edge = which(edge.length == 1)
    #print(ind.fg.edge)
    nodeIds.fg.edge = edge[ind.fg.edge,] # all foreground edges in the temp tree
    #print(nodeIds.fg.edge)

    tip.sisters = vector("integer",length=0)
    for (i in 1:length(sisters_list)){
      sisters = sisters_list[[i]]
      #print(sisters)
      nodeId.sisters = which(useSpecies %in% sisters)
      if (length(nodeId.sisters)>0){
        tip.sisters = c(tip.sisters,nodeId.sisters)
      }
    }
    #print(tip.sisters)

    # Find and correct the pairs
    fg.sisters.table = matrix(nrow=0,ncol=2)
    colnames(fg.sisters.table) = c("species1","species2")
    if (length(as.vector(nodeIds.fg.edge)) > 2){
      all.nodeId.ca = sort(nodeIds.fg.edge[,1])
      count_all_nodeId_ca = table(all.nodeId.ca)
      unq.nodeId.ca = unique(all.nodeId.ca)
      fg_ca = vector("integer",length=0) # node IDs of the common ancestor foregrounds
      nodes_addressed = NULL
      #print("Starting while loop")
      #print(unq.nodeId.ca)
      while (length(unq.nodeId.ca) != length(nodes_addressed)){
        nodeId.ca = sort(all.nodeId.ca[which(!(all.nodeId.ca %in% nodes_addressed))])
        #print(nodeId.ca)
        if (length(nodeId.ca) == 1){
          nodes_addressed = c(nodes_addressed, nodeId.ca)
        } else {
          for (nn in 1:(length(nodeId.ca)-1)){

            if (nodeId.ca[nn] == nodeId.ca[nn+1]){
              nodeId.desc = nodeIds.fg.edge[which(nodeIds.fg.edge[,1]==nodeId.ca[nn]),2]
              #print(nodeId.desc)  
            if (length(which(nodeId.desc %in% tip.sisters)) > 0){
                fg_ca = c(fg_ca,nodeId.ca[nn])
                #print(fg_ca)
                fg.sisters.table = rbind(fg.sisters.table, nodeIds.fg.edge[which(nodeIds.fg.edge[,1]==nodeId.ca[nn]),2])
                #print(fg.sisters.table)
                #print("here1")
                nodes_addressed = c(nodes_addressed, nodeId.ca[nn])
              } else {
                if (length(which(fg_tree$tip.label[nodeId.desc] %in% fg_vec)) == 2){
                  fg_tree$edge.length[which(edge[,2]==nodeId.ca[nn])] = 0
                  #print("here2")
                  nodes_addressed = c(nodes_addressed, nodeId.ca[nn])
                } else {
                  if (length(which(nodeId.desc %in% nodes_addressed)) == 2){
                    fg_ca = c(fg_ca,nodeId.ca[nn])
                    fg.sisters.table = rbind(fg.sisters.table, nodeIds.fg.edge[which(nodeIds.fg.edge[,1]==nodeId.ca[nn]),2])
                    #print("here3")
                    nodes_addressed = c(nodes_addressed, nodeId.ca[nn])
                  }
                  #print("no else here...")
                  #print(length(which(nodeId.desc %in% nodes_addressed)))
                }
              }
            } else {
              nodeId.desc = nodeIds.fg.edge[which(nodeIds.fg.edge[,1]==nodeId.ca[nn]),2]
              if (length(nodeId.desc) == 2){
                if (nodeId.ca[nn] != nodeId.ca[nn-1]){
                  fg_tree$edge.length[which(edge[,2] == nodeId.ca[nn])] = 0
                  #print("here4") 
                  nodes_addressed = c(nodes_addressed, nodeId.ca[nn])
                  nodes_addressed = unique(nodes_addressed)
                }
              } else {
                #print("here5")
                nodes_addressed = c(nodes_addressed, nodeId.ca[nn])
              }
            }
          }
        }
        #print(nodes_addressed)
      }
      rownames(fg.sisters.table) = fg_ca
    }
  }
  if (plotTree==T){
    plot(fg_tree)
  }
  output.list = list("fg.sisters.table"=fg.sisters.table,"tree"=fg_tree)
  output.list
}


#NOTE: Below are some functions I didn't change but copied because they are marked as "internal" so the functions that call them can't find them unless they're in the same file

#' @keywords internal
getSpeciesMembershipStats = function(tree,masterTree,foregrounds){
  master.tips = masterTree$tip.label
  tips = tree$tip.label
  spec_membership = which(master.tips %in% tips)
  fg_membership = which(foregrounds %in% tips)

  num_spec = length(spec_membership)
  num_fg = length(fg_membership)

  spec.members = rep(0,length(master.tips))
  spec.members[spec_membership] = 1
  spec.members = toString(spec.members)

  fg.members = rep(0,length(foregrounds))
  fg.members[fg_membership] = 1
  fg.members = toString(fg.members)

  df = data.frame("num.fg"=as.character(num_fg), "num.spec"=as.character(num_spec), "spec.members"=spec.members, "fg.members"=fg.members)

  return(df)
}

#' @keywords internal
groupTrees = function(spec.members){
  unique.trees = unique(spec.members)
  ind.tree.groups = lapply(unique.trees,findGroupedTrees,spec.members=spec.members)
  ind.unique.trees = lapply(ind.tree.groups, function(i) min(i))
  output.list = list()
  output.list$ind.unique.trees = ind.unique.trees
  output.list$ind.tree.groups = ind.tree.groups
  return(output.list)
}

#' @keywords internal
findGroupedTrees = function(unique.tree,spec.members){
  ind.grouped.trees = which(spec.members == unique.tree)
  return(ind.grouped.trees)
}

#'Calculates the clade mappings between the gene tree and the master tree (with the complete topology)
#' @param gene_tree A binary phenotype tree of a gene
#' @param treesObj treesObj from \code{\link{readTrees}}
#' @return output.map A list containing a dataframe of clades mapping
#' @export
matchAllNodesClades=function(gene_tree, treesObj){
  foregrounds = getForegroundsFromBinaryTree(gene_tree)
  tree1 = foreground2Tree(foregrounds,treesObj,clade="all",plotTree = F)

  map=matchNodesInject_mod(tree1,treesObj$masterTree)
  map=map[order(map[,1]),]
  #map

  output.map = list()
  output.map[[1]]=map
  output.map
}

#' @keywords internal
compareClades=function(clade2index,desc.tree2,clade1){
  output=NA
  clade2 = desc.tree2[[clade2index]]
  if (all(clade1 %in% clade2)){
    output = as.numeric(names(desc.tree2)[clade2index])
  }
  return(output)
}

#' @keywords internal
findMappedClade2Node=function(desc.tr1.index,desc.tr2.index.list,desc.tree1,desc.tree2){
  clade1 = desc.tree1[[desc.tr1.index]]
  mapped.clades.list = lapply(desc.tr2.index.list,compareClades,desc.tree2=desc.tree2,clade1=clade1)
  mapped.clades = unlist(mapped.clades.list)
  mapped.clades.nonNA = mapped.clades[!is.na(mapped.clades)]
  return(max(mapped.clades.nonNA))
}

#' @keywords  internal
matchNodesInject_mod=function (tr1, tr2){
  if(length(tmpsp<-setdiff(tr1$tip.label, tr2$tip.label))>0){
    #stop(paste(paste(tmpsp, ","), "in tree1 do not exist in tree2"))
    stop(c("The following species in tree1 do not exist in tree2: ",paste(tmpsp, ", ")))
  }
  commontiplabels <- intersect(tr1$tip,tr2$tip)
  if(RF.dist(pruneTree(tr1,commontiplabels),pruneTree(tr2,commontiplabels))>0){
    stop("Discordant tree topology detected - gene/trait tree and treesObj$masterTree have irreconcilable topologies")
  }
  #if(RF.dist(tr1,tr2)>0){
  #  stop("Discordant tree topology detected - trait tree and treesObj$masterTree have irreconcilable topologies")
  #}

  toRm=setdiff(tr2$tip.label, tr1$tip.label)
  desc.tr1 <- lapply(1:tr1$Nnode + length(tr1$tip), function(x) extract.clade(tr1,
                                                                              x)$tip.label)
  names(desc.tr1) <- 1:tr1$Nnode + length(tr1$tip)
  desc.tr2 <- lapply(1:tr2$Nnode + length(tr2$tip), function(x) extract.clade(tr2,
                                                                              x)$tip.label)
  names(desc.tr2) <- 1:tr2$Nnode + length(tr2$tip)
  Nodes <- matrix(NA, length(desc.tr1), 2, dimnames = list(NULL,
                                                           c("tr1", "tr2")))
  Nodes[,1] = as.numeric(names(desc.tr1))
  desc.tr1.index.list = as.list(1:length(desc.tr1))
  desc.tr2.index.list = as.list(1:length(desc.tr2))


  mapped.clade.list = lapply(desc.tr1.index.list,findMappedClade2Node,desc.tr2.index.list,desc.tree1=desc.tr1,desc.tree2=desc.tr2)
  Nodes[,2] = unlist(mapped.clade.list)

  iim=match(tr1$tip.label, tr2$tip.label)
  Nodes=rbind(cbind(1:length(tr1$tip.label),iim),Nodes)
  if(any(table(Nodes[,2])>1)){
    stop("Incorrect pseudorooting detected - use fixPseudoroot() function to correct trait tree topology")
  }

  Nodes
}

#' @keywords internal
calculatePermulatedPaths=function(permulated.trees,map,treesObj){
  permulated.paths=lapply(permulated.trees,tree2PathsCladesQuiet,trees=treesObj)
  output.list = list()
  output.list[[1]] = permulated.paths
  output.list
}


#EK - Edited to NOT print "Species from tree not present in master tree or useSpecies:"
#' @keywords internal
tree2Paths_mapQuiet=function(tree, map, treesObj, binarize=NULL, useSpecies=NULL){
  if (class(tree)[1]=="phylo"){
    stopifnot(class(tree)[1]=="phylo")
    stopifnot(class(treesObj)[2]=="treesObj")

    if (is.null(tree$tip.label)){
      vals=as.double(rep(NA,length(treesObj$ap$dist)))
    } else {
      foregrounds = getForegroundsFromBinaryTree(tree)
      tree = foreground2Tree(foregrounds,treesObj,clade="all",plotTree = F)


      isbinarypheno <- sum(tree$edge.length %in% c(0,1)) == length(tree$edge.length) #Is the phenotype tree binary or continuous?
      if (is.null(binarize)) { #unless specified, determine default for binarize based on type of phenotype tree
        if (isbinarypheno) {
          binarize = T #default for binary phenotype trees: set all positive paths = 1
        } else {
          binarize = F #default for continuous phenotype trees: do not convert to binary
        }
      }

      #unroot if rooted
      if (is.rooted(tree)) {
        tree = unroot(tree)
      }

      #reduce tree to species in master tree and useSpecies
      sp.miss = setdiff(tree$tip.label, union(treesObj$masterTree$tip.label, useSpecies))
      #if (length(sp.miss) > 0) {
      #  message(paste0("Species from tree not present in master tree or useSpecies: ", paste(sp.miss,
      #                                                                                       collapse = ",")))
      #}

      if (!is.null(useSpecies)) {
        tree = pruneTree(tree, intersect(intersect(tree$tip.label, treesObj$masterTree$tip.label), useSpecies))
      } else {
        tree = pruneTree(tree, intersect(tree$tip.label, treesObj$masterTree$tip.label))
      }
      treePaths=allPaths(tree)

      #remap the nodes
      treePaths$nodeId[,1]=map[treePaths$nodeId[,1],2 ]
      treePaths$nodeId[,2]=map[treePaths$nodeId[,2],2 ]

      #indices for which paths to return
      ii=treesObj$ap$matIndex[(treePaths$nodeId[,2]-1)*nrow(treesObj$ap$matIndex)+treePaths$nodeId[,1]]

      vals=double(length(treesObj$ap$dist))
      vals[]=NA
      vals[ii]=treePaths$dist
      if(binarize){
        if(isbinarypheno) {
          vals[vals>0]=1
        } else {
          mm=mean(vals)
          vals[vals>mm]=1
          vals[vals<=mm]=0
        }
      }
    }
  } else {
    vals=as.double(rep(NA,length(treesObj$ap$dist)))
  }

  vals
}

#EK - Edited to NOT print "Species from tree not present in master tree or useSpecies:"
#' Generate a phenotype paths vector from a phenotype tree
#'
#' \code{tree2Paths} generates a phenotype paths vector matching the treesObject
#'     from a tree where branches specify phenotypes.
#'
#' The tree topology of the phenotype tree must match that of the master tree within the treesObject.
#'
#' @param tree A phenotype tree, with branch length encoding a phenotype.
#' @param treesObj A treesObject created by \code{\link{readTrees}}
#' @param binarize Force binary path representation. Default action depends upon the type of data within the phenotype tree
#'     (binary or continuous).
#'     \itemize{
#'    \item If binary (all branch lengths == 0 or 1): Sets all positive path values to 1. Useful if the tree has non-zero branch lengths
#'        for an internal branch or branches; otherwise, values are simply added along branches when calculating paths.
#'        Default behavior: binarize = TRUE.
#'    \item If continuous (not all branch lengths == 0 or 1): Sets all path values > the mean to 1 and all those <= the mean to 0.
#'        Converts a continuous phenotype to a binary phenotype, with state determined by comparison to the mean across all paths.
#'        Default behavior: binarize = FALSE.
#'        }
#' @param useSpecies Give only a subset of the species to use for ancestral state reconstruction
#' (e.g., only those species for which the trait can be reliably determined).
#' @return A vector of length equal to the number of paths in treesObj
#' @export
tree2Paths=function(tree, treesObj, binarize=NULL, useSpecies=NULL, categorical = F){
  stopifnot(class(tree)[1]=="phylo")
  stopifnot(class(treesObj)[2]=="treesObj")

  isbinarypheno <- sum(tree$edge.length %in% c(0,1)) == length(tree$edge.length) #Is the phenotype tree binary or continuous?
  if (is.null(binarize)) { #unless specified, determine default for binarize based on type of phenotype tree
    if (isbinarypheno) {
      binarize = T #default for binary phenotype trees: set all positive paths = 1
    } else {
      binarize = F #default for continuous phenotype trees: do not convert to binary
    }
  }
  #unroot if rooted
  if (is.rooted(tree)) {
    tree = unroot(tree)
  }

  #reduce tree to species in master tree and useSpecies
  sp.miss = setdiff(tree$tip.label, union(treesObj$masterTree$tip.label, useSpecies))
  if (length(sp.miss) > 0) {
    #message(paste0("Species from tree not present in master tree or useSpecies: ", paste(sp.miss,
    #                                                                                     collapse = ",")))
  }
  if (!is.null(useSpecies)) {
    tree = pruneTree(tree, intersect(intersect(tree$tip.label, treesObj$masterTree$tip.label), useSpecies))
  } else {
    tree = pruneTree(tree, intersect(tree$tip.label, treesObj$masterTree$tip.label))
  }

  treePaths=allPaths(tree, categorical = categorical)
  map=matchAllNodes(tree,treesObj$masterTree)

  #remap the nodes
  treePaths$nodeId[,1]=map[treePaths$nodeId[,1],2 ]
  treePaths$nodeId[,2]=map[treePaths$nodeId[,2],2 ]

  #indices for which paths to return
  ii=treesObj$ap$matIndex[(treePaths$nodeId[,2]-1)*nrow(treesObj$ap$matIndex)+treePaths$nodeId[,1]]

  vals=double(length(treesObj$ap$dist))
  vals[]=NA
  vals[ii]=treePaths$dist
  if(binarize){
    if(isbinarypheno) {
      vals[vals>0]=1
    } else {
      mm=mean(vals)
      vals[vals>mm]=1
      vals[vals<=mm]=0
    }
  }
  vals
}



#EK - Edited to NOT print "Species from tree not present in master tree or useSpecies:"
#' Creates a binary trait tree from a set of foreground species.
#' @param foreground. A character vector containing the foreground species
#' @param treesObj A treesObj created by \code{\link{readTrees}}
#' @param collapse2anc Put all the weight on the ancestral branch when the trait appears on a whole clade
#' (redundant to "clade", kept for backwards compatibility)
#' @param plotTree Plot a tree representation of the result
#' @param wholeClade Whether to implement the weighted edge option across
#' all members of a foreground clade (redundant to "clade", kept for backwards compatibility)
#' @param clade A character string indicating which branches within the clade
#' containing the foreground species should be set to foreground. Must be one
#' of the strings "ancestral", "terminal", "all".
#' @param transition A character string indicating whether transitions between background and foreground branches
#' are "bidirectional" or "unidirectional" (no foreground to background transitions, the default)
#' @param weighted if set to TRUE weights foreground edges belonging to the same clade such that their branch lengths sum up to 1 (only done for clade options "all" and "terminal").
#' @param useSpecies Give only a subset of the species to use for ancestral state reconstruction
#' (e.g., only those species for which the trait can be reliably determined).
#' @return A tree with edge.lengths representing phenotypic states
#' @export
foreground2Tree = function(foreground,treesObj, plotTree=T, clade=c("ancestral","terminal","all"), weighted = F, transition = "unidirectional", useSpecies=NULL){
  clade <- match.arg(clade)
  res = treesObj$masterTree
  if (!is.null(useSpecies)) {
    sp.miss = setdiff(res$tip.label, useSpecies)
    if (length(sp.miss) > 0) {
      #message(paste0("Species from master tree not present in useSpecies: ", paste(sp.miss,
      #                                                                            collapse = ",")))
    }
    useSpecies = intersect(useSpecies, res$tip.label)
    res = pruneTree(res, useSpecies)
  } else {
    useSpecies = res$tip.label
  }
  foreground = intersect(foreground, useSpecies)
  res$edge.length <- rep(0,length(res$edge.length))
  if(clade == "terminal"){
    res$edge.length[nameEdges(res) %in% foreground] = 1
    names(res$edge.length) = nameEdges(res)
  }else if(clade == 'ancestral'){
    weighted = F
    if(transition == 'bidirectional'){
      res <- inferBidirectionalForegroundClades(res,foreground,ancestralOnly =T)
    }else{
      res <- inferUnidirectionalForegroundClades(res,foreground,ancestralOnly =T)
    }
  }else{
    if(transition == 'bidirectional'){
      res <- inferBidirectionalForegroundClades(res,foreground,ancestralOnly = F)
    }else{
      res <- inferUnidirectionalForegroundClades(res,foreground,ancestralOnly = F)
    }
  }
  if(weighted){
    if(clade == 'all'){
      tobeweighted <- rep(TRUE,length(res$edge.length))
      tobeweighted[res$edge.length == 0] <- FALSE
      while(sum(tobeweighted)>0){
        edgetodo <- which(tobeweighted == T)[1]
        clade.down.edges = getAllCladeEdges(res, edgetodo)
        if(length(clade.down.edges) > 1){
          clade.edges = c(clade.down.edges, edgetodo)
          clade.edges.toweight <- clade.edges[res$edge.length[clade.edges] == 1]
          res$edge.length[clade.edges.toweight] <- 1.0/(length(clade.edges.toweight))
          tobeweighted[clade.edges] <- FALSE
        } else{
          tobeweighted[clade.down.edges] <- FALSE
        }
      }
    } else if(clade == 'terminal'){
      tobeweightededgeterminalnode <- unique(res$edge[(res$edge[,2] %in% c(1:length(res$tip.label))),1])
      tobeweighted <- setdiff(match(tobeweightededgeterminalnode, res$edge[,2]), NA)
      for(edgetodo in tobeweighted){
        clade.down.edges = getAllCladeEdges(res, edgetodo)
        if(all(res$edge.length[clade.down.edges]==1)){
          res$edge.length[clade.down.edges] <- 0.5
        }
      }
    }
  }
  if(plotTree){
    res2=res
    mm=min(res2$edge.length[res2$edge.length>0])
    res2$edge.length[res2$edge.length==0]=max(0.02,mm/20)
    plot(res2, main = paste0("Clade: ",clade,'\nTransition: ',transition,'\nWeighted: ',weighted), cex = 0.5)
    if(weighted){
      labs <- round(res$edge.length,3)
      labs[labs == 0] <- NA
      edgelabels(labs, col = 'black', bg = 'transparent', adj = c(0.5,-0.5),cex = 0.4,frame='n')
    }
  }
  res
}

#EK - Edited to NOT print "Species from tree not present in master tree or useSpecies:"
#' main RER computation function
#' @param treesObj A treesObj created by \code{\link{readTrees}}
#' @param a cutoff value for branch lengths bellow which the branch lengths will be discarded, very data dependent but should roughly correspond to 0 or 1 sequence change on that branch. If left NULL this whill be set to the bottom 0.05 quantile. Set to 0 for no cutoff.
#' @param transform The transformation to apply to the trees branch values before computing relative rates. Available options are sqrt and log, sqrt is recommended.
#' @param weighted Use weighted regression to compute relative rates, meant to correct for the non-constant mean-variance relationship in evolutionary rate data.
#' @param useSpecies Give only a subset of the species to use for RER calculation. Some times excluding unusually long branches can provide more stable results
#' @param min.sp The minimum number of species needed to compute RER
#' @param scale Scale relative rates internally for each species subset. Increases computation time with little apparent benefit. Better to scale the final matrix.
#' @param doOnly The index of a specific tree in the treesObj to calculate RER for. Useful if a single result is needed quickly.
#' @param maxT The maximum number of trees to compute results for. Since this function takes some time this is useful for debugging.
#' @param plot Whether to plot the output of the correction for mean-variance relationship.
#' @return A numer of trees by number of paths matrix of relative evolutionary rates. Only an independent set of paths has non-NA values for each tree.
#' @export
getAllResiduals=function(treesObj, cutoff=NULL, transform="sqrt", weighted=T,  useSpecies=NULL,  min.sp=10, scale=T,  doOnly=NULL, maxT=NULL, scaleForPproj=F, mean.trim=0.05, plot=T){

  if(is.null(cutoff)){
    cutoff=quantile(treesObj$paths, 0.05, na.rm=T)
    message(paste("cutoff is set to", cutoff))
  }
  if (weighted){
    weights=computeWeightsAllVar(treesObj$paths, transform=transform, plot=plot)
    residfunc=fastLmResidMatWeighted
  }
  else{
    residfunc=fastLmResidMat
  }
  # residfunc=naresid

  if (is.null(useSpecies)){
    useSpecies=treesObj$masterTree$tip.label
    #mappedEdges=trees$mappedEdges
  }
  if(is.null(maxT)){
    maxT=treesObj$numTrees
  }
  if(transform!="none"){
    transform=match.arg(transform,c("sqrt", "log"))
    transform=get(transform)
  }
  else{
    transform=NULL
  }



  #cm is the names of species that are included in useSpecies and the master tree
  cm=intersect(treesObj$masterTree$tip.label, useSpecies)
  sp.miss = setdiff(treesObj$masterTree$tip.label, useSpecies)
  if (length(sp.miss) > 0) {
    #message(paste0("Species from master tree not present in useSpecies: ", paste(sp.miss,
    #                                                                             collapse = ",")))

  }

  rr=matrix(nrow=nrow(treesObj$paths), ncol=ncol(treesObj$paths))

  #maximum number of present species
  maxn=rowSums(treesObj$report[,cm])

  if(is.null(doOnly)){
    doOnly=1
  }
  else{
    maxT=1
  }
  skipped=double(nrow(rr))
  skipped[]=0

  for (i in doOnly:(doOnly+maxT-1)){

    if(sum(!is.na(rr[i,]))==0&&!skipped[i]==1){


      #get the ith tree
      tree1=treesObj$trees[[i]]

      #get the common species, prune and unroot
      both=intersect(tree1$tip.label, cm)
      if(length(both)<min.sp){
        next
      }
      tree1=unroot(pruneTree(tree1,both))

      #do the same for the refTree


      #find all the genes that contain all of the species in tree1
      allreport=treesObj$report[,both]
      ss=rowSums(allreport)
      iiboth=which(ss==length(both)) #this needs to be >1
      if (length(iiboth) < 2) {
        message(paste("Skipping i =",i,"(no other genes with same species set)"))
        next
      }

      nb=length(both)
      ai=which(maxn[iiboth]==nb)


      message(paste("i=", i))


      if(T){

        ee=edgeIndexRelativeMaster(tree1, treesObj$masterTree)

        ii= treesObj$matIndex[ee[, c(2,1)]]

        allbranch=treesObj$paths[iiboth,ii]
        if (is.null(dim(allbranch))) {
          message(paste("Issue with gettiing paths for genes with same species as tree",i))
          return(list("iiboth"=iiboth,"ii"=ii))
        }

        if(weighted){
          allbranchw=weights[iiboth,ii]
        }
        if(scaleForPproj){
          nv=apply(scaleMatMean(allbranch), 2, mean, na.rm=T, trim=mean.trim)
        }
        else{
          nv=apply(allbranch, 2, mean, na.rm=T, trim=mean.trim)
        }

        iibad=which(allbranch<cutoff)
        #don't scale
        #allbranch=scaleMat(allbranch)
        if(!is.null(transform)){
          nv=transform(nv)
          allbranch=transform(allbranch)
        }
        allbranch[iibad]=NA




        if(!scale){
          if(!weighted){
            proj=residfunc(allbranch[ai, ,drop=F], model.matrix(~1+nv))

          }
          else{

            proj=residfunc(allbranch[ai, ,drop=F], model.matrix(~1+nv), allbranchw[ai, ,drop=F])

          }
        }

        else{

          if(!weighted){
            proj=residfunc(allbranch[, ,drop=F], model.matrix(~1+nv))
          }
          else{

            proj=residfunc(allbranch[, ,drop=F], model.matrix(~1+nv),allbranchw)
          }

          proj=scale(proj, center = F)[ai, , drop=F]

        }


        #we have the projection



        rr[iiboth[ai],ii]=proj

      }

    }}
  message("Naming rows and columns of RER matrix")
  rownames(rr)=names(treesObj$trees)
  colnames(rr)=namePathsWSpecies(treesObj$masterTree)
  rr
}


#EK - Edited to NOT print "Species from tree not present in master tree or useSpecies:"
#' Creates a categorical trait tree from a set of tip species.
#'@param tipvals the trait/phenotype/character value at the tip, \code{names(tip.vals)} should match some of the \code{mastertree$tip.label}, though a perfect match is not required
#'@param treesObj A treesObj created by \code{\link{readTrees}}
#'@param useSpecies Give only a subset of the species to use for ancestral state reconstruction
#' (e.g., only those species for which the trait can be reliably determined).
#'@param model Specifies what rate model to use
#'@param plot Plots a phenotype tree
#'@param anctrait The trait to use for all ancestral species instead of inferring ancestral states if not NULL. The default is NULL.
#'@return A tree with edge.length representing phenotype states
#'@export
char2TreeCategorical <- function (tipvals, treesObj, useSpecies = NULL,
                                     model = "ER", root_prior = "auto",
                                     plot = FALSE, anctrait = NULL)
{
  # get the master tree and prune to include useSpecies/species with phenotype data
  mastertree = treesObj$masterTree
  if (!is.null(useSpecies)) {
    sp.miss = setdiff(mastertree$tip.label, useSpecies)
    if (length(sp.miss) > 0) {
      #message(paste0("Species from master tree not present in useSpecies: ",
      #               paste(sp.miss, collapse = ",")))
    }
    useSpecies = intersect(mastertree$tip.label, useSpecies)
    mastertree = pruneTree(mastertree, useSpecies)
  }
  else {
    mastertree = pruneTree(mastertree, intersect(mastertree$tip.label,
                                                 names(tipvals)))
  }
  # use ASR to infer phenotype tree
  if (is.null(anctrait)) {

    tipvals <- tipvals[mastertree$tip.label]
    intlabels <- map_to_state_space(tipvals)
    print("The integer labels corresponding to each category are:")
    print(intlabels$name2index)

    ancliks = getAncLiks(mastertree, intlabels$mapped_states, rate_model = model,
                         root_prior = root_prior)

    states = rep(0, nrow(ancliks))
    for (i in 1:length(states)) {
      states[i] = which.max(ancliks[i,])
    }
    states = c(intlabels$mapped_states, states)
    tree = mastertree
    tree$edge.length = states[tree$edge[, 2]]

    # convert to binary tree if necessary, plot, & return tree
    if(length(unique(tipvals)) == 2) {
      if(sum(! unique(tipvals) %in% c(TRUE,FALSE)) > 0) { # check that the two categories are TRUE/FALSE
        message("Returning categorical tree for binary phenotype because phenotype values are not TRUE/FALSE")
      } else {
        tree$edge.length = ifelse(tree$edge.length == 2, 1, 0)
        print("There are only 2 categories: returning a binary phenotype tree.")
        if (plot) {
          plotTree(tree)
        }
        return(tree)
      }
    }
    if (plot) {
      plotTreeCategorical(tree, category_names = intlabels$state_names,
                          master = mastertree, node_states = states)
    }
    return(tree)
  }
  else {
    if (length(unique(tipvals)) <= 2) {
      fgspecs <- names(tipvals)[tipvals != anctrait]
      res <- foreground2Tree(fgspecs, treesObj, plotTree = plot,
                             clade = "terminal", useSpecies = useSpecies)
      print("There are only 2 categories: returning a binary phenotype tree.")
      if(plot) {
        plotTree(res)
      }
      return(res)
    }
    else {
      tipvals <- tipvals[mastertree$tip.label]
      intlabels <- map_to_state_space(tipvals)
      j <- which(intlabels$state_names == anctrait)
      if (length(j) < 1) {
        warning("The ancestral trait provided must match one of the traits in the phenotype vector.")
      }
      res = mastertree
      res$edge.length <- rep(j, length(res$edge.length))
      traits <- intlabels$state_names
      for (trait in traits) {
        if (trait == anctrait) {
          next
        }
        i <- which(intlabels$state_names == trait)
        res$edge.length[nameEdges(res) %in% names(tipvals)[tipvals ==
                                                             trait]] = i
      }
      names(res$edge.length) = nameEdges(res)
      if (plot) {
        # get states for plotting
        states = res$edge.length[order(res$edge[,2])]
        states = c(j, states) # add root since it's not included in res$edge.length (no edge leading to the root)
        plotTreeCategorical(res, category_names = traits,
                            master = treesObj$masterTree,
                            node_states = states)
      }
      print("Category names are mapped to integers as follows:")
      print(intlabels$name2index)
      return(res)
    }
  }
}

#' @keywords internal
inferUnidirectionalForegroundClades <- function(tree, fgd = NULL, ancestralOnly = F){
  finaltree <- tree
  finaltree$edge.length <- rep(0, length(tree$edge.length))
    finaltree$edge.length[nameEdges(finaltree) %in% fgd] <- 1
  #figure out node depth - terminal nodes have depth of 1; higher numbers indicate ancestral nodes;
  nodedepths <- node.depth(finaltree)
  edgeterminalnodedepths <- nodedepths[finaltree$edge[,2]]
  #going from 1-away from terminal ancestral branch to the base of the tree, figure out branches where all downstream clades are foreground
  for(inode in sort(unique(edgeterminalnodedepths))[-1]){
    edgesToDo <- which(edgeterminalnodedepths == inode)
    for(edgeindex in edgesToDo){
      clade.edges = getAllCladeEdges(finaltree, edgeindex)
      if(all(finaltree$edge.length[clade.edges]==1)){
        finaltree$edge.length[edgeindex] <- 1
      }
    }
  }
  if(ancestralOnly){
    for(edgeii in 1:length(finaltree$edge.length)){
      if(finaltree$edge.length[edgeii] == 1){
        if(nameEdges(finaltree)[edgeii]==""){
          clade.edges = setdiff(getAllCladeEdges(finaltree, edgeii), edgeii)
          finaltree$edge.length[clade.edges] <- 0
        }
      }
    }
  }
  finaltree
}

#' @keywords internal
nameEdges=function(tree){
  nn=character(nrow(tree$edge))
  iim=match(1:length(tree$tip.label), tree$edge[,2])
  nn[iim]=tree$tip.label
  nn
}

getAllCladeEdges=function(tree, AncEdge){
  node=tree$edge[AncEdge,2]
  #get descendants
  iid=getDescendants(tree, node)
  #find their edges
  iim=match(iid, tree$edge[,2])
  iim
}

#' @keywords internal
getDepthOrder=function(fgTree){
  unq_edge_lengths = unique(fgTree$edge.length)
  if (length(which(!(unq_edge_lengths %in% c(0,1)))) > 0){
    stop('Phenotype must be binary.')
  }
  all_edges = fgTree$edge
  num_tip_species = length(fgTree$tip.label)

  idx_fg_branches = which(fgTree$edge.length == 1)
  if (length(idx_fg_branches)==1){
	fg_edges = fgTree$edge[idx_fg_branches,]
	fg_edges = t(as.data.frame(fg_edges))
  } else {
	fg_edges = fgTree$edge[idx_fg_branches,]
	tip_fg_edges = fg_edges[which(fg_edges[,2] <= num_tip_species),]
	tip_foregrounds = fgTree$tip.label[tip_fg_edges[,2]]
	node_fg_edges = fg_edges[which(fg_edges[,2] > num_tip_species),]
  }

  idx_node_edges = which(fg_edges[,2] > num_tip_species)
  if (length(idx_node_edges) == 1){
    node_fg_edges = fg_edges[which(fg_edges[,2] > num_tip_species),]
    node_fg_edges = t(as.data.frame(node_fg_edges))
  }
  if (length(idx_node_edges) == 0) {
    sisters_list = NULL
    depth_order = NULL
  } else {
    #node_fg_edges = fg_edges[which(fg_edges[,2] > num_tip_species),]
    daughters_info_list = list()
    parents = NULL
    for (i in 1:nrow(node_fg_edges)){
      edge_i = node_fg_edges[i,]
      # find the daughters of this node
      idx_daughters_i = which(all_edges[,1] == edge_i[2])
      daughter_edges = all_edges[idx_daughters_i,]
      daughters_info_list[[i]] = daughter_edges[,2]
      parents = c(parents, edge_i[2])
    }
    names(daughters_info_list) = parents
    ### write something to order the branches based on depth
    tip_fg_ids = tip_fg_edges[,2]
    depth_order = rep(NA, length(daughters_info_list))
    names(depth_order) = names(daughters_info_list)
    order_assigned = NULL
    while(length(which(is.na(depth_order))) > 0){
      idx_na = which(is.na(depth_order))
      if (length(idx_na) > 0){
        for (j in 1:length(idx_na)){
          idx_na_j = idx_na[j]
          parent_j = parents[idx_na_j]
          daughters_j = daughters_info_list[[idx_na_j]]
          num_tip_daughters = length(which(daughters_j %in% tip_fg_ids))
          if (num_tip_daughters == 2){
            depth_order[idx_na_j] = 1
            order_assigned = c(order_assigned, parent_j)
          } else if (num_tip_daughters==1){
            node_daughter = daughters_j[which(daughters_j > length(fgTree$tip.label))]
            if (node_daughter %in% order_assigned){
              depth_order[idx_na_j] = depth_order[as.character(node_daughter)] + 1
              order_assigned = c(order_assigned, parent_j)
            }
          } else if (num_tip_daughters==0){
            node_daughters = daughters_j[which(daughters_j > length(fgTree$tip.label))]
            if (length(which(node_daughters %in% order_assigned)) == 2){
              node_daughters_depths = depth_order[as.character(node_daughters)]
              depth_order[idx_na_j] = max(node_daughters_depths) + 1
              order_assigned = c(order_assigned, parent_j)
            }
          }
        }
      }
    }
  }
  depth_order
}

#' @keywords internal
extractCorResults=function(corResultsList,numperms,mode="Rho"){
  output = matrix(NA,nrow=length(corResultsList),ncol=length(corResultsList[[1]]))
  for (i in 1:length(corResultsList)){
    gene = corResultsList[[i]]
    output[i,] = extractPermulationResults(gene,numperms,mode)
  }
  return(output)
}

#' @keywords internal
extractPermulationResults=function(gene,numperms,mode="Rho"){
  table_perm = lapply(gene,linearizeCorResults)
  df_perm = do.call(rbind,table_perm)
  if (mode=="Rho"){
    output=df_perm[,1]
  } else if (mode=="P"){
    output=df_perm[,3]
  }
  return(output)
}

#' @keywords internal
linearizeCorResults=function(cor_result){
  vec.cor = unlist(cor_result)
  return(vec.cor)
}

#' @keywords internal
allPaths=function(tree, categorical = F){
  dd=dist.nodes(tree) # returns a matrix with col_names and row_names denoting the numbers of the tips and the nodes, containing distances between nodes
  allD=double() # this is the 'path' vector that will be outputted in the end
  nn=matrix(nrow=0, ncol=2) # initialize an empty matrix that will store nodeIds of each node and its ancestors corresponding to each position in the path vector
  nA=length(tree$tip.label)+tree$Nnode # nA is the total number of nodes -- tree$Nnode is the number of internal nodes, length(tree$tip.label) is the number of tip nodes
  matIndex=matrix(nrow=nA, ncol=nA)
  index=1
  # Below here is where the vector of 'paths' are generated. Starting from the first to last tip labels, and then the internal nodes from root to last.
  for ( i in 1:nA){ # for each node i
    ia=getAncestors(tree,i) # find the ancestor nodes of node i
    if(length(ia)>0){ # If node i has ancestors at all
      allD=c(allD, dd[i, ia]) # extends allD per node, where each extension is the distances between node i and its ancestors
      nn=rbind(nn,cbind(rep(i, length(ia)), ia)) # extend the xxx-by-2 matrix by nodeId pairs of node i and its ancestors
      for (j in ia){
        matIndex[i,j]=index
        index=index+1
      }
    }
  }
  return(list(dist=allD, nodeId=nn, matIndex=matIndex))
}

#' @keywords internal
getAncestors=function(tree, nodeN){
  if(is.character(nodeN)){
    nodeN=which(tree$tip.label==nodeN)
  }
  im=which(tree$edge[,2]==nodeN)
  if(length(im)==0){
    return()
  }
  else{
    anc=tree$edge[im,1]
    return(c(anc, getAncestors(tree, anc)))
  }
}

#' @keywords internal
findPairs=function(binary.tree){
  tip.labels = binary.tree$tip.label
  edge = binary.tree$edge
  edge.length = binary.tree$edge.length
  ind.fg.edge = which(edge.length == 1)
  fg.edges = edge[ind.fg.edge,]

  # Find the pairs
  fg.pairs.table = matrix(nrow=0,ncol=2)
  colnames(fg.pairs.table) = c("species1","species2")

  if (length(as.vector(fg.edges)) > 2){
    nodeId.start = sort(fg.edges[,1])
    fg_ca = vector("integer",length=0) # node IDs of the common ancestor foregrounds
    for (nn in 1:(length(nodeId.start)-1)){
      if (nodeId.start[nn] == nodeId.start[nn+1]){
        fg_ca = c(fg_ca,nodeId.start[nn])
        fg.pairs.table = rbind(fg.pairs.table, fg.edges[which(fg.edges[,1]==nodeId.start[nn]),2])
      }
    }
    rownames(fg.pairs.table) = fg_ca
  }
  fg.pairs.table
}

#'Calculate permulation correlation statistics
#' @param permulated.paths A nested list of permulated paths (e.g., output of \code{\link{calculatePermulatedPaths_apply}}
#' @param RERmat An RER matrix calculated using \code{\link{getAllResiduals}}.
#' @param min.sp Minimum number of species that must be present for a gene
#' @param min.pos Minimum number of species that must be present in the foreground (non-zero phenotype values)
#' @param method Method used to compute correlations. Accepts the same arguments as \code{\link{cor}}. Set to "auto" to select automatically based on the number of unique values in charP. This will also auto set the winsorization for Pearson correlation. Set winsorizetrait=some number and winsorizeRER=some number to override
#' @return A nested list containing the correlation statistics for the permulations
#' @export
calculateCorPermuted=function(permulated.paths,RERmat,min.sp=10,min.pos=2,method="k"){
  corMatList = lapply(permulated.paths,getAllCorSSM,RERmat,min.sp=min.sp,min.pos=min.pos,method=method)
  output.list <- list()
  output.list[[1]] <- corMatList
  return(output.list)
}

#' @keywords internal
getAllCorSSM=function(charP, RERmat, method="auto",min.sp=10, min.pos=2, winsorizeRER=NULL, winsorizetrait=NULL, weighted=F){
  if (method=="auto"){
    lu=length(unique(charP))
    if(lu==2){
      method="k"
      message("Setting method to Kendall")
    }
    else if (lu<=5){
      method="s"
      message("Setting method to Spearman")
    }
    else{
      method="p"
      message("Setting method to Pearson")
      if(is.null(winsorizeRER)){
        message("Setting winsorizeRER=3")
        winsorizeRER=3
      }
      if(is.null(winsorizetrait)){
        message("Setting winsorizetrait=3")
        winsorizetrait=3
      }
    }
  }
  win=function(x,w){
    xs=sort(x[!is.na(x)], decreasing = T)
    xmax=xs[w]
    xmin=xs[length(xs)-w+1]

    x[x>xmax]=xmax
    x[x<xmin]=xmin
    x
  }
  dim(RERmat) <- c(1,length(RERmat))
  corout=matrix(nrow=nrow(RERmat), ncol=3)
  rownames(corout)=rownames(RERmat)

  colnames(corout)=c("Rho", "N", "P")

  for( i in 1:nrow(corout)){

    if(((nb<-sum(ii<-(!is.na(charP)&!is.na(RERmat[i,]))))>=min.sp)){
      if (method!="p"&&sum(charP[ii]!=0)<min.pos){
        next
      }

      if(!weighted){

        x=RERmat[i,]

        #winsorize
        indstouse=which(!is.na(x) & !is.na(charP))
        if(!is.null(winsorizeRER)){
          x=win(x[indstouse], winsorizeRER)
        }else{
          x=x[indstouse]
        }
        if(!is.null(winsorizetrait)){
          y=win(charP[indstouse], winsorizetrait)
        }else{
          y=charP[indstouse]
        }

        cres=cor.test(x, y, method=method, exact=F)
        corout[i,1:3]=c(cres$estimate, nb, cres$p.value)
      }
      else{
        charPb=(charP[ii]>0)+1-1

        weights=charP[ii]
        weights[weights==0]=1

        cres=wtd.cor(RERmat[i,ii], charPb, weight = weights, mean1 = F)
        corout[i, 1:3]=c(cres[1], nb, cres[4])
      }
    }
    else{
      #show(i)
      #show(c(nb, charP[ii]))
    }

  }

  corout=as.data.frame(corout)
  corout$p.adj=p.adjust(corout$P, method="BH")
  corout
}
