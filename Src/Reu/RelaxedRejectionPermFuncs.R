# Functions for Permulations for Categorical Traits

# Include fixes in improveTree that NEED TO GET UPDATED ON GITHUB!!!

# internal function that generates a set of N null tips with the number of species in each category matching the actual phenotype data
getNullTips <- function(tree, Q, N, intlabels, root_prob = "stationary", percent_relax) {
  
  # GET TRUE TIP COUNTS
  true_counts = table(intlabels$mapped_states)
  
  # MAKE MATRIX TO STORE THE SETS OF NULL TIPS AND SETS OF INTERNAL NODES
  tips = matrix(nrow = N, ncol = length(tree$tip.label), dimnames = list(NULL, tree$tip.label))
  nodes = matrix(nrow = N, ncol = tree$Nnode)
  
  cnt = 0
  while(cnt < N) {
    # SIMULATE STATES
    sim = simulate_mk_model(tree, Q, root_probabilities = root_prob)
    sim_counts = table(sim$tip_states)
    
    # CHECK THAT ALL STATES GET SIMULATED IN THE TIPS
    if(length(unique(sim$tip_states)) < length(true_counts)) { 
      next
    }
    
    # IF THE TIP COUNTS MATCH THE TIP COUNTS IN THE REAL DATA, ADD TO THE LIST
    # sum(true_counts == sim_counts) == length(true_counts)
    if(sum(abs(sim_counts - true_counts) <= true_counts*percent_relax) == length(true_counts)) {
      cnt = cnt + 1
      
      print(cnt)
      
      tips[cnt,] = sim$tip_states
      nodes[cnt,] = sim$node_states
    }
  }
  return(list(tips = tips, nodes = nodes))
}

# internal function that shuffles the categories around the tree based on ancestral likelihoods
# serves as a starting point for the function improveTree
shuffleInternalNodes <- function(shuffled_states, available_nodes, ancliks, Nnode, ntips) {
  internal_states = vector(mode = "numeric", length = Nnode)
  
  if(length(shuffled_states) != nrow(ancliks)){
    stop("number of shuffled states and number of nodes with ancestral likelihoods do not match")
  }
  
  for(state in shuffled_states){
    if(length(available_nodes) > 1) {
      liks = ancliks[,state]
      node = sample(available_nodes, size = 1, prob = liks)
      available_nodes = available_nodes[- which(available_nodes == node)]
      ancliks = ancliks[- which(rownames(ancliks) == as.character(node)),]
      internal_states[node - ntips] = state
    }
    else { # only one node left
      internal_states[available_nodes - ntips] = state
    }
  }
  return(internal_states)
}

# internal function that for a set of null tips, shuffles the correct number of each category around the internal nodes
getNullTrees <- function(node_states, null_tips, tree, Q) {
  
  nullTrees = list()
  
  for(i in 1:nrow(null_tips)) {
    print(i)
    tips = null_tips[i,]
    
    ancliks = getAncLiks(tree, tipvals = tips, Q = Q)
    
    ntips = length(tree$tip.label)
    available_nodes = (ntips + 1):(tree$Nnode + ntips)
    rownames(ancliks) = available_nodes
    
    shuffled_states = sample(node_states)
    internal_states = shuffleInternalNodes(shuffled_states, 
                                           available_nodes = available_nodes,
                                           ancliks = ancliks, 
                                           Nnode = tree$Nnode, ntips = ntips)
    tr = list(tips = tips, nodes = internal_states)
    nullTrees = append(nullTrees, list(tr))
  }
  return(nullTrees)
}

# internal function that rearranges the shuffled internal nodes to improve the likelihoods of the permulated trees
improveTree <- function(tree, Q, P, nodes, tips, T0, Nk, cycles, alpha) {
  
  # get ancliks and max_states
  ancliks = getAncLiks(tree, tips, Q)
  
  states = c(tips, nodes)
  curr_lik = 1
  for(i in 1:nrow(tree$edge)){
    a = states[tree$edge[i,1]]
    d = states[tree$edge[i,2]]
    curr_lik = curr_lik * P[[i]][a, d]
  }
  
  # calculate initial ratios
  nstates = nrow(Q)
  ratios = c() # list of ratios
  ratio_info = matrix(nrow = (nstates - 1) * tree$Nnode, ncol = 3, dimnames = list(NULL, c("node", "state", "other.state"))) # info for each ratio
  
  # ns aren't the node numbers in the tree - they are the index of the internal node in nodes, node number in tree is n + ntips
  for(n in 1:tree$Nnode) {
    # calculate ratios
    pie = ancliks[n,]
    rr = pie[-nodes[n]] / pie[nodes[n]] # other states / state
    ratios = c(ratios, rr)
    # fill in ratio_info 
    # rows = c((n-1)*3 + 1, (n-1)*3 + 2, (n-1)*3 + 3)
    rows = ((n-1)*(nstates-1) + 1):((n-1)*(nstates-1) + (nstates-1))
    ratio_info[rows,"node"] = rep(n, nstates - 1)
    ratio_info[rows,"state"] = rep(nodes[n], nstates - 1)
    ratio_info[rows,"other.state"] = (1:nstates)[-nodes[n]]
  }
  
  # pre-calculate and store edge numbers for each node
  ntips = length(tree$tip.label)
  edg_nums = lapply(seq_along(vector(mode = "list", length = tree$Nnode + ntips)), function(x){
    c(which(tree$edge[,1] == x),(which(tree$edge[,2] == x)))
  })
  
  j = 1 # iteration counter
  k = 1 # cycle counter 
  Tk = T0 
  
  while(k <= cycles) { 
    
    # get 2 nodes to swap
    nn = nodes
    
    # 1: pick a node randomly, weighted by the ratios
    r1 = sample(1:length(ratios), 1, prob = ratios)
    n1 = ratio_info[r1, "node"] # node 1
    s1 = ratio_info[r1, "state"] # state1
    s2 = ratio_info[r1, "other.state"] # state2
    
    # 2: pick a node to swap it with 
    ii = intersect(which(ratio_info[,"state"] == s2), which(ratio_info[,"other.state"] == s1))
    if(length(ii) > 1) {
      n2 = sample(ratio_info[ii,"node"], 1, prob = ratios[ii]) # node2
    } else { # only one node with state2
      n2 = ratio_info[ii,"node"]
    }
    
    # make the swap
    nn[n1] = s2
    nn[n2] = s1
    
    # calculate new likelihood
    states_new = c(tips, nn)
    states_old = c(tips, nodes)
    
    r = 1
    for(i in unique(c(edg_nums[[n1 + ntips]], edg_nums[[n2 + ntips]]))){ # check this over many cases including when n1 and n2 effect the same edge
      ao = states_old[tree$edge[i,1]]
      do = states_old[tree$edge[i,2]]
      
      an = states_new[tree$edge[i,1]]
      dn = states_new[tree$edge[i,2]]
      r = r * (P[[i]][an,dn] / P[[i]][ao, do])
    }
    
    if(r >= 1) { # if the swap increases likelihood, commit to the swap
      
      nodes = nn
      
      curr_lik = curr_lik * r # this should do the same thing, BUT CHECK THIS GETS THE SAME RESULT IN MULTIPLE CASES!
      
      # update ratios
      # rows1 = c((n1-1)*3 + 1, (n1-1)*3 + 2, (n1-1)*3 + 3) # rows to update ratios for n1
      rows1 = ((n1-1)*(nstates-1) + 1):((n1-1)*(nstates-1) + (nstates-1))
      
      pie = ancliks[n1,]
      rr = pie[-s2] / pie[s2] # other states / state
      ratios[rows1] = rr
      
      # fill in ratio_info 
      ratio_info[rows1,"state"] = rep(s2, nstates - 1) 
      ratio_info[rows1,"other.state"] = (1:nstates)[-s2]
      
      # rows2 = c((n2-1)*3 + 1, (n2-1)*3 + 2, (n2-1)*3 + 3) # rows to update ratios for n2
      rows2 = ((n2-1)*(nstates-1) + 1):((n2-1)*(nstates-1) + (nstates-1))
      
      pie = ancliks[n2,]
      rr = pie[-s1] / pie[s1] # other states / state
      ratios[rows2] = rr
      
      # fill in ratio_info 
      ratio_info[rows2,"state"] = rep(s1, nstates - 1) 
      ratio_info[rows2,"other.state"] = (1:nstates)[-s1]
      
    } 
    else { # make jump with probability u
      
      # calculate u which includes dividing by tmp
      dh = -log(curr_lik * r) + log(curr_lik)
      u = exp(-dh/Tk)
      if(u == 0) warning("u is zero")
      
      # if(u == 0) stop(paste("temp is", Tk))
      
      if(runif(1) <= u) {
        
        nodes = nn
        
        curr_lik = curr_lik * r # CHECK THIS GETS THE SAME RESULT
        
        # update ratios
        # rows1 = c((n1-1)*3 + 1, (n1-1)*3 + 2, (n1-1)*3 + 3) # rows to update ratios for n1
        rows1 = ((n1-1)*(nstates-1) + 1):((n1-1)*(nstates-1) + (nstates-1))
        
        pie = ancliks[n1,]
        rr = pie[-s2] / pie[s2] # other states / state
        ratios[rows1] = rr
        # fill in ratio_info
        
        ratio_info[rows1,"state"] = rep(s2, nstates - 1)
        ratio_info[rows1,"other.state"] = (1:nstates)[-s2]
        
        # rows2 = c((n2-1)*3 + 1, (n2-1)*3 + 2, (n2-1)*3 + 3) # rows to update ratios for n2
        rows2 = ((n2-1)*(nstates-1) + 1):((n2-1)*(nstates-1) + (nstates-1))
        
        pie = ancliks[n2,]
        rr = pie[-s1] / pie[s1] # other states / state
        ratios[rows2] = rr
        # fill in ratio_info
        ratio_info[rows2,"state"] = rep(s1, nstates - 1)
        ratio_info[rows2,"other.state"] = (1:nstates)[-s1]
      }
    }
    
    # increment j
    j = j + 1
    
    # print(curr_lik)
    
    # move to next cycle if necessary
    if(j >= Nk) {
      j = 1 # reset j
      # Tk = T0 * alpha^k
      # Tk = T0 / (1 + alpha * log(k))
      Tk = T0 / (1 + alpha*k)
      k = k + 1
    }
    
  }
  end = Sys.time()
  return(list(nodes = nodes, lik = log10(curr_lik)))
}

# returns null trees!!! 
# treesObj - trees object returned by readTrees
# phenvals - the named phenotype vector
# rm - rate model (should be the same as the one used to reconstruct ancestral history of the tree)
# rp - root prior (default is auto)
# ntrees - number of null trees to generate
# percent_relax is the percentage of the category size (between 0 and 1) that the category counts can be off by
# e.g. the counts can be true_counts +/- percent_relax * true_counts
  # percent_relax can be a single value (use same percent for each category)
  # percent_relax can be a vector of percentages (between 0 and 1) for each category
    # should be in the same order as the numerical values assigned to each category (printed out by char2TreeCategorical)
categoricalPermulations <- function(treesObj, phenvals, rm, rp = "auto", ntrees, percent_relax = 0){
  
  # check percent_relax is one value or a vector of length = # traits
  if(!(length(percent_relax) == 1 || length(percent_relax) == length(unique(phenvals)))) {
    stop("percent_relax is the wrong length")
  }
  
  # PRUNE TREE, ORDER PHENVALS, MAP TO STATE SPACE
  tree = treesObj$masterTree
  keep = intersect(names(phenvals), tree$tip.label)
  tree = pruneTree(tree, keep)
  phenvals = phenvals[tree$tip.label]
  intlabels = map_to_state_space(phenvals)
  
  # FIT A TRANSITION MATRIX ON THE DATA
  message("Fitting transition matrix")
  Q = fit_mk(tree, intlabels$Nstates, intlabels$mapped_states,
             rate_model = rm, root_prior = rp)$transition_matrix
  
  # GET NULL TIPS (AND STORE INTERNAL NODES FROM SIMULATIONS TOO)
  message("Simulating trees")
  simulations = getNullTips(tree, Q, ntrees, intlabels, 
                            percent_relax = percent_relax)
  
  ancliks = getAncLiks(tree, intlabels$mapped_states, Q = Q)
  node_states = getStatesAtNodes(ancliks)
  
  # GET SHUFFLED STARTING-POINT TREES
  message("Shuffling internal states")
  nullTrees = getNullTrees(node_states, simulations$tips, tree, Q)
  
  P = lapply(tree$edge.length, function(x){expm(Q * x)})
  
  # IMPROVE LIKELIHOOD OF EACH NULL TREE
  message("Improving tree likelihoods")
  improvedNullTrees = lapply(nullTrees, function(x){
    list(tips = x$tips, nodes = improveTree(tree, Q, P, x$nodes, x$tips, 10, 10, 100, 0.9)$nodes)
  })
  
  # RETURN
  message("Done")
  return(list(sims = simulations, trees = improvedNullTrees, startingTrees = nullTrees))
}

# calculates empirical p-values from the null trees
# realCors - output of correlateWithCategoricalPhenotype
# nullPhens - categoricalPermulations(...)$trees
# phenvals - the named phenotype vector
# treesObj - trees object returned by readTrees
# RERmat - from getAllResiduals
# method - "kw" or "aov"
getPermPvalsCategorical <- function(realCors, nullPhens, phenvals, treesObj, RERmat, method = "kw") {
  # CHECK IF TRAIT IS BINARY
  binary = FALSE
  if(method != "kw" & method != "aov"){
    binary = TRUE
  }
  
  # PRUNE TREE
  tree = treesObj$masterTree
  keep = intersect(names(phenvals), tree$tip.label)
  tree = pruneTree(tree, keep)
  
  # UNROOT THE TREE IF IT IS ROOTED
  if (is.rooted(tree)) {
    tree = unroot(tree)
  }
  
  # generate the paths 
  if(!binary) {
    message("Generating null paths")
    nullPaths = lapply(nullPhens, function(x){
      tr = tree # make a copy of the tree
      tr$edge.length = c(x$tips, x$nodes)[tr$edge[,2]] # assign states to edges
      tree2Paths(tr, treesObj, categorical = TRUE, useSpecies = names(phenvals)) # calculate paths
    })
  } else {
    message("Generating null paths")
    nullPaths = lapply(nullPhens, function(x){
      tr = tree # make a copy of the tree
      tr$edge.length = c(x$tips, x$nodes)[tr$edge[,2]] # assign states to edges
      # subtract 1 to convert 1 - FALSE, 2 - TRUE to 0 and 1
      # THIS IS A QUICK FIX, WON'T WORK IF THE FOREGROUND IS 1 AND BACKGROUND IS 2!!!
      tree2Paths(tr, treesObj, categorical = TRUE, useSpecies = names(phenvals)) - 1 # calculate paths
    })
  }
  
  # calculate correlation statistics
  message("Calculating correlation statistics")
  # make matrices to store the results
  corsMatPvals = matrix(nrow = nrow(RERmat), ncol = length(nullPhens),
                        dimnames = list(rownames(RERmat), NULL))
  corsMatEffSize = matrix(nrow = nrow(RERmat), ncol = length(nullPhens),
                          dimnames = list(rownames(RERmat), NULL))
  
  if(!binary){
    # make matrices for the pairwise testing
    Ppvals = lapply(1:length(realCors[[2]]), matrix, data = NA, nrow = nrow(RERmat),
                    ncol = length(nullPhens), dimnames = list(rownames(RERmat), NULL))
    names(Ppvals) = names(realCors[[2]])
    
    Peffsize = lapply(1:length(realCors[[2]]), matrix, data = NA, nrow = nrow(RERmat),
                      ncol = length(nullPhens), dimnames = list(rownames(RERmat), NULL))
    names(Peffsize) = names(realCors[[2]])
  }
  
  # run getAllCor on every null phenotype and store the pval/effect size for every gene
  if(!binary) {
    for(i in 1:length(nullPaths)) {
      cors = getAllCor(RERmat, nullPaths[[i]], method = method)
      corsMatPvals[,i] = cors[[1]]$P # store p values
      corsMatEffSize[,i] = cors[[1]]$Rho # store effect size
      
      # add results of pairwise tests
      for(j in 1:length(cors[[2]])){ # loop through each table in cors[[2]]
        Ppvals[[names(cors[[2]])[j]]][,i] = cors[[2]][[j]]$P
        Peffsize[[names(cors[[2]])[j]]][,i] = cors[[2]][[j]]$Rho
      }
    }
  } else {
    for(i in 1:length(nullPaths)) {
      cors = getAllCor(RERmat, nullPaths[[i]], method = method)
      corsMatPvals[,i] = cors$P # store p values
      corsMatEffSize[,i] = cors$Rho # store effect size
    }
  }
  
  message("Obtaining permulations p-values")
  if(!binary){
    # calculate empirical pvals
    N = nrow(realCors[[1]])
    
    realCors[[1]]$permP = rep(NA, N) # add a column to the real results for the empirical pvals
    
    for(j in 1:length(realCors[[2]])){ # loop through each table in realCors[[2]]
      realCors[[2]][[j]]$permP = rep(NA, N)
    }
    
    for(gene in 1:N) {
      # check whether the gene is NA in realCors and if so, set p = NA
      if(is.na(realCors[[1]]$Rho[gene])) {
        p = NA
      } else {
        # count number of times null is more extreme than observed
        # effect size for Kruskal-Wallis is epsilon squared, only assumes non-negative values
        # ANOVA is still using eta2 (as of right now)
        p = sum(corsMatEffSize[gene,] > realCors[[1]]$Rho[gene], na.rm = TRUE) / sum(!is.na(corsMatEffSize[gene,]))
      }
      realCors[[1]]$permP[gene] = p
      
      for(j in 1:length(realCors[[2]])) {
        if(is.na(realCors[[2]][[j]]$Rho[gene])) {
          p = NA
        } else {
          # I think the effect size for Tukey and Dunn are bidirectional - using absolute value, NEED TO CHECK WITH AMANDA
          p = sum(abs(Peffsize[[names(realCors[[2]][j])]][gene,]) > abs(realCors[[2]][[j]]$Rho[gene]), na.rm = TRUE) / sum(!is.na(Peffsize[[names(realCors[[2]][j])]][gene,]))
        }
        realCors[[2]][[j]]$permP[gene] = p
      }
    }
  } else {
    # calculate empirical pvals
    N = nrow(realCors)
    
    realCors$permP = rep(NA, N) # add a column to the real results for the empircal pvals
    
    for(gene in 1:N) {
      if(is.na(realCors$Rho[gene])){
        p = NA
      } else {
        # count number of times null is more extreme than observed
        # use abs because stat can be positive or negative
        p = sum(abs(corsMatEffSize[gene,]) > abs(realCors$Rho[gene]), na.rm = TRUE) / sum(!is.na(corsMatEffSize[gene,]))
      }
      realCors$permP[gene] = p
    }
  }
  
  message("Done")
  # return results
  if(!binary){
    return(list(res = realCors, pvals = list(corsMatPvals,Ppvals), effsize = list(corsMatEffSize,Peffsize)))
  }
  else {
    return(list(res = realCors, pvals = corsMatPvals, effsize = corsMatEffSize))
  }
}

plotTreeWStates <- function(tree, tips, internal_states) {
  plot(tree, cex = 0.2)
  tiplabels(pie = to.matrix(tips, sort(unique(tips))),cex = 0.2)
  nodelabels(pie = to.matrix(internal_states, sort(unique(internal_states))), cex = 0.2)
}
