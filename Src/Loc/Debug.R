ERTimes = readRDS("Output/CategoricalPermulationsTimesER.rds")
SYMTimes = readRDS("Output/CategoricalPermulationsTimesSYM.rds")


newMainPhen = phenotypeVectorMain
oldMainPhen = readRDS("Data/CategoricalPermulationsTimingHillerPhenotypes.rds")

all.equal(newMainPhen, oldMainPhen)

oldMainPhen[!names(oldMainPhen) %in% names(newMainPhen)]


all.equal(ERTimes, SYMTimes)

ERTimesNumericMinutes = unlist(ERTimes[-c(1,2,3,4,7,8,9,10,11,22)])
ERTimesSeconds = unlist(ERTimes[c(1,2,3,4,7,8,9,10,11)])
ERTimesHours = unlist(ERTimes[c(22)])
ERTimeTotalHours = ((sum(ERTimesSeconds) + (sum(ERTimesNumericMinutes)*60) + (sum(ERTimesHours)*60*60))/(60*60))

SYMTimesMinutes = unlist(SYMTimes[-c(1,2,3,4,7,9,10,11,22,34)])
SYMTimesSeconds = unlist(SYMTimes[c(1,2,3,4,7,9,10,11)])
SYMTimesHours = unlist(SYMTimes[c(22,34)])
SYMTimeTotalHours = ((sum(SYMTimesSeconds) + (sum(SYMTimesMinutes)*60) + (sum(SYMTimesHours)*60*60))/(60*60))
timeDifferenceMinutes = ((SYMTimeTotal - ERTimeTotal)*60)


all.equal(ERTimes, SYMTimes)
ERTimeTotalHours
SYMTimeTotalHours
timeDifferenceMinutes


ERTimes

attributes(ERTimes[[1]])
ERTimesdf = unlist(ERTimes)


phenotypeVectorMain = readRDS("Data/CategoricalPermulationsTimingHillerPhenotypes.rds")
numberofEachPhen = table(phenotypeVectorMain)



function(treesObj, phenvals, rm, rp = "auto", ntrees, percent_relax = 0){
  
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





function(treesObj, phenvals, rm, rp = "auto", ntrees, percent_relax = 0){
  
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

function(tree, Q, P, nodes, tips, T0, Nk, cycles, alpha) {
  
  # get ancliks and max_states
  ancliks = getAncLiks(tree, tips, Q)
  max_states = getStatesAtNodes(ancliks)
  
  # calculate tree likelihoods
  # states = c(tips, max_states)
  # max_lik = 1 
  # for(i in 1:nrow(tree$edge)){
  #   a = states[tree$edge[i,1]]
  #   d = states[tree$edge[i,2]]
  #   max_lik = max_lik * P[[i]][a, d]
  # }
  
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