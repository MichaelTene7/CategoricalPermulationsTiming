library(RERconverge)
#My custom permulation function
source("simBinPhenoSSM_fromMasterTree.r") 
#Additional functions needed for running permulations - loading RERconverge doesn't work for some reason when I run permulations from a source file
source("neededPermulationFunctions.r")

#Number of permulations to run for testing
numPerms<-1000
#Which gene tree to use for testing
#Tree 3 has 72 out of 80 species present, including the two outgroups (Lophuromys_woosnami_LSUMZ37793, Lophiomys_imhausi_UM5152)
#Tree 4 has 65 out of 80 species present, NOT the two outgroups
useTree<-4

#Load RERconverge output from murine rodent analysis
load("RERconverge_output.logRTM_binary.OUmodel.RTMspeciesOnly.rds")

#Setup foreground species vector
fgspec<-c("Pseudomys_novaehollandiae_ABTC08140", "Pseudomys_delicatulus_U1509", "Notomys_alexis_U1308", "Notomys_fuscus_M22830", "Notomys_mitchellii_M21518", "Zyzomys_pedunculatus_Z34925", "Bandicota_indica_ABTC119185", "Nesokia_indica_ABTC117074", "Hyomys_goliath_ABTC42697", "Pseudomys_shortridgei_Z25113", "Paruromys_dominator_JAE4870", "Eropeplus_canus_NMVZ21733")

#Run one permulation with custom script
permTree<-simBinPhenoSSM_fromMasterTree(tree=myTrees$trees[[useTree]], trees=myTrees, fg_vec=fgspec, sisters_list=NULL, pathvec=pathvec)

#Get simulated foregrounds from permulation
fgEdges<-permTree$edge[which(permTree$edge.length==1),2]
permFgs<-permTree$tip.label[fgEdges]

#Run n permulations; get simulated foreground taxa for each permulation
#Also keep track of how long it takes to run n permulations
print(paste("Running", numPerms, "test permulations..."))
#List of permulated trees
permTrees<-list()
#List of vectors of permulated foreground species
permFG_list<-list()
#Data frame to keep track of whether or not each species is in the foreground for each permulation
inFG<-c()
Sys.time()
for(i in 1:numPerms){
	permTree<-simBinPhenoSSM_fromMasterTree(tree=myTrees$trees[[useTree]], trees=myTrees, fg_vec=fgspec, pathvec=pathvec)
	permTrees[[i]]<-permTree
	fgEdges<-permTree$edge[which(permTree$edge.length==1),2]
	permFgs<-permTree$tip.label[fgEdges]
	permFG_list[[i]]<-permFgs
	inFG<-cbind(inFG, unlist(sapply(myTrees$masterTree$tip.label, function(x) if(x %in% permFgs){1} else{0})))
}
Sys.time()
inFG_sums<-rowSums(inFG)
inFG_props<-inFG_sums/numPerms
shortnames<-unlist(sapply(names(inFG_sums), function(x) paste(strsplit(x, "_")[[1]][1:2], collapse="_")))
names(inFG_sums)<-shortnames
names(inFG_props)<-shortnames
pdf("permulations_test_barplot.pdf", onefile=TRUE, height=8.5, width=11)
barplot(height=inFG_sums, main=paste("Number of times each species appeared in the foreground out of", numPerms, "permulations"), las=2, cex.names=0.5)
barplot(height=inFG_props, main=paste("Proportion of times each species appeared in the fourground out of", numPerms, "permulations"), las=2, cex.names=0.5)
dev.off()

inFG_noMidpointRoot<-c()
Sys.time()
for(i in 1:numPerms){
	permTree<-simBinPhenoSSM_fromMasterTree_noMidpointRoot(tree=myTrees$trees[[useTree]], trees=myTrees, fg_vec=fgspec, pathvec=pathvec)
	permTrees[[i]]<-permTree
	fgEdges<-permTree$edge[which(permTree$edge.length==1),2]
	permFgs<-permTree$tip.label[fgEdges]
	inFG_noMidpointRoot<-cbind(inFG, unlist(sapply(myTrees$masterTree$tip.label, function(x) if(x %in% permFgs){1} else{0})))
}
Sys.time()
inFG_noMidpointRoot_sums<-rowSums(inFG_noMidpointRoot)
inFG_noMidpointRoot_props<-inFG_noMidpointRoot_sums/numPerms
names(inFG_noMidpointRoot_sums)<-shortnames
names(inFG_noMidpointRoot_props)<-shortnames
pdf("permulations_test_barplot.noMidpointRoot.pdf", onefile=TRUE, height=8.5, width=11)
barplot(height=inFG_noMidpointRoot_sums, main=paste("Number of times each species appeared in the foreground out of", numPerms, "permulations; no midpoint root"), las=2, cex.names=0.5)
barplot(height=inFG_noMidpointRoot_props, main=paste("Proportion of times each species appeared in the fourground out of", numPerms, "permulations; no midpoint root"), las=2, cex.names=0.5)
dev.off()

inFG_notMasterTree_noMidpointRoot<-c()
Sys.time()
for(i in 1:numPerms){
	permTree<-simBinPhenoSSM_noMidpointRoot(tree=myTrees$trees[[useTree]], trees=myTrees, fg_vec=fgspec, pathvec=pathvec)
	permTrees[[i]]<-permTree
	fgEdges<-permTree$edge[which(permTree$edge.length==1),2]
	permFgs<-permTree$tip.label[fgEdges]
	inFG_notMasterTree_noMidpointRoot<-cbind(inFG, unlist(sapply(myTrees$masterTree$tip.label, function(x) if(x %in% permFgs){1} else{0})))
}
Sys.time()
inFG_notMasterTree_noMidpointRoot_sums<-rowSums(inFG_notMasterTree_noMidpointRoot)
inFG_notMasterTree_noMidpointRoot_props<-inFG_notMasterTree_noMidpointRoot_sums/numPerms
names(inFG_notMasterTree_noMidpointRoot_sums)<-shortnames
names(inFG_notMasterTree_noMidpointRoot_props)<-shortnames
pdf("permulations_test_barplot.notMasterTree_noMidpointRoot.pdf", onefile=TRUE, height=8.5, width=11)
barplot(height=inFG_notMasterTree_noMidpointRoot_sums, main=paste("Number of times each species appeared in the foreground out of", numPerms, "permulations\nnot using master tree and no midpoint root"), las=2, cex.names=0.5)
barplot(height=inFG_notMasterTree_noMidpointRoot_props, main=paste("Proportion of times each species appeared in the fourground out of", numPerms, "permulations\nnot using master tree and no midpoint root"), las=2, cex.names=0.5)
dev.off()
