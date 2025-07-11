convertBinaryToCategoricalCorrelationsFormat = function(correlation){
  noNCorrelation = correlation[,-2]
  fakePairwiseCorrelation = list(noNCorrelation)
  names(fakePairwiseCorrelation) = "1 - 2"
  fakeCombinedCorrelation = list(correlation, fakePairwiseCorrelation)
  fakeCombinedCorrelation
}
