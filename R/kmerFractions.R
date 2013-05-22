#This function extracts and prints over-repreented 5-mers e.g typically AAAAA or TTTTT
kmerFractions= function(fastqc, frac=3)
{
  
  # create a vector of all possible 5mers (excluding N) 
  acgt=c('A','C', 'G', 'T')
  acgtCombn = expand.grid(acgt,acgt,acgt,acgt,acgt)
  acgtCombn = apply(acgtCombn,1, function(x)paste0(x,collapse=''))
  
  # the c++ workhorse function
  kmerData=fastqc$kmers
  nucf= fastqc$nucFreq
  
  nucfCombn = expand.grid(nucf,nucf,nucf,nucf,nucf)
  expec = apply(nucfCombn,1, prod)*(fastqc$Cycles-4)
  
  # the c++ function returns kmer data in lexical order so we reorder
  o = order(acgtCombn)
  df=data.frame(kmers=acgtCombn[o], fraction=kmerData, expected=expec[o])
  
  return(df[which((df$fraction/df$expected)>frac),])
  
}