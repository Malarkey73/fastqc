# Here we are printing a table of duplicated sequences
seqDupl= function(fastqc)
{
  dups.df=data.frame(Bins=1:10, Counts=fastqc$DupBins/fastqc$DupBins[1])
  dupFraction = sum(fastqc$DupBins[2:10])/sum(fastqc$DupBins)
  
  seqs.df= data.frame(Duplicates=fastqc$DupSeqs, Count=fastqc$DupSeqsCount, 
                      Fraction= fastqc$DupSeqsCount/fastqc$Reads)
  
  print(seqs.df)
  qplot(x=Bins, y=Counts, data=dups.df, geom='line', main = paste('Duplicate Fraction =', 
                                                                  format(dupFraction, digits=3)))+theme_bw()
  
  
}