# Here we are printing a table of duplicated sequences
seqDupl= function(fastqc)
{
  dups.df=data.frame(Bins=1:10, Counts=fastqc$DupBins/fastqc$DupBins[1])
  dupFraction = sum(fastqc$DupBins[2:10])/sum(fastqc$DupBins)
  
  seqs.df= data.frame(Duplicates=fastqc$DupSeqs, Fraction=fastqc$DupSeqsCount, 
                      Fraction= fastqc$DupSeqsCount/fastqc$Reads)
  
  print(seqs.df)
  qplot(x=Bins, y=Counts, data=dups.df[2:10,], geom=c('path','bar'), stat="identity", 
    main = paste('Duplicate Fraction =', format(dupFraction, digits=3)))+theme_bw()+ 
    scale_x_continuous(breaks=2:10, labels=c(2:9, "10+"))
  
  
}