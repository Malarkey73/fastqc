# The histogram of mean Phred Scores of reads
qualBySeq = function(fastqc)
{
  
  # tidy into a data.frame of bins and read counts
  qualSeqFrame = data.frame(MeanSeqQual= (fastqc$minQual+1):fastqc$maxQual, Reads = fastqc$qualBySeq)
  
  g <- ggplot(qualSeqFrame, aes(x = MeanSeqQual, y = Reads))
  g + geom_line()+theme_bw()+labs(list(x='Mean Sequence Quality (Phred Scores)', y= 'Reads'))
  
}
