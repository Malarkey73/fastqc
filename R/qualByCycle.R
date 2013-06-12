# The qualByCycle boxplots the phred quality score from beginning to end of the reads cycle
qualByCycle = function(fastqc)
{
  boxStatsFrame = fastqc$Fivenum
  boxtitles = paste('no of reads = ', fastqc$Reads)
  boxStatsFrame= t(boxStatsFrame)
  colnames(boxStatsFrame)=c('min', 'Q1', 'median', 'Q3', 'max')
  boxStatsFrame=as.data.frame(boxStatsFrame)
  boxStatsFrame$Cycle=factor(1:fastqc$Cycles)
  
  
  # make a plot
  g <- ggplot(boxStatsFrame, aes(x = Cycle, ymin = min, lower = Q1, middle = median, upper = Q3, ymax = max))
  g + geom_boxplot(stat = "identity")+theme_bw()+labs(list(title=boxtitles, y= 'Phred Score'))  
  
}