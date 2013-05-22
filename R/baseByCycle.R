# The fraction of each base from beginning to end of the reads cycle
baseByCycle = function(fastqc)
{
  
  bbc= t(fastqc$baseByCycle/fastqc$Reads)
  baseByCycleFrame = as.data.frame(bbc)
  colnames(baseByCycleFrame)= c('A', 'C', 'G', 'T', 'N')
  baseByCycleFrame$Cycle = 1:nrow(bbc)
  baseByCycleFrame= melt(baseByCycleFrame, id.vars= 'Cycle', variable.name='Base', value.name='Fraction')  
  g <- ggplot(baseByCycleFrame, aes(x = Cycle, y = Fraction, color = Base))
  g + geom_line()+theme_bw()+ labs(title='Base Fraction per Cycle')
  
}