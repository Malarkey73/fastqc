\name{baseByCycle}
\alias{baseByCycle}
\title{Plotting bases across the read cycle.}
\description{
  Particularly at the 5' end of reads there is often a GC bias. This function plots the fraction of bases (A/T/G/C/N) from    the beginning to the end of the read.

  This example shows (in the corresponding C++ code) how to access,
  modify and create a data frame.
}
\details{
  The \code{baseByCycle} function accepts the return result from the \code{fastqc} function and produces a plot using ggplot2.

  
  \preformatted{%
      
  }
}
\author{Stephen Henderson

s.henderson@ucl.ac.uk}
\examples{
  \dontrun{
  test2=fastqc("test2.fastq")
  baseByCycle(test2)
  }
}
\keyword{fastq}
\keyword{genomics}
