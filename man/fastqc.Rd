\name{fastqc}
\alias{fastqc}
\title{fastq QC workhorse}
\description{
  This function calls the C++ workhorse code that does all the numerical calculations for the fastq QC. All the other functions in the \code{fastqc} package take the \code{fastqc} result as input.
}
\details{
  The \code{fastqc} function depends upon Rcpp, RcppArmadillo and especially the seqan C++ library for DNA sequences.

  
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
