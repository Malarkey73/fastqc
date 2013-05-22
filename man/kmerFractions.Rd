\name{kmerFractions}
\alias{kmerFractions}
\title{Over Represented kmers(5)}
\description{
  Sometimes short read data can contain biased or artefactual repetitive DNA sequence - sometimes also called low complexity. You can spot such repetitive atrefacts by screening for over-represented kmers (k=5).
}
\details{
  The expected fraction of a kmer of 5 bases can be estimated from the overall base (A/T/G/C) composition. Ignoring N there are 1024 possible kmers. The fastqc function calculates the actual and expected base fraction for each. Then for this function there is a second argument \code{frac}. If the raito of actual kmers to expected kmers > \code{frac} then it is returned in the results table. The results table is sorted lexicographically i.e. AAAAA to TTTTT.

  
  \preformatted{%
      
  }
}
\author{Stephen Henderson

s.henderson@ucl.ac.uk}
\examples{
  \dontrun{
  test2=fastqc("test2.fastq")
  kmerfractions(test2)
  }
}
\keyword{fastq}
\keyword{genomics}
