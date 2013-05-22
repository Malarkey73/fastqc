\name{seqDupl}
\alias{seqDupl}
\title{Detecting duplicated artefacts}
\description{
  When you sequence DNA (or cDNA) you expect to get some depth - or overlapping reads -at least some of which will be identical. However when you see too many identical reads this is likely due to experimental artefacts.
}
\details{
  This function plots both the distribution of occurence of sequences from 1 (singletons) to 10 copies and finally >10. It also prints a table showing highly repeated sequences. Often the top repeated sequence is a poly-A artefact.


  \preformatted{%
      
  }
}
\references{ http://en.wikipedia.org/wiki/Phred_quality_score}
\author{Stephen Henderson

s.henderson@ucl.ac.uk}
\examples{
  \dontrun{
  test2=fastqc("test2.fastq")
  seqDupl(test2)
  }
}
\keyword{fastq}
\keyword{genomics}
