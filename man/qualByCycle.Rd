\name{qualByCycle}
\alias{qualByCycle}
\title{Phred quality across the reads}
\description{
  There is a tendency for short read DNA sequencing to deteriorate towards the end of the read. The most crucial QC step is thus to find the poin tto trim the ends.
}
\details{
  Phred Quality scores are a character based quality score: \eqn{Q=-log10P}
  
  Where P is the estimated error rate (or probability). The \code{qualByCycle} function takes \code{fastqc} function results as input and plots this quality score from first to last base of the read (a.k.a. the cycles).

  
  \preformatted{%
      
  }
}
\references{ http://en.wikipedia.org/wiki/Phred_quality_score}
\author{Stephen Henderson

s.henderson@ucl.ac.uk}
\examples{
  \dontrun{
  test2=fastqc("test2.fastq")
  qualByCycle(test2)
  }
}
\keyword{fastq}
\keyword{genomics}
