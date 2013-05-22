\name{qualBySeq}
\alias{qualBySeq}
\title{Average Phred Quality of reads}
\description{
  This function plots an overview of the average read quality (across each entire read).
}
\details{
  Phred Quality scores are a character based quality score: \eqn{Q=-log10P}
  
  Where P is the estimated error rate (or probability). The \code{qualBySeq} function takes \code{fastqc} function results as input and plots the distribution of average quality scores of the reads.


  \preformatted{%
      
  }
}
\references{ http://en.wikipedia.org/wiki/Phred_quality_score}
\author{Stephen Henderson

s.henderson@ucl.ac.uk}
\examples{
  \dontrun{
  test2=fastqc("test2.fastq")
  qualBySeq(test2)
  }
}
\keyword{fastq}
\keyword{genomics}
