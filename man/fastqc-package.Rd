\name{fastqc-package}
\alias{fastqc-package}
\alias{fastqc}
\docType{package}
\title{
Quality control of short read fastq files.
}
\description{
This package reimplements the fastqc quality control steps in R. The original fastqc is a standalone java program. The fastqc package uses Rcpp and RcppArmadillo to interface with the seqan C++ library. It's pretty fast and uses efficient DNA containers and unsigned integer matrices to minimise the memory footprint. On my laptop (4GB RAM) it will process 2 million reads by 50 cycles in ~ 30secs. The original fastqc program restricted itself to sampling the first ~200k reads - which is linearly ~ 3 secs for this implementation. /cr

Future updates will be handling multiple files, and multithreading different tasks.

}
\details{
\tabular{ll}{
Package: \tab fastqc\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-05-21\cr
License: \tab GPL?\cr
}
The workhorse function is \code{fastqc} which calls the C++ code to calculate all the necessary numerical data into a list of named results. The other functions are simple printing and plotting of these results using R routines.

}
\author{
Stephen Henderson\cr
UCL Cancer Institute

Maintainer: Stephen Henderson <s.henderson@ucl.ac.uk>

}
\references{
The original implementation of fastqc is below. This R version gives near identical results. \cr
http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ package }

\examples{
fastqcResults = fastqc("test7.fastq")
}
