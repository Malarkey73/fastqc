# This is the C++ workhorse using Rcpp, Armadillo and especially seqan
# The resulting object contains all the necessary data for other functions.
fastqc= function(file)
{
  fastqCpp(file)
  #.Call("fastqCpp", file, package="fastqc")  
}















