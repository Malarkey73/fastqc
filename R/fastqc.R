# This is the C++ workhorse using Rcpp, Armadillo and especially seqan
# The resulting object contains all the necessary data for other functions.
# If mc=TRUE the the function will attempt parallel processing
# If sampled = FALSE all the fastq file wil be read, if sampled= TRUE the first 200K
fastqc= function(files, mc=TRUE, numreads= 200000)
{
  # for a single file
  if (length(files)==1)
    res= fastqCpp(files, numreads)
  
  # multi files multi core
  if (length(files) > 1 & mc==TRUE)
  {
    require(parallel)
    res= mclapply(X=files, FUN=fastqCpp, numreads=numreads, mc.preschedule=TRUE, mc.cores=detectCores()+1)  
  }
  
  # multi files single core
  if (length(files) > 1 & mc==FALSE)
  {
    res= lapply(X=files, FUN=fastqCpp, numreads=numreads)
  }
  
  return(res)
  #.Call("fastqCpp", file, package="fastqc")  
}















