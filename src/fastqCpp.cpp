#include "fastqCpp.h"
using namespace Rcpp;
using namespace seqan;

//// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List fastqCpp(std::string argv)
{    
    MultiSeqFile multiSeqFile;
    if (!open(multiSeqFile.concat, argv.c_str(), OPEN_RDONLY))
        return 1;
        
    AutoSeqFormat format;
    guessFormat(multiSeqFile.concat, format);
    split(multiSeqFile, format);    
    unsigned numreads = length(multiSeqFile); // it is read once througn
    seqan::StringSet<seqan::Dna5String> seqs;
    reserve(seqs, numreads, Exact());
    seqan::Dna5String seq, seqLast;
    seqan::CharString qual;
    assignSeq(seq, multiSeqFile[0], format); // I read the first seq and assume all are the same length
    unsigned numcycles = length(seq);
    
    
    // armadillo Matrices
    arma::umat qualMatrix(numreads, numcycles);
    qualMatrix.zeros();   
    arma::umat baseCounterMatrix(5, numcycles); //there are 5 bases ACGTN hence 5 rows.
    baseCounterMatrix.zeros();    
    arma::vec perReadQuality(numreads);
    perReadQuality.zeros();
    //arma::umat kmerMatrix(kmers.size(), numcycles-4); // presume 5 mers
    //kmerMatrix.zeros();
    
    // seqan kmer
    seqan::String<unsigned> kmerCounts;
    seqan::String<double> nucleotideFrequencies; 
    unsigned k = 5;  // Count all 5-mers
    NumericVector kmers(1024); // this is 4^5 i.e. (A/T/G/C)^5
    NumericVector nucs(4);
    
    
    for (unsigned i = 0; i < numreads; ++i) // outer read loop
    {
        assignSeq(seq, multiSeqFile[i], format);// read sequence
        appendValue(seqs, seq, Exact()); // write to the StringSet
        assignQual(qual, multiSeqFile[i], format);  // ascii quality values
        
        //
        countKmers(kmerCounts, nucleotideFrequencies, seq, k);
        for(unsigned km =0; km < 1024; ++km)
        {
          kmers[km] += kmerCounts[km];
        }
        for(unsigned base =0; base<4; ++base)
        {
          nucs[base] += nucleotideFrequencies[base];
        }
  
                     
        for (unsigned cycle = 0; cycle < length(qual) && cycle < length(seq); ++cycle) // inner cycle loop
        {
            qualMatrix(i, cycle) = (int)ordValue(qual[cycle])-33;   //record each qual value
            perReadQuality[i] += (int)ordValue(qual[cycle])-33; //sum each row
        
        // a matrix of baseCounting
          if(seq[cycle] =='A') 
            baseCounterMatrix(0,cycle) += 1;           
          if(seq[cycle] =='C') 
            baseCounterMatrix(1,cycle) +=1;          
          if(seq[cycle] =='G') 
            baseCounterMatrix(2,cycle) +=1;         
          if(seq[cycle] =='T') 
            baseCounterMatrix(3,cycle) +=1;          
          if(seq[cycle] =='N') 
            baseCounterMatrix(4,cycle) +=1;

        } // inner cycles loop
    
    }// end of outer reads loop
    
    
    // counting the kmer frequencies
    kmers= kmers/numreads;
    nucs= nucs/numreads;
    
    
    // check for duplicates
    typedef Iterator<StringSet<Dna5String> >::Type TStringSetIter;    
    std::sort(begin(seqs), end(seqs)); // LEXICOGRAPHICAl SORT
    TStringSetIter seqIt1 = begin(seqs);
    TStringSetIter seqIt2 = seqIt1 + 1 ; // this is 1 in front of the other iterator
    unsigned dupCount = 0;
    NumericVector bins(10);
    unsigned dupThreshold = numreads/10000;
    std::vector<int > dupSeqsCount;
    std::string dupSeqs(numcycles, 'N');
    std::vector<std::string> dupSeqsVec;
    seqan::Dna5String seq1, seq2;
    
    
    for( ; seqIt2 != end(seqs); ++seqIt1, ++seqIt2)
    {
      seq1= value(seqIt1);
      seq2=value(seqIt2);
      // if first and next are the same we start counting how many
      if (seq1 == seq2)
      {        
        ++dupCount;
      }
      // till we come to a different read sequence
      if (seq1 != seq2)
      {
        // bin 1-10 increment
        if(dupCount <= 9)  
          bins[dupCount] += 1;
        // if >10 increment 10th bin
        if(dupCount > 9)
          bins[9] += 1;
          
        // if a threshold is met we record the sequence and how many there were  
        if(dupCount > dupThreshold)
        {
          dupSeqsCount.push_back(dupCount);
          for(unsigned i=0; i< numcycles; ++i)
              dupSeqs[i] = seq1[i];
            
          dupSeqsVec.push_back(dupSeqs);
        }
      // then reset the counter after bin increments and recording is done  
      dupCount = 0;    
      }
    
    }
    
    
    // qualities by read
    perReadQuality = perReadQuality/numcycles; // make the sum phred score to a mean phred score    
    int minQual = arma::min(perReadQuality);
    int maxQual = arma::max(perReadQuality);
    int qualRange = maxQual-minQual;
    arma::uvec qualityBins = hist(perReadQuality, qualRange);
    
    
    // qualities by cycle
    qualMatrix= sort(qualMatrix,0,0);
    int q1 = numreads/4; // integer division
    int mid = numreads/ 2;
    int q3 = q1*3;
    NumericMatrix fivenum(5,numcycles);    
    for (unsigned cycle = 0; cycle < numcycles; ++cycle) 
    { 
        
        fivenum(0,cycle) = qualMatrix(0,cycle);
        fivenum(1,cycle) = qualMatrix(q1,cycle);
        fivenum(2,cycle) = qualMatrix(mid,cycle);
        fivenum(3,cycle) = qualMatrix(q3,cycle);
        fivenum(4,cycle) = qualMatrix(numreads-1,cycle);        
    }

    
    
    
    // return final List of Results
    return List::create(
                    Named("File")           = argv,
                    Named("Reads")          = numreads,
                    Named("Cycles")         = numcycles,
                    Named("Fivenum")        = fivenum,
                    Named("baseByCycle")    = baseCounterMatrix,
                    Named("DupBins")        = bins,
                    Named("DupSeqs")        = dupSeqsVec,
                    Named("DupSeqsCount")   = dupSeqsCount,
                    Named("minQual")        = minQual,
                    Named("maxQual")        = maxQual,
                    Named("qualBySeq")      = qualityBins,
                    Named("kmers")          = kmers,
                    Named("nucFreq")        = nucs);
                      
}
 