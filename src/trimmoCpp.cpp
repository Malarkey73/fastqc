//ILLUMINACLIP: Cut adapter and other illumina-specific sequences from the read.
//SLIDINGWINDOW: Perform a sliding window trimming, cutting once the average quality within the window falls below a threshold.
//LEADING: Cut bases off the start of a read, if below a threshold quality
//TRAILING: Cut bases off the end of a read, if below a threshold quality
//CROP: Cut the read to a specified length
//HEADCROP: Cut the specified number of bases from the start of the read
//MINLEN: Drop the read if it is below a specified length
//TOPHRED33: Convert quality scores to Phred-33
//TOPHRED64: Convert quality scores to Phred-64


#include "fastqCpp.h"
using namespace Rcpp;
using namespace seqan;

// [[Rcpp::export]]

int trimmoCpp(std::string inp1, 
              int minlen, 
              int crop, 
              int headcrop, 
              bool slidingwindow, 
              int slidingwindowlen)

{    
    std::fstream in(inp1.c_str(), std::ios::binary | std::ios::in);
    seqan::RecordReader<std::fstream, seqan::SinglePass<> > reader(in);

    // Read file record-wise.
    seqan::CharString id;
    seqan::Dna5String seq;
    seqan::CharString qual;
    
    unsigned drop = 0; // write unless some drop criterium is met
    
    while (!atEnd(reader))
    {
        if (readRecord(id, seq, qual, reader, seqan::Fastq()) != 0)
            return 1;  // Could not record from file.
        
        //CROP and HEADCROP  

        if(headcrop > 0 && headcrop < length(seq))
          erase(seq, 0, headcrop);
        if(crop > 0 && crop < length(seq))
          erase(seq, length(seq)-crop, length(seq));
        if(headcrop > 0 && headcrop < length(qual))  
          erase(qual, 0, headcrop);
        if(crop > 0 && crop < length(qual))
          erase(qual,length(qual)-crop, length(qual));
          
        //SLIDINGWINDOW
        if(slidingwindow == TRUE)  
        {
          for(unsigned i = 0; i< length(qual); ++i)
            
            (int)ordValue(qual[i])-33;
 
        
        }
        
        if(length(seq) < minlen) // if it is shorter than MINLEN then just leave it.
          drop=1;
          
        if(drop != 1)
        {
          std::cout << id << "\n" << seq << "\n" << qual << "\n" << "\n";
        }
    
    }
        
    
    
    return 1;

}