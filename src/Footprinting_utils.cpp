#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
IntegerVector rleSumsStrandedChr(S4 rle, IntegerVector x, IntegerVector strand, int width){
  
  //Adapted from IRanges viewSums

  int ranges_length, *rle_lengths, upper_run, lower_run, lower_bound, upper_bound;
  int strand_i;
  
  // Stuff from RLE
  IntegerVector values =  rle.slot("values");
  IntegerVector lengths =  rle.slot("lengths");
  int rle_length = sum(lengths);
  
  // Stuff from Ranges
  ranges_length = x.size();
  
  // Initialize
  rle_lengths = INTEGER(lengths);
  upper_run = *rle_lengths;
  int index = 0;
  int position, i, y;
  int max_index = lengths.size() - 1;
  int start_i;
  int stretch;
  IntegerVector ans_pos = IntegerVector(width);
  IntegerVector ans_minus = IntegerVector(width);
  IntegerMatrix tmp = IntegerMatrix(3,width);
  
  for (i = 0; i < ranges_length; i++) {
    
    position = 0;
    start_i = x[i];
    strand_i = strand[i];
    
    if(start_i > 0 && start_i < rle_length){
      
      while (index > 0 && upper_run > start_i) {
        upper_run -= *rle_lengths;
        rle_lengths--;
        index--;
      }
      
      while (upper_run < start_i) {
        rle_lengths++;
        index++;
        upper_run += *rle_lengths;
      }
      
      lower_run = upper_run - *rle_lengths + 1;
      upper_bound = start_i + width - 1;
      lower_bound = start_i;
      
      while (lower_run <= upper_bound) {
        
        stretch = (1 + (upper_bound < upper_run ? upper_bound : upper_run) -
                  (lower_bound > lower_run ? lower_bound : lower_run));
        
        if (INTEGER(values)[index] == NA_INTEGER){

          for(y = 0; y < stretch; y++){
            position += 1;
          } 

        }else{

          for(y = 0; y < stretch; y++){
            tmp(strand_i - 1, position) += INTEGER(values)[index];
            position += 1;
          } 

        }

        if (index < max_index) {
          rle_lengths++;
          index++;
          lower_run = upper_run + 1;
          lower_bound = lower_run;
          upper_run += *rle_lengths;
        } else {
          break;
        }

      }
    }
  }
  
  // 0 = +, 1 = -, 2 = *.

  IntegerVector out = tmp(0, _ ) + tmp(2, _ ) + rev(tmp(1, _ ));
  
  return (out);

}

// [[Rcpp::export]]
IntegerVector rleSumsStranded(List rleList, List grList, int width, Function as_integer){
  
  //This will iterate over a coverage object

  IntegerVector strand, debug, start;
  IntegerVector out = IntegerVector(width);
  
  // Clone grList
  grList = Rcpp::clone(grList);

  int n = grList.size();
  int shift = floor(width/2);

  for(int i = 0; i < n; i++){
    //rle
    S4 rle = rleList[i];
    //gr
    S4 gr = grList[i];
    //strand
    S4 gr_strand = gr.slot("strand");
    strand = as_integer(gr_strand);
    //start
    S4 ranges = gr.slot("ranges");
    start = ranges.slot("start");
    start = start - shift;
    out += rleSumsStrandedChr(rle, start, strand, width);
  }
  
  debug = IntegerVector::create(grList.size());
  
  return (out);

}

