#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
IntegerVector kmerIdxCpp(const std::string& str, const int window, const int n, CharacterVector &kmer){
  CharacterVector result( window );
  for ( int j = 0; j < window; j++ ){ 
    result[j] = str.substr( j, n );
  }
  IntegerVector out = match( result , kmer );
  return out;
}

// [[Rcpp::export]]
IntegerMatrix kmerPositionFrequencyCpp(StringVector &string_vector, IntegerVector &strand_vector, const int window, const int w, CharacterVector &kmer){
  
  // Initialize Matrix
  IntegerMatrix out = IntegerMatrix(kmer.size(),window);
  rownames(out) = kmer;
  
  // Get Constants
  int n = string_vector.size();
  std::string str_i;

  //Simple Vector for Storing matches
  IntegerVector m(window);

  for(int i=0; i<n; i++){

    str_i = string_vector[i];

    // Match Kmer over window
    m = kmerIdxCpp(str_i,window,w,kmer);

    for(int j = 0; j < window; j++){

      if(!IntegerVector::is_na(m[j])){

        if(strand_vector[i] == 2){ // Minus Stranded

          out( m[j] - 1, window - j - 1) = out( m[j] - 1, window - j - 1) + 1; 

        }else{ // Other / Plus Stranded
          
          out( m[j] - 1, j) = out( m[j] - 1, j) + 1; 

        }
      }
    }
  }

  return out;

}

// [[Rcpp::export]]
IntegerMatrix kmerIDFrequencyCpp(StringVector &string_vector, IntegerVector &id_vector, const int n_id, const int window, const int w, CharacterVector &kmer){
  
  // Initialize Matrix
  IntegerMatrix out = IntegerMatrix(kmer.size() , n_id);
  rownames(out) = kmer;
  
  // Get Constants
  int n = string_vector.size();
  std::string str_i;
  int id_i;

  //Simple Vector for Storing matches
  IntegerVector m(window);

  for(int i = 0; i < n; i++){

    str_i = string_vector[i];
    id_i = id_vector[i];

    // Match Kmer over window
    m = kmerIdxCpp(str_i, window, w, kmer);

    // Add Matched Value if not NA ie containing an N
    for(int j = 0; j < window; j++){

      if(!IntegerVector::is_na(m[j])){

        out(m[j] - 1, id_i - 1) = out(m[j] - 1, id_i - 1) + 1; 

      }

    }

  }

  return out;
}




