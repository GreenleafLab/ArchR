#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
Rcpp::IntegerVector determineOverlapCpp(IntegerMatrix m, int overlapCut){

  int k2 = 2 * m.ncol();
  int nr = m.nrow();
  int nUnion;
  int maxOverlap;
  IntegerVector unionVector;
  IntegerVector testVector = IntegerVector(nr);
  IntegerVector nOverlap = IntegerVector(nr);
  NumericVector maxOverlapVector = NumericVector(nr);
  IntegerVector vi;
  IntegerVector vj;

  for (int i = 1; i < nr; i++){
   
    if (i % 500 == 0) Rcpp::Rcout << "Completed Computing KNN Overlap " << i << " of " << nr << endl;
    
    for(int j = 0; j < i; j++){
      
      if(testVector(j) == 0){
        vi = m(i, _);
        vj = m(j, _);
        unionVector = union_( vi , vj );
        nUnion = unionVector.size();
        nOverlap(j) = k2 - nUnion;
      }else{
        nOverlap(j) = 0;
      }
    }

    maxOverlap = max( nOverlap );
    maxOverlapVector(i) = maxOverlap;
    if(maxOverlap > overlapCut){
      testVector(i) = -1;
    }

  }

  return testVector;

}