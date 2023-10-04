#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

// Pearson Correlation, Adapted from https://github.com/AEBilgrau/correlateR/blob/master/src/auxiliary_functions.cpp
// [[Rcpp::export]]
Rcpp::NumericVector rowCorCpp(IntegerVector idxX, IntegerVector idxY, Rcpp::NumericMatrix X, Rcpp::NumericMatrix Y) {
  
  if(X.ncol() != Y.ncol()){
    stop("Columns of Matrix X and Y must be equal length!");
  }

  if(max(idxX) > X.nrow()){
    stop("Idx X greater than nrow of Matrix X");
  }

  if(max(idxY) > Y.nrow()){
    stop("Idx Y greater than nrow of Matrix Y");
  }
    
  // Transpose Matrices
  X = transpose(X);
  Y = transpose(Y);
  
  const int nx = X.ncol();
  const int ny = Y.ncol();

  // Centering the matrices
  for (int j = 0; j < nx; ++j) {
    X(Rcpp::_, j) = X(Rcpp::_, j) - Rcpp::mean(X(Rcpp::_, j));
  }

  for (int j = 0; j < ny; ++j) {
    Y(Rcpp::_, j) = Y(Rcpp::_, j) - Rcpp::mean(Y(Rcpp::_, j));
  }

  // Compute 1 over the sample standard deviation
  Rcpp::NumericVector inv_sqrt_ss_X(nx);
  for (int i = 0; i < nx; ++i) {
    inv_sqrt_ss_X(i) = 1/sqrt(Rcpp::sum( X(Rcpp::_, i) * X(Rcpp::_, i) ));
  }

  Rcpp::NumericVector inv_sqrt_ss_Y(ny);
  for (int i = 0; i < ny; ++i) {
    inv_sqrt_ss_Y(i) = 1/sqrt(Rcpp::sum( Y(Rcpp::_, i) * Y(Rcpp::_, i) ));
  }

  //Calculate Correlations
  const int n = idxX.size();
  Rcpp::NumericVector cor(n);
  for(int k = 0; k < n; k++){
    cor[k] = Rcpp::sum( X(Rcpp::_, idxX[k] - 1) * Y(Rcpp::_, idxY[k] - 1) ) * inv_sqrt_ss_X( idxX[k] - 1) * inv_sqrt_ss_Y( idxY[k] - 1);
  } 

  return(cor);

}
