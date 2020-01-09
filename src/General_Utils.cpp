#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

// // [[Rcpp::export]]
// IntegerMatrix tabulate1dCpp(IntegerVector x1, int xmin, int xmax){
//   IntegerVector x = clone(x1);
//   IntegerVector r = seq(xmin,xmax);
//   IntegerMatrix out(r.size(),2);
//   out(_, 0) = r;
//   int n = x.size();
//   int xi;
//   for(int i = 0; i < n; i++){
//     xi = (x[i] - xmin);
//     if(xi >= 0 && xi <= r.size()){
//       out( xi , 1 ) = out( xi , 1 ) + 1;
//     }
//   }
//   return out;
// }

// [[Rcpp::export]]
IntegerMatrix tabulate2dCpp(IntegerVector &x, int &xmin, int &xmax, IntegerVector &y, int &ymin, int &ymax){
  if(x.size() != y.size()){
    stop("width must equal size!");
  }
  int n = x.size();
  IntegerVector rx = seq(xmin, xmax);
  IntegerVector ry = seq(ymin, ymax);
  IntegerMatrix mat( ry.size() , rx.size() );
  int rys = ry.size();
  int rxs = rx.size();
  int xi,yi;
  for(int i = 0; i < n; i++){
    xi = (x[i] - xmin);
    yi = (y[i] - ymin);
    if(yi >= 0 && yi < rys){
      if(xi >= 0 && xi < rxs){
        mat( yi , xi ) = mat( yi , xi ) + 1; 
      }
    }
  }
  return mat;
}

// // [[Rcpp::export]]
// IntegerMatrix tabulate2dCpp(IntegerVector x1, int xmin, int xmax, IntegerVector y1, int ymin, int ymax){
//   if(x1.size() != y1.size()){
//     stop("width must equal size!");
//   }
//   IntegerVector x = clone(x1);
//   IntegerVector y = clone(y1);
//   int n = x.size();
//   IntegerVector rx = seq(xmin,xmax);
//   IntegerVector ry = seq(ymin,ymax);
//   IntegerMatrix mat( ry.size() , rx.size() );
//   int xi,yi;
//   for(int i = 0; i < n; i++){
//     xi = (x[i] - xmin);
//     yi = (y[i] - ymin);
//     if(yi >= 0 && yi < ry.size()){
//       if(xi >= 0 && xi < rx.size()){
//         mat( yi , xi ) = mat( yi , xi ) + 1; 
//       }
//     }
//   }
//   return mat;
// }

// [[Rcpp::export]]
Rcpp::NumericVector computeSparseRowVariances(IntegerVector j, NumericVector val, NumericVector rm, int n) {
  const int nv = j.size();
  const int nm = rm.size();
  Rcpp::NumericVector rv(nm);
  Rcpp::NumericVector rit(nm);
  int current;
  // Calculate RowVars Initial
  for (int i = 0; i < nv; ++i) {
    current = j(i) - 1;
    rv(current) = rv(current) + (val(i) - rm(current)) * (val(i) - rm(current));
    rit(current) = rit(current) + 1;
  }
  // Calculate Remainder Variance
  for (int i = 0; i < nm; ++i) {
    rv(i) = rv(i) + (n - rit(i))*rm(i)*rm(i);
  }
  rv = rv / (n - 1);
  return(rv);
}

