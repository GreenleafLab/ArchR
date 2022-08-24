#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::depends('RcppArmadillo')]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::uvec find_non_zero_indices(arma::sp_mat &X, int col) {
    IntegerVector tmp;
    sp_mat::const_col_iterator start = X.begin_col(col);
    sp_mat::const_col_iterator end = X.end_col(col);
    for ( sp_mat::const_col_iterator i = start; i != end; ++i ){
        tmp.push_back(i.row());
    }
    arma::uvec out = as<uvec>(tmp);
    return out;
}

// [[Rcpp::export]]
Rcpp::List deviations_cpp(
   const arma::sp_mat &X,
   const arma::umat &B,
   arma::sp_mat &anno_mat,
   const arma::mat &expect,
   const arma::mat &CpS,
   const bool verbose = true,
   const std::string &prefix = "",
   const int print_every = 100
   ){
   
    // dev, Z
    arma::mat out_dev(X.n_rows, anno_mat.n_cols);
    arma::mat out_z(X.n_rows, anno_mat.n_cols);

    for (int i = 0; i < anno_mat.n_cols; ++i){

        if(verbose){
            if(i > 0){
                if(i % print_every == 0){
                    Rcout << prefix << " : Deviations for Annotation " << i + 1<< " of " << anno_mat.n_cols << std::endl;
                }
            }           
        }

        //Anno idx
        arma::uvec idx_anno = find_non_zero_indices(anno_mat, i);

        // Foreground Deviations
        arma::mat observed = X.cols(idx_anno) * ones(idx_anno.size(), 1);
        arma::mat expected = accu(expect.rows(idx_anno)) * CpS;
        arma::mat dev = (observed - expected) / expected;

        // Background Deviations
        arma::umat bg_mat = B.rows(idx_anno);
        arma::mat observed_bg(X.n_rows, B.n_cols);
        arma::mat expected_bg(X.n_rows, B.n_cols);

        for (int j = 0; j < B.n_cols; ++j) {

            arma::uvec idx_j = bg_mat.col(j) - 1;
            observed_bg.col(j) = X.cols(idx_j) * ones(idx_anno.size(), 1);
            expected_bg.col(j) = accu(expect.rows(idx_j)) * CpS;

        }
        arma::mat dev_bg = (observed_bg - expected_bg) / expected_bg;

        //Compute Background Deviations Stats
        arma::mat mean_dev_bg = mean(dev_bg, 1);
        arma::mat sd_dev_bg = stddev(dev_bg, 0, 1);

        //Compute Corrected Deviations
        out_dev.col(i) = dev.col(0) - mean_dev_bg.col(0);
        out_z.col(i) = out_dev.col(i) / sd_dev_bg.col(0);

    }

    return(List::create(out_dev, out_z));

}
