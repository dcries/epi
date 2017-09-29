#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
// // [[Rcpp::export]]
// arma::uvec sub2(arma::mat x, arma::vec index, int val){
//   arma::uvec ids = find(index == val);
//   return (ids);
// }
// 
// // [[Rcpp::export]]
// arma::mat sub(arma::mat x, arma::vec index, int val){
//   arma::uvec ids = find(index == val);
//   return x.rows(ids);
// }

// // [[Rcpp::export]]
// arma::vec sub2(arma::vec x, arma::vec index, int val){
//   arma::uvec ids = find(index == val);
//   return x(ids);
// }
// [[Rcpp::export]]
arma::ivec pick(IntegerVector x, arma::vec probs, int k){
  //Rcpp::IntegerVector frame = Rcpp::Range(0,k-1);
  arma::ivec out;
  out = Rcpp::RcppArmadillo::sample(x, k, true, probs);
  return out;
}
