// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
#include <RcppArmadillo.h>
#include <RcppDist.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec dmvnorm_rcpp(const arma::mat& x, const arma::vec& mean, const arma::mat& sigma, bool logd = false) {

  return dmvnorm(x, mean, sigma, logd);
}

// [[Rcpp::export]]
arma::mat rmvnorm_rcpp(int n, const arma::vec& mean, const arma::mat& sigma) {
  return rmvnorm(n, mean, sigma);
}