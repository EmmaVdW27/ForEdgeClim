
#include <RcppArmadillo.h>
// Function that searches for the fastest way to solve the matrix equation.
// The method most commonly used for square and non-symmetrical matrices is the LU-decomposition (Lower-Upper decomposition).
// This method decomposes the matrix S into 2 triangular matrices that are more easily inverted than square matrices.
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec solve_rtm(const arma::mat& S, const arma::vec& y) {
  return arma::solve(S, y, arma::solve_opts::fast);
}
