#include <Rcpp.h>
using namespace Rcpp;

// signC2: helper-function: the sign-function
// input: double x: the number to determine the sign of
// output: int: 1, if x is positive, 0 if x is 0, -1, if x is negative
int signC2(double x) {
  if (x > 0) {
    return 1;
  } else if (x == 0) {
    return 0;
  } else {
    return -1;
  }
}

//' @title Simplified K-sign depth
//'
//' @description
//' \code{calcSimpDepth} calculates the simplified K-sign depth of a given vector of residuals.
//'
//' @param res [\code{numeric}]\cr
//'   numeric vector of residuals
//' @param K [\code{integer(1)}]\cr
//'   the parameter K of the simplified K-sign depth
//' @return [\code{numeric(1)}] the calculated simplified K-sign depth of res.
//'
//' @examples calcSimpDepth(rnorm(10), 3)
//' @export
// [[Rcpp::export]]
double calcSimpDepth(NumericVector res, int K) {
  int n = res.size();
  if(n < K)
    stop("Can't calculate Depth when k is larger than the numbers of residuals");

  double num = 0;
  double denom = n - K + 1;

  int counter, tmp;
  for(int i = 0; i < denom; i++) {
    counter = 1;
    tmp = signC2(res[i]);
    for(int j = 1; j < K; j++) {
      if(signC2(res[i+j]) == tmp) {
        counter = 0;
        break;
      } else {
        tmp *= -1;
      }
    }
    num += counter;
  }

  return num / denom;
}
