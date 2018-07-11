#include <Rcpp.h>
using namespace Rcpp;

// signC: helper-function: the sign-function
// input: double x: the number to determine the sign of
// output: int: 1, if x is positive, 0 if x is 0, -1, if x is negative
int signC(double x) {
  if (x > 0) return 1;
  if (x < 0) return -1;
  return 0;
}


// calc: helper-function: counting the number of sign changes
int calc(NumericVector x, int k, int ind, int signOld) {
  int erg = 0;
  if (signC(x[ind]) == signOld) return 0;
  if (k == 1) return 1;
  int n = x.size();
  for (int i = ind + 1; i < n - k + 2; i++)
    erg += calc(x, k - 1, i, signC(x[ind]));
  return erg;
}


// THE MAIN FUNCTION

//' @title k Sign Depth
//'
//' @description
//' \code{calcDepth} calculates the k Sign Depth of a given vector of residuals.
//'
//' @param res [\code{numeric}]\cr
//'   numeric vector of residuals
//' @param k [\code{integer(1)}]\cr
//'   the k of the depth-formula
//' @return [\code{numeric(1)}] the calculated k-Depth of res.
//'
//' @examples calcDepth(rnorm(10), 3)
//' @export
// [[Rcpp::export]]
double calcDepth(NumericVector res, int k) {
  int n = res.size();
  if(n < k)
    stop("Can't calculate Depth when k is larger than the numbers of residuals");
  int num = 0;
  for(int i = 0; i < n - k + 1; i++) {
    num += calc(res, k, i, 0);
  }
  return num / Rf_choose(n, k);
}
