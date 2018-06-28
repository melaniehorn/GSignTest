#include <Rcpp.h>
using namespace Rcpp;

// signC: helper-function: the sign-function
// input: double x: the number to determine the sign of
// output: int: 1, if x is positive, 0 if x is 0, -1, if x is negative
int signC(double x) {
  if (x > 0) {
    return 1;
  } else if (x == 0) {
    return 0;
  } else {
    return -1;
  }
}


// calcSwitchSign: helper-function: determines if a vector has alternating signes
// input: NumericVector res: The residuals
//        IntegerVector ind: the indices you want to look at (ind has length k)
// output: double: 1 if res[ind] has alternating signes, 0 otherwise
double calcSwitchSign(NumericVector res, IntegerVector ind) {
  int k = ind.size();
  double counter = 1;
  int tmp = signC(res[ind[0]]);
  for(int j = 1; j < k; j++) {
    if(signC(res[ind[j]]) == tmp) {
      counter = 0;
      break;
    } else {
      tmp *= -1;
    }
  }
  return counter;
}


// calcNextVec: helper-function: calculates the "next" combination of
//              residuals given an old combination (Ordering is like in the
//              R-function combn())
// input: IntegerVector oldV: an "old" index-vector
//        int n: the length of the residuals
//        int k: the length of oldV
// output: IntegerVector: output has length k. The "next" index-vector
IntegerVector calcNextVec(IntegerVector oldV, int n, int k) {
  IntegerVector newV(k);
  int pos = k - 1;
  while(oldV[pos] == n - k + pos)
    pos--;
  for(int j = 0; j < pos; j++)
    newV[j] = oldV[j];
  newV[pos] = oldV[pos] + 1;
  for(int j = pos + 1; j < k; j++)
    newV[j] = newV[j - 1] + 1;
  return newV;
}


// THE MAIN FUNCTION

//' @title Simplicial k-Depth
//'
//' @description
//' \code{calcDepth} calculates the k-Depth of a given vector of residuals.
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

  double num = 0;
  double denom = Rf_choose(n, k);

  IntegerVector oldV(k);
  for(int i = 0; i < k; i++)
    oldV[i] = i;
  num = calcSwitchSign(res, oldV);

  IntegerVector newV(k);
  for(int i = 1; i < denom; i++) {
    newV = calcNextVec(oldV, n, k);
    num += calcSwitchSign(res, newV);
    oldV = newV;
  }

  return num / denom;
}
