#include <Rcpp.h>
using namespace Rcpp;


// signC: helper-function: the sign-function
// input: double x: the number to determine the sign of
// output: int: 1, if x is positive, 0 if x is 0, -1, if x is negative
// [[Rcpp::export]]
int signC(double x) {
  if (x > 0) return 1;
  if (x < 0) return -1;
  return 0;
}


// calc: helper-function: counting the number of sign changes
// [[Rcpp::export]]
int calc(NumericVector x, int k, int ind, int signOld) {
  int erg = 0;
  if (signC(x[ind]) == signOld) return 0;
  if (k == 1) return 1;
  int n = x.size();
  for (int i = ind + 1; i < n - k + 2; i++)
    erg += calc(x, k - 1, i, signC(x[ind]));
  return erg;
}

// calcDepth_def: helper-function: calculates the K Sign Depth via definition
// input: res: Vector of residuals; k: Parameter of the K Sign Depth.
// output: The K Sign Depth of res.
// [[Rcpp::export]]
double calcDepth_def(NumericVector res, int k) {
  int n = res.size();
  if(n < k)
    stop("Can't calculate Depth when k is larger than the numbers of residuals");
  int num = 0;
  for(int i = 0; i < n - k + 1; i++) {
    num += calc(res, k, i, 0);
  }
  return num / Rf_choose(n, k);
}


//RcppCalcDepthBlock: helper-function: calculates the K-sign depth via block
//                                     implementation
// input: res: Vector of residuals; k: Parameter of the K Sign Depth.
// output: The K Sign Depth of res.
// [[Rcpp::export]]
double RcppCalcKDepthBlock(NumericVector resSigns, int K) {

  int N = resSigns.size();

  NumericVector temp(N);
  int z = 0; // length of block
  int y = 0; // block number for correct index
  for (int i = 0; i < N - 1; i++) {
    z = z + 1;
    if (resSigns[i+1] + resSigns[i] == 0) {
      temp(y) = z;
      z = 0;
      y = y + 1;
    }
  }
  z = z + 1;
  temp(y) = z;

  NumericVector blocks(temp.size());
  blocks = temp[temp > 0];
  int B = blocks.size();


  if (B < K) {
    return 0;
  } else if (B == K) {
    NumericVector value = cumprod(blocks);
    return value[B-1];
  } else if (B == K + 1) {
    blocks[0] = blocks[0] + blocks[B-1];
    NumericVector value = cumprod(blocks);
    return value[B-2];
  }

  if ((B + K) % 2 == 1) {

    // can be reduced to the even case by merging first and last block

    blocks[0] = blocks[0] + blocks[B-1];
    B = B - 1;

  }

  double value = 0;
  NumericMatrix Q1(K-1,(B-K)/2+1); // for odd index
  NumericMatrix Q2(K-1,(B-K)/2); // for even index

  // odd index

  Q1(0, 0) = blocks[B-1];

  for (int i = 1; i < (B - K) / 2 + 1; i++) {
    Q1(0, i) = Q1(0, i - 1) + blocks[B-1-2*i];
  }

  for (int k = 1; k < K - 1; k++) {
    Q1(k, 0) = Q1(k-1, 0) * blocks[B-1-k];

    for (int i = 1; i < (B - K) / 2 + 1; i++) {
      Q1(k, i) = Q1(k, i-1) + Q1(k-1, i) * blocks[B-1-2*i-k];
    }
  }

  // even index

  Q2(0, 0) = blocks[B-1-1];

  for (int i = 1; i < (B - K) / 2; i++) {
    Q2(0, i) = Q2(0, i-1) + blocks[B-1-1-2*i];
  }

  for (int k = 1; k < K-1; k++) {
    Q2(k, 0) = Q2(k-1, 0) * blocks[B-1-1-k];

    for (int i = 1; i < (B - K) / 2; i++) {
      Q2(k, i) = Q2(k, i-1) + Q2(k-1, i) * blocks[B-1-1-2*i-k];
    }
  }

  // do the summation

  // odd and even in one loop:

  for (int b = 0; b < (B - K) / 2; b++) {
    value = value + blocks[2*b] * Q1(K-2, (B-K)/2-b);
    value = value + blocks[2*b+1] * Q2(K-2, (B-K)/2-1-b);
  }

  // Note that Q1 has one column more than Q2!

  // last case for the odd index!

  value = value + blocks[B-K] * Q1(K-2, 0);

  return value;

}
