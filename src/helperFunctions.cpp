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


// prod_one_factor: helper-function: Will be needed for calculation of the
//                  exact depth in linear runtime
// input: residuals: Vector of residuals; K: Parameter of the K Sign Depth;
//        L: Parameter
// output: Needed result
// [[Rcpp::export]]
double prod_one_factor(NumericVector residuals, int K, int L) {
  int N = residuals.size();
  IntegerVector S = sign(residuals);
  IntegerVector M(N + 1);
  M[0] = 0;
  for(int i = 1; i <= N; i++) {
    int temp = 0;
    if(residuals[i-1] < 0) temp = 1;
    M[i] = M[i-1] + temp;
  }
  double erg = 0;
  for(int j = L; j <= N - K + L; j++) {
    double p1 = 0;
    for(int k1 = 0; k1 <= L - 1; k1++) {
      p1 += pow(-1, k1) * Rf_choose(M[j-1], k1) * Rf_choose(j - M[j-1] - 1, L - 1 - k1);
    }
    double p2 = 0;
    for(int k2 = 0; k2 <= K - L; k2++) {
      p2 += pow(-1, k2) * Rf_choose(M[N] - M[j], k2) * Rf_choose(N - j - M[N] + M[j], K - L - k2);
    }
    erg += j * S[j-1] * p1 * p2;
  }
  return erg;
}


// asym_K_depth: helper-function: Calculates the transformed approximative
//               K Sign Depth
// input: residuals: Vector of residuals; K: Parameter of the K Sign Depth
// output: approximative transformed K Sign Depth
// [[Rcpp::export]]
double asymp_K_depth(NumericVector residuals, int K) {
  int N = residuals.size();
  IntegerVector res = sign(residuals);
  NumericMatrix S(N+1, K-1);
  for(int i = 0; i < K - 1; i++) S(0, i) = 0;
  for(int j = 0; j <= K - 2; j++) {
    for(double k = 1; k <= N; k++)
      S(k, j) = S(k-1, j) + pow(k / N, j) * res[k-1];
  }
  double erg = 0;
  for(int j = 0; j <= K - 2; j++) {
    double tmp = 0;
    for(double n = 2; n <= N; n++) {
      tmp += pow(0.5 - n/N, K - 2 - j) * res[n-1] * S(n-1, j);
    }
    erg += Rf_choose(K - 2, j) * tmp;
  }
  double tmp2 = 1;
  for(double i = N - K + 1; i <= N; i++) tmp2 *= i;
  erg = - K * (K - 1) * pow(N, K - 1) / (2 * tmp2) * erg;
  return erg;
}
