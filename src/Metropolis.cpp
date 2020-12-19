#include <cmath>
#include <Rcpp.h>
using namespace Rcpp;
//' @title The probability density of standard Laplace distribution using Rcpp
//' @description The probability density of standard Laplace distribution using Rcpp
//' @param x the value of the variable
//' @return the probability density of the value of x
//' @examples
//' \dontrun{
//'   LAPvalue <-f(4)
//' }
//' @export
// [[Rcpp::export]]
double f(double x) {
  return exp(-abs(x));
}

//' @title A generation of standard Laplace distribution using Rcpp
//' @description A generation of standard Laplace distribution using Rcpp
//' @param sigma the parameter of proposal distribution
//' @param x0 the initial values
//' @param N the number of samples
//' @return a random sample of size \code{N}
//' @examples
//' \dontrun{
//' LAPsamples <- Metropolis(.2,5,3000)
//' }
//' @export
// [[Rcpp::export]]
NumericVector Metropolis (double sigma, double x0, int N) {
  NumericVector y(N);
  y[0] = x0; 
  NumericVector u = runif(N);
  for (int i = 1; i < N;i++ ) {
    NumericVector z = rnorm(1, y[i-1], sigma);
    if (u[i] <= (f(z[0]) / f(y[i-1]))){
      y[i] = z[0];
    }
    else { 
      y[i] = y[i-1]; 
    }
  }
  return(y);
} 

