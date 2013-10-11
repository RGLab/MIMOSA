#include <Rcpp.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
//#include <math.h>
#include <stdlib.h>
//#include <time.h>

#include <R_ext/Utils.h>
using std::vector;

void updatemus(vector<double>& xmust, vector<double>& xmuut, vector<int>& xn_s, vector<int>& xn_u, int xI, int xK, int xM, int K1,
               vector<double>& xp_var, Rcpp::NumericVector& sqrt_var1, Rcpp::NumericVector& sqrt_var2, int xtt, vector<int>& xgammat, Rcpp::IntegerMatrix& xd, 
               Rcpp::NumericMatrix& xybar_s, Rcpp::NumericMatrix& xybar_u, Rcpp::NumericMatrix& pt1,Rcpp::NumericMatrix& pt3,
               double xlambda, double xbeta, double xalpha, vector<double>& xms,vector<double>& xSigs ,vector<int>& xAmus);
