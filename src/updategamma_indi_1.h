#include <Rcpp.h>

#include <R_ext/Utils.h>
#include <stdlib.h>
using std::vector;

void updategamma_indi_1(vector<int>& xn_s,vector<int>& xn_u, vector<int>& gamma_tt, int xI, int xK, int xM, int K1,
                     vector<double>& xalphau, vector<double>& xalphas, double xa, double xb, Rcpp::IntegerVector& xmk, int& xIstar, int& xmKstar, 
                     Rcpp::NumericVector& xpp, double xpb1, double xpb2, double xlambda, double xbeta, double xalpha, vector<double>& xmu_s, 
                     vector<double>& xmu_u, Rcpp::IntegerMatrix& xd, Rcpp::NumericMatrix& xybar_s, Rcpp::NumericMatrix& xybar_u, 
                     Rcpp::NumericMatrix& pt1, Rcpp::NumericMatrix& pt2,Rcpp::NumericMatrix& pt3, Rcpp::IntegerMatrix& xindicator, vector<int>& xAg,
                     double ab, double abI, double bI);