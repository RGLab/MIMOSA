#include <Rcpp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdlib.h>


#include <R_ext/Utils.h>
using std::vector;

void updatealphas(vector<double>& xalphast,vector<int>& xn_s, int xK, int xI, vector<double>& xlambda_s, vector<int>& xgammat, 
                  Rcpp::NumericVector& sqrt_var1,Rcpp::NumericVector& sqrt_var2, vector<double>& xp_var, int xtt, vector<int>& xAalphas);
