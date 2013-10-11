#include <Rcpp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
//#include <math.h>
#include <stdlib.h>
//#include <time.h>


#include <R_ext/Utils.h>
using std::vector;
void updatealphau(vector<double>& xalphaut, vector<int>& xn_s, vector<int>& xn_u, int xI, int xK, vector<double>& xlambda_u, vector<double>& sqrt_var, int xtt, vector<int>& xgammat,vector<int>& xAalphau);
