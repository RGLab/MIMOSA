#include <Rcpp.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <stdlib.h>
//#include <time.h>

#include <R_ext/Utils.h>
#include <boost/math/special_functions/digamma.hpp>

using namespace std;
using namespace Rcpp;


void updatealphas_RW(vector<double>& xalphast,vector<int>& xn_s, int xK, int xI, vector<double>& xlambda_s, vector<int>& xgammat, 
                  Rcpp::NumericVector& sqrt_var1, Rcpp::NumericVector& sqrt_var2, vector<double>& xp_var, int xtt, vector<int>& xAalphas)
{

    double log1 = 0.0;
    double log2 = 0.0;  int flag1=0;   int flagkk=0; int lp1=0;
    double sum_alp_ns = 0.0;
    double sum_alp = 0.0;
    double sum_gl_alp = 0.0;
    double sum_gl_alp_ns = 0.0;
    double sum_alp_ns1 = 0.0;
    double sum_alp1 = 0.0;
    double sum_gl_alp1 = 0.0;
    double sum_gl_alp_ns1 = 0.0;
    for (int kk = 0; kk < xK; kk++) {
       Rcpp::NumericVector alpha_s_p= Rcpp::rnorm(1, xalphast[kk], sqrt_var1[kk]);
       if (Rcpp::as<double>(Rcpp::rbinom(1,1,xp_var[kk])) == 1) {
           alpha_s_p= Rcpp::rnorm(1,xalphast[kk], sqrt_var1[kk]);}
       else { alpha_s_p = Rcpp::rnorm(1, xalphast[kk], sqrt_var2[kk]);}
      if (alpha_s_p[0]>0.0 && alpha_s_p[0]<=xlambda_s[kk]) {
        log1 = 0.0;
        log2 = 0.0;
        double alp[xK];
           for (int i = 0; i<xK; i++) {
               alp[i] = xalphast[i];
           }
        alp[kk] = alpha_s_p[0];
           
        for (int i = 0; i < xI; i++) {
            lp1=0;
            for (int k = 0; k < xK; k++) {
              if (xgammat[i+xI*k] == 1) { lp1 +=1;}
            }
            int p1[lp1]; flag1=0; flagkk=0;
            for (int k = 0; k < xK; k++) {
              if (xgammat[i+xI*k] == 1) { 
                p1[flag1] = k;
                flag1 += 1;
                if (k == kk) {flagkk = 1;} 
              }       
            }
            sum_alp_ns = 0.0;
            sum_alp = 0.0;
            sum_gl_alp = 0.0;
            sum_gl_alp_ns = 0.0;
            sum_alp_ns1 = 0.0;
            sum_alp1 = 0.0;
            sum_gl_alp1 = 0.0;
            sum_gl_alp_ns1 = 0.0;
            for (int k = 0; k<lp1; k++){
                double sums = xalphast[p1[k]] + xn_s[i+xI*p1[k]];
                double sums1 =alp[p1[k]] + xn_s[i+xI*p1[k]];
                sum_alp_ns += sums; sum_alp_ns1 += sums1;
                sum_alp += xalphast[p1[k]]; sum_alp1 += alp[p1[k]];
                sum_gl_alp += lgamma(xalphast[p1[k]]); sum_gl_alp1 += lgamma(alp[p1[k]]);
                sum_gl_alp_ns += lgamma(sums); sum_gl_alp_ns1 += lgamma(sums1);
            }
            if (lp1>0) {
               log2 += -(sum_gl_alp-lgamma(sum_alp)) + (sum_gl_alp_ns - lgamma(sum_alp_ns));
               log1 += -(sum_gl_alp1-lgamma(sum_alp1)) + (sum_gl_alp_ns1 - lgamma(sum_alp_ns1));
            }
        }
             
           if (log(Rcpp::as<double>(Rcpp::runif(1))) <= (log1 - log2)) {
                xalphast[kk] = alp[kk];
                xAalphas[kk] = 1;
           } 
        }
    
    }

}
