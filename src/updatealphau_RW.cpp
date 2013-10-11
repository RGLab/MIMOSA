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


void updatealphau_RW(vector<double>& xalphaut, vector<int>& xn_s, vector<int>& xn_u, int xI, int xK, vector<double>& xlambda_u, Rcpp::NumericVector& sqrt_var, int xtt, vector<int>& xgammat,vector<int>& xAalphau)
{        
    double log1 = 0.0; double log2 = 0.0; double sum_alphau = 0.0; int flag1 = 0; int flag0 = 0;int flagkk = 0; double sum_alphau1 = 0.0;
    int temp=0; int lp0=0; int lp1 = 0;
    for (int kk = 0; kk < xK; kk++) {
        Rcpp::NumericVector alpha_u_p = Rcpp::rnorm(1,xalphaut[kk], sqrt_var[kk]);       
        if (alpha_u_p[0]>0.0 && alpha_u_p[0]<=xlambda_u[kk]) {
        log1 = 0.0;
        log2 = 0.0;
        double alp[xK];
        for (int i = 0; i<xK; i++) {
            alp[i] = xalphaut[i];
          }
        alp[kk] = alpha_u_p[0];
        sum_alphau = 0.0; sum_alphau1 = 0.0;
        for (int s = 0; s < xK; s++) {
            sum_alphau += xalphaut[s]; 
            sum_alphau1 += alp[s];
        }
        log2 -= xI*lgamma(xalphaut[kk]);  
        log2 += xI*lgamma(sum_alphau);       
        log1 -= xI*lgamma(alp[kk]);
        log1 += xI*lgamma(sum_alphau1);
        
        for (int i = 0; i < xI; i++) {
            lp1 = 0; 
            for (int k = 0; k < xK; k++) {
              if (xgammat[xI*k+i] == 1) { lp1 +=1;}       
            }
            lp0 = xK-lp1;
            int p1[lp1];  flag1 = 0;
            int p0[lp0];  flag0 = 0;
            flagkk = 0; // whether gamma_k = 1
           
            for (int k= 0; k < xK; k++) {
               if (xgammat[xI*k+i] == 1) {
                  p1[flag1] = k;
                  flag1 += 1;
                  if (k == kk) {flagkk = 1;}
               } else {
                  p0[flag0] = k;
                  flag0 +=1;
               }
            }
            if (flagkk==1) {
               log2 += lgamma(xn_u[i+xI*kk]+xalphaut[kk]);
               log1 += lgamma(xn_u[i+xI*kk]+alp[kk]);
               double sum_nualphau = 0.0;
               double sum_nusalphau = 0.0;
               double sum_nualphau1 = 0.0;
               double sum_nusalphau1 = 0.0;
               for (int k = 0; k<lp1; k++) {
                 temp = i+xI*p1[k];
                   double sum = xn_u[temp]+xalphaut[p1[k]];
                   double sum1 = xn_u[temp]+alp[p1[k]];
                   sum_nualphau += sum;
                   sum_nusalphau += (sum+xn_s[temp]);
                    sum_nualphau1 += sum1;
                    sum_nusalphau1 += (sum1+xn_s[temp]);
               }
               log2 -=lgamma(sum_nualphau);
               log2 += lgamma(sum_nusalphau+1);
               log1 -=lgamma(sum_nualphau1);
               log1 += lgamma(sum_nusalphau1+1);
               for (int k= 0; k<lp0; k++) {
                    temp = i+xI*p0[k];
                    sum_nusalphau +=(xn_u[temp]+xalphaut[p0[k]]+xn_s[temp]);
                    sum_nusalphau1 +=(xn_u[temp]+alp[p0[k]]+xn_s[temp]);
               }
             
               log2 -= lgamma(sum_nusalphau+1);
                log1 -= lgamma(sum_nusalphau1+1);
            } else {
               log2 += lgamma(xn_u[i+xI*kk]+xalphaut[kk]+xn_s[i+xI*kk]);
               log1 += lgamma(xn_u[i+xI*kk]+alp[kk]+xn_s[i+xI*kk]);
               double sum_nusalphau = 0.0; double sum_nusalphau1 = 0.0;
               for ( int k = 0; k<xK; k++) {
                   temp = i+xI*k;
                   sum_nusalphau +=xn_u[temp]+xalphaut[k]+xn_s[temp];
                   sum_nusalphau1 +=xn_u[temp]+alp[k]+xn_s[temp];
               }
               log2 -= lgamma(sum_nusalphau+1);
               log1 -= lgamma(sum_nusalphau1+1);
               
           }
 
        }
            if (log(Rcpp::as<double>(Rcpp::runif(1)) ) <= (log1 - log2)) {
                xalphaut[kk] = alp[kk];
                xAalphau[kk] = 1;
            } 
        }
    }
   
}

