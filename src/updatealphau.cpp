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


void updatealphau(vector<double>& xalphaut, vector<int>& xn_s, vector<int>& xn_u, int xI, int xK, vector<double>& xlambda_u, vector<double>& sqrt_var, int xtt, vector<int>& xgammat,vector<int>& xAalphau)
{    
    
    double delF = 0.0; double log1 = 0.0; double log2 = 0.0; double sum_alphau = 0.0; int flag1 = 0; int flag0 = 0;int flagkk = 0; 
    int temp=0;
    for (int kk = 0; kk < xK; kk++) {
        delF = 0.0;
        log1 = 0.0;
        log2 = 0.0;
        sum_alphau = 0.0;
        for (int s = 0; s < xK; s++) {
            sum_alphau += xalphaut[s];  
        }
        log2 -= xI*lgamma(xalphaut[kk]);
        delF += xI*(boost::math::digamma(sum_alphau)- boost::math::digamma(xalphaut[kk]));
        log2 += xI*lgamma(sum_alphau);
        for (int i = 0; i < xI; i++) {
            int lp1 = 0; 
            for (int k = 0; k < xK; k++) {
              if (xgammat[xI*k+i] == 1) { lp1 +=1;}       
            }
            int lp0 = xK-lp1;
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
               delF +=boost::math::digamma(xn_u[i+xI*kk]+xalphaut[kk]);
               double sum_nualphau = 0.0;
               double sum_nusalphau = 0.0;
               for (int k = 0; k<lp1; k++) {
                 temp = i+xI*p1[k];
                   double sum = xn_u[temp]+xalphaut[p1[k]];
                   sum_nualphau += sum;
                   sum_nusalphau += (sum+xn_s[temp]);
               }
               log2 -=lgamma(sum_nualphau);
               log2 += lgamma(sum_nusalphau+1);
               delF -=boost::math::digamma(sum_nualphau);
               delF += boost::math::digamma(sum_nusalphau+1);
              
               for (int k= 0; k<lp0; k++) {
                    temp = i+xI*p0[k];
                    sum_nusalphau +=(xn_u[temp]+xalphaut[p0[k]]+xn_s[temp]);
               }
               delF -= boost::math::digamma(sum_nusalphau+1);
               log2 -= lgamma(sum_nusalphau+1);
            } else {
               log2 += lgamma(xn_u[i+xI*kk]+xalphaut[kk]+xn_s[i+xI*kk]);
               delF += boost::math::digamma(xn_u[i+xI*kk]+xalphaut[kk]+xn_s[i+kk*xI]);
               double sum_nusalphau = 0.0;
               for ( int k = 0; k<xK; k++) {
                   sum_nusalphau +=xn_u[i+xI*k]+xalphaut[k]+xn_s[i+xI*k];
               }
               log2 -= lgamma(sum_nusalphau+1);
               delF -= boost::math::digamma(sum_nusalphau+1);
           }
 
        }
        double mean_p = std::max(0.01, xalphaut[kk]+delF/xtt);
        Rcpp::NumericVector alpha_u_p = Rcpp::rnorm(1, mean_p, sqrt_var[kk]);
        
        if (alpha_u_p[0]>0.0 && alpha_u_p[0]<=xlambda_u[kk]) {
            double alp[xK];
            for (int i = 0; i<xK; i++) {
               alp[i] = xalphaut[i];
            }
            alp[kk] = alpha_u_p[0];
            log2 += log(gsl_ran_gaussian_pdf(alp[kk]-mean_p, sqrt_var[kk]));
            delF = 0.0; sum_alphau = 0.0;
            for (int s = 0; s < xK; s++) {
                sum_alphau += alp[s];
            }
            log1 -= xI*lgamma(alp[kk]);
            delF += xI*(boost::math::digamma(sum_alphau)- boost::math::digamma(alp[kk]));
            log1 += xI*lgamma(sum_alphau);
            for (int i = 0; i < xI; i++ ){
                int lp1 = 0; 
                for (int k = 0; k < xK; k++) {
                    if (xgammat[xI*k+i] == 1) { lp1 +=1;}       
                 }
                 int lp0 = xK-lp1;
                 int p1[lp1]; flag1 = 0;
                 int p0[lp0]; flag0 = 0;
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
                   log1 += lgamma(xn_u[i+xI*kk]+alp[kk]);
                   delF +=boost::math::digamma(xn_u[i+xI*kk]+alp[kk]);
                   double sum_nualphau = 0.0;
                   double sum_nusalphau = 0.0;
                   for (int k = 0; k<lp1; k++) {
                       temp = i+xI*p1[k];
                       double sum = xn_u[temp]+alp[p1[k]];
                       sum_nualphau += sum;
                       sum_nusalphau += (sum+xn_s[temp]);
                   }
                   log1 -=lgamma(sum_nualphau);
                   log1 += lgamma(sum_nusalphau+1);
                   delF -=boost::math::digamma(sum_nualphau);
                   delF += boost::math::digamma(sum_nusalphau+1);
              
                   for (int k= 0; k<lp0; k++) {
                       sum_nusalphau +=(xn_u[i+xI*p0[k]]+alp[p0[k]]+xn_s[i+xI*p0[k]]);
                   }
                   delF -= boost::math::digamma(sum_nusalphau+1);
                   log1 -= lgamma(sum_nusalphau+1);
                 } else {
                   log1 += lgamma(xn_u[i+xI*kk]+alp[kk]+xn_s[i+xI*kk]);
                   delF += boost::math::digamma(xn_u[i+xI*kk]+alp[kk]+xn_s[i+xI*kk]);
                   double sum_nusalphau = 0.0;
                   for ( int k = 0; k<xK; k++) {
                     temp = i+xI*k;
                      sum_nusalphau +=xn_u[temp]+alp[k]+xn_s[temp];
                   }
                   log1 -= lgamma(sum_nusalphau+1);
                   delF -= boost::math::digamma(sum_nusalphau+1);
                }
                
            }
            mean_p = std::max(0.01, alp[kk] + delF/xtt);
            log1 +=log(gsl_ran_gaussian_pdf(xalphaut[kk]-mean_p, sqrt_var[kk]));
            if (log(Rcpp::as<double>(Rcpp::runif(1)) ) <= (log1 - log2)) {
                xalphaut[kk] = alp[kk];
                xAalphau[kk] = 1;
            } 
        }
    }
   
}


