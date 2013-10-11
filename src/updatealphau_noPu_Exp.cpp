#include <Rcpp.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <stdlib.h>
//#include <time.h>

#include <R_ext/Utils.h>
#include <boost/math/special_functions/digamma.hpp>

RcppExport SEXP updatealphau_noPu_Exp(SEXP alphaut, SEXP n_s, SEXP n_u, SEXP I, SEXP K, SEXP lambda_u, SEXP var_p, SEXP ttt, SEXP gammat)
{
    BEGIN_RCPP
    Rcpp::IntegerMatrix xgammat(gammat);
    Rcpp::NumericVector xalphaut(alphaut);
    Rcpp::IntegerMatrix xn_s(n_s);
    Rcpp::IntegerMatrix xn_u(n_u);
    int xI = Rcpp::as<int>(I);
    int xK = Rcpp::as<int>(K);  
    Rcpp::NumericVector sqrt_var(var_p);
    int xtt = Rcpp::as<int>(ttt);
    Rcpp::NumericVector xlambda_u(lambda_u);
    Rcpp::IntegerVector xAalphau(xK);

    Rcpp::RNGScope scope;

    double delF = 0.0;
    double log1 = 0.0;
    double log2 = 0.0;
    double sum_alphau = 0.0;
    int flag1 = 0; int flag0 = 0; int flagkk = 0;
    int lp0 = 0; int lp1 = 0; 
    double sum_nusalphau = 0.0;
    double sum_nualphau = 0.0;
    double sums = 0.;
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
            lp1 = 0; 
            for (int k = 0; k < xK; k++) {
              if (xgammat(i,k) == 1) { lp1 +=1;}       
            }
            lp0 = xK-lp1;
            int p1[lp1]; flag1 = 0;
            int p0[lp0]; flag0 = 0;
            flagkk = 0; // whether gamma_k = 1
           
            for (int k= 0; k < xK; k++) {
               if (xgammat(i,k) == 1) {
                  p1[flag1] = k;
                  flag1 += 1;
                  if (k == kk) {flagkk = 1;}
               } else {
                  p0[flag0] = k;
                  flag0 +=1;
               }
            }
            if (flagkk==1) {
               log2 += lgamma(xn_u(i,kk)+xalphaut[kk]);
               delF +=boost::math::digamma(xn_u(i,kk)+xalphaut[kk]);
               sum_nualphau = 0.0;
               sum_nusalphau = 0.0;
               for (int k = 0; k<lp1; k++) {
                   sums = xn_u(i,p1[k])+xalphaut[p1[k]];
                   sum_nualphau += sums;
                   sum_nusalphau += (sums+xn_s(i,p1[k]));
               }
               log2 -=lgamma(sum_nualphau);
               log2 += lgamma(sum_nusalphau+1);
               delF -=boost::math::digamma(sum_nualphau);
               delF += boost::math::digamma(sum_nusalphau+1);
              
               for (int k= 0; k<lp0; k++) {
                    sum_nusalphau +=(xn_u(i,p0[k])+xalphaut[p0[k]]+xn_s(i,p0[k]));
               }
               delF -= boost::math::digamma(sum_nusalphau+1);
               log2 -= lgamma(sum_nusalphau+1);
            } else {
               log2 += lgamma(xn_u(i,kk)+xalphaut[kk]+xn_s(i,kk));
               delF += boost::math::digamma(xn_u(i,kk)+xalphaut[kk]+xn_s(i,kk));
               sum_nusalphau = 0.0;
               for ( int k = 0; k<xK; k++) {
                   sum_nusalphau +=xn_u(i,k)+xalphaut[k]+xn_s(i,k);
               }
               log2 -= lgamma(sum_nusalphau+1);
               delF -= boost::math::digamma(sum_nusalphau+1);
           }
 
        }
        double mean_p = std::max(0.01, xalphaut[kk]+delF/xtt);
        Rcpp::NumericVector alpha_u_p = Rcpp::rnorm(1, mean_p, sqrt_var[kk]);
        if (alpha_u_p[0]>0.0) {
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
                lp1 = 0; 
                for (int k = 0; k < xK; k++) {
                    if (xgammat(i,k) == 1) { lp1 +=1;}       
                 }
                 lp0 = xK-lp1;
                 int p1[lp1];  flag1 = 0;
                 int p0[lp0];  flag0 = 0;
                 flagkk = 0; // whether gamma_k = 1
           
                 for (int k= 0; k < xK; k++) {
                     if (xgammat(i,k) == 1) {
                      p1[flag1] = k;
                      flag1 += 1;
                     if (k == kk) {flagkk = 1;}
                     } else {
                       p0[flag0] = k;
                       flag0 +=1;
                     }
                 }
                 if (flagkk==1) {
                   log1 += lgamma(xn_u(i,kk)+alp[kk]);
                   delF +=boost::math::digamma(xn_u(i,kk)+alp[kk]);
                   sum_nualphau = 0.0;
                   sum_nusalphau = 0.0;
                   for (int k = 0; k<lp1; k++) {
                       sums = xn_u(i,p1[k])+alp[p1[k]];
                       sum_nualphau += sums;
                       sum_nusalphau += (sums+xn_s(i,p1[k]));
                   }
                   log1 -=lgamma(sum_nualphau);
                   log1 += lgamma(sum_nusalphau+1);
                   delF -=boost::math::digamma(sum_nualphau);
                   delF += boost::math::digamma(sum_nusalphau+1);
              
                   for (int k= 0; k<lp0; k++) {
                       sum_nusalphau +=(xn_u(i,p0[k])+alp[p0[k]]+xn_s(i,p0[k]));
                   }
                   delF -= boost::math::digamma(sum_nusalphau+1);
                   log1 -= lgamma(sum_nusalphau+1);
                 } else {
                   log1 += lgamma(xn_u(i,kk)+alp[kk]+xn_s(i,kk));
                   delF += boost::math::digamma(xn_u(i,kk)+alp[kk]+xn_s(i,kk));
                   sum_nusalphau = 0.0;
                   for ( int k = 0; k<xK; k++) {
                      sum_nusalphau +=xn_u(i,k)+alp[k]+xn_s(i,k);
                   }
                   log1 -= lgamma(sum_nusalphau+1);
                   delF -=boost::math::digamma(sum_nusalphau+1);
                }
                
            }
            mean_p = std::max(0.01, alp[kk] + delF/xtt);
            log1 +=log(gsl_ran_gaussian_pdf(xalphaut[kk]-mean_p, sqrt_var[kk]));
            log1 += log(gsl_ran_exponential_pdf(alp[kk],xlambda_u[kk])); //exponential prior
            log2 += log(gsl_ran_exponential_pdf(xalphaut[kk],xlambda_u[kk])); //exponential prior
            //if (alp[kk]<0 || alp[kk]>xlambda_u[kk]) {log1+=log(0);} //Uniform prior
            //if (xalphaut[kk]<0 || xalphaut[kk]>xlambda_u[kk]) {log2+=log(0);} //Uniform prior
            
            if (log(Rcpp::as<double>(Rcpp::runif(1)) ) <= (log1 - log2)) {
                xalphaut[kk] = alp[kk];
                xAalphau[kk] = 1;
            } else{
                xAalphau[kk] = 0;
            }
        }
    }



    return Rcpp::List::create(Rcpp::Named("alphau_tt") = xalphaut, Rcpp::Named("Aalphau") = xAalphau);

    END_RCPP
}


