#include <Rcpp.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <stdlib.h>
#include <R_ext/Utils.h>
#include <boost/math/special_functions/digamma.hpp>


RcppExport SEXP updatealphas_Exp(SEXP alphast,SEXP n_s, SEXP K, SEXP I, SEXP lambda_s, SEXP gammat, SEXP var_1,SEXP var_2, SEXP p_var, SEXP ttt)
{
    BEGIN_RCPP
    Rcpp::NumericVector xalphast(alphast);

    Rcpp::IntegerMatrix xn_s(n_s);
    Rcpp::IntegerMatrix xgammat(gammat);
    int xI = Rcpp::as<int>(I);
    int xK = Rcpp::as<int>(K);
    Rcpp::NumericVector sqrt_var1(var_1);
    Rcpp::NumericVector sqrt_var2(var_2);
    int xtt = Rcpp::as<int>(ttt);
    Rcpp::NumericVector xlambda_s(lambda_s);
    Rcpp::IntegerVector xAalphas(xK);

    Rcpp::RNGScope scope;
    Rcpp::NumericVector xp_var(p_var); //proposal mixture

    double delF = 0.0;
    double psik =0.;
    double log1 = 0.0;
    double log2 = 0.0;
    double sums = 0.;
    double sum_alp_ns = 0.0;
    double sum_alp = 0.0;
    double sum_gl_alp = 0.0;
    double sum_gl_alp_ns = 0.0;
    int flag1 =0; int flagkk = 0; int lp1 = 0; 
    for (int kk = 0; kk < xK; kk++) {
        delF = 0.0;
        psik =boost::math::digamma(xalphast[kk]);
        log1 = 0.0;
        log2 = 0.0;
        for (int i = 0; i < xI; i++) {
            lp1 = 0; 
            for (int k = 0; k < xK; k++) {
              if (xgammat(i,k) == 1) { lp1 +=1;}       
            }

            int p1[lp1]; flag1 =0;
            flagkk = 0;
            for (int k= 0; k < xK; k++) {
               if (xgammat(i,k) == 1) {
                  p1[flag1] = k;
                  flag1 += 1;
                  if (k == kk) {flagkk = 1;}
               } 
            }
            sum_alp_ns = 0.0;
            sum_alp = 0.0;
            sum_gl_alp = 0.0;
            sum_gl_alp_ns = 0.0;
            for (int k = 0; k<lp1; k++){
                sums = xalphast[p1[k]] + xn_s(i,p1[k]);
                sum_alp_ns += sums;
                sum_alp += xalphast[p1[k]];
                sum_gl_alp += lgamma(xalphast[p1[k]]);
                sum_gl_alp_ns += lgamma(sums);
            }
            if (flagkk > 0) {
               delF += boost::math::digamma(xn_s(i,kk)+xalphast[kk]) - psik - boost::math::digamma(sum_alp_ns) + boost::math::digamma(sum_alp);
            }
            if (lp1 >0) {
               log2 += -(sum_gl_alp-lgamma(sum_alp)) + (sum_gl_alp_ns - lgamma(sum_alp_ns));
            }
        }
        double mean_p = std::max(0.01, xalphast[kk]+delF/xtt);
        Rcpp::NumericVector alpha_s_p= Rcpp::rnorm(1, mean_p, sqrt_var1[kk]);

        if (Rcpp::as<double>(Rcpp::rbinom(1,1,xp_var[kk])) == 1) {
            alpha_s_p= Rcpp::rnorm(1, mean_p, sqrt_var1[kk]);}
        else { alpha_s_p = Rcpp::rnorm(1, mean_p, sqrt_var2[kk]);}    

        if (alpha_s_p[0]>0.0) {
           double alp[xK];
        
           for (int i = 0; i<xK; i++) {
               alp[i] = xalphast[i];
           }
           alp[kk] = alpha_s_p[0];
           log2 += log(xp_var[kk]*gsl_ran_gaussian_pdf(alp[kk]-mean_p, sqrt_var1[kk])+(1-xp_var[kk])*gsl_ran_gaussian_pdf(alp[kk]-mean_p, sqrt_var2[kk])); 
           delF = 0.0; psik = boost::math::digamma(alp[kk]);
           for ( int i = 0; i < xI; i++) {
               lp1 = 0; 
               for (int k = 0; k < xK; k++) {
                  if (xgammat(i,k) == 1) { lp1 +=1;}       
               }

               int p1[lp1]; flag1 =0;
               flagkk = 0;
               for (int k= 0; k < xK; k++) {
                  if (xgammat(i,k) == 1) {
                     p1[flag1] = k;
                     flag1 += 1;
                     if (k == kk) {flagkk = 1;}
                   } 
               }

               sum_alp_ns = 0.0;
               sum_alp = 0.0;
               sum_gl_alp = 0.0;
               sum_gl_alp_ns = 0.0;
               for (int k = 0; k<lp1; k++){
                   sums = alp[p1[k]] + xn_s(i,p1[k]);
                   sum_alp_ns += sums;
                   sum_alp += alp[p1[k]];
                   sum_gl_alp += lgamma(alp[p1[k]]);
                   sum_gl_alp_ns += lgamma(sums);
               }
               if (flagkk > 0) {
                   delF += boost::math::digamma(xn_s(i,kk)+xalphast[kk]) - psik - boost::math::digamma(sum_alp_ns) +boost::math::digamma(sum_alp);
               }
               if (lp1 >0) {
                   log1 += -(sum_gl_alp-lgamma(sum_alp)) + (sum_gl_alp_ns - lgamma(sum_alp_ns));
               }
           }
           mean_p = std::max(0.01, alp[kk] + delF/xtt);
           log1 +=log(xp_var[kk]*gsl_ran_gaussian_pdf(xalphast[kk]-mean_p, sqrt_var1[kk])+(1-xp_var[kk])*gsl_ran_gaussian_pdf(xalphast[kk]-mean_p, sqrt_var2[kk])); 
           log1 += log(gsl_ran_exponential_pdf(alp[kk],xlambda_s[kk])); //exponential prior
           log2 += log(gsl_ran_exponential_pdf(xalphast[kk],xlambda_s[kk]));//exponential prior
           //if (alp[kk]<0 || alp[kk]>xlambda_s[kk]) {log1+=log(0);} //Uniform prior
           //if (xalphast[kk]<0 || xalphast[kk]>xlambda_s[kk]) {log2+=log(0);} //Uniform prior

           if (log(Rcpp::as<double>(Rcpp::runif(1)) ) <= (log1 - log2)) {
                xalphast[kk] = alp[kk];
                xAalphas[kk] = 1;
           } else{
                xAalphas[kk] = 0;
           }
        }
    
    }

    return Rcpp::List::create(Rcpp::Named("alphas_tt") = xalphast, Rcpp::Named("Aalphas") = xAalphas);

    END_RCPP
}


