#include <Rcpp.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <R_ext/Utils.h>
#include <stdlib.h>
using namespace std;
using namespace Rcpp;



void updatemus(vector<double>& xmust, vector<double>& xmuut, vector<int>& xn_s, vector<int>& xn_u, int xI, int xK, int xM, int K1,
               vector<double>& xp_var, Rcpp::NumericVector& sqrt_var1, Rcpp::NumericVector& sqrt_var2, int xtt, vector<int>& xgammat, Rcpp::IntegerMatrix& xd, 
               Rcpp::NumericMatrix& xybar_s, Rcpp::NumericMatrix& xybar_u, Rcpp::NumericMatrix& pt1,Rcpp::NumericMatrix& pt3,
               double xlambda, double xbeta, double xalpha, vector<double>& xms,vector<double>& xSigs ,vector<int>& xAmus)
{  
    int temp=0;  double delF=0.; double log1=0.; double log2=0.; int nik=0; int xIK=xI*xK;
    double alpha_np=0.;;  double beta_np=0.;
    
    for ( int p=0; p<xM; p++) {
    delF=0.; log1=0.; log2=0.; 
    nik=0; 
    for (int i=0; i<xI;i++) {
         for(int k=0; k<K1; k++) {
            temp = xI*k+i;
             if(xd(k,p)==1) {
                nik = xn_u[temp]+xn_s[temp];
                if(xgammat[temp]==1) {
                    if (xn_u[temp] >0 && xn_s[temp]>0) {
                           beta_np = pt3(temp,p)+xbeta+pow((xmuut[p]-xybar_u(temp,p)),2)/(2*(xlambda+1/xn_u[temp]))+
                                        pow((xmust[p]-xybar_s(temp,p)),2)/(2*(xlambda+1/xn_s[temp]));
                           alpha_np = nik/2+xalpha;
                           log2+=-alpha_np*log(beta_np);
                           delF+=(-alpha_np/beta_np)*(xmust[p]-xybar_s(temp,p))/(xlambda+1/xn_s[temp]);    
                     }else if (xn_u[temp]==0 && xn_s[temp]>0) {
                           beta_np = pt1(temp,p)+xbeta+pow((xmust[p]-xybar_s(temp,p)),2)/(2*(xlambda+1/xn_s[temp]));
                           alpha_np= nik/2+xalpha;
                           log2+=-alpha_np*log(beta_np);    
                           delF+=(-alpha_np/beta_np)*(xmust[p]-xybar_s(temp,p))/(xlambda+1/xn_s[temp]);                 
                     }
                 }
              }
          }
    }
    double mean_p = std::max(0.01, xmust[p]+delF/xtt);
    Rcpp::NumericVector mus_p= Rcpp::rnorm(1, mean_p, sqrt_var1[p]);
    if (Rcpp::as<double>(Rcpp::rbinom(1,1,xp_var[p])) == 1) {
         mus_p= Rcpp::rnorm(1, mean_p, sqrt_var1[p]);}
    else {mus_p = Rcpp::rnorm(1, mean_p, sqrt_var2[p]);}
    
    if(mus_p[0]>0.0) {
         log2+=log(xp_var[p]*gsl_ran_gaussian_pdf(mus_p[0]-mean_p, sqrt_var1[p])+(1-xp_var[p])*gsl_ran_gaussian_pdf(mus_p[0]-mean_p, sqrt_var2[p])); 
         delF = 0.;
         for (int i=0; i<xI;i++) {
                for(int k=0; k<K1; k++) {
                  temp = xI*k+i;
                  if(xd(k,p)==1) {
                     nik = xn_u[temp]+xn_s[temp];
                     if(xgammat[temp]==1) {
                        if (xn_u[temp] >0 && xn_s[temp]>0) {
                                 beta_np = pt3(temp,p)+xbeta+pow((xmuut[p]-xybar_u(temp,p)),2)/(2*(xlambda+1/xn_u[temp]))+
                                        pow((mus_p[0]-xybar_s(temp,p)),2)/(2*(xlambda+1/xn_s[temp]));
                                 alpha_np= nik/2+xalpha;
                                 log1+=-alpha_np*log(beta_np);
                                 delF+=(-alpha_np/beta_np)*(mus_p[0]-xybar_s(temp,p))/(xlambda+1/xn_s[temp]);    
                 
                         }else if (xn_u[temp]==0 && xn_s[temp]>0) {
                                 beta_np = pt1(temp,p)+xbeta+pow((mus_p[0]-xybar_s(temp,p)),2)/(2*(xlambda+1/xn_s[temp]));
                                 alpha_np= nik/2+xalpha;
                                log1+=-alpha_np*log(beta_np);    
                                delF+=(-alpha_np/beta_np)*(mus_p[0]-xybar_s(temp,p))/(xlambda+1/xn_s[temp]);                 
                         }
                    }
                  }
                }
          }
          mean_p = std::max(0.01, mus_p[0]+delF/xtt);
          log1+=log(xp_var[p]*gsl_ran_gaussian_pdf(xmust[p]-mean_p, sqrt_var1[p])+(1-xp_var[p])*gsl_ran_gaussian_pdf(xmust[p]-mean_p, sqrt_var2[p]));
          log2+=log(gsl_ran_gaussian_pdf(xmust[p]-xms[p], xSigs[p]));
          log1+=log(gsl_ran_gaussian_pdf(mus_p[0]-xms[p], xSigs[p]));
          //std::cout <<"log1= "<< log1 << " log2 = "<<log2 <<" p="<<p<<std::endl;
          if (log(Rcpp::as<double>(Rcpp::runif(1)) ) <= (log1 - log2)) {
               xmust[p] = mus_p[0];
               xAmus[p] = 1;
          } 
    }
    }
}
