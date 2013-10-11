#include <Rcpp.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <R_ext/Utils.h>
#include <stdlib.h>
using namespace std;
using namespace Rcpp;


void updatebeta_scalar(vector<double>& xmust, vector<double>& xmuut, vector<int>& xn_s, vector<int>& xn_u, int xI, int xK, int xM,
                      int K1, int t1, double xp_var, double sqrt_var1, double sqrt_var2, int xtt, vector<int>& xgammat, 
                      Rcpp::IntegerMatrix& xd, Rcpp::NumericMatrix& xybar_s, Rcpp::NumericMatrix& xybar_u, Rcpp::NumericMatrix& pt1, Rcpp::NumericMatrix& pt2, 
                      Rcpp::NumericMatrix& pt3, double xlambda, Rcpp::NumericVector& beta, double xalpha, double xsig_beta1, 
                      double xbeta1, Rcpp::IntegerVector& xAbeta)
{
    double xbeta = beta[t1]; int tt=t1+1; int temp=0; beta[tt] = xbeta; 
    double delF=0.; double log1=0.; double log2=0.; int nik=0; double beta_np=0.; double alpha_np=0.;
    double alb= xalpha*log(xbeta); double adb = xalpha/xbeta;
    
    for (int p=0; p<xM; p++) { 
    for (int i=0; i<xI;i++) {
         for(int k=0; k<K1; k++) {
           temp = xI*k+i;
             if(xd(k,p)==1) {
                nik = xn_u[temp]+xn_s[temp];
                if(xgammat[temp]==1) {
                    if (xn_u[temp]>0 && xn_s[temp]>0) {
                           beta_np = pt3(temp,p)+xbeta+pow((xmuut[p]-xybar_u(temp,p)),2)/(2*(xlambda+1/xn_u[temp]))+
                                        pow((xmust[p]-xybar_s(temp,p)),2)/(2*(xlambda+1/xn_s[temp]));
                           alpha_np = nik/2+xalpha;
                           log2 += alb-alpha_np*log(beta_np);
                           delF+=adb-alpha_np/beta_np;
                    }else if (xn_u[temp]==0 && xn_s[temp]>0) {
                           beta_np = pt1(temp,p)+xbeta+pow((xmust[p]-xybar_s(temp,p)),2)/(2*(xlambda+1/xn_s[temp]));
                           alpha_np = nik/2+xalpha;
                           log2 += alb-alpha_np*log(beta_np);
                           delF+=adb-alpha_np/beta_np;                 
                    } else if (xn_u[temp]>0 && xn_s[temp]==0) {
                           beta_np = pt2(temp,p)+xbeta+pow((xmuut[p]-xybar_u(temp,p)),2)/(2*(xlambda+1/xn_u[temp]));
                           alpha_np = nik/2+xalpha;
                           log2 += alb-alpha_np*log(beta_np);
                           delF+=adb-alpha_np/beta_np;
                    }         
                }else{
                     if (nik>0) {
                           beta_np =pt3(temp,p)+xbeta+pow((xmuut[p]-(xn_s[temp]*xybar_s(temp,p)+xn_u[temp]*xybar_u(temp,p))/nik),2)/(2*(xlambda+1/nik))+
                                           0.5*xn_s[temp]*pow(xybar_s(temp,p),2)+ 0.5*xn_u[temp]*pow(xybar_u(temp,p),2)-0.5*pow((xn_s[temp]*xybar_s(temp,p)+xn_u[temp]*xybar_u(temp,p)),2)/nik;
                          alpha_np = nik/2+xalpha;
                           log2 += alb-alpha_np*log(beta_np);
                           delF+=adb-alpha_np/beta_np;
                     }       
                }
             }
         }
    }}
    double mean_p = std::max(0.01, xbeta+delF/xtt);
    Rcpp::NumericVector betap= Rcpp::rnorm(1, mean_p, sqrt_var1);
    if (Rcpp::as<double>(Rcpp::rbinom(1,1,xp_var)) == 1) {
         betap= Rcpp::rnorm(1, mean_p, sqrt_var1);}
    else {betap = Rcpp::rnorm(1, mean_p, sqrt_var2);}
    //if(betap[0][0]>0.0 && betap[0][0]<=xlambda_beta[p]) {
//std::cout <<"log2 = "<< log2 <<" betap[0] = "<< betap[0]<< std::endl;
    if(betap[0]>0.0) {
         log2+=log(xp_var*gsl_ran_gaussian_pdf(betap[0]-mean_p, sqrt_var1)+(1-xp_var)*gsl_ran_gaussian_pdf(betap[0]-mean_p, sqrt_var2)); 
         delF = 0.; alb= xalpha*log(betap[0]); adb = xalpha/betap[0];
         for (int p=0; p<xM; p++) {
         for (int i=0; i<xI;i++) {
         for(int k=0; k<K1; k++) {
           temp=xI*k+i;
             if(xd(k,p)==1) {
                nik = xn_u[temp]+xn_s[temp];
                if(xgammat[temp]==1) {
                    if (xn_u[temp] >0 && xn_s[temp]>0) {
                           beta_np = pt3(temp,p)+betap[0]+pow((xmuut[p]-xybar_u(temp,p)),2)/(2*(xlambda+1/xn_u[temp]))+
                                        pow((xmust[p]-xybar_s(temp,p)),2)/(2*(xlambda+1/xn_s[temp]));
                           alpha_np = nik/2+xalpha;
                           log1 += alb-alpha_np*log(beta_np);
                           delF+=adb-alpha_np/beta_np;
                    }else if (xn_u[temp]==0 && xn_s[temp]>0) {
                           beta_np = pt1(temp,p)+betap[0]+pow((xmust[p]-xybar_s(temp,p)),2)/(2*(xlambda+1/xn_s[temp]));
                            alpha_np = nik/2+xalpha;
                           log1 += alb-alpha_np*log(beta_np);
                           delF+=adb-alpha_np/beta_np;                 
                    } else if (xn_u[temp]>0 && xn_s[temp]==0) {
                           beta_np =pt2(temp,p)+betap[0]+pow((xmuut[p]-xybar_u(temp,p)),2)/(2*(xlambda+1/xn_u[temp]));
                           alpha_np = nik/2+xalpha;
                           log1 += alb-alpha_np*log(beta_np);
                           delF+=adb-alpha_np/beta_np;
                    }         
                }else{
                     if (nik>0) {
                           beta_np = pt3(temp,p)+betap[0]+pow((xmuut[p]-(xn_s[temp]*xybar_s(temp,p)+xn_u[temp]*xybar_u(temp,p))/nik),2)/(2*(xlambda+1/nik))+
                                          0.5*xn_s[temp]*pow(xybar_s(temp,p),2)+ 0.5*xn_u[temp]*pow(xybar_u(temp,p),2)-0.5*pow((xn_s[temp]*xybar_s(temp,p)+xn_u[temp]*xybar_u(temp,p)),2)/nik;
                           alpha_np = nik/2+xalpha;
                           log1 += alb-alpha_np*log(beta_np);
                           delF+=adb-alpha_np/beta_np;
                     }       
                }
             }
         }
        }}
        mean_p = std::max(0.01, betap[0]+delF/xtt);
        log1+=log(xp_var*gsl_ran_gaussian_pdf(xbeta-mean_p, sqrt_var1)+(1-xp_var)*gsl_ran_gaussian_pdf(xbeta-mean_p, sqrt_var2));
        log2+=log(gsl_ran_gaussian_pdf(xbeta-xbeta1, xsig_beta1));
        log1+=log(gsl_ran_gaussian_pdf(betap[0]-xbeta1, xsig_beta1));
        //std::cout <<"log1 = "<< log1 <<" log2 = "<< log2 << std::endl;
        if (log(Rcpp::as<double>(Rcpp::runif(1)) ) <= (log1 - log2)) {
              beta[tt] = betap[0];
              xAbeta[tt] = 1;
        } 
    }
}
