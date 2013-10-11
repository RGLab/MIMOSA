#include <Rcpp.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <R_ext/Utils.h>
#include <stdlib.h>
using namespace std;
using namespace Rcpp;


void updatemuu(vector<double>& xmust, vector<double>& xmuut, vector<int>& xn_s, vector<int>& xn_u, int xI, int xK, int xM, int K1,
               vector<double>& xp_var, Rcpp::NumericVector& sqrt_var1, Rcpp::NumericVector& sqrt_var2, int xtt, vector<int>& xgammat, Rcpp::IntegerMatrix& xd, 
                Rcpp::NumericMatrix& xybar_s,  Rcpp::NumericMatrix& xybar_u,  Rcpp::NumericMatrix& pt2,  Rcpp::NumericMatrix& pt3,
                double xlambda,  double xbeta, double xalpha, vector<double>& xmu, vector<double>& xSigu, vector<int>& xAmuu)
{
    int temp=0;
     double beta_np=0.; double alpha_np=0.;
    double delF=0.; double log1=0.; double log2=0.;  int nik=0;  
    
for (int p=0; p<xM; p++) {
    delF=0.; log1=0.; log2=0.;  nik=0;  
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
                           delF+=(-alpha_np/beta_np)*(xmuut[p]-xybar_u(temp,p))/(xlambda+1/xn_u[temp]);    
                     }else if (xn_u[temp]>0 && xn_s[temp]==0) {
                           beta_np = pt2(temp,p)+xbeta+pow((xmuut[p]-xybar_u(temp,p)),2)/(2*(xlambda+1/xn_u[temp]));
                           alpha_np = nik/2+xalpha;
                           log2 += -alpha_np*log(beta_np);    
                           delF+=(-alpha_np/beta_np)*(xmuut[p]-xybar_u(temp,p))/(xlambda+1/xn_u[temp]);    
                     }                  
                }else{
                     if (nik>0) {
                            beta_np =pt3(temp,p)+xbeta+pow((xmuut[p]-(xn_s[temp]*xybar_s(temp,p)+xn_u[temp]*xybar_u(temp,p))/nik),2)/(2*(xlambda+1/nik))+
                                             0.5*xn_s[temp]*pow(xybar_s(temp,p),2)+ 0.5*xn_u[temp]*pow(xybar_u(temp,p),2)-0.5*pow((xn_s[temp]*xybar_s(temp,p)+xn_u[temp]*xybar_u(temp,p)),2)/nik;
                            alpha_np = nik/2+xalpha;
                            log2+= -alpha_np*log(beta_np); 
                            delF +=(-alpha_np/beta_np)*(xmuut[p]-(xn_s[temp]*xybar_s(temp,p)+xn_u[temp]*xybar_u(temp,p))/nik)/(xlambda+1/nik);
                     }       
                }
             }

         }
    }
    double mean_p = std::max(0.01, xmuut[p]+delF/xtt);
    Rcpp::NumericVector muu_p= Rcpp::rnorm(1, mean_p, sqrt_var1[p]);
    if (Rcpp::as<double>(Rcpp::rbinom(1,1,xp_var[p])) == 1) {
         muu_p= Rcpp::rnorm(1, mean_p, sqrt_var1[p]);}
    else {muu_p = Rcpp::rnorm(1, mean_p, sqrt_var2[p]);}
    
    if(muu_p[0]>0.0) {
         log2+=log(xp_var[p]*gsl_ran_gaussian_pdf(muu_p[0]-mean_p, sqrt_var1[p])+(1-xp_var[p])*gsl_ran_gaussian_pdf(muu_p[0]-mean_p, sqrt_var2[p])); 
         delF = 0.;
         for (int i=0; i<xI;i++) {
            for(int k=0; k<K1; k++) {
              temp=xI*k+i;
                if(xd(k,p)==1) {
                    nik = xn_u[temp]+xn_s[temp];
                    if(xgammat[temp]==1) {
                       if (xn_u[temp] >0 && xn_s[temp]>0) {
                           beta_np = pt3(temp,p)+xbeta+pow((muu_p[0]-xybar_u(temp,p)),2)/(2*(xlambda+1/xn_u[temp]))+
                                      pow((xmust[p]-xybar_s(temp,p)),2)/(2*(xlambda+1/xn_s[temp]));
                           alpha_np = nik/2+xalpha;
                           log1+=-alpha_np*log(beta_np);
                           delF+=(-alpha_np/beta_np)*(muu_p[0]-xybar_u(temp,p))/(xlambda+1/xn_u[temp]);    
                        }else if (xn_u[temp]>0 && xn_s[temp]==0) {
                          beta_np = pt2(temp,p)+xbeta+pow((muu_p[0]-xybar_u(temp,p)),2)/(2*(xlambda+1/xn_u[temp]));
                          alpha_np = nik/2+xalpha;
                           log1 += -alpha_np*log(beta_np);    
                           delF+=(-alpha_np/beta_np)*(muu_p[0]-xybar_u(temp,p))/(xlambda+1/xn_u[temp]);    
                        }                  
                    }else{
                      if (nik>0) {
                            beta_np = pt3(temp,p)+xbeta+pow((muu_p[0]-(xn_s[temp]*xybar_s(temp,p)+xn_u[temp]*xybar_u(temp,p))/nik),2)
                                    /(2*(xlambda+1/nik))+0.5*xn_s[temp]*pow(xybar_s(temp,p),2)+ 0.5*xn_u[temp]*pow(xybar_u(temp,p),2)-0.5*pow((xn_s[temp]*xybar_s(temp,p)+xn_u[temp]*xybar_u(temp,p)),2)/nik;
                            alpha_np = nik/2+xalpha;
                            log1+= -alpha_np*log(beta_np); 
                            delF +=(-alpha_np/beta_np)*(muu_p[0]-(xn_s[temp]*xybar_s(temp,p)+xn_u[temp]*xybar_u(temp,p))/nik)/(xlambda+1/nik);
                       }       
                   }
              }
           }
        }
        mean_p = std::max(0.01, muu_p[0]+delF/xtt);
        log1+=log(xp_var[p]*gsl_ran_gaussian_pdf(xmuut[p]-mean_p, sqrt_var1[p])+(1-xp_var[p])*gsl_ran_gaussian_pdf(xmuut[p]-mean_p, sqrt_var2[p]));
        log2+=log(gsl_ran_gaussian_pdf(xmuut[p]-xmu[p], xSigu[p]));
        log1+=log(gsl_ran_gaussian_pdf(muu_p[0]-xmu[p], xSigu[p]));
       //std::cout <<"log1= "<< log1 << " log2 = "<<log2 <<" p="<<p<<std::endl;

          if (log(Rcpp::as<double>(Rcpp::runif(1)) ) <= (log1 - log2)) {
               xmuut[p] = muu_p[0];
               xAmuu[p] = 1;
          } 
    }
    }
}
