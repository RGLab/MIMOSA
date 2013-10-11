#include <Rcpp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdlib.h>
#include "model.h"
#include "updatealphau_RW.h"
#include "updategamma_indi.h"
#include "updategamma_indi_1.h"
#include "updatealphas_RW.h"
#include "updatemus.h"
#include "updatemuu.h"
#include "updatealpha_RW.h"
#include "updatebeta_RW.h"
#include <R_ext/Utils.h>

using namespace Rcpp;
using namespace std;

RcppExport SEXP model(SEXP T, SEXP I, SEXP K, SEXP M, SEXP ttt, SEXP SS, SEXP alpha_u, SEXP alpha_s, SEXP mu_u, SEXP mu_s, SEXP alpha, 
                      SEXP beta, SEXP gamma, SEXP n_s, SEXP n_u, SEXP varp_u, SEXP lambda_u,SEXP indi, SEXP d, SEXP ybar_s, SEXP ybar_u,
                      SEXP ys2_s, SEXP ys2_u, SEXP a, SEXP b, SEXP lambda, SEXP mk, SEXP Istar, SEXP mKstar, SEXP pp, SEXP pb1,
                      SEXP pb2, SEXP lambda_s,SEXP var_1, SEXP var_2, SEXP p_var, SEXP p_vars, SEXP var_1s, SEXP var_2s, SEXP m_s,SEXP Sigma_s,
                      SEXP p_varu, SEXP var_1u,SEXP var_2u, SEXP m_u, SEXP Sigma_u, SEXP p_vara, SEXP var_1a, SEXP var_2a, SEXP sig_alpha1, SEXP alpha1,
                      SEXP p_varb, SEXP var_1b, SEXP var_2b, SEXP sig_beta1, SEXP beta1, SEXP A_alphau,  SEXP A_alphas, SEXP A_gm, SEXP A_mus, SEXP A_muu, 
                      SEXP A_alpha, SEXP A_beta, SEXP Tune, SEXP pgamma)
{
  BEGIN_RCPP
  int xT = Rcpp::as<int>(T); 
  int xI = Rcpp::as<int>(I); 
  int xK = Rcpp::as<int>(K); 
  int xtt = Rcpp::as<int>(ttt);
  int xM = Rcpp::as<int>(M); 
  int K1 = xK-1; 
  int xSS = Rcpp::as<int>(SS);
  int xIK = xI*xK;
  int xTune = Rcpp::as<int>(Tune);
  
  vector<int> xn_s = Rcpp::as<vector<int> >(n_s); 
  vector<int> xn_u = Rcpp::as<vector<int> >(n_u); 
  Rcpp::IntegerMatrix xd(d);
  Rcpp::IntegerVector xmk(mk);
  int xIstar = Rcpp::as<int>(Istar);
  int xmKstar = Rcpp::as<int>(mKstar);    
  double xpgamma = Rcpp::as<double>(pgamma);
  
  Rcpp::NumericMatrix xybar_s(ybar_s);
  Rcpp::NumericMatrix xybar_u(ybar_u);
  Rcpp::NumericMatrix xys2_s(ys2_s);
  Rcpp::NumericMatrix xys2_u(ys2_u);
  Rcpp::IntegerMatrix xindicator(indi);
  
  // mcmc tuning parameters
  Rcpp::NumericVector sqrt_var(varp_u); //for alpha_u 
  vector<double> xlambda_u =Rcpp::as<vector<double> >(lambda_u);
  double xa = Rcpp::as<double>(a); //for gamma
  double xb = Rcpp::as<double>(b);
  double ab = xa+xb; double abI = ab+xI; double bI = xI+xb; 
  double xlambda = Rcpp::as<double>(lambda); 
  Rcpp::NumericVector xpp(pp); 
  double xpb1 = Rcpp::as<double>(pb1);
  double xpb2 = Rcpp::as<double>(pb2);
  vector<double> xlambda_s = Rcpp::as<vector<double> >(lambda_s);
  Rcpp::NumericVector sqrt_var1(var_1); //for alpha_s
  Rcpp::NumericVector sqrt_var2(var_2);
  vector<double> xp_var = Rcpp::as<vector<double> >(p_var); 
  Rcpp::NumericVector sqrt_var1s(var_1s); // for mu_s
  Rcpp::NumericVector sqrt_var2s(var_2s); 
  vector<double> xp_vars = Rcpp::as<vector<double> >(p_vars);
  vector<double> xms = Rcpp::as<vector<double> >(m_s); 
  vector<double> xSigs = Rcpp::as<vector<double> >(Sigma_s);
  Rcpp::NumericVector sqrt_var1u(var_1u); // for mu_u
  Rcpp::NumericVector sqrt_var2u(var_2u);
  vector<double> xp_varu = Rcpp::as<vector<double> >(p_varu);
  vector<double> xmu = Rcpp::as<vector<double> >(m_u); 
  vector<double> xSigu = Rcpp::as<vector<double> >(Sigma_u);
  double xalpha1 = Rcpp::as<double>(alpha1); //for alpha
  double xsig_alpha1 = Rcpp::as<double>(sig_alpha1);
  double sqrt_var1a = Rcpp::as<double>(var_1a);
  double sqrt_var2a = Rcpp::as<double>(var_2a); double xp_vara = Rcpp::as<double>(p_vara);
  double sqrt_var1b = Rcpp::as<double>(var_1b); // for beta
  double sqrt_var2b = Rcpp::as<double>(var_2b);
  double xp_varb = Rcpp::as<double>(p_varb);
  double xbeta1 = Rcpp::as<double>(beta1);
  double xsig_beta1 = Rcpp::as<double>(sig_beta1);
  

  // Initialize model parameters to store mcmc traces
  Rcpp::NumericMatrix alphaut(alpha_u); 
  Rcpp::IntegerMatrix Aalphau(A_alphau);
  Rcpp::NumericMatrix alphast(alpha_s);
  Rcpp::IntegerMatrix Aalphas(A_alphas);
  Rcpp::NumericMatrix muut(mu_u);  
  Rcpp::IntegerMatrix Amuu(A_muu);
  Rcpp::NumericMatrix must(mu_s); 
  Rcpp::IntegerMatrix Amus(A_mus);
  Rcpp::NumericVector alphat(alpha); 
  Rcpp::IntegerVector Aalpha(A_alpha);
  Rcpp::NumericVector betat(beta); 
  Rcpp::IntegerVector Abeta(A_beta);
  Rcpp::IntegerMatrix gammat(gamma);
  Rcpp::IntegerMatrix Agamma(A_gm);
  
  
  vector<double> alphau_t1(xK);
  vector<double> alphas_t1(xK);
  vector<int> gamma_t1(xIK);
  vector<int> Agm_t1(xI);
  vector<int> Aalp_t1(xK);
  vector<int> Aalps_t1(xK);
  vector<int> Amus_t1(xM); vector<int> Amuu_t1(xM);
  vector<double> mus_t1(xM); vector<double> muu_t1(xM); 

  Rcpp::NumericMatrix pt1(xIK,xM);
  Rcpp::NumericMatrix pt2(xIK,xM);
  Rcpp::NumericMatrix pt3(xIK,xM);
  for ( int p = 0; p<xM; p++) {
      for (int ik=0; ik<xIK; ik++){
        pt1(ik,p) = 0.5*xn_s[ik]*xys2_s(ik,p);
        pt2(ik,p) = 0.5*xn_u[ik]*xys2_u(ik,p);
        pt3(ik,p) = pt1(ik,p) + pt2(ik,p);
      }
  }
  
  
  int t1 = 0;
  for ( int k=0; k<xK; k++){
      alphau_t1[k] = alphaut(t1,k); 
      alphas_t1[k] = alphast(t1,k);
  }
  for (int m=0; m<xM; m++){
      mus_t1[m] = must(t1,m);
      muu_t1[m] = muut(t1,m);
   }
  for ( int ik=0; ik<(xIK); ik++) {gamma_t1[ik] = gammat(ik,t1);}
  
  Rcpp::RNGScope scope;

  for (int tt=1; tt<xT; tt++){
    t1=tt-1; 
      
     // update alpha_u
    updatealphau_RW(alphau_t1, xn_s, xn_u, xI, xK, xlambda_u, sqrt_var, xtt, gamma_t1, Aalp_t1);
    for (int k=0; k<xK; k++){
      alphaut(tt,k) = alphau_t1[k]; //alphau_t1 is updated
      Aalphau(k,tt) = Aalp_t1[k];
      Aalp_t1[k] = 0;
    }
     
    // update gamma
    //if (Rcpp::as<double>(Rcpp::runif(1)) <=xpgamma){
      updategamma_indi(xn_s,xn_u, gamma_t1, xI,  xK,  xM,  xSS, K1, alphau_t1, alphas_t1, xa, xb,  xmk,  xIstar, xmKstar, 
                     xpp,  xpb1,  xpb2, xlambda, betat[t1], alphat[t1], mus_t1, muu_t1, xd,xybar_s, xybar_u, 
                     pt1,pt2,pt3, xindicator,  Agm_t1,ab, abI, bI);
    //}else{
    //  updategamma_indi_1(xn_s,xn_u, gamma_t1, xI,  xK,  xM, K1, alphau_t1, alphas_t1, xa, xb,  xmk,  xIstar, xmKstar, 
    //                 xpp,  xpb1,  xpb2, xlambda, betat[t1], alphat[t1], mus_t1, muu_t1, xd,xybar_s, xybar_u, 
    //                 pt1,pt2,pt3, xindicator,  Agm_t1,ab, abI, bI);
    //}

    for (int ik=0; ik<(xI*xK); ik++) {
      gammat(ik,tt) = gamma_t1[ik]; //gamma_t1 is updated
    }
    for (int i=0; i<xI; i++){Agamma(i,tt)=Agm_t1[i]; Agm_t1[i]=0;}
    
   // update alpha_s
    updatealphas_RW(alphas_t1, xn_s, xK, xI,  xlambda_s, gamma_t1, sqrt_var1, sqrt_var2, xp_var, xtt, Aalps_t1);
    for (int k=0; k<xK; k++){
      alphast(tt,k) = alphas_t1[k]; //alphas_t1 is updated
      Aalphas(k,tt) = Aalps_t1[k];
      Aalps_t1[k] = 0;
    }
    
    //update mu_s
    updatemus(mus_t1, muu_t1, xn_s, xn_u, xI, xK,  xM, K1, xp_vars, sqrt_var1s, sqrt_var2s, xtt, gamma_t1, xd, xybar_s, 
              xybar_u, pt1, pt3, xlambda, betat[t1], alphat[t1], xms,xSigs , Amus_t1);
   for (int m=0; m<xM; m++){
     must(tt,m) = mus_t1[m];
     Amus(m,tt) = Amus_t1[m];
     Amus_t1[m]=0;
   }
   
   // update mu_u
   updatemuu(mus_t1,muu_t1, xn_s, xn_u, xI, xK, xM, K1, xp_varu, sqrt_var1u, sqrt_var2u, xtt, gamma_t1, xd, 
             xybar_s, xybar_u, pt2, pt3, xlambda, betat[t1], alphat[t1], xmu, xSigu, Amuu_t1);
   for (int m=0; m<xM; m++){
     muut(tt,m) = muu_t1[m];
     Amuu(m,tt) = Amuu_t1[m];
     Amuu_t1[m] = 0;
   }
    // update alpha 
   updatealpha_RW(mus_t1, muu_t1, xn_s, xn_u, xI, xK, xM,  K1,  t1, xp_vara, sqrt_var1a, sqrt_var2a, xtt, gamma_t1,
                      xd, xybar_s, xybar_u, pt1, pt2, pt3, xlambda, betat[t1], alphat, xsig_alpha1, xalpha1, Aalpha);
                      
   // update beta
   updatebeta_RW(mus_t1, muu_t1, xn_s, xn_u, xI, xK, xM, K1, t1, xp_varb, sqrt_var1b, sqrt_var2b, xtt, gamma_t1, 
                    xd, xybar_s, xybar_u, pt1, pt2, pt3, xlambda, betat, alphat[(t1+1)], xsig_beta1, xbeta1, Abeta);
                
   if (xTune==1){
     // parameter tuning
    if ((tt+1)>4000 & ((tt+1) % 4000)==0) {
      int sec1 = tt-4000+1; int sec2 = tt;
      for (int kk = 0; kk<xK; kk++){
        double Au =0.; double As = 0.;
        for (int e = sec1; e<=sec2; e++){
          Au += Aalphau(kk,e);
          As += Aalphas(kk,e);
        }
        Au = Au/4000; As = As/4000;
        if (Au < 0.001) {sqrt_var[kk] = sqrt_var[kk]*sqrt(0.1);}
        else if (Au <0.05) {sqrt_var[kk] = sqrt_var[kk]*sqrt(0.5);}
        else if (Au < 0.2) {sqrt_var[kk] = sqrt_var[kk]*sqrt(0.9);}
        else if (Au>0.5) {sqrt_var[kk] = sqrt_var[kk]*sqrt(1.1);}
        else if (Au>0.75) {sqrt_var[kk] = sqrt_var[kk]*sqrt(2);}
        else if (Au>0.95) {sqrt_var[kk] = sqrt_var[kk]*sqrt(10);}
        if (As<0.001) {sqrt_var1[kk] = sqrt_var1[kk]*sqrt(0.1);}
        else if (As <0.05) {sqrt_var1[kk] = sqrt_var1[kk]*sqrt(0.5);}
        else if (As < 0.2) {sqrt_var1[kk] = sqrt_var1[kk]*sqrt(0.9);}
        else if (As>0.5) {sqrt_var2[kk] = sqrt_var2[kk]*sqrt(1.1);}
        else if (As>0.75) {sqrt_var2[kk] = sqrt_var2[kk]*sqrt(2);}
        else if (As>0.95) {sqrt_var2[kk] = sqrt_var2[kk]*sqrt(10);}          
      }
      for ( int ii=0; ii<xI; ii++) {  
        double Ag =0.;
        for (int e = sec1; e<=sec2; e++){
          Ag += Agamma(ii,e);
        }
        Ag = Ag/4000;
        if (Ag<0.001) {xpp[ii] = min(0.9,xpp[ii]*1.4);}
        else if (Ag<0.05) {xpp[ii] = min(0.9,xpp[ii]*1.2);}
        else if (Ag<0.2) {xpp[ii] = min(0.9,xpp[ii]*1.1);}
        else if (Ag>0.6) {xpp[ii] = max(0.1,xpp[ii]*0.8);}
        else if (Ag>0.75) {xpp[ii] = max(0.1,xpp[ii]*0.5);}
        else if (Ag>0.95) {xpp[ii] = max(0.1,xpp[ii]*0.2);}
      }
      
      for ( int p=0; p<xM; p++){
        double Ams = 0.; double Amu=0.;
        for (int e=sec1; e<=sec2; e++){
          Ams += Amus(p,e);
          Amu += Amuu(p,e);
        }
        Ams=Ams/4000; Amu=Amu/4000;
        if (Ams<0.001) {sqrt_var1s[p] = sqrt_var1s[p]*sqrt(0.1);}
        else if (Ams <0.05) {sqrt_var1s[p] = sqrt_var1s[p]*sqrt(0.5);}
        else if (Ams < 0.2) {sqrt_var1s[p] = sqrt_var1s[p]*sqrt(0.9);}
        else if (Ams>0.5) {sqrt_var2s[p] = sqrt_var2s[p]*sqrt(1.1);}
        else if (Ams>0.75) {sqrt_var2s[p] = sqrt_var2s[p]*sqrt(2);}
        else if (Ams>0.95) {sqrt_var2s[p] = sqrt_var2s[p]*sqrt(10);}
        if (Amu<0.001) {sqrt_var1u[p] = sqrt_var1u[p]*sqrt(0.1);}
        else if (Amu <0.05) {sqrt_var1u[p] = sqrt_var1u[p]*sqrt(0.5);}
        else if (Amu < 0.2) {sqrt_var1u[p] =sqrt_var1u[p]*sqrt(0.9);}
        else if (Amu>0.5) {sqrt_var2u[p] = sqrt_var2u[p]*sqrt(1.1);}
        else if (Amu>0.75) {sqrt_var2u[p] =sqrt_var2u[p]*sqrt(2);}
        else if (Amu>0.95) {sqrt_var2u[p] =sqrt_var2u[p]*sqrt(10);}
      }
      
      double Aal = 0.; double Abe = 0.;
      for (int e=sec1; e<=sec2; e++){
        Aal += Aalpha[e];
        Abe += Abeta[e];
      }
      Aal = Aal/4000; Abe = Abe/4000;
      if (Aal<0.001) {sqrt_var1a = sqrt_var1a*sqrt(0.1);}
      else if (Aal <0.05) {sqrt_var1a = sqrt_var1a*sqrt(0.5);}
      else if (Aal< 0.2) {sqrt_var1a =sqrt_var1a*sqrt(0.9);}
      else if (Aal>0.5) {sqrt_var2a = sqrt_var2a*sqrt(1.1);}
      else if (Aal>0.75) {sqrt_var2a=sqrt_var2a*sqrt(2);}
      else if (Aal>0.95) {sqrt_var2a =sqrt_var2a*sqrt(10);}
      if (Abe<0.001) {sqrt_var1b = sqrt_var1b*sqrt(0.1);}
      else if (Abe <0.05) {sqrt_var1b= sqrt_var1b*sqrt(0.5);}
      else if (Abe< 0.2) {sqrt_var1b =sqrt_var1b*sqrt(0.9);}
      else if (Abe>0.5) {sqrt_var2b = sqrt_var2b*sqrt(1.1);}
      else if (Abe>0.75) {sqrt_var2b=sqrt_var2b*sqrt(2);}
      else if (Abe>0.95) {sqrt_var2b=sqrt_var2b*sqrt(10);}  
     }
   }
   
  }
  
  return Rcpp::List::create(Rcpp::Named("mk") = xmk, Rcpp::Named("Istar") = xIstar, Rcpp::Named("mKstar") = xmKstar, Rcpp::Named("var_1a") = sqrt_var1a,
                            Rcpp::Named("var_2a") = sqrt_var2a, Rcpp::Named("var_1b") = sqrt_var1b, Rcpp::Named("var_2b") = sqrt_var2b);
END_RCPP
  
}
  
  
  
  
