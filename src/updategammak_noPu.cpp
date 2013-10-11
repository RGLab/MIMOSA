#include <Rcpp.h>

#include <R_ext/Utils.h>
#include <stdlib.h>

RcppExport SEXP updategammak_noPu(SEXP n_s,SEXP n_u, SEXP gammat, SEXP I, SEXP K, SEXP SS, SEXP alphau, SEXP alphas, SEXP alpha, SEXP mk, SEXP Istar, 
                                  SEXP mKstar, SEXP pp, SEXP pb1, SEXP pb2, SEXP indi)
{
    BEGIN_RCPP

    Rcpp::IntegerMatrix xn_s(n_s);
    Rcpp::IntegerMatrix xn_u(n_u);
    Rcpp::NumericVector xalphas(alphas); 
    Rcpp::IntegerMatrix xindicator(indi); 
    Rcpp::NumericVector xalphau(alphau);
    double xalpha = Rcpp::as<double>(alpha);
    int xK = Rcpp::as<int>(K);
    int xI = Rcpp::as<int>(I); 
    int xSS = Rcpp::as<int>(SS); 
    Rcpp::IntegerVector xAg(xI);
    Rcpp::IntegerVector gamma_old(xK);
    Rcpp::IntegerMatrix gamma_tt(gammat);
    Rcpp::IntegerVector xmk(mk);
    int xIstar = Rcpp::as<int>(Istar);
    int xmKstar = Rcpp::as<int>(mKstar);    
    Rcpp::NumericVector xpp(pp); 
    double xpb1 = Rcpp::as<double>(pb1);
    double xpb2 = Rcpp::as<double>(pb2);
    Rcpp::RNGScope scope;
    Rcpp::NumericVector m(1); 
    int K1 = xK-1;
    int array1K[K1];
    int mik[K1];
    double alphaK = xalpha/xK;
    int xindiK1 = 0;
    
    int miKstar = 0;
    int Iistar = 0;
    int flag =0; int flag_old = 0;
    double log_likeli_deno = 0.; double log_likeli_nume = 0.0;
    int ct=0;
    double log_probK = 0.;
    double log_1probK = 0.;
    int vec = 0; int arr = 0;
    int lplace0_new = 0; int lplace0 = 0;
    int flag_prop = 0;
    int lplace1_new = 0;
    int lplace1 = 0;
    for (int i = 0; i < xI; i++) {
        xindiK1 = xindicator(i,K1);
        int array1K[xindiK1];
        for (int j = 0; j < xK; j++) {
            gamma_old[j] = gamma_tt(i,j);
            
        }
        for (int ss = 0; ss < xSS; ss++) {
          if(xindiK1>0) {
            miKstar = xmKstar;
            Iistar = xIstar;
            flag =0; flag_old = 0;
            log_likeli_deno = 0.;
            ct=0;
            for (int j = 0; j < K1; j++) {
                gamma_tt(i,j) = gamma_old[j];
                mik[j] = xmk[j]-gamma_old[j];
                flag +=gamma_tt(i,j);
                if (xindicator(i,j) ==1) {array1K[ct] = j; ct +=1;}
                log_likeli_deno+=gamma_old[j]*(log(mik[j]+alphaK)-log(xI+alphaK))+(1-gamma_old[j])*(log(xI-mik[j])-log(xI+alphaK));
            }
            if (flag >1 ) {//*
                Iistar = xIstar-1; 
                miKstar = xmKstar-gamma_old[K1];
                flag_old = 1;
            }
            gamma_tt(i,K1) = gamma_old[K1];
            log_probK =log(miKstar+alphaK)-log(Iistar+alphaK+1); 
            log_1probK = log(Iistar-miKstar+1)-log(Iistar+alphaK+1);
            if(flag_old == 1) {
               log_likeli_deno+=gamma_old[K1]*(log_probK)+(1-gamma_old[K1])*(log_1probK);
            }
            if(xindiK1<=2) {m=1;}
            else{
              if (Rcpp::as<double>(Rcpp::runif(1)) <=xpp[i]) 
                 { m = Rcpp::rbinom(1,xindiK1,xpb1); }//m[0] +=1;}
              else {m =Rcpp::rbinom(1,xindiK1,xpb2); }//m[0] +=1;}
            } 
            /* random sample */
            if (m[0]>0) {
             for (int k = 0; k<m[0]; k++) {             
               Rcpp::NumericVector tmp = Rcpp::runif(1,k,xindiK1);
               vec = (int) tmp[0];
               arr = array1K[k];
               array1K[k] = array1K[vec];
               array1K[vec] = arr;
               gamma_tt(i,array1K[k]) = 1-gamma_tt(i,array1K[k]);
               //if (Rcpp::as<double>(Rcpp::runif(1)) < 0.5) {gamma_tt(i,array1K[k]) = 1;}
               //else {gamma_tt(i,array1K[k]) = 0;}
              }
            }
            lplace0_new = 0; lplace0 = 0;
            log_likeli_nume = 0.0;
            for (int k = 0; k < K1; k++) {
               log_likeli_nume += gamma_tt(i,k)*(log(mik[k]+alphaK)-log(xI+alphaK))+(1-gamma_tt(i,k))*(log(xI-mik[k])-log(xI+alphaK));
               if (gamma_tt(i,k) == 0) { lplace0_new +=1;}
               if (gamma_old[k] == 0){ lplace0 += 1;}            
            }
            lplace1_new = K1-lplace0_new;
            lplace1 = K1-lplace0;
            flag_prop = 0;
            if (lplace1_new == 1) { //*
               gamma_tt(i,K1) = 1; 
               lplace1_new +=1;
             } else if (lplace1_new == 0) {
               gamma_tt(i,K1) = 0;
               lplace0_new += 1;
            } else {
               flag_prop = 1;
               log_likeli_deno += log(0.5); // proposal
               if (Rcpp::as<double>(Rcpp::runif(1)) <= 0.5) {
                    gamma_tt(i,K1) = 1;
               }else{
                    gamma_tt(i,K1) = 0;
               }
               if (gamma_tt(i,K1) == 1) {
                   lplace1_new +=1;
                   log_likeli_nume +=log_probK; 
               } else {
                   lplace0_new += 1;
                   log_likeli_nume += log_1probK; 
               }
            
            } 
             if ( lplace1 == 1) {   //*          
               lplace1 +=1;
            } else if (lplace1 == 0) {
               lplace0+= 1;
            } else {
               log_likeli_nume += log(0.5); // proposal 
               if (gamma_old[K1] == 1) {
                   lplace1 +=1;
               } else {
                   lplace0 += 1;
               }
       
            } 
                    
            int place0_new[lplace0_new];  
            int place0[lplace0];   
            int place1[lplace1];   
            int flag0 = 0; int flag1 = 0;
            int flag0_old = 0; int flag1_old = 0;
            int place1_new[lplace1_new];
            for (int k= 0; k < xK; k++) {
               if (gamma_tt(i,k) == 0) {
                  place0_new[flag0] = k;
                  flag0 += 1;
               } else {
                  place1_new[flag1] = k;
                  flag1 += 1;
               }   
               if (gamma_old[k] == 0) {
                  place0[flag0_old] = k;
                  flag0_old += 1;
               } else{
                  place1[flag1_old] = k;
                  flag1_old +=1;
               }          
            }
            
            if (lplace1_new == 0) {
                double sum_alphau_nus = 0.0;
                double sum_gammaln_nu = 0.0;
               
                for ( int k = 0; k<xK; k++) {
                    double summ = xalphau[k]+xn_u(i,k)+xn_s(i,k);
                    sum_alphau_nus += summ; 
                    sum_gammaln_nu +=lgamma(summ); 
                }
                log_likeli_nume +=(sum_gammaln_nu-lgamma(sum_alphau_nus));

            } else if (lplace0_new == 0) {
                double sum_alphau_nus = 0.0;
                double sum_gammaln_nu = 0.0;
                double sum_gammaln = 0.0;
                double sum_alphas = 0.0;
                double sum_gammaln_ns = 0.0;
                double sum_alphas_ns = 0.0;

                for ( int k=0; k<xK; k++) {
                    double summ = xalphau[k]+xn_u(i,k);
                    sum_alphau_nus += summ; 
                    sum_gammaln_nu +=lgamma(summ); 

                    sum_gammaln +=lgamma(xalphas[k]);
                    sum_alphas +=xalphas[k];
                    sum_gammaln_ns += (lgamma(xalphas[k] + xn_s(i,k)));
                    sum_alphas_ns += (xalphas[k]+xn_s(i,k));
                }
                log_likeli_nume +=(sum_gammaln_nu-lgamma(sum_alphau_nus));
                log_likeli_nume -= (sum_gammaln - lgamma(sum_alphas));
                log_likeli_nume +=(sum_gammaln_ns -lgamma(sum_alphas_ns));
            } else {
                double sum_gammaln = 0.0;
                double sum_alphas = 0.0;
                double sum_gammaln_ns = 0.0;
                double sum_alphas_ns = 0.0;
                double sum_alphau_nu = 0.0;
                double sum_alphau_nus = 0.0;
                double sum_gammaln_nu = 0.0;
                for (int k = 0; k<lplace1_new; k++) {
                  sum_alphau_nu +=(xalphau[place1_new[k]]+xn_u(i,place1_new[k]));
                  sum_alphau_nus+=(xalphau[place1_new[k]]+xn_u(i,place1_new[k])+xn_s(i,place1_new[k]));
                  sum_gammaln_nu +=lgamma(xalphau[place1_new[k]]+xn_u(i,place1_new[k]));
                  sum_gammaln +=lgamma(xalphas[place1_new[k]]);
                  sum_alphas +=xalphas[place1_new[k]];
                  sum_gammaln_ns += (lgamma(xalphas[place1_new[k]] + xn_s(i,place1_new[k])));
                  sum_alphas_ns += (xalphas[place1_new[k]]+xn_s(i,place1_new[k]));
                }
                log_likeli_nume -= (sum_gammaln - lgamma(sum_alphas));
                log_likeli_nume +=(sum_gammaln_ns -lgamma(sum_alphas_ns));
                log_likeli_nume += (sum_gammaln_nu-lgamma(sum_alphau_nu));

                log_likeli_nume +=lgamma(sum_alphau_nus+1);
                sum_gammaln_nu = 0.0;
                for (int k = 0; k<lplace0_new; k++) {
                    sum_alphau_nus+=(xalphau[place0_new[k]]+xn_u(i,place0_new[k])+xn_s(i,place0_new[k]));
                    sum_gammaln_nu +=lgamma(xalphau[place0_new[k]]+xn_u(i,place0_new[k])+xn_s(i,place0_new[k]));
                } 
                log_likeli_nume +=sum_gammaln_nu;   
                log_likeli_nume -=lgamma(sum_alphau_nus+1);
            }
   
        //
            if (lplace1 == 0) {
                double sum_alphau_nus = 0.0;
                double sum_gammaln_nu = 0.0;
                for ( int k = 0; k<xK; k ++) {
                    double summ = xalphau[k]+xn_u(i,k)+xn_s(i,k);
                    sum_alphau_nus += summ; 
                    sum_gammaln_nu +=lgamma(summ); 
                }
                log_likeli_deno +=(sum_gammaln_nu-lgamma(sum_alphau_nus));
            } else if (lplace0 == 0) {
                double sum_alphau_nus = 0.0;
                double sum_gammaln_nu = 0.0;
                double sum_gammaln = 0.0;
                double sum_alphas = 0.0;
                double sum_gammaln_ns = 0.0;
                double sum_alphas_ns = 0.0;

                for ( int k=0; k<xK; k++) {
                    double summ = xalphau[k]+xn_u(i,k);
                    sum_alphau_nus += summ; 
                    sum_gammaln_nu +=lgamma(summ); 

                    sum_gammaln +=lgamma(xalphas[k]);
                    sum_alphas +=xalphas[k];
                    sum_gammaln_ns += (lgamma(xalphas[k] + xn_s(i,k)));
                    sum_alphas_ns += (xalphas[k]+xn_s(i,k));
                }
                log_likeli_deno +=(sum_gammaln_nu-lgamma(sum_alphau_nus));
                log_likeli_deno -= (sum_gammaln - lgamma(sum_alphas));
                log_likeli_deno +=(sum_gammaln_ns -lgamma(sum_alphas_ns));
            } else {
                double sum_gammaln = 0.0;
                double sum_alphas = 0.0;
                double sum_gammaln_ns = 0.0;
                double sum_alphas_ns = 0.0;
                double sum_alphau_nu = 0.0;
                double sum_alphau_nus = 0.0;
                double sum_gammaln_nu = 0.0;
                for (int k = 0; k<lplace1; k++) {
                  sum_alphau_nu +=(xalphau[place1[k]]+xn_u(i,place1[k]));
                  sum_alphau_nus+=(xalphau[place1[k]]+xn_u(i,place1[k])+xn_s(i,place1[k]));
                  sum_gammaln_nu +=lgamma(xalphau[place1[k]]+xn_u(i,place1[k]));
                  sum_gammaln +=lgamma(xalphas[place1[k]]);
                  sum_alphas +=xalphas[place1[k]];
                  sum_gammaln_ns += (lgamma(xalphas[place1[k]] + xn_s(i,place1[k])));
                  sum_alphas_ns += (xalphas[place1[k]]+xn_s(i,place1[k]));
                }
                log_likeli_deno -= (sum_gammaln - lgamma(sum_alphas));
                log_likeli_deno +=(sum_gammaln_ns -lgamma(sum_alphas_ns));
                log_likeli_deno += (sum_gammaln_nu-lgamma(sum_alphau_nu));

                log_likeli_deno +=lgamma(sum_alphau_nus+1);
                sum_gammaln_nu = 0.0;
                for (int k = 0; k<lplace0; k++) {
                    sum_alphau_nus+=(xalphau[place0[k]]+xn_u(i,place0[k])+xn_s(i,place0[k]));
                    sum_gammaln_nu +=lgamma(xalphau[place0[k]]+xn_u(i,place0[k])+xn_s(i,place0[k]));
                } 
                log_likeli_deno +=sum_gammaln_nu;   
                log_likeli_deno -=lgamma(sum_alphau_nus+1);
            }
      
           //std::cout <<"log_nume = "<< log_likeli_nume << std::endl;
           //std::cout <<"log_deo = "<< log_likeli_deno << std::endl;
            if (log(Rcpp::as<double>(Rcpp::runif(1)) ) <= (log_likeli_nume - log_likeli_deno)) {
               for (int j = 0; j < K1; j++) {
                   xmk[j] = xmk[j]-gamma_old[j]+gamma_tt(i,j);
                   gamma_old[j] = gamma_tt(i,j);
                  
                 }
               gamma_old[K1] = gamma_tt(i,K1);
               xIstar = Iistar+flag_prop;
               xmKstar = miKstar+flag_prop*gamma_tt(i,K1);
               xAg[i] = 1;
            } else {
               for (int j = 0; j < xK; j++) {
                  gamma_tt(i,j) = gamma_old[j];
               }
                
            }
        }
       }
    }

    return Rcpp::List::create(Rcpp::Named("gamma_tt") = gamma_tt,Rcpp::Named("Ag") = xAg, Rcpp::Named("mKstar") = xmKstar, Rcpp::Named("xmk") = xmk, Rcpp::Named("xIstar")=xIstar);

    END_RCPP
}

