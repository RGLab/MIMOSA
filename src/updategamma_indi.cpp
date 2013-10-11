#include <Rcpp.h>
#include <R_ext/Utils.h>
#include <stdlib.h>
using namespace std;
using namespace Rcpp;


void updategamma_indi(vector<int>& xn_s,vector<int>& xn_u, vector<int>& gamma_tt, int xI, int xK, int xM, int xSS, int K1,
                     vector<double>& xalphau, vector<double>& xalphas, double xa, double xb, Rcpp::IntegerVector& xmk, int& xIstar, int& xmKstar, 
                     Rcpp::NumericVector& xpp, double xpb1, double xpb2, double xlambda, double xbeta, double xalpha, vector<double>& xmu_s, 
                     vector<double>& xmu_u, Rcpp::IntegerMatrix& xd, Rcpp::NumericMatrix& xybar_s, Rcpp::NumericMatrix& xybar_u, 
                     Rcpp::NumericMatrix& pt1, Rcpp::NumericMatrix& pt2,Rcpp::NumericMatrix& pt3, Rcpp::IntegerMatrix& xindicator, vector<int>& xAg,
                     double ab, double abI, double bI)
{  
    vector<int> gamma_old(xK);
    int xindiK1 = 0;
    Rcpp::NumericVector m(1); 
    int mik[K1]; int temp=0;
    int nik = 0; double test=0.; 
    for (int i = 0; i < xI; i++) {
        xindiK1 = xindicator(i,K1);
        int array1K[xindiK1];
        for (int j = 0; j < xK; j++) {
            gamma_old[j] = gamma_tt[xI*j+i];         
        }
        for (int ss = 0; ss < xSS; ss++) {
            int miKstar = xmKstar;
            int Iistar = xIstar;
            int flag =0; int flag_old = 0;
            double log_likeli_deno = 0.;
            int ct = 0;
            for (int j = 0; j < K1; j++) {
                gamma_tt[xI*j+i] = gamma_old[j];
                mik[j] = xmk[j]-gamma_old[j];
                flag +=gamma_tt[xI*j+i];
                if (xindicator(i,j) ==1) {array1K[ct] = j; ct +=1;}
                log_likeli_deno+=gamma_old[j]*(log(mik[j]+xa)-log(abI-1))+(1-gamma_old[j])*(log(bI-mik[j]-1)-log(abI-1));
            } 
            if (flag >1 ) {//*
                Iistar = xIstar-1; 
                miKstar = xmKstar-gamma_old[K1];
                flag_old = 1;
            }
            gamma_tt[i+xI*K1] = gamma_old[K1];
            double log_probK =log(miKstar+xa)-log(Iistar+ab); 
            double log_1probK = log(Iistar-miKstar+xb)-log(Iistar+ab);
            if(flag_old == 1) {
               log_likeli_deno+=gamma_old[K1]*(log_probK)+(1-gamma_old[K1])*(log_1probK);
            }
            // test
            //std::cout <<"llike_deno_counts = "<< log_likeli_deno << std::endl;
            if(xindiK1<=2) {m=1;}
            else{
              if (Rcpp::as<double>(Rcpp::runif(1)) <=xpp[i]) 
                 { m = Rcpp::rbinom(1,xindiK1,xpb1); }//m[0] +=1;}
              else {m =Rcpp::rbinom(1,xindiK1,xpb2); }//m[0] +=1;}
            } 
            //std::cout <<"m = "<< m[0] << std::endl;
            /* random sample */ 
           if (m[0]>0) {
              for (int k = 0; k<m[0]; k++) {             
                Rcpp::NumericVector tmp = Rcpp::runif(1,k,xindiK1);
                int vec = (int) tmp[0];
                int arr = array1K[k];
                array1K[k] = array1K[vec];
                array1K[vec] = arr;
                gamma_tt[xI*array1K[k]+i] = 1-gamma_tt[xI*array1K[k]+i];
              }
            } 
            //if (i == 1) {std::cout <<"m = "<< m[0] << std::endl; }
            int lplace0_new = 0; int lplace0 = 0;
            double log_likeli_nume = 0.0;
            for (int k = 0; k < K1; k++) {
               temp = xI*k+i;
               log_likeli_nume += gamma_tt[temp]*(log(mik[k]+xa)-log(abI-1))+(1-gamma_tt[temp])*(log(bI-mik[0]-1)-log(abI-1));
               if (gamma_tt[temp] == 0) { lplace0_new +=1;}
               if (gamma_old[k] == 0){ lplace0 += 1;}            
            }
            int lplace1_new = K1-lplace0_new;
            int lplace1 = K1-lplace0;
            int flag_prop = 0; temp = xI*K1+i;
            if (lplace1_new == 1) { //*
               gamma_tt[i+xI*K1] = 1; 
               lplace1_new +=1;
             } else if (lplace1_new == 0) {
               gamma_tt[i+xI*K1] = 0;
               lplace0_new += 1;
            } else {
               flag_prop = 1;
               log_likeli_deno += log(0.5); // proposal
               //test
               //std::cout <<"llike_deno_prop = "<< log_likeli_deno << std::endl;
               if (Rcpp::as<double>(Rcpp::runif(1)) <= 0.5) {
                    gamma_tt[i+xI*K1] = 1;
               }else{
                    gamma_tt[xI*K1+i] = 0;
               }
               if (gamma_tt[xI*K1+i] == 1) {
                   lplace1_new +=1;
                   log_likeli_nume +=log_probK; 
               } else {
                   lplace0_new += 1;
                   log_likeli_nume += log_1probK;                   
               }           
            } 
            //test
            //std::cout <<"llike_nume_counts = "<< log_likeli_nume << std::endl;
             if ( lplace1 == 1) {   //*          
               lplace1 +=1;
            } else if (lplace1 == 0) {
               lplace0+= 1;
            } else {
               log_likeli_nume += log(0.5); // proposal 
               //test
               //std::cout <<"llike_nume_prop = "<< log_likeli_nume << std::endl;
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
               if (gamma_tt[xI*k+i] == 0) {
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
                   temp = xI*k+i;
                    double summ = xalphau[k]+xn_u[temp]+xn_s[temp];
                    sum_alphau_nus += summ; 
                    sum_gammaln_nu +=lgamma(summ); 
                }
                log_likeli_nume +=(sum_gammaln_nu-lgamma(sum_alphau_nus));

                //L0ik
                for ( int k = 0; k<K1; k++) {
                  temp = xI*k+i;
                   nik = xn_u[temp]+xn_s[temp];
                   if (nik>0) {
                     for ( int p = 0; p<xM; p++) {
                        if (xd(k,p)==1){
                        double beta_np = pt3(temp,p)+xbeta+pow((xmu_u[p]-(xn_s[temp]*xybar_s(temp,p)+xn_u[temp]*xybar_u(temp,p))/nik),2)/(2*(xlambda+1/nik))+
                                 0.5*xn_s[temp]*pow(xybar_s(temp,p),2)+ 0.5*xn_u[temp]*pow(xybar_u(temp,p),2)-0.5*pow((xn_s[temp]*xybar_s(temp,p)+xn_u[temp]*xybar_u(temp,p)),2)/nik;
                          double alpha_np = nik/2+xalpha;
                          test = xalpha*log(xbeta)+lgamma(alpha_np)-0.5*log(xlambda*nik+1)-lgamma(xalpha)-alpha_np*log(beta_np);
                          //std::cout <<"test1 = "<< test << std::endl; //
                          log_likeli_nume += test; 
                        }                     
                      }
                    }         
               }            
            } else if (lplace0_new == 0) {
                double sum_alphau_nus = 0.0;
                double sum_gammaln_nu = 0.0;
                double sum_gammaln = 0.0;
                double sum_alphas = 0.0;
                double sum_gammaln_ns = 0.0;
                double sum_alphas_ns = 0.0;

                for ( int k=0; k<xK; k++) {
                    temp = xI*k+i;
                    double summ = xalphau[k]+xn_u[temp];
                    sum_alphau_nus += summ; 
                    sum_gammaln_nu +=lgamma(summ); 

                    sum_gammaln +=lgamma(xalphas[k]);
                    sum_alphas +=xalphas[k];
                    sum_gammaln_ns += (lgamma(xalphas[k] + xn_s[temp]));
                    sum_alphas_ns += (xalphas[k]+xn_s[temp]);
                }
                log_likeli_nume +=(sum_gammaln_nu-lgamma(sum_alphau_nus));
                log_likeli_nume -= (sum_gammaln - lgamma(sum_alphas));
                log_likeli_nume +=(sum_gammaln_ns -lgamma(sum_alphas_ns));
  
                 for ( int k = 0; k<K1; k++) {
                   temp = xI*k+i;
                   nik = xn_u[temp]+xn_s[temp];
                   if (xn_u[temp] >0 && xn_s[temp]>0) {
                     for ( int p = 0; p<xM; p++) {
                       if (xd(k,p)==1){
                         double beta_np = pt3(temp,p)+xbeta+pow((xmu_u[p]-xybar_u(temp,p)),2)/(2*(xlambda+1/xn_u[temp]))+
                            pow((xmu_s[p]-xybar_s(temp,p)),2)/(2*(xlambda+1/xn_s[temp]));
                         double alpha_np = nik/2+xalpha;
                         test = xalpha*log(xbeta)+lgamma(alpha_np)-0.5*log(xlambda*xn_u[temp]+1) -0.5*log(xlambda*xn_s[temp]+1) 
                             -lgamma(xalpha)-alpha_np*log(beta_np); 
                         //std::cout <<"test2 = "<< test << std::endl;//
                         log_likeli_nume += test;
                        }                  
                      }
                    }else if (xn_u[temp]==0 && xn_s[temp]>0) {
                       for ( int p = 0; p<xM; p++) {
                         if (xd(k,p)==1){
                            double beta_np = pt1(temp,p)+xbeta+pow((xmu_s[p]-xybar_s(temp,p)),2)/(2*(xlambda+1/xn_s[temp]));
                            double alpha_np = nik/2+xalpha;
                            test=xalpha*log(xbeta)+lgamma(alpha_np)-0.5*log(xlambda*xn_s[temp]+1) 
                                 -lgamma(xalpha)-alpha_np*log(beta_np);    
                            //std::cout <<"test3 = "<< test << std::endl;//
                            log_likeli_nume += test;
                          }                  
                        }
                    } else if (xn_u[temp]>0 && xn_s[temp]==0) {
                       for ( int p = 0; p<xM; p++) {
                          if (xd(k,p)==1){
                             double beta_np =pt2(temp,p)+xbeta+pow((xmu_u[p]-xybar_u(temp,p)),2)/(2*(xlambda+1/xn_u[temp]));
                             double alpha_np = nik/2+xalpha;
                             test = xalpha*log(xbeta)+lgamma(alpha_np)-0.5*log(xlambda*xn_u[temp]+1) -lgamma(xalpha)-alpha_np*log(beta_np);  
                             //std::cout <<"test4 = "<< test << std::endl;//
                             log_likeli_nume += test;
                           }                  
                        }
                    }
               }             
  
            } else {
                double sum_gammaln = 0.0;
                double sum_alphas = 0.0;
                double sum_gammaln_ns = 0.0;
                double sum_alphas_ns = 0.0;
                double sum_alphau_nu = 0.0;
                double sum_alphau_nus = 0.0;
                double sum_gammaln_nu = 0.0;
                for (int k = 0; k<lplace1_new; k++) {
                  sum_alphau_nu +=(xalphau[place1_new[k]]+xn_u[i+xI*place1_new[k]]);
                  sum_alphau_nus+=(xalphau[place1_new[k]]+xn_u[i+xI*place1_new[k]]+xn_s[i+xI*place1_new[k]]);
                  sum_gammaln_nu +=lgamma(xalphau[place1_new[k]]+xn_u[i+xI*place1_new[k]]);
                  sum_gammaln +=lgamma(xalphas[place1_new[k]]);
                  sum_alphas +=xalphas[place1_new[k]];
                  sum_gammaln_ns += (lgamma(xalphas[place1_new[k]] + xn_s[i+xI*place1_new[k]]));
                  sum_alphas_ns += (xalphas[place1_new[k]]+xn_s[i+xI*place1_new[k]]);
                }
                log_likeli_nume -= (sum_gammaln - lgamma(sum_alphas));
                log_likeli_nume +=(sum_gammaln_ns -lgamma(sum_alphas_ns));
                log_likeli_nume += (sum_gammaln_nu-lgamma(sum_alphau_nu));

                log_likeli_nume +=lgamma(sum_alphau_nus+1);
                sum_gammaln_nu = 0.0;

                for (int k = 0; k<lplace0_new; k++) {
                    sum_alphau_nus+=(xalphau[place0_new[k]]+xn_u[i+xI*place0_new[k]]+xn_s[i+xI*place0_new[k]]);
                    sum_gammaln_nu +=lgamma(xalphau[place0_new[k]]+xn_u[i+xI*place0_new[k]]+xn_s[i+xI*place0_new[k]]);
                } 
                log_likeli_nume +=sum_gammaln_nu;   
                log_likeli_nume -=lgamma(sum_alphau_nus+1);
                for (int k = 0; k<K1; k++) {
                   temp = xI*k+i;
                    nik = xn_u[temp]+xn_s[temp];
                    if(gamma_tt[xI*k+i]==1) {
                        if (xn_u[temp] >0 && xn_s[temp]>0) {
                           for ( int p = 0; p<xM; p++) {
                               if (xd(k,p)==1){
                                 double beta_np = pt3(temp,p)+xbeta+pow((xmu_u[p]-xybar_u(temp,p)),2)/(2*(xlambda+1/xn_u[temp]))+
                                        pow((xmu_s[p]-xybar_s(temp,p)),2)/(2*(xlambda+1/xn_s[temp]));
                                 double alpha_np = nik/2+xalpha;
                                 test = xalpha*log(xbeta)+lgamma(alpha_np)-0.5*log(xlambda*xn_u[temp]+1) -0.5*log(xlambda*xn_s[temp]+1) 
                                       -lgamma(xalpha)-alpha_np*log(beta_np);
                                 //std::cout <<"test5 = "<< test << std::endl;//  
                                 log_likeli_nume += test;
                               }                  
                            }
                        }else if (xn_u[temp]==0 && xn_s[temp]>0) {
                           for ( int p = 0; p<xM; p++) {
                               if (xd(k,p)==1){
                                  double beta_np = pt1(temp,p)+xbeta+pow((xmu_s[p]-xybar_s(temp,p)),2)/(2*(xlambda+1/xn_s[temp]));
                                  double alpha_np = nik/2+xalpha;
                                  test = xalpha*log(xbeta)+lgamma(alpha_np)-0.5*log(xlambda*xn_s[temp]+1) 
                                       -lgamma(xalpha)-alpha_np*log(beta_np); 
                                  //std::cout <<"test6 = "<< test << std::endl;//    
                                  log_likeli_nume += test; 
                                }                   
                           }
                        } else if (xn_u[temp]>0 && xn_s[temp]==0) {
                           for ( int p = 0; p<xM; p++) {
                               if (xd(k,p)==1){
                                 double beta_np = pt2(temp,p)+xbeta+pow((xmu_u[p]-xybar_u(temp,p)),2)/(2*(xlambda+1/xn_u[temp]));
                                 double alpha_np = nik/2+xalpha;
                                 test =  xalpha*log(xbeta)+lgamma(alpha_np)-0.5*log(xlambda*xn_u[temp]+1) -lgamma(xalpha)-alpha_np*log(beta_np); 
                                 //std::cout <<"test7 = "<< test << std::endl;// 
                                 log_likeli_nume += test;  
                               }                  
                           }
                        }
                     }else{
                          nik = xn_u[temp]+xn_s[temp];
                          if (nik>0) {
                            for ( int p = 0; p<xM; p++) {
                               if (xd(k,p)==1){
                                  double beta_np = pt3(temp,p)+xbeta+pow((xmu_u[p]-(xn_s[temp]*xybar_s(temp,p)+xn_u[temp]*xybar_u(temp,p))/nik),2)
                                            /(2*(xlambda+1/nik))+0.5*xn_s[temp]*pow(xybar_s(temp,p),2)+ 0.5*xn_u[temp]*pow(xybar_u(temp,p),2)-
                                            0.5*pow((xn_s[temp]*xybar_s(temp,p)+xn_u[temp]*xybar_u(temp,p)),2)/nik;
                                  double alpha_np = nik/2+xalpha;
                                  test = xalpha*log(xbeta)+lgamma(alpha_np)-0.5*log(xlambda*nik+1)-lgamma(xalpha)-alpha_np*log(beta_np);
                                  //std::cout <<"test8 = "<< test << std::endl;//  
                                  log_likeli_nume +=test;
                                }                     
                            }
                          }       
                     }
                }
            }
        //
            if (lplace1 == 0) {
                double sum_alphau_nus = 0.0;
                double sum_gammaln_nu = 0.0;
                for ( int k = 0; k<xK; k ++) {
                    temp = xI*k+i;
                    double summ = xalphau[k]+xn_u[temp]+xn_s[temp];
                    sum_alphau_nus += summ; 
                    sum_gammaln_nu +=lgamma(summ); 
                }
                log_likeli_deno +=(sum_gammaln_nu-lgamma(sum_alphau_nus));

                for ( int k = 0; k<K1; k++) {
                  temp = xI*k+i;
                   nik = xn_u[temp]+xn_s[temp];
                   if (nik>0) {
                     for ( int p = 0; p<xM; p++) {
                        if (xd(k,p)==1){
                          double beta_np = pt3(temp,p)+xbeta+pow((xmu_u[p]-(xn_s[temp]*xybar_s(temp,p)+xn_u[temp]*xybar_u(temp,p))/nik),2)/
                                         (2*(xlambda+1/nik))+0.5*xn_s[temp]*pow(xybar_s(temp,p),2)+ 0.5*xn_u[temp]*pow(xybar_u(temp,p),2)
                                         -0.5*pow((xn_s[temp]*xybar_s(temp,p)+xn_u[temp]*xybar_u(temp,p)),2)/nik;
                          double alpha_np = nik/2+xalpha;
                          test = xalpha*log(xbeta)+lgamma(alpha_np)-0.5*log(xlambda*nik+1)-lgamma(xalpha)-alpha_np*log(beta_np); 
                          //std::cout <<"test9 = "<< test << std::endl;//
                          log_likeli_deno += test;
                        }                     
                      }
                    }         
                }            
            } else if (lplace0 == 0) {
                double sum_alphau_nus = 0.0;
                double sum_gammaln_nu = 0.0;
                double sum_gammaln = 0.0;
                double sum_alphas = 0.0;
                double sum_gammaln_ns = 0.0;
                double sum_alphas_ns = 0.0;

                for ( int k=0; k<xK; k++) {
                   temp = xI*k+i;
                    double summ = xalphau[k]+xn_u[temp];
                    sum_alphau_nus += summ; 
                    sum_gammaln_nu +=lgamma(summ); 

                    sum_gammaln +=lgamma(xalphas[k]);
                    sum_alphas +=xalphas[k];
                    sum_gammaln_ns += (lgamma(xalphas[k] + xn_s[temp]));
                    sum_alphas_ns += (xalphas[k]+xn_s[temp]);
                }
                log_likeli_deno +=(sum_gammaln_nu-lgamma(sum_alphau_nus));
                log_likeli_deno -= (sum_gammaln - lgamma(sum_alphas));
                log_likeli_deno +=(sum_gammaln_ns -lgamma(sum_alphas_ns));

                for ( int k = 0; k<K1; k++) {
                   temp = xI*k+i;
                   nik = xn_u[temp]+xn_s[temp];
                   if (xn_u[temp] >0 && xn_s[temp]>0) {
                     for ( int p = 0; p<xM; p++) {
                       if (xd(k,p)==1){
                         double beta_np = pt3(temp,p)+xbeta+pow((xmu_u[p]-xybar_u(temp,p)),2)/(2*(xlambda+1/xn_u[temp]))+
                            pow((xmu_s[p]-xybar_s(temp,p)),2)/(2*(xlambda+1/xn_s[temp]));
                         double alpha_np = nik/2+xalpha;
                         test= xalpha*log(xbeta)+lgamma(alpha_np)-0.5*log(xlambda*xn_u[temp]+1) -0.5*log(xlambda*xn_s[temp]+1) 
                             -lgamma(xalpha)-alpha_np*log(beta_np);    
                         //std::cout <<"test10 = "<< test << std::endl;//
                         log_likeli_deno +=test;
                        }                  
                      }
                    }else if (xn_u[temp]==0 && xn_s[temp]>0) {
                       for ( int p = 0; p<xM; p++) {
                         if (xd(k,p)==1){
                            double beta_np = pt1(temp,p)+xbeta+pow((xmu_s[p]-xybar_s(temp,p)),2)/(2*(xlambda+1/xn_s[temp]));
                            double alpha_np = nik/2+xalpha; 
                            test = xalpha*log(xbeta)+lgamma(alpha_np)-0.5*log(xlambda*xn_s[temp]+1) 
                                 -lgamma(xalpha)-alpha_np*log(beta_np);  
                            //std::cout <<"test11 = "<< test << std::endl;//
                            log_likeli_deno += test;  
                          }                  
                        }
                    } else if (xn_u[temp]>0 && xn_s[temp]==0) {
                       for ( int p = 0; p<xM; p++) {
                          if (xd(k,p)==1){
                             double beta_np =pt2(temp,p)+xbeta+pow((xmu_u[p]-xybar_u(temp,p)),2)/(2*(xlambda+1/xn_u[temp]));
                             double alpha_np = nik/2+xalpha;
                             test=xalpha*log(xbeta)+lgamma(alpha_np)-0.5*log(xlambda*xn_u[temp]+1) -lgamma(xalpha)-alpha_np*log(beta_np);
                             //std::cout <<"test12 = "<< test << std::endl;// 
                             log_likeli_deno += test;
                           }                  
                        }
                    }
               }            
            } else {
                double sum_gammaln = 0.0;
                double sum_alphas = 0.0;
                double sum_gammaln_ns = 0.0;
                double sum_alphas_ns = 0.0;
                double sum_alphau_nu = 0.0;
                double sum_alphau_nus = 0.0;
                double sum_gammaln_nu = 0.0;
                for (int k = 0; k<lplace1; k++) {
                  sum_alphau_nu +=(xalphau[place1[k]]+xn_u[i+xI*place1[k]]);
                  sum_alphau_nus+=(xalphau[place1[k]]+xn_u[i+xI*place1[k]]+xn_s[i+xI*place1[k]]);
                  sum_gammaln_nu +=lgamma(xalphau[place1[k]]+xn_u[i+xI*place1[k]]);
                  sum_gammaln +=lgamma(xalphas[place1[k]]);
                  sum_alphas +=xalphas[place1[k]];
                  sum_gammaln_ns += (lgamma(xalphas[place1[k]] + xn_s[i+xI*place1[k]]));
                  sum_alphas_ns += (xalphas[place1[k]]+xn_s[i+xI*place1[k]]);
                }
                log_likeli_deno -= (sum_gammaln - lgamma(sum_alphas));
                log_likeli_deno +=(sum_gammaln_ns -lgamma(sum_alphas_ns));
                log_likeli_deno += (sum_gammaln_nu-lgamma(sum_alphau_nu));

                log_likeli_deno +=lgamma(sum_alphau_nus+1);
                sum_gammaln_nu = 0.0;
                for (int k = 0; k<lplace0; k++) {
                    sum_alphau_nus+=(xalphau[place0[k]]+xn_u[i+xI*place0[k]]+xn_s[i+xI*place0[k]]);
                    sum_gammaln_nu +=lgamma(xalphau[place0[k]]+xn_u[i+xI*place0[k]]+xn_s[i+xI*place0[k]]);
                } 
                log_likeli_deno +=sum_gammaln_nu;   
                log_likeli_deno -=lgamma(sum_alphau_nus+1);

                for (int k = 0; k<K1; k++) {
                    temp = xI*k+i;
                    nik = xn_u[temp]+xn_s[temp];
                    if(gamma_old[k]==1) {
                        if (xn_u[temp] >0 && xn_s[temp]>0) {
                           for ( int p = 0; p<xM; p++) {
                               if (xd(k,p)==1){
                                 double beta_np = pt3(temp,p)+xbeta+pow((xmu_u[p]-xybar_u(temp,p)),2)/(2*(xlambda+1/xn_u[temp]))+
                                        pow((xmu_s[p]-xybar_s(temp,p)),2)/(2*(xlambda+1/xn_s[temp]));
                                 double alpha_np = nik/2+xalpha;
                                 test =  xalpha*log(xbeta)+lgamma(alpha_np)-0.5*log(xlambda*xn_u[temp]+1) -0.5*log(xlambda*xn_s[temp]+1) 
                                       -lgamma(xalpha)-alpha_np*log(beta_np);    
                                 //std::cout <<"test13 = "<< test << std::endl;//
                                 log_likeli_deno +=test;
                               }                  
                            }
                        }else if (xn_u[temp]==0 && xn_s[temp]>0) {
                           for ( int p = 0; p<xM; p++) {
                               if (xd(k,p)==1){
                                  double beta_np =pt1(temp,p)+xbeta+pow((xmu_s[p]-xybar_s(temp,p)),2)/(2*(xlambda+1/xn_s[temp]));
                                  double alpha_np = nik/2+xalpha;
                                  test= xalpha*log(xbeta)+lgamma(alpha_np)-0.5*log(xlambda*xn_s[temp]+1) 
                                       -lgamma(xalpha)-alpha_np*log(beta_np);  
                                  //std::cout <<"test14 = "<< test << std::endl;//
                                  log_likeli_deno += test;
                                }                   
                           }
                        } else if (xn_u[temp]>0 && xn_s[temp]==0) {
                           for ( int p = 0; p<xM; p++) {
                               if (xd(k,p)==1){
                                 double beta_np = pt2(temp,p)+xbeta+pow((xmu_u[p]-xybar_u(temp,p)),2)/(2*(xlambda+1/xn_u[temp]));
                                 double alpha_np = nik/2+xalpha;
                                 test= xalpha*log(xbeta)+lgamma(alpha_np)-0.5*log(xlambda*xn_u[temp]+1) -lgamma(xalpha)-alpha_np*log(beta_np);    
                                 //std::cout <<"test15 = "<< test << std::endl;//
                                 log_likeli_deno +=test;
                               }                  
                           }
                        }
                     }else{
                          nik = xn_u[temp]+xn_s[temp];
                          if (nik>0) {
                            for ( int p = 0; p<xM; p++) {
                               if (xd(k,p)==1){
                               double beta_np = pt3(temp,p)+xbeta+pow((xmu_u[p]-(xn_s[temp]*xybar_s(temp,p)+xn_u[temp]*xybar_u(temp,p))/nik),2)
                                               /(2*(xlambda+1/nik))+0.5*xn_s[temp]*pow(xybar_s(temp,p),2)+ 0.5*xn_u[temp]*pow(xybar_u(temp,p),2)
                                               -0.5*pow((xn_s[temp]*xybar_s(temp,p)+xn_u[temp]*xybar_u(temp,p)),2)/nik;
                                  double alpha_np = nik/2+xalpha;
                                  test=xalpha*log(xbeta)+lgamma(alpha_np)-0.5*log(xlambda*nik+1)-lgamma(xalpha)-alpha_np*log(beta_np); 
                                  //std::cout <<"test16 = "<< test << std::endl;//
                                  log_likeli_deno += test;
                                }                     
                            }
                          }       
                     }
                }
            }
            //test
            //std::cout <<"llike_deno_cts = "<< log_likeli_deno << std::endl;//
            //std::cout <<"llike_nume_cts = "<< log_likeli_nume << std::endl;//
            if (log(Rcpp::as<double>(Rcpp::runif(1)) ) <= (log_likeli_nume - log_likeli_deno)) {
               for (int j = 0; j < K1; j++) {
                   xmk[j] = xmk[j]-gamma_old[j]+gamma_tt[xI*j+i];
                   gamma_old[j] = gamma_tt[xI*j+i];
                  
                 }
               gamma_old[K1] = gamma_tt[xI*K1+i];
               xIstar = Iistar+flag_prop;
               xmKstar = miKstar+flag_prop*gamma_tt[xI*K1+i];
               xAg[i] = 1;
            } else {
              xAg[i] = 0;
              for (int j = 0; j < xK; j++) {
                gamma_tt[xI*j+i] = gamma_old[j];
              }
                
            }
        } 
    }

}
