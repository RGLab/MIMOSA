#ifndef MIMOSA_H
#define MIMOSA_H
#include <Rcpp.h>


#define RATE 2.4
#define DEFAULT_RATE  0.4


using namespace Rcpp;


/*
 * 10 parameters
 */
RcppExport SEXP fitMCMC(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
//double lkbeta(const NumericMatrix::Row& ,NumericVector&);
//double lkbeta(const NumericVector &,NumericVector&);
double lkbeta(const NumericMatrix::Row& );
double lkbeta(const NumericVector &);

double op_lgamma(double);
void loglikenull(const NumericMatrix&  ,const NumericVector&  ,NumericVector& ,NumericMatrix& );
void loglikeresp(const NumericMatrix &stim,const NumericVector &alphas, const NumericVector  &unstim, const NumericVector  &alphau,NumericVector &output, NumericMatrix &sum_dat_alphas,NumericMatrix &sum_dat_alphau);
double alphaProposal(const NumericVector&,double, int);
void completeLL(NumericMatrix &z,NumericVector &lnull, NumericVector &lresp,NumericVector &cll, LogicalVector &filter);
void simZ(NumericVector &,NumericVector&, NumericVector&,NumericMatrix&,NumericVector&,LogicalVector &filter);
double simQ(NumericMatrix &z);
bool FILTER = false,FAST=false;
double LOWER=0.15, UPPER=0.5;
#endif
