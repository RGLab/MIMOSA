#include <Rcpp.h>
using namespace Rcpp;
/*
 * 10 parameters
 */
RcppExport SEXP fitMCMC(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
double lkbeta(const NumericMatrix::Row& ,NumericVector&);
double lkbeta(const NumericVector &,NumericVector&);

double op_lgamma(double);
void loglikenull(const NumericMatrix&  ,const NumericVector&  ,NumericVector& ,NumericMatrix& ,NumericVector&);
void loglikeresp(const NumericMatrix &stim,const NumericVector &alphas, const NumericVector  &unstim, const NumericVector  &alphau,NumericVector &output, NumericMatrix &sum_dat_alphas,NumericMatrix &sum_dat_alphau, NumericVector&);
double alphaProposal(const NumericVector&,double, int);
void completeLL(NumericMatrix &z,NumericVector &lnull, NumericVector &lresp,NumericVector &cll);
void simZ(NumericVector &,NumericVector&, NumericVector&,NumericMatrix&,NumericVector&);
double simQ(NumericMatrix &z);
