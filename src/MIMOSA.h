#ifndef MIMOSA_H
#define MIMOSA_H
#include <Rcpp.h>
#include <Rcpp/routines.h>


#define RATE 2.4
#define DEFAULT_RATE  0.4


using namespace Rcpp;


/*
 * 18 parameters
 */
RcppExport SEXP fitMCMC(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
//RcppExport SEXP fitMCMCMultiStim(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);


//double lkbeta(const NumericMatrix::Row& ,NumericVector&);
//double lkbeta(const NumericVector &,NumericVector&);
double lkbeta(const std::vector<double>&,int,int,int);
double lkbeta(const std::vector<double>&);

double op_lgamma(double);
void loglikenull(const std::vector<double>&  ,const std::vector<double>&  ,std::vector<double>& ,std::vector<double>&,int ,int );
void loglikeresp(const std::vector<double> &stim,const std::vector<double> &alphas, const std::vector<double>  &unstim, const std::vector<double>  &alphau,std::vector<double> &output, std::vector<double> &sum_dat_alphas,std::vector<double> &sum_dat_alphau,int,int);
double alphaProposal(const std::vector<double>&,double, int);
void completeLL(std::vector<double> &z,std::vector<double> &lnull, std::vector<double> &lresp,std::vector<double> &cll, std::vector<bool> &filter,int, int);
void simZ(double &,std::vector<double>&, std::vector<double>&,std::vector<double>&,std::vector<double>&,std::vector<bool> &filter,int, int);
double simQ(std::vector<double> &z,int ,int,NumericVector);
void normalizingConstant(std::vector<double> &stim,std::vector<double> &unstim,std::vector<double> &alphas,std::vector<double> &alphau,std::vector<double> &normconst, int,int);
bool FILTER, FAST;
double EXPRATE=1000;
int ntune=0;
//double normconstIBeta(double as, double bs, double au, double bu);
double normconstIBeta(double au, double bu, double as, double bs);
void sampleP(std::vector<double>&,std::vector<double>&,std::vector<double>&,std::vector<double>&,std::vector<double>&,std::vector<double>&,std::vector<double>&,std::vector<double>&,int,int);
double nc(double as, double bs, double au,double bu,double B);
double normconstMC(double as, double bs,double au, double bu);
double dgeom(int k,double p);
double alphaDiscreteProposal(const std::vector<double> &alpha, double d, int i);

#endif
