#include <Rcpp.h>
#include <armadillo>
#include <Rmath.h>
using namespace Rcpp;
using namespace arma;

double op_sum(double, double);
double op_lbeta(double,double);

NumericVector betaintegralRcpp(NumericVector&,NumericVector&, NumericVector&, NumericVector&, IntegerVector&, IntegerVector&, IntegerVector&, IntegerVector&,IntegerVector&);
RcppExport SEXP betaintegral(SEXP ,SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP );
NumericVector normalizingConstant(NumericVector&, NumericVector&, NumericVector&, NumericVector&);

/*Compute the integral via Monte-Carlo integration*/
NumericVector betaintegralRcpp(NumericVector &alphaS,NumericVector &betaS, NumericVector &alpha0, NumericVector &beta0, IntegerVector &Nu, IntegerVector &nu, IntegerVector &Ns, IntegerVector &ns, IntegerVector &MCITER){
	NumericVector samps(MCITER(0));
	NumericVector SAMPS(MCITER(0));
	NumericVector S(MCITER(0));
	NumericVector SS(MCITER(0));
	NumericVector result(Nu.length());
	NumericVector temp(Nu.length());
	NumericVector nc(0);
    RNGScope scope ;


	double init=0;
	arma::vec a(alpha0.length()), b(beta0.length()),A(alphaS.length()), B(betaS.length());
	NumericVector aa(1);
	NumericVector bb(1);
	NumericVector AA(1);
	NumericVector BB(1);

	std::transform(alpha0.begin(),alpha0.end(),nu.begin(),a.begin(),::op_sum);

	std::transform(beta0.begin(),beta0.end(),Nu.begin(),b.begin(),::op_sum);

	std::transform(alphaS.begin(),alphaS.end(),ns.begin(),A.begin(),::op_sum);

	std::transform(betaS.begin(),betaS.end(),Ns.begin(),B.begin(),::op_sum);


	std::transform(a.begin(),a.end(),b.begin(),result.begin(),::op_lbeta);
//
	std::transform(A.begin(),A.end(),B.begin(),temp.begin(),::op_lbeta);
//
	std::transform(result.begin(),result.end(),temp.begin(),result.begin(),::op_sum);
	aa(0)=alpha0(0);
	bb(0)=beta0(0);
	AA(0)=alphaS(0);
	BB(0)=betaS(0);
	for(int i=0;i<(int)alpha0.length();i++){
		samps = rbeta(MCITER(0),a(i),b(i));
		SAMPS = rbeta(MCITER(0),aa(0),bb(0));

		//pbeta lower tail, log=F
		S=pbeta(samps,A(i),B(i),false,false);
		SS=pbeta(SAMPS,AA(0),BB(0),false,false);
		double NUM=std::accumulate(SS.begin(),SS.end(),init,::op_sum)/SS.length();
		double num=std::accumulate(S.begin(),S.end(),init,::op_sum)/S.length();

		result(i)=result(i)+log(num)-log(NUM);
	}
	return result;
}

RcppExport SEXP betaintegral(SEXP _alphaS,SEXP _betaS, SEXP _alpha0, SEXP _beta0, SEXP _Nu, SEXP _nu, SEXP _Ns, SEXP _ns, SEXP _MCITER){
 Rcpp::NumericVector alphaS(_alphaS);
 Rcpp::NumericVector betaS(_betaS);
 Rcpp::NumericVector alpha0(_alpha0);
 Rcpp::NumericVector beta0(_beta0);
 Rcpp::IntegerVector Nu(_Nu);
 Rcpp::IntegerVector nu(_nu);
 Rcpp::IntegerVector Ns(_Ns);
 Rcpp::IntegerVector ns(_ns);
 Rcpp::IntegerVector MCITER(_MCITER);
 Rcpp::NumericVector result(Ns.length());
 result=betaintegralRcpp(alphaS,betaS,alpha0,beta0,Nu,nu,Ns,ns,MCITER);
 return(wrap(result));
}
