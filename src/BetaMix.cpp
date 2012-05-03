#include <Rcpp.h>;
#include <cmath>;
#define NSAMPS 100
using namespace Rcpp;

RcppExport SEXP CompleteDataLLRcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
RcppExport SEXP MarginalNULL(SEXP ,SEXP , SEXP , SEXP , SEXP ,SEXP ,SEXP , SEXP,SEXP ,SEXP  );
RcppExport SEXP MarginalGT(SEXP ,SEXP , SEXP , SEXP , SEXP ,SEXP ,SEXP ,SEXP ,SEXP ,SEXP,SEXP  );
RcppExport SEXP MarginalNE(SEXP ,SEXP , SEXP , SEXP , SEXP ,SEXP ,SEXP ,SEXP ,SEXP ,SEXP  );


NumericVector MarginalNULLRcpp(NumericVector&, NumericVector&, IntegerVector&, IntegerVector&, NumericVector&, IntegerVector&, IntegerVector&, NumericVector&, NumericVector&, NumericVector&, NumericVector&);
NumericVector MarginalGTRcpp(NumericVector&, NumericVector&, IntegerVector&, IntegerVector&, NumericVector&, IntegerVector&, IntegerVector&, NumericVector&, NumericVector&, NumericVector&, NumericVector&,IntegerVector&);
NumericVector MarginalNERcpp(NumericVector&, NumericVector&, IntegerVector&, IntegerVector&, NumericVector&, IntegerVector&, IntegerVector&, NumericVector&, NumericVector&, NumericVector&, NumericVector&);
NumericVector betaintegralRcpp(NumericVector&,NumericVector&, NumericVector&, NumericVector&, IntegerVector&, IntegerVector&, IntegerVector&, IntegerVector&,IntegerVector&);
//NumericVector normalizingConstant(NumericVector&, NumericVector&, NumericVector&, NumericVector&);

double op_lchoose(int i, int j){return (double) ::Rf_lchoose(i+j,j);}
double op_sum(double i, double j){return i+j;}
double op_diff(double i, double j){return i-j;}
double op_lbeta(double a,double b){return ::Rf_lbeta(a,b);}
double op_fpclassify(double i){if(i<DBL_MAX&&i>-DBL_MAX){return i;}else{return -DBL_MAX;}}
double op_prod(double i, double j){return i*j;}





//R wrapper for MarginalNULL
RcppExport SEXP MarginalNULL(SEXP ns,SEXP Ns, SEXP nu, SEXP Nu, SEXP alpha0,SEXP beta0,SEXP alphaS,SEXP betaS,SEXP w,SEXP log){
	Rcpp::IntegerVector _ns(ns);
	NumericVector result(_ns.length());
	Rcpp::IntegerVector _Ns(Ns);
	Rcpp::IntegerVector _nu(nu);
	Rcpp::IntegerVector _Nu(Nu);
	Rcpp::NumericVector _alpha0(alpha0);
	Rcpp::NumericVector _beta0(beta0);
	Rcpp::NumericVector _alphaS(alphaS);
	Rcpp::NumericVector _betaS(betaS);
	Rcpp::NumericVector _w(w);
	Rcpp::NumericVector M(_Nu.length());
	Rcpp::NumericVector retval(_Ns.length());
	Rcpp::LogicalVector Log(log);

	result=MarginalNULLRcpp(M,_w,_Ns,_ns,retval,_Nu,_nu,_alpha0,_beta0, _alphaS,_betaS);
	if(!Log(0)){
		std::transform(result.begin(),result.end(),result.begin(),::exp);
	}
	return(wrap(result));
}

//R wrapper for MarginalGT
RcppExport SEXP MarginalGT(SEXP ns,SEXP Ns, SEXP nu, SEXP Nu, SEXP alpha0,SEXP beta0,SEXP alphaS,SEXP betaS,SEXP w,SEXP log, SEXP _MCINTS){
	Rcpp::IntegerVector MCINTS(_MCINTS); // number of monte carlo iterations
	Rcpp::IntegerVector _ns(ns);//vectors of data
	Rcpp::IntegerVector _Ns(Ns);
	Rcpp::IntegerVector _nu(nu);
	Rcpp::IntegerVector _Nu(Nu);
	Rcpp::NumericVector _alpha0(alpha0); //vectors of hyperparameters
	Rcpp::NumericVector _beta0(beta0);
	Rcpp::NumericVector _alphaS(alphaS);
	Rcpp::NumericVector _betaS(betaS);
	Rcpp::NumericVector _w(w); //w's
	Rcpp::LogicalVector Log(log); //return log or exponentiated values

	Rcpp::NumericVector M(_Ns.length()); //Deprecated, not used anymore
	Rcpp::NumericVector retval(_Ns.length());
	NumericVector result(_Ns.length());
	result=MarginalGTRcpp(M,_w,_Ns,_ns,retval,_Nu,_nu,_alpha0,_beta0, _alphaS,_betaS,MCINTS);
	if(!Log(0)){
	std::transform(result.begin(),result.end(),result.begin(),::exp);
	}
	return(wrap(result));
}
//R wrapper for MarginalNE
RcppExport SEXP MarginalNE(SEXP ns,SEXP Ns, SEXP nu, SEXP Nu, SEXP alpha0,SEXP beta0,SEXP alphaS,SEXP betaS,SEXP w,SEXP log){
	Rcpp::IntegerVector _ns(ns);
	NumericVector result(_ns.length());
	Rcpp::IntegerVector _Ns(Ns);
	Rcpp::IntegerVector _nu(nu);
	Rcpp::IntegerVector _Nu(Nu);
	Rcpp::NumericVector _alpha0(alpha0);
	Rcpp::NumericVector _beta0(beta0);
	Rcpp::NumericVector _alphaS(alphaS);
	Rcpp::NumericVector _betaS(betaS);
	Rcpp::NumericVector _w(w);
	Rcpp::NumericVector M(_Nu.length());
	Rcpp::NumericVector retval(_Ns.length());
	Rcpp::LogicalVector Log(log);
	result=MarginalNERcpp(M,_w,_Ns,_ns,retval,_Nu,_nu,_alpha0,_beta0, _alphaS,_betaS);
	if(!Log(0)){
	std::transform(result.begin(),result.end(),result.begin(),::exp);
	}
	return(wrap(result));
}

//Compute the Complete Data Log Likelihood
RcppExport SEXP CompleteDataLLRcpp(SEXP ns,SEXP Ns, SEXP nu, SEXP Nu, SEXP alpha0,SEXP beta0,SEXP alphaS,SEXP betaS,SEXP z,SEXP w,SEXP alternative,SEXP _MCMCINTS){
	Rcpp::IntegerVector MCMCINTS(_MCMCINTS);
	std::string ALT_GT("greater");
	std::string ALT_NE("not equal");
	double ll=0;
	Rcpp::IntegerVector _ns(ns);
	Rcpp::IntegerVector _Ns(Ns);
	Rcpp::IntegerVector _nu(nu);
	Rcpp::IntegerVector _Nu(Nu);
	Rcpp::NumericVector _alpha0(alpha0);
	Rcpp::NumericVector _beta0(beta0);
	Rcpp::NumericVector _alphaS(alphaS);
	Rcpp::NumericVector _betaS(betaS);
	Rcpp::NumericMatrix _z(z);
	Rcpp::NumericVector _w(w);

	Rcpp::NumericVector M(_z.nrow());
	Rcpp::NumericVector retval(_Ns.length());
	Rcpp::NumericVector LL(_Ns.length());
	Rcpp::NumericVector LL2(_Ns.length());


	std::string Alt = Rcpp::as<std::string>(alternative);

	retval=MarginalNULLRcpp(M,_w,_Ns,_ns,retval,_Nu,_nu,_alpha0,_beta0, _alphaS,_betaS);

	std::transform(retval.begin(),retval.end(),_z(_,0).begin(),LL.begin(),op_prod);


	if(Alt == ALT_GT){
		retval=MarginalGTRcpp(M,_w,_Ns,_ns,retval,_Nu,_nu,_alpha0,_beta0,_alphaS, _betaS,MCMCINTS);

		std::transform(retval.begin(),retval.end(),_z(_,1).begin(),LL2.begin(),op_prod);

	}else if(Alt == ALT_NE){
		retval=MarginalNERcpp(M,_w,_Ns,_ns,retval,_Nu,_nu,_alpha0,_beta0, _alphaS,_betaS);

		std::transform(retval.begin(),retval.end(),_z(_,1).begin(),LL2.begin(),op_prod);

	}
	for(int i=0;i<LL.length();i++){

		LL(i)=LL(i)+LL2(i);
	}

	ll=std::accumulate(LL.begin(),LL.end(),0.0,op_sum);

	//should be finite, but just in case..
	if(!std::isfinite(ll)){
		ll=(-DBL_MAX);
	}
	return(wrap(ll));
}

//Compute the marginal log likelihood for the alternative model where ps>pu
NumericVector MarginalGTRcpp(NumericVector &M, NumericVector &_w,IntegerVector &_Ns,IntegerVector &_ns, NumericVector &retval, IntegerVector &_Nu, IntegerVector &_nu, NumericVector &_alpha0, NumericVector &_beta0, NumericVector &_alphaS, NumericVector &_betaS,IntegerVector& MCMCINTS){
	NumericVector temp(_Ns.length());
	NumericVector nc(1);
	retval.fill(0);temp.fill(0);

	std::transform(_Ns.begin(),_Ns.end(),_ns.begin(),retval.begin(),op_lchoose);


	std::transform(_Nu.begin(),_Nu.end(),_nu.begin(),temp.begin(),op_lchoose);

	std::transform(temp.begin(),temp.end(),retval.begin(),retval.begin(),op_sum);

	std::transform(_alphaS.begin(),_alphaS.end(),_betaS.begin(),temp.begin(),op_lbeta);


	std::transform(retval.begin(),retval.end(),temp.begin(),retval.begin(),op_diff);

	std::transform(_alpha0.begin(),_alpha0.end(),_beta0.begin(),temp.begin(),op_lbeta);


	std::transform(retval.begin(),retval.end(),temp.begin(),retval.begin(),op_diff);
	temp = (betaintegralRcpp(_alphaS,_betaS,_alpha0,_beta0,_Nu, _nu, _Ns, _ns,MCMCINTS));
	std::transform(retval.begin(),retval.end(),temp.begin(),retval.begin(),op_sum);

	return retval;
}

//Compute the marginal log likelihood for the null model where ps=pu
NumericVector MarginalNULLRcpp(NumericVector &M, NumericVector &_w, IntegerVector &_Ns, IntegerVector &_ns, NumericVector &retval, IntegerVector &_Nu, IntegerVector &_nu, NumericVector &_alpha0, NumericVector &_beta0, NumericVector &_alphaS, NumericVector &_betaS){
	NumericVector temp(_Ns.length());
	NumericVector tempB(_Ns.length());
	retval.fill(0);temp.fill(0);tempB.fill(0);
	std::transform(_Ns.begin(),_Ns.end(),_ns.begin(),retval.begin(),op_lchoose);
	std::transform(_Nu.begin(),_Nu.end(),_nu.begin(),temp.begin(),op_lchoose);
	std::transform(retval.begin(),retval.end(),temp.begin(),retval.begin(),op_sum);
	std::transform(_alpha0.begin(),_alpha0.end(),_beta0.begin(),temp.begin(),op_lbeta);
	std::transform(retval.begin(),retval.end(),temp.begin(),retval.begin(),op_diff);
	std::transform(_ns.begin(),_ns.end(),_nu.begin(),temp.begin(),op_sum);
	std::transform(temp.begin(),temp.end(),_alpha0.begin(),temp.begin(),op_sum);
	std::transform(_Ns.begin(),_Ns.end(),_Nu.begin(),tempB.begin(),op_sum);
	std::transform(tempB.begin(),tempB.end(),_beta0.begin(),tempB.begin(),op_sum);
	std::transform(temp.begin(),temp.end(),tempB.begin(),temp.begin(),op_lbeta);
	std::transform(retval.begin(),retval.end(),temp.begin(),retval.begin(),op_sum);

	return(retval);
}

//Compute the Marginal log likelihood for the alternative model where ps != pu.
NumericVector MarginalNERcpp(NumericVector &M, NumericVector &_w, IntegerVector &_Ns, IntegerVector &_ns, NumericVector &retval, IntegerVector &_Nu, IntegerVector &_nu, NumericVector &_alpha0, NumericVector &_beta0, NumericVector &_alphaS, NumericVector &_betaS){
	NumericVector tempSTIM(_Ns.length());
	NumericVector tempUNSTIM(_Ns.length());
	retval.fill(0);tempSTIM.fill(0);tempUNSTIM.fill(0);
	std::transform(_Ns.begin(),_Ns.end(),_ns.begin(),tempSTIM.begin(),op_lchoose);
	std::transform(_Nu.begin(),_Nu.end(),_nu.begin(),tempUNSTIM.begin(),op_lchoose);
	std::transform(tempUNSTIM.begin(),tempUNSTIM.end(),tempSTIM.begin(),retval.begin(),op_sum);

	std::transform(_alpha0.begin(),_alpha0.end(),_beta0.begin(),tempUNSTIM.begin(),op_lbeta);
	std::transform(retval.begin(),retval.end(),tempUNSTIM.begin(),retval.begin(),op_diff);

	std::transform(_alphaS.begin(),_alphaS.end(),_betaS.begin(),tempUNSTIM.begin(),op_lbeta);
	std::transform(retval.begin(),retval.end(),tempUNSTIM.begin(),retval.begin(),op_diff);

	std::transform(_nu.begin(),_nu.end(),_alpha0.begin(),tempUNSTIM.begin(),op_sum);
	std::transform(_Nu.begin(),_Nu.end(),_beta0.begin(),tempSTIM.begin(),op_sum);

	std::transform(tempUNSTIM.begin(),tempUNSTIM.end(),tempSTIM.begin(),tempUNSTIM.begin(),op_lbeta);
	std::transform(retval.begin(),retval.end(),tempUNSTIM.begin(),retval.begin(),op_sum);

	std::transform(_ns.begin(),_ns.end(),_alphaS.begin(),tempUNSTIM.begin(),op_sum);
	std::transform(_Ns.begin(),_Ns.end(),_betaS.begin(),tempSTIM.begin(),op_sum);
	std::transform(tempUNSTIM.begin(),tempUNSTIM.end(),tempSTIM.begin(),tempUNSTIM.begin(),op_lbeta);

	std::transform(retval.begin(),retval.end(),tempUNSTIM.begin(),retval.begin(),op_sum);

	return(retval);
}



