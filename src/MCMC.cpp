#include <Rcpp.h>
#include <armadillo>
#include <RcppArmadillo.h>
#include <assert.h>
#include <omp.h>
#include <Rmath.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <R_ext/Utils.h>
#include "MIMOSA.h"
//#define NDEBUG
/*
 * 16 parameters
 */
RcppExport SEXP fitMCMC(SEXP _stim, SEXP _unstim, SEXP _alphas, SEXP _alphau, SEXP _q, SEXP _z,SEXP _iter, SEXP _burn, SEXP _thin, SEXP _tune,SEXP _outfile, SEXP _filter, SEXP _UPPER, SEXP _LOWER,SEXP _FILTER, SEXP _FAST,SEXP _EXPRATE,SEXP _fixedNULL){
	BEGIN_RCPP
	//TODO add argument to pass the complete list of unstimulated samples.
	using namespace Rcpp;
	using namespace arma;
	using namespace std;
	bool fixed = false;

	Rcpp::RNGScope globalscope;
	//	printf("%f\n",normconstIBeta(10,3000,1.5,3000));
	//	exit(0);
	/*
	 * Copy R variables to standard vectors
	 */
	Rcpp::LogicalVector rfilter(_filter);
	std::vector<bool> filter(rfilter.length(),0);
	filter.resize(rfilter.length());
	copy(rfilter.begin(),rfilter.end(),filter.begin());

	FILTER = Rcpp::as<bool> (_FILTER);
	FAST = Rcpp::as<bool> (_FAST);
	bool fixedNULL = Rcpp::as<bool>(_fixedNULL);
	double UPPER = Rcpp::as<double>(_UPPER);
	double LOWER = Rcpp::as<double>(_LOWER);
	EXPRATE = Rcpp::as<double>(_EXPRATE);

	std::vector< double > stdalphas = Rcpp::as<std::vector<double> >(_alphas);
	std::vector< double > stdalphau = Rcpp::as<std::vector<double> 	>(_alphau);
	double q = Rcpp::as<double>(_q);

	std::vector< double > z = Rcpp::as<vector <double> >(_z);
	std::vector< double > stdstim = Rcpp::as < vector < double > > (_stim);
	std::vector< double > stdunstim = Rcpp::as < vector < double > > (_unstim);


	/*
	 * Wrap R variables in Rcpp objects
	 */
	Rcpp::NumericMatrix const stim(_stim);
	Rcpp::NumericMatrix const unstim(_unstim);
	Rcpp::NumericVector const iter(_iter);
	Rcpp::NumericVector const burn(_burn);
	Rcpp::NumericVector const thin(_thin);
	Rcpp::NumericVector const tune(_tune);

	std::string outfile = Rcpp::as<std::string>(_outfile);

	std::string outfilep(outfile.data());
	outfilep.append("P");

	printf("Creating %s\n",outfile.data());
	FILE* file = fopen(outfile.data(),"w");
	FILE* fileP = fopen(outfilep.data(),"w");
	if(file==NULL|fileP==NULL){
		return(wrap("Can't open file!"));
	}
	/*
	 * Parameters
	 */
	const double NITERS=iter[0];
	const double BURNIN=burn[0];
	const double THINNING=thin[0];
	const double TUNING=tune[0];
	double realitcounter=0;

	/*
	 * Assert that dimensions match
	 */
	assert(stdalphas.size()==stim.ncol());
	assert(stim.nrow()==unstim.nrow());
	assert(stim.ncol()==unstim.ncol());
	assert(stdalphas.size()==stdalphau.size());
	/*
	 * Dimensions of the problem
	 */
	const int k = stim.ncol();
	const int P = stim.nrow();

	/*
	 * Output variables
	 */
	std::vector<double> ps(P,0),pu(P,0);
	std::vector <double> stdllnullRes(P,0);
	std::vector <double> stdllrespRes(P,0);
	std::vector <double> stdllnullResNew(P,0);
	std::vector <double> stdllrespResNew(P,0);


	std::vector <double> stdsum_data_alpha(P*k,0);
	std::vector <double> stdsum_data_alphau(P*k,0);

	arma::vec asi, aui;
	/*
	 * Armadillo objects will store the simulations from the tuning phase
	 * so we can easily compute variances and so forth.
	 */
	arma::mat Ms(int(TUNING),stdalphas.size());
	arma::mat Mu(int(TUNING),stdalphau.size());

	/*
	 * Precompute and preallocate a few things we'll need
	 */
	std::vector <double> stdsum_stim_unstim(P*k,0);
	std::transform(stdstim.begin(),stdstim.end(),stdunstim.begin(),stdsum_stim_unstim.begin(),plus<double>());


	std::vector<double> stdnextalphavec(stdalphas.size(),0);
	stdnextalphavec.resize(stdalphas.size());

	NumericVector sigmas(stdalphas.size(),10.0);
	NumericVector sigmau(stdalphau.size(),10.0);

	arma::vec rateS(stdalphas.size());
	arma::vec rateU(stdalphau.size());

	rateS.fill(RATE);
	rateU.fill(RATE);
	NumericVector accepts(stdalphas.size(),0.0);
	NumericVector acceptu(stdalphas.size(),0.0);
	std::vector<double> cll(P,0.0);
	std::vector<double> p(P,0.0);
	std::vector<double> cz(P,0.0);


	double prior=0,priornext=0;
	double oldll=0,newll=0;


	/*
	 * Construct the header string
	 *
	 */
	std::stringstream headers(stringstream::in|stringstream::out);
	std::stringstream headersP(stringstream::in|stringstream::out);

	//write the z's and p's to the second file
	for(int i=0;i<P;i++){
		headersP<<"z."<<i<<"\t";
	}
	for(int i=0;i<P-1;i++){
		headersP<<"ps."<<i<<"\t"<<"pu."<<i<<"\t";
	}
	headersP<<"ps."<<P-1<<"\t"<<"pu."<<P-1<<std::endl;



	for(int i=0;i<stdalphas.size();i++){
		headers<<"alphas."<<i<<"\t";
	}
	for(int i=0;i<stdalphau.size();i++){
		headers<<"alphau."<<i<<"\t";
	}
	headers<<"q"<<std::endl;
	//counters
	int iteration=0,j=0;
	/*
	 * Run the MCMC algorithm
	 */

	/*
	 * Initialize the normalizing constants
	 */


	for(iteration = 0; iteration < NITERS; iteration++){
		for(j=0;j<k;j++){
			/*
			 * prior for the current alphas_j
			 */
			prior=::Rf_dexp(stdalphas[j],EXPRATE,true);
			std::copy(stdalphas.begin(),stdalphas.end(),stdnextalphavec.begin());

			//current null marginal log likelihood
			loglikenull(stdsum_stim_unstim,stdalphau,stdllnullRes,stdsum_data_alphau,P,k);

			//If alternative is greater compute the ratio of normalizing constants for the alternative marginal log likelihood
			if(FILTER&k==2&!FAST){
				normalizingConstant(stdstim,stdunstim,stdalphas,stdalphau,stdllrespRes,P,k);
			}else{
				//otherwise the two sided marginal log likelihood
				loglikeresp(stdstim,stdalphas,stdunstim,stdalphau,stdllrespRes,stdsum_data_alpha,stdsum_data_alphau,P,k);
			}

			//compute z1*lnull+z2*lresp+prior
			completeLL(z,stdllnullRes,stdllrespRes,cll,filter,P,k);
			oldll=std::accumulate(cll.begin(),cll.end(),0.0)+prior;


			//simulate alphas_j
			stdnextalphavec[j]=alphaProposal(stdalphas,sigmas[j]*rateS[j],j);
			priornext=::Rf_dexp(stdnextalphavec[j],EXPRATE,true);

			//don't need to recompute the null marginal log likelihood since it doesn't depend on alphas_j
			if(FILTER&k==2&!FAST){
				//compute one sided marginal log likelihood
				normalizingConstant(stdstim,stdunstim,stdnextalphavec,stdalphau,stdllrespResNew,P,k);
			}else{
				//two sided
				loglikeresp(stdstim,stdnextalphavec,stdunstim,stdalphau,stdllrespResNew,stdsum_data_alpha,stdsum_data_alphau,P,k);
			}
			//compute z1*lnull+z2*lresp+prior
			completeLL(z,stdllnullRes,stdllrespResNew,cll,filter,P,k);
			newll=std::accumulate(cll.begin(),cll.end(),0.0)+priornext;




			if(stdnextalphavec[j]>0&&(::log(Rf_runif(0.0,1.0)) <= (newll-oldll) )&&(!ISNAN(newll-oldll))){
#ifdef NDEBUG
				printf("ACCEPTED alphas_%d %f prob: %f newll %f oldll %f\n",j,stdnextalphavec[j],::exp(newll-oldll), newll, oldll);
#endif
				accepts[j]=accepts[j]+1;//increment acceptance count
				stdalphas[j]=stdnextalphavec[j]; //new alphas_j is accepted
				oldll=newll-priornext; //save the new complete data log likelihood (minus the prior for alphas_j) so we don't recompute it for the next step
				std::copy(stdllnullResNew.begin(),stdllnullResNew.end(),stdllnullRes.begin()); //ditto for the null and alternative marginal log likelihood
				std::copy(stdllrespResNew.begin(),stdllrespResNew.end(),stdllrespRes.begin());
			}else{
				oldll=oldll-prior; //reject so just subtract the alpha-specific prior
#ifdef NDEBUG
				printf("REJECTED alphas_%d %f, deltall: %f newll %f oldll %f\n",j,stdnextalphavec[j],(newll-oldll),newll, oldll);
#endif
			}



			/*
			 * Simulate alphau_j
			 */
			//prior for the current alphau_j
			prior=Rf_dexp(stdalphau[j],EXPRATE,true);
			oldll=oldll+prior;
			//copy the alphau vector to the proposal vector.
			std::copy(stdalphau.begin(),stdalphau.end(),stdnextalphavec.begin());//copy the current alpha vector to the new alpha vector prior to drawing a sample

			//simulate alphau)j
			stdnextalphavec[j]=alphaProposal(stdalphau,sigmau[j]*rateU[j],j);
			priornext=Rf_dexp(stdnextalphavec[j],EXPRATE,true);

			//don't need to recompute the current complete data log likelihood. It's the same as before, just differs by the prior.

			//compute z1*lnull+z2*lresp+prior
			loglikenull(stdsum_stim_unstim,stdnextalphavec,stdllnullResNew,stdsum_data_alphau,P,k); //new null marginal likelihood.
			if(FILTER&k==2&!FAST){
				normalizingConstant(stdstim,stdunstim,stdalphas,stdnextalphavec,stdllrespResNew,P,k); //new responder marginal LL - one sided
			}else{
				//two sided
				loglikeresp(stdstim,stdalphas,stdunstim,stdnextalphavec,stdllrespResNew,stdsum_data_alpha,stdsum_data_alphau,P,k);
			}
			//compute z1*lnull+z2*lresp+prior
			completeLL(z,stdllnullResNew,stdllrespResNew,cll,filter,P,k);
			newll=std::accumulate(cll.begin(),cll.end(),0.0)+priornext;



			if(stdnextalphavec[j]&&(log(Rf_runif(0.0,1.0)) <= (newll-oldll) )&&(!ISNAN(newll-oldll))){
				//increment acceptance count for alphauj
#ifdef NDEBUG
				printf("ACCEPTED alphau_%d %f, prob ratio %f newll %f oldll %f\n",j,stdnextalphavec[j],::exp(newll-oldll),newll, oldll);
#endif
				acceptu[j]=acceptu[j]+1;
				stdalphau[j]=stdnextalphavec[j];//new alphau_j is accepted
				oldll=newll-priornext; //complete data log likelihood (minus the prior)
				//marginal null and alternative log likelihoods for the accepted parameter are saved so we don't have to recompute them
				std::copy(stdllnullResNew.begin(),stdllnullResNew.end(),stdllnullRes.begin());
				std::copy(stdllrespResNew.begin(),stdllrespResNew.end(),stdllrespRes.begin());
			}else{
				oldll=oldll-prior;
#ifdef NDEBUG
				printf("REJECTED alphau_%d %f, deltall: %f newll %f oldll %f\n",j,stdnextalphavec[j],(newll-oldll),newll, oldll);
#endif
			}
			//simulate q (w)
			q=simQ(z,P,k);
			//simulate z
			simZ(q,stdllnullRes,stdllrespRes,z,p,filter,P,k); //overwrites the current z. A running average is stored in cz

		}

		/*
		 * If we haven't fixed the step sizes yet..
		 */
		if(!fixed){
			/*
			 * Tuning phase
			 */
			//Fill in the Ms and Mu matrices
			asi = conv_to<arma::vec>::from(stdalphas);
			aui = conv_to<arma::vec>::from(stdalphau);

			Ms.row(iteration%(int)TUNING) = asi;
			Mu.row(iteration%(int)TUNING) = aui;
			//Tuning
			if(((iteration+1) % (int)TUNING)==0){
				//Compute the covariances
				arma::mat covarianceS = arma::cov(Ms);
				arma::mat covarianceU = arma::cov(Mu);
				arma::vec dS = arma::sqrt(covarianceS.diag());
				arma::vec dU = arma::sqrt(covarianceU.diag());
				//Weight and assign
				//Tweak the acceptance rates
				accepts=accepts/TUNING;
				acceptu=acceptu/TUNING;

				if((Rcpp::any(accepts > UPPER).is_true() || Rcpp::any(acceptu > UPPER).is_true() || Rcpp::any(accepts < LOWER).is_true() || Rcpp::any(acceptu < LOWER).is_true())){
					for(j=0;j<accepts.length();j++){
						//stimulated
						if((accepts[j] ) > UPPER || (accepts[j] ) < LOWER){
							//rateS[j]=accepts[j]/DEFAULT_RATE;
						}
						//unstimulated
						if((acceptu[j] ) > UPPER || (acceptu[j] ) < LOWER){
							//rateU[j]=acceptu[j]/DEFAULT_RATE;
						}
					}
					for(j=0;j<sigmas.length();j++){
						sigmas[j] = dS(j);
						sigmau[j] = dU(j);
					}

					for(j=0;j<k;j++){
						if(rateS[j]==0){
							R_CheckUserInterrupt();
							rateS[j]=1;
							accepts[j]=0;
							iteration=0;
						}
						if(ISNAN(sigmas[j])){
							R_CheckUserInterrupt();
							sigmas[j]=1;
							accepts[j]=0;
							iteration=0;
						}
						if(rateU[j]==0){
							R_CheckUserInterrupt();
							rateU[j]=1;
							iteration=0;
							acceptu[j]=0;
						}if(ISNAN(sigmau[j])){
							R_CheckUserInterrupt();
							sigmau[j]=1;
							iteration=0;
							acceptu[j]=0;
						}
					}

				}else{
					printf("Fixed step size\n");
					printf("sigmas: ");
					for(int j=0;j<sigmas.length();j++){
						printf("%f %f ",sigmas[j],sigmau[j]);
					}
					printf("\n Acceptance rates:");
					for(int j = 0; j < accepts.length(); j++){
						printf("%f %f",accepts[j],acceptu[j]);
					}
					printf("\n");
					fixed=true;
					iteration = 0;
					//write out the headers for the two data files
					fprintf(file,"%s",headers.str().data());
					fprintf(fileP,"%s",headersP.str().data());
				}
				//reset the acceptance counts
				accepts.fill(0);
				acceptu.fill(0);
			}
		}
		/*
		 * Running average of z1's after the burn-in
		 */

		if((iteration+1)%10000==0){
			R_CheckUserInterrupt();
			if(iteration==0){
				std::copy(p.begin(),p.end(),cz.begin());
			}
			printf("--- Done %i iterations ---\n",(int)iteration+1);
			printf("\n Acceptance rates:");
			for(int j = 0; j < accepts.length(); j++){
				printf("%f %f",accepts[j],acceptu[j]);
			}
			printf("\n");
		}

		if(iteration>=BURNIN&&iteration%(int)THINNING==0&&fixed){
			realitcounter++;
			for(int j=0;j<p.size();j++){
				double f = realitcounter/(realitcounter+1.0);
				double foo = p[j]/realitcounter;
				cz[j]=(cz[j]+foo)*f;
			}
			//write chain to file
			//write out the z's
			//when the model is two sided and two dimensional..
			if(k==2&!FILTER){
				//and the proportions sampled from each model
				sampleP(stdsum_stim_unstim,stdstim,stdunstim,stdalphas,stdalphau,z,ps,pu,P,k);
				for(int obs=0;obs<P;obs++){
					fprintf(fileP,"%f\t", z[obs]);
				}
				for(int obs=0;obs<P-1;obs++){
					fprintf(fileP,"%f\t%f\t",ps[obs],pu[obs]);
				}
				fprintf(fileP,"%f\t%f\n",ps[P-1],pu[P-1]);
			}
			for(int obs=0;obs<stdalphas.size();obs++){
				fprintf(file,"%f\t", stdalphas[obs]);
			}
			for(int obs=0;obs<stdalphau.size();obs++){
				fprintf(file,"%f\t", stdalphau[obs]);
			}
			fprintf(file,"%f\n",q);
		}
#ifdef NDEBUG
		printf("alphas: %f %f alphau: %f %f\n",stdalphas[0],stdalphas[1],stdalphau[0],stdalphau[1]);
#endif
	}
	if(!fixed){
		cout<<"Failed to set step size. Run a longer chain.\n";
	}

	/*
	 * Close the file and return some stuff to R.
	 */
	fflush(file);fflush(fileP);
	fclose(file);fclose(fileP);
	return Rcpp::List::create(
			Rcpp::Named("z") = cz,
			Rcpp::Named("stepsizeS") = sigmas,
			Rcpp::Named("stepsizeU") = sigmau);
	END_RCPP
}

/*
 * Null component log likelihood
 * data is the data, alpha is the parameters, output is the result, sum_dat_alphau is data+alpha
 */
void loglikenull(const std::vector<double> &data,const std::vector<double>  &alpha,std::vector<double> &output, std::vector<double> &sum_dat_alphau,int P, int k){
	int i=0,j=0;
	double da=0,a=0;
		a=lkbeta(alpha);
		for(i=0;i<P;i++){
			for(j=0;j<k;j++){
				sum_dat_alphau[j*P+i]=data[j*P+i]+alpha[j];
			}
			da=lkbeta(sum_dat_alphau,i,k,P);
			output[i]=da-a;
		}

}
/*
 * Responder component log-likelihood
 */
void loglikeresp(const std::vector<double>  &stim, const std::vector<double>  &alphas,const  std::vector<double>  &unstim, const std::vector<double>  &alphau,std::vector<double> &output, std::vector<double> &sum_dat_alphas,std::vector<double> &sum_dat_alphau,int P, int k){
	int i=0,j=0;
	double da,db,a,b;
	b=lkbeta(alphau);
	a=lkbeta(alphas);
	for(i=0;i<P;i++){
		for(j=0;j<k;j++){
			sum_dat_alphas[j*P+i]=stim[j*P+i]+alphas[j];
			sum_dat_alphau[j*P+i]=unstim[j*P+i]+alphau[j];
		}
		da=lkbeta(sum_dat_alphas,i,k,P);
		db=lkbeta(sum_dat_alphau,i,k,P);
		output[i]=da+db-a-b;
	}
}


/*
 * K-dimensional Beta function
 */
inline double lkbeta(const std::vector<double>& alpha,int I,int k,int P){
	double sum_alpha=0;
	double sum_log_gamma_alpha=0;
	double log_gamma_sum_alpha=0;
	for(int j = 0;j<k;j++){
		sum_alpha=sum_alpha+alpha[I+j*P];
		sum_log_gamma_alpha = sum_log_gamma_alpha+lgamma(alpha[j*P+I]);
	}
	log_gamma_sum_alpha = lgamma(sum_alpha);
	return sum_log_gamma_alpha-log_gamma_sum_alpha;
}

inline double lkbeta(const std::vector<double> &alpha){
	double sum_alpha = std::accumulate(alpha.begin(),alpha.end(),0.0);
	double log_gamma_sum_alpha = lgamma(sum_alpha);
	double sum_log_gamma_alpha=0;
	for(int i=0;i<alpha.size();i++){
		sum_log_gamma_alpha=sum_log_gamma_alpha+lgamma(alpha[i]);
	}
	return sum_log_gamma_alpha-log_gamma_sum_alpha;
}


/*
 * log-gamma function for use with std::transform
 */
double op_lgamma(double i){
	return lgamma(i);
}

/*
 * Draw a proposal for the ith component of an alpha vector
 */
double alphaProposal(const std::vector<double> &alpha, double sigma, int i){
	double na;
	na = ::Rf_rnorm(alpha[i],sigma);
	return na;
}

/*
 * Compute the complete data log-likelihood
 * If FILTER is true, we check each index against the value of filter, and set the posterior probability to zero for FILTER_j = true
 */
void completeLL(std::vector<double> &z,std::vector<double> &lnull, std::vector<double> &lresp,std::vector<double> &cll,std::vector<bool> &filter,int P, int k){
	int i;
	for(i=0;i< P;i++){
		if(FAST&filter[i]){
			z[i+P]=0.0;z[i]=1.0;
		}
		cll[i] = z[i]*lnull[i]+z[i+P]*(lresp[i]);
	}
}
void simZ(double &q,std::vector<double> &lnull, std::vector<double> &lresp,std::vector<double>& z,std::vector<double> &p,std::vector<bool> &filter,int P, int k){
	int i;
	double lq = ::log(q);
	double mlq = ::log(1.0-q);
	for(i=0;i < lnull.size(); i++){
		lnull[i]=lnull[i]+lq;
		lresp[i]=lresp[i]+mlq;
		double mx=std::max(lnull[i],lresp[i]);
		//printf("null:%f\talternative: %f\n",lnull[i],lresp[i]);
		if(FAST&filter[i]){
			p[i]=1;
		}else{
			p[i] = ::exp(lnull[i]-::log(::exp(lnull[i]-mx)+::exp(lresp[i]-mx))-mx);
			z[i] = ::Rf_rbinom(1.0,p[i]);
			z[i+P] = 1.0-z[i];
		}
	}
}
inline double simQ(std::vector<double> &z, int P,int k){
	std::vector<double> ab(2,0);
	double q;
	for(int j=0;j<2;j++){
		for(int i=0;i<P;i++){
			ab[j]=ab[j]+z[j*P+i];
		}
	}
	q = 1.0-::Rf_rbeta(ab[1]+1,ab[0]+1);
	return q;
}
void normalizingConstant(std::vector<double> &stim,std::vector<double> &unstim,std::vector<double> &alphas,std::vector<double> &alphau,std::vector<double> &llresp, int P,int k){
	assert(k==2);
	double numerator=0,denominator=0,nummc,denommc;
	double C=1,CC=1;
	//If any alphas are <= 0 fill with nan;
	std::vector<double> u(2,0), s(2,0);
	int i=0,j=0;
	for(i=0;i<P;i++){
		for(j=0;j<2;j++){
			s[j]=stim[i+j*P]+alphas[j];
			u[j]=unstim[i+j*P]+alphau[j];
		}
		//alphas[1] is alpha, alphas[0] is beta
		//data+hyperparameters
		numerator=normconstIBeta((double)s[0],(double)s[1],(double)u[0],(double)u[1]);
		//hyperparameters only
		denominator=normconstIBeta((double)alphas[0],(double)alphas[1],(double)alphau[0],(double)alphau[1]);
//		nummc = normconstMC((double)s[1],(double)s[0],(double)u[1],(double)u[0]);
//		denommc = normconstMC((double)alphas[1],(double)alphas[0],(double)alphau[1],(double)alphau[1]);

		CC=numerator-denominator;
		CC=log(1-exp(CC))+1;
//		C=nummc-denommc;
		//		printf("C=%f\n",::exp(C));
		double K=lgamma((double)s[1])+lgamma((double)s[0])-lgamma(s[1]+s[0])+lgamma((double)u[1])+lgamma((double)u[0])-lgamma(u[1]+u[0])-lgamma((double)alphas[1])-lgamma((double)alphas[0])+lgamma(alphas[1]+alphas[0])-lgamma((double)alphau[1])-lgamma((double)alphau[0])+lgamma(alphau[1]+alphau[0]);
//		double K = ::Rf_lbeta((double) s[1],(double) s[0])+::Rf_lbeta((double) u[1], (double) u[0])-::Rf_lbeta((double) alphas[1],(double)alphas[0])-::Rf_lbeta((double)alphau[1],(double)alphau[0]);
		if(ISNAN(CC)){
//			printf("s0: %f s1: %f u0: %f u1: %f as0: %f as1: %f au0: %f au1: %f \n",s[0],s[1],u[0],u[1],alphas[0],alphas[1],alphau[0],alphau[1]);
//			printf("numerator: %f  denominator %f  C: %f\n",numerator,denominator,numerator-denominator);
//			printf("log(1-exp(C))=%f\n",log(1-exp(numerator-denominator)));
			nummc = normconstMC((double)s[1],(double)s[0],(double)u[1],(double)u[0]);
			denommc = normconstMC((double)alphas[1],(double)alphas[0],(double)alphau[1],(double)alphau[1]);
			C=nummc-denommc;
			CC=C;
		}
//		printf("mc: %f  approx: %f, diff: %f\n",K+CC,K+C,CC-C);
		llresp[i]=K+C;
		//printf("%f\n",llresp[i]);
	}
}

double normconstMC(double as, double bs, double au, double bu){
	double res;
	NumericVector r = rbeta(50,au,bu);
	r=pbeta(r,as,bs,false,false);
	res = std::accumulate(r.begin(),r.end(),0.0)/r.length();
	return res;
}

double normconstIBeta(double as, double bs, double au, double bu){
	double alphas=(double)ceil(as);
	double betas = (double)ceil(bs);
	double alphau = (double)ceil(au);
	double betau = (double)ceil(bu);
	if(alphas<=0||betas<=0||alphau<=0||betau<=0){
		return 0.0/0.0; //if any parameters are negative return nan.
	}
	double sum=0,mx=0;
	double upper = (double) (alphau+betau);
	std::vector<double> res((int)betau,0.0);
	//	printf("as=%f bs=%f au=%f bu=%f\n",as,bs,au,bu);
	//	printf("INTS: as=%f bs=%f au=%f bu=%f\n",alphas,betas,alphau,betau);
	//	printf("upper: %d, alphau: %d size: %d\n",upper, alphau,res.size());
	double K = -::Rf_lbeta(alphas,betas)+Rf_lgammafn(alphau+betau)-Rf_lgammafn(alphas+betas+alphau+betau-1);
#ifdef FOO
	printf("upper=%f\n",upper);
	printf("K=%f\n",K);
	printf("-::Rf_lbeta(alphas,betas)=%f\n",-::Rf_lbeta(alphas,betas));
	printf("Rf_lgammafn(alphau+betau)=%f\n",Rf_lgammafn(alphau+betau));
	printf("-Rf_lgammafn(alphas+betas+alphau-1)=%f\n",-Rf_lgammafn(alphas+betas+alphau-1));
#endif
	for(int j = (int)alphau;j< ((int)upper);j++){
#ifdef FOO
		printf("j=%d\n",j);
		printf("Rf_lgammafn(alphas+j)=%f\n",Rf_lgammafn(alphas+j));
		printf("Rf_lgammafn(alphau+betau+betas-1)=%f\n",Rf_lgammafn(alphau+betau+betas-j-1));
		printf("-Rf_lgammafn(j+1)=%f\n",-Rf_lgammafn(j+1));
		printf("-Rf_lgammafn(alphau+betau-j)=%f\n",-Rf_lgammafn(alphau+betau-j));
		printf("sm = %f\n",K+Rf_lgammafn(alphas+j)+Rf_lgammafn(alphau+betau+betas-j-1)-Rf_lgammafn(j+1)-Rf_lgammafn(alphau+betau-j));
#endif
		res[j-(int)alphau]=K+Rf_lgammafn(alphas+j)+Rf_lgammafn(alphau+betau+betas-j-1)-Rf_lgammafn(j+1)-Rf_lgammafn(alphau+betau-j);
	}

	//todo normalized sum of exponentials
	std::vector<double>::iterator where=std::max_element(res.begin(),res.end());
	mx = *where;
	for(int i=0;i<res.size();i++){
		sum=sum+::exp(res[i]-mx);
	}
	sum=log(sum)+mx;
	//printf("%f \n",sum);
	return(sum);
}


//samples P's for the 2-d case only.
void sampleP(std::vector<double>& sumdata,std::vector<double>& stim,std::vector<double>& unstim,std::vector<double>& alphas,std::vector<double>& alphau,std::vector<double>& z, std::vector<double> &ps, std::vector<double> &pu, int P,int k){
	for(int i=0;i<P;i++){
		if(z[i+P]==0){
			//sample from the null model
			ps[i]=Rf_rbeta(sumdata[i+1*P]+alphas[1]+alphau[1],sumdata[i+0*P]+alphas[0]+alphau[0]);
			pu[i]=ps[i];
		}else{
			//otherwise sample from the responder model
			ps[i]=Rf_rbeta(stim[i+1*P]+alphas[1],stim[i+0*P]+alphas[0]);
			pu[i]=Rf_rbeta(unstim[i+1*P]+alphau[1],unstim[i+0*P]+alphau[0]);
		}
	}
}

double nc(double as, double bs, double au,double bu,double B){
	double K,mx,sm=0;
	std::vector<double> s(B+1,0);
	K=::Rf_lbeta(au+as,bu+bs)-::log(au);
	s[0]=K;
	for(int i=0;i<(s.size()-1);i++){
		s[i+1]=(Rf_lbeta(au+1,i+1)+Rf_lbeta(au+as+i+1,bu+bs)-Rf_lbeta(au+bu,i+1)-log(au));
	}
	std::vector<double>::iterator where=std::max_element(s.begin(),s.end());
	mx=*where;
	for(int i=0;i<s.size();i++){
		s[i]=::exp(s[i]-mx);
		sm=sm+s[i];
	}
	sm=log(sm)+mx;
	return(sm);
}


