#include <Rcpp.h>
#include <armadillo>
#include <RcppArmadillo.h>
#include <assert.h>
#include <omp.h>
#include <Rmath.h>
#include <iostream>
#include <fstream>
#include <R_ext/Utils.h>
#include "MIMOSA.h"
#undef NDEBUG

/*
 * 10 parameters
 */
RcppExport SEXP fitMCMC(SEXP _stim, SEXP _unstim, SEXP _alphas, SEXP _alphau, SEXP _q, SEXP _z,SEXP _iter, SEXP _burn, SEXP _thin, SEXP _tune,SEXP _outfile, SEXP _filter){
	BEGIN_RCPP
	using namespace Rcpp;
	using namespace arma;
	using namespace std;

	bool fixed = false;

	Rcpp::RNGScope globalscope;

	/*
	 * Copy R variables to standard vectors
	 */
//	std::vector<bool> stdfilter = Rcpp::as<std::vector<bool> >(_filter);
//	std::vector<double> stdalphas = Rcpp::as<std::vector<double> >(_alphas);
//	std::vector<double> stdalphau = Rcpp::as<std::vector<double> 	>(_alphau);
//	double stdq = Rcpp::as<double>(_q);
//	double stditer = Rcpp::as<double>(_z);
//	double stdburn = Rcpp::as<double>(_burn);
//	double stdthin = Rcpp::as<double>(_thin);
//	double stdtune = Rcpp::as<double>(_tune);


	/*
	 * Wrap R variables in Rcpp objects
	 */
	Rcpp::NumericMatrix const stim(_stim);
	Rcpp::NumericMatrix const unstim(_unstim);
	Rcpp::NumericVector alphas(_alphas);
	Rcpp::NumericVector alphau(_alphau);
	Rcpp::NumericVector q(_q);
	Rcpp::NumericMatrix z(_z);
	Rcpp::NumericVector const iter(_iter);
	Rcpp::NumericVector const burn(_burn);
	Rcpp::NumericVector const thin(_thin);
	Rcpp::NumericVector const tune(_tune);
	Rcpp::LogicalVector filter(_filter);
	std::string outfile = Rcpp::as<std::string>(_outfile);
	//normalizing constant for the alternative in the one-sided 2x2 case.
	std::vector<double> normconst(0,z.nrows());

	printf("Creating %s\n",outfile.data());
	FILE* file = fopen(outfile.data(),"w");

	if(file==NULL){
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
	assert(alphas.length()==stim.ncol());
	assert(stim.nrow()==unstim.nrow());
	assert(stim.ncol()==unstim.ncol());
	assert(alphas.length()==alphau.length());

	if(any(filter).is_true()){
		FILTER=true;
	}

	/*
	 * Dimensions of the problem
	 */
	const int k = stim.ncol();
	const int P = stim.nrow();

	/*
	 * Output variables
	 */
	NumericVector llnullRes(P); //null component log likelihood
	NumericVector llrespRes(P);// responder component log likelihood
	NumericMatrix sum_data_alpha(P,k); //set aside memory for sum(data+alpha)
	NumericMatrix sum_data_alphau(P,k);
	arma::vec asi, aui;
	/*
	 * Armadillo objects will store the simulations from the tuning phase
	 * so we can easily compute variances and so forth.
	 */
	arma::mat Ms(int(TUNING),alphas.length());
	arma::mat Mu(int(TUNING),alphau.length());

	/*
	 * Precompute and preallocate a few things we'll need
	 */
	NumericMatrix sum_stim_unstim(P,k,(stim+unstim).begin()); //memory for sum of stim and unstim data
	NumericVector nextalphavec(alphas.length());
	NumericVector sigmas(alphas.length(),10.0);
	NumericVector sigmau(alphau.length(),10.0);
	arma::vec rateS(alphas.length());
	arma::vec rateU(alphau.length());
	rateS.fill(RATE);
	rateU.fill(RATE);
	NumericVector accepts(alphas.length(),0.0);
	NumericVector acceptu(alphas.length(),0.0);
	NumericVector cll(stim.nrow(),0.0);
	NumericVector p(z.nrow(),0.0);
	NumericVector cz(z.nrow(),0.0);
	//std::copy((z.column(0)).begin(),(z.column(0)).end(),cz.begin());


	double prior=0,priornext=0;
	double oldll=0,newll=0;


	/*
	 * Construct the header string
	 *
	 */
	std::stringstream headers(stringstream::in|stringstream::out);
	for(int i=0;i<z.nrow();i++){
		headers<<"z."<<i<<"\t";
	}
	for(int i=0;i<alphas.length();i++){
		headers<<"alphas."<<i<<"\t";
	}
	for(int i=0;i<alphau.length();i++){
		headers<<"alphau."<<i<<"\t";
	}
	headers<<"q"<<std::endl;
	//counters
	int iteration=0,j=0;
	/*
	 * Run the MCMC algorithm
	 */
	for(iteration = 0; iteration < NITERS; iteration++){
		for(j=0;j<k;j++){
			/*
			 * Simulate one alphas step
			 */
			//prior logdexp(0.0001)
			prior=::Rf_dexp(alphas[j],10000,true);

			std::copy(alphas.begin(),alphas.end(),nextalphavec.begin());


			nextalphavec[j]=alphaProposal(alphas,sigmas[j]*rateS[j],j);
			priornext=::Rf_dexp(nextalphavec[j],10000,true);

			loglikenull(sum_stim_unstim,alphau,llnullRes,sum_data_alphau);
			loglikeresp(stim,alphas,unstim,alphau,llrespRes,sum_data_alpha,sum_data_alphau);

			//compute z1*lnull+z2*lresp+prior
			completeLL(z,llnullRes,llrespRes,cll,filter);
			oldll=std::accumulate(cll.begin(),cll.end(),0.0)+prior;

			loglikenull(sum_stim_unstim,alphau,llnullRes,sum_data_alphau);
			loglikeresp(stim,nextalphavec,unstim,alphau,llrespRes,sum_data_alpha,sum_data_alphau);

			//compute z1*lnull+z2*lresp+prior
			completeLL(z,llnullRes,llrespRes,cll,filter);
			newll=std::accumulate(cll.begin(),cll.end(),0.0)+priornext;

			if(all(nextalphavec>0).is_true()&&(::log(Rf_runif(0.0,1.0)) <= (newll-oldll) )&&(!ISNAN(newll-oldll))){
				//increment acceptance count for alphasj
#ifdef NDEBUG
				printf("ACCEPTED alphas_%d %f prob: %f newll %f oldll %f\n",j,nextalphavec[j],::exp(newll-oldll), newll, oldll);
#endif
				accepts[j]=accepts[j]+1;
				alphas[j]=nextalphavec[j];
			}else{
#ifdef NDEBUG
				printf("REJECTED alphas_%d %f prob %f newll %f oldll %f\n",j,nextalphavec[j],::exp(newll-oldll),newll,oldll);
#endif
			}



			/*
			 * Simulate one alphau step
			 */
			//prior logexp(0.0001)
			prior=Rf_dexp(alphau[j],10000,true);
			std::copy(alphau.begin(),alphau.end(),nextalphavec.begin());//copy the current alpha vector to the new alpha vector prior to drawing a sample

			//TODO test for valid sigma and valid acceptance rate
			nextalphavec[j]=alphaProposal(alphau,sigmau[j]*rateU[j],j);
			priornext=Rf_dexp(nextalphavec[j],10000,true);

			/*
			 * Compute log likelihood
			 */
			loglikenull(sum_stim_unstim,alphau,llnullRes,sum_data_alphau);
			loglikeresp(stim,alphas,unstim,alphau,llrespRes,sum_data_alpha,sum_data_alphau);
			//compute z1*lnull+z2*lresp+prior
			completeLL(z,llnullRes,llrespRes,cll,filter);
			oldll=std::accumulate(cll.begin(),cll.end(),0.0)+prior;
			loglikenull(sum_stim_unstim,nextalphavec,llnullRes,sum_data_alphau);
			loglikeresp(stim,alphas,unstim,nextalphavec,llrespRes,sum_data_alpha,sum_data_alphau);

			//compute z1*lnull+z2*lresp+prior
			completeLL(z,llnullRes,llrespRes,cll,filter);
			newll=std::accumulate(cll.begin(),cll.end(),0.0)+priornext;

			if(all(nextalphavec>0).is_true()&&(log(Rf_runif(0.0,1.0)) <= (newll-oldll) )&&(!ISNAN(newll-oldll))){
				//increment acceptance count for alphauj
#ifdef NDEBUG
				printf("ACCEPTED alphau_%d %f, prob ratio %f newll %f oldll %f\n",j,nextalphavec[j],::exp(newll-oldll),newll, oldll);
#endif
				acceptu[j]=acceptu[j]+1;
				alphau[j]=nextalphavec[j];
			}else{
#ifdef NDEBUG
				printf("REJECTED alphau_%d %f, deltall: %f newll %f oldll %f\n",j,nextalphavec[j],(newll-oldll),newll, oldll);
#endif
			}

			//simulate q
			q[0]=simQ(z);

			//simulate z
			loglikenull(sum_stim_unstim,alphau,llnullRes,sum_data_alphau);
			loglikeresp(stim,alphas,unstim,alphau,llrespRes,sum_data_alpha,sum_data_alphau);
			simZ(q,llnullRes,llrespRes,z,p,filter); //overwrites the current z
		}

		/*
		 * If we haven't fixed the step sizes yet..
		 */
		if(!fixed){
			/*
			 * Tuning phase
			 */
			//Fill in the Ms and Mu matrices
			asi = Rcpp::as<arma::vec>(alphas.asSexp());
			aui = Rcpp::as<arma::vec>(alphau.asSexp());
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
				//#ifdef DNDEBUG
				printf("acceptance rates  %f %f %f %f\t%f %f %f %f\n",accepts[0],accepts[1],accepts[2],accepts[3],acceptu[0],acceptu[1],acceptu[2],acceptu[3]);
				printf("acceptance ratios %f %f %f %f\t%f %f %f %f\n",rateS[0],rateS[1],rateS[2],rateS[3],rateU[0],rateU[1],rateU[2],rateU[3]);
				printf("sigmas            %f %f %f %f\t%f %f %f %f\n",sigmas[0],sigmas[1],sigmas[2],sigmas[3],sigmau[0],sigmau[1],sigmau[2],sigmau[3]);
				printf("alphas            %f %f %f %f\t%f %f %f %f\n\n\n",alphas[0],alphas[1],alphas[2],alphas[3],alphau[0],alphau[1],alphau[2],alphau[3]);
				//#endif
				if((Rcpp::any(accepts > UPPER).is_true() || Rcpp::any(acceptu > UPPER).is_true() || Rcpp::any(accepts < LOWER).is_true() || Rcpp::any(acceptu < LOWER).is_true())){
					for(j=0;j<accepts.length();j++){
						//stimulated
						if((accepts[j] ) > UPPER || (accepts[j] ) < LOWER){
							rateS[j]=accepts[j]/DEFAULT_RATE;
						}
						//unstimulated
						if((acceptu[j] ) > UPPER || (acceptu[j] ) < LOWER){
							rateU[j]=acceptu[j]/DEFAULT_RATE;
						}
					}
					for(j=0;j<sigmas.length();j++){
						sigmas[j] = sigmas[j]*0.5+0.5*dS(j);
						sigmau[j] = sigmau[j]*0.5+0.5*dU(j);
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
					fprintf(file,"%s",headers.str().data());
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
		}

		if(iteration>=BURNIN&&iteration%(int)THINNING==0&&fixed){
			realitcounter++;
			for(int j=0;j<p.length();j++){
				double f = realitcounter/(realitcounter+1.0);
				double foo = p[j]/realitcounter;
				cz[j]=(cz[j]+foo)*f;
			}
			//write chain to file
			for(int obs=0;obs<z.nrow();obs++){
				fprintf(file,"%f\t", z(obs,0));
			}
			for(int obs=0;obs<alphas.length();obs++){
				fprintf(file,"%f\t", alphas(obs));
			}
			for(int obs=0;obs<alphau.length();obs++){
				fprintf(file,"%f\t", alphau(obs));
			}
			fprintf(file,"%f\n",q[0]);
		}
#ifdef NDEBUG
		printf("alphas: %f %f alphau: %f %f\n",alphas[0],alphas[1],alphau[0],alphau[1]);
#endif
	}
	if(!fixed){
		cout<<"Failed to set step size. Run a longer chain.\n";
	}

	/*
	 * Close the file and return some stuff to R.
	 */
	fflush(file);
	fclose(file);
	//TODO don't return these here.. either read the mcmc results and compute the mean or compute the mean inside the C code and return that.. these are just the last samples taken. Not what we want.
	return Rcpp::List::create(Rcpp::Named("q") = q,
			Rcpp::Named("z") = cz,
			Rcpp::Named("alpha.u") = alphau,
			Rcpp::Named("alpha.s") = alphas,
			Rcpp::Named("stepsizeS") = sigmas,
			Rcpp::Named("stepsizeU") = sigmau);
	END_RCPP
}

/*
 * Null component log likelihood
 * data is the data, alpha is the parameters, output is the result, sum_dat_alphau is data+alpha
 */
void loglikenull(const NumericMatrix &data,const NumericVector  &alpha,NumericVector &output, NumericMatrix &sum_dat_alphau){
	int i=0,j=0;
	double da=0,a=0;
	a=lkbeta(alpha);
	//#pragma omp parallel for firstprivate(j) private(da) lastprivate(i)
	for(i=0;i<data.nrow();i++){
		for(j=0;j<data.ncol();j++){
			sum_dat_alphau(i,j)=data(i,j)+alpha[j];
		}
		da=lkbeta(sum_dat_alphau.row(i));
		output[i]=da-a;
	}
	//printf("Thread started on obeservation %i",i);
}
/*
 * Responder component log-likelihood
 */
void loglikeresp(const NumericMatrix  &stim, const NumericVector  &alphas,const  NumericVector  &unstim, const NumericVector  &alphau,NumericVector &output, NumericMatrix &sum_dat_alphas,NumericMatrix &sum_dat_alphau){
	int i,j;
	double da,db,a,b;
	b=lkbeta(alphau);
	a=lkbeta(alphas);
	//#pragma omp parallel for firstprivate(j) private(da,db) lastprivate(i)
	for(i=0;i<stim.nrow();i++){
		for(j=0;j<stim.ncol();j++){
			sum_dat_alphas(i,j)=stim(i,j)+alphas[j];
			sum_dat_alphau(i,j)=unstim(i,j)+alphau[j];
		}
		da=lkbeta(sum_dat_alphas.row(i));
		db=lkbeta(sum_dat_alphau.row(i));
		output(i)=da+db-a-b;
	}

}

/*
 * K-dimensional beta function
 */
//inline double lkbeta(const NumericVector  &alpha,NumericVector &work){
//	double sum_alpha = std::accumulate(alpha.begin(),alpha.end(),0.0);
//	double log_gamma_sum_alpha=lgamma(sum_alpha);
//	std::transform(alpha.begin(),alpha.end(),work.begin(),op_lgamma);
//	double sum_log_gamma_alpha = std::accumulate(work.begin(),work.end(),0.0);
//	return sum_log_gamma_alpha - log_gamma_sum_alpha;
//}

/*
 * K-dimensional Beta function
 */
inline double lkbeta(const NumericVector& alpha){
	double sum_alpha = std::accumulate(alpha.begin(),alpha.end(),0.0);
	double log_gamma_sum_alpha = lgamma(sum_alpha);
	double sum_log_gamma_alpha=0;
	for(int i=0;i<alpha.length();i++){
		sum_log_gamma_alpha=sum_log_gamma_alpha+lgamma(alpha(i));
	}
	return sum_log_gamma_alpha-log_gamma_sum_alpha;
}

inline double lkbeta(const NumericMatrix::Row &alpha){
	double sum_alpha = std::accumulate(alpha.begin(),alpha.end(),0.0);
	double log_gamma_sum_alpha = lgamma(sum_alpha);
	double sum_log_gamma_alpha=0;
	for(int i=0;i<alpha.size();i++){
		sum_log_gamma_alpha=sum_log_gamma_alpha+lgamma(alpha[i]);
	}
	return sum_log_gamma_alpha-log_gamma_sum_alpha;
}

//inline double lkbeta(const NumericMatrix::Row  &alpha,NumericVector &work){
//	double sum_alpha = std::accumulate(alpha.begin(),alpha.end(),0.0);
//	double log_gamma_sum_alpha=lgamma(sum_alpha);
//	std::transform(alpha.begin(),alpha.end(),work.begin(),op_lgamma);
//	double sum_log_gamma_alpha = std::accumulate(work.begin(),work.end(),0.0);
//	return sum_log_gamma_alpha-log_gamma_sum_alpha;
//}

/*
 * log-gamma function for use with std::transform
 */
double op_lgamma(double i){
	return lgamma(i);
}

/*
 * Draw a proposal for the ith component of an alpha vector
 */
double alphaProposal(const NumericVector &alpha, double sigma, int i){
	double na;
	na = ::Rf_rnorm(alpha[i],sigma);
	return na;
}

/*
 * Compute the complete data log-likelihood
 * If FILTER is true, we check each index against the value of filter, and set the posterior probability to zero for FILTER_j = true
 */
void completeLL(NumericMatrix &z,NumericVector &lnull, NumericVector &lresp,NumericVector &cll,LogicalVector &filter){
	int i;
	for(i=0;i< z.nrow();i++){
		if(filter[i]){
			z(i,1)=0.0;z(i,0)=1.0;
		}
		cll[i] = z(i,0)*lnull(i)+z(i,1)*lresp(i);
#ifdef NDEBUG
		printf("%f = %f*%f + %f*%f\n",cll[i],z(i,0),lnull(i),z(i,1),lresp(i));
#endif
	}
}
void simZ(NumericVector &q,NumericVector &lnull, NumericVector &lresp,NumericMatrix& z,NumericVector &p,LogicalVector &filter){
	int i;
	double lq = ::log(q[0]);
	double mlq = ::log(1.0-q[0]);

#ifdef NDEBUG
	printf("Prelim setup in simZ okay\n");
	printf("lq=%f mlq=%f\n",lq,mlq);
#endif

	for(i=0;i < lnull.length(); i++){
		lnull[i]=lnull[i]+lq;
		lresp[i]=lresp[i]+mlq;
		double mx=std::max(lnull[i],lresp[i]);
		p[i] = ::exp(lnull[i]-::log(::exp(lnull[i]-mx)+::exp(lresp[i]-mx))-mx);
		if(!filter[i]){
			z(i,0) = ::Rf_rbinom(1.0,p[i]);
			z(i,1) = 1.0-z(i,0);
		}else{
			p[i]=1;
		}

#ifdef NDEBUG
		printf("z1: %f z2: %f\n",z(i,0),z(i,1));
		printf("Loop %d in simZ okay\n",i);
#endif
	}
}
inline double simQ(NumericMatrix &z){
	double alpha=std::accumulate((z.column(0)).begin(),(z.column(0)).end(),0.0);
	double beta=std::accumulate((z.column(1)).begin(),(z.column(1)).end(),0.0);
	double q = ::Rf_rbeta(alpha+1,beta+1);
	return q;
}
