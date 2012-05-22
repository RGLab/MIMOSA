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
	std::vector<bool> filter = Rcpp::as<std::vector<bool> >(_filter);
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
	//normalizing constant for the alternative in the one-sided 2x2 case.
	std::vector<double> normconst(0,(int)z.size()/2);

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
	assert(stdalphas.size()==stim.ncol());
	assert(stim.nrow()==unstim.nrow());
	assert(stim.ncol()==unstim.ncol());
	assert(stdalphas.size()==stdalphau.size());

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

	std::vector <double> stdllnullRes(0,P);
	std::vector <double> stdllrespRes(0,P);

	std::vector <double> stdsum_data_alpha(0,P*k);
	std::vector <double> stdsum_data_alphau(0,P*k);

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
	std::vector <double> stdsum_stim_unstim(0,P*k);
	std::transform(stdstim.begin(),stdstim.end(),stdunstim.begin(),stdsum_stim_unstim.begin(),plus<double>());



	std::vector stdnextalphavec(0,stdalphas.size());

	NumericVector sigmas(stdalphas.size(),10.0);
	NumericVector sigmau(stdalphau.size(),10.0);

	arma::vec rateS(stdalphas.size());
	arma::vec rateU(stdalphau.size());

	rateS.fill(RATE);
	rateU.fill(RATE);
	NumericVector accepts(stdalphas.size(),0.0);
	NumericVector acceptu(stdalphas.size(),0.0);
	NumericVector cll(P,0.0);
	NumericVector p(P,0.0);
	NumericVector cz(P,0.0);


	double prior=0,priornext=0;
	double oldll=0,newll=0;


	/*
	 * Construct the header string
	 *
	 */
	std::stringstream headers(stringstream::in|stringstream::out);
	//don't write the z's
	//	for(int i=0;i<z.nrow();i++){
	//		headers<<"z."<<i<<"\t";
	//	}
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
	for(iteration = 0; iteration < NITERS; iteration++){
		for(j=0;j<k;j++){
			/*
			 * Simulate one alphas step
			 */
			//prior logdexp(1/10000)
			prior=::Rf_dexp(stdalphas[j],10000,true);

			std::copy(stdalphas.begin(),stdalphas.end(),stdnextalphavec.begin());


			stdnextalphavec[j]=alphaProposal(stdalphas,sigmas[j]*rateS[j],j);
			priornext=::Rf_dexp(stdnextalphavec[j],10000,true);

			loglikenull(stdsum_stim_unstim,stdalphau,stdllnullRes,stdsum_data_alphau,P,k);
			loglikeresp(stdstim,stdalphas,stdunstim,stdalphau,stdllrespRes,stdsum_data_alpha,stdsum_data_alphau,P,k);

			//compute z1*lnull+z2*lresp+prior
			completeLL(z,stdllnullRes,stdllrespRes,cll,stdfilter);
			oldll=std::accumulate(cll.begin(),cll.end(),0.0)+prior;

			loglikenull(stdstdsum_stim_unstim,stdalphau,stdllnullRes,stdsum_data_alphau,P,k);
			loglikeresp(stdstim,stdnextalphavec,stdunstim,stdalphau,stdllrespRes,stdsum_data_alpha,stdsum_data_alphau,P,k);

			//compute z1*lnull+z2*lresp+prior
			completeLL(z,stdllnullRes,stdllrespRes,cll,stdfilter);
			newll=std::accumulate(cll.begin(),cll.end(),0.0)+priornext;

			if(all(stdnextalphavec>0).is_true()&&(::log(Rf_runif(0.0,1.0)) <= (newll-oldll) )&&(!ISNAN(newll-oldll))){
				//increment acceptance count for alphasj
#ifdef NDEBUG
				printf("ACCEPTED alphas_%d %f prob: %f newll %f oldll %f\n",j,stdnextalphavec[j],::exp(newll-oldll), newll, oldll);
#endif
				accepts[j]=accepts[j]+1;
				stdalphas[j]=stdnextalphavec[j];
			}else{
#ifdef NDEBUG
				printf("REJECTED alphas_%d %f prob %f newll %f oldll %f\n",j,stdnextalphavec[j],::exp(newll-oldll),newll,oldll);
#endif
			}



			/*
			 * Simulate one alphau step
			 */
			//prior logexp(0.0001)
			prior=Rf_dexp(stdalphau[j],10000,true);
			std::copy(stdalphau.begin(),stdalphau.end(),stdnextalphavec.begin());//copy the current alpha vector to the new alpha vector prior to drawing a sample

			//TODO test for valid sigma and valid acceptance rate
			stdnextalphavec[j]=alphaProposal(stdalphau,sigmau[j]*rateU[j],j);
			priornext=Rf_dexp(stdnextalphavec[j],10000,true);

			/*
			 * Compute log likelihood
			 */
			loglikenull(stdsum_stim_unstim,stdalphau,stdllnullRes,stdsum_data_alphau,P,k);
			loglikeresp(stdstim,stdalphas,stdunstim,stdalphau,stdllrespRes,stdsum_data_alpha,stdsum_data_alphau,P,k);
			//compute z1*lnull+z2*lresp+prior
			completeLL(z,stdllnullRes,stdllrespRes,cll,stdfilter);
			oldll=std::accumulate(cll.begin(),cll.end(),0.0)+prior;
			loglikenull(stdsum_stim_unstim,stdnextalphavec,stdllnullRes,stdsum_data_alphau,P,k);
			loglikeresp(stdstim,stdalphas,stdunstim,stdnextalphavec,stdllrespRes,stdsum_data_alpha,stdsum_data_alphau,P,k);

			//compute z1*lnull+z2*lresp+prior
			completeLL(z,stdllnullRes,stdllrespRes,cll,stdfilter);
			newll=std::accumulate(cll.begin(),cll.end(),0.0)+priornext;

			if(all(stdnextalphavec>0).is_true()&&(log(Rf_runif(0.0,1.0)) <= (newll-oldll) )&&(!ISNAN(newll-oldll))){
				//increment acceptance count for alphauj
#ifdef NDEBUG
				printf("ACCEPTED alphau_%d %f, prob ratio %f newll %f oldll %f\n",j,stdnextalphavec[j],::exp(newll-oldll),newll, oldll);
#endif
				acceptu[j]=acceptu[j]+1;
				stdalphau[j]=stdnextalphavec[j];
			}else{
#ifdef NDEBUG
				printf("REJECTED alphau_%d %f, deltall: %f newll %f oldll %f\n",j,stdnextalphavec[j],(newll-oldll),newll, oldll);
#endif
			}

			//simulate q
			q[0]=simQ(z);

			//simulate z
			loglikenull(stdsum_stim_unstim,stdalphau,stdllnullRes,stdsum_data_alphau,P,k);
			loglikeresp(stdstim,stdalphas,stdunstim,stdalphau,stdllrespRes,stdsum_data_alpha,stdsum_data_alphau,P,k);
			simZ(q,stdllnullRes,stdllrespRes,z,p,stdfilter); //overwrites the current z
		}

		/*
		 * If we haven't fixed the step sizes yet..
		 */
		if(!fixed){
			/*
			 * Tuning phase
			 */
			//Fill in the Ms and Mu matrices
//			asi = Rcpp::as<arma::vec>(alphas.asSexp());
//			aui = Rcpp::as<arma::vec>(alphau.asSexp());
			asi = conv_to<arma::vec>::from(stdalphas);
			aui = conv_to<arma::vec>>::from(stdalphau);
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
#ifdef DNDEBUG
				printf("acceptance rates  %f %f %f %f\t%f %f %f %f\n",accepts[0],accepts[1],accepts[2],accepts[3],acceptu[0],acceptu[1],acceptu[2],acceptu[3]);
				printf("acceptance ratios %f %f %f %f\t%f %f %f %f\n",rateS[0],rateS[1],rateS[2],rateS[3],rateU[0],rateU[1],rateU[2],rateU[3]);
				printf("sigmas            %f %f %f %f\t%f %f %f %f\n",sigmas[0],sigmas[1],sigmas[2],sigmas[3],sigmau[0],sigmau[1],sigmau[2],sigmau[3]);
				printf("alphas            %f %f %f %f\t%f %f %f %f\n\n\n",alphas[0],alphas[1],alphas[2],alphas[3],alphau[0],alphau[1],alphau[2],alphau[3]);
#endif
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
			//don't write out the z's
			//			for(int obs=0;obs<z.nrow();obs++){
			//				fprintf(file,"%f\t", z(obs,0));
			//			}
			for(int obs=0;obs<stdalphas.size();obs++){
				fprintf(file,"%f\t", stdalphas[obs]);
			}
			for(int obs=0;obs<stdalphau.size();obs++){
				fprintf(file,"%f\t", stdalphau[obs]);
			}
			fprintf(file,"%f\n",q[0]);
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
	fflush(file);
	fclose(file);
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
	//#pragma omp parallel for firstprivate(j) private(da) lastprivate(i)
	for(i=0;i<P;i++){
		for(j=0;j<k;j++){
			sum_dat_alphau[j*P+i]=data[j*P+i]+alpha[j];
		}
		da=lkbeta(sum_dat_alphau,i,k,P);
		output[i]=da-a;
	}
	//printf("Thread started on obeservation %i",i);
}
/*
 * Responder component log-likelihood
 */
void loglikeresp(const std::vector<double>  &stim, const std::vector<double>  &alphas,const  std::vector<double>  &unstim, const std::vector<double>  &alphau,std::vector<double> &output, std::vector<double> &sum_dat_alphas,std::vector<double> &sum_dat_alphau){
	int i=0,j=0;
	double da,db,a,b;
	b=lkbeta(alphau);
	a=lkbeta(alphas);
	//#pragma omp parallel for firstprivate(j) private(da,db) lastprivate(i)
	for(i=0;i<P;i++){
		for(j=0;j<k;j++){
			ssum_dat_alphas[j*P+i]=stim[j*P+i]+alphas[j];
			ssum_dat_alphau[j*P+i]=unstim[j*P+i]+alphau[j];
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
		if(filter[i]){
			z[i+P]=0.0;z[i]=1.0;
		}
		cll[i] = z[i]*lnull[i]+z[i+P]*lresp[i];
#ifdef NDEBUG
		printf("%f = %f*%f + %f*%f\n",cll[i],z[i],lnull[i],z[i+P],lresp[i]);
#endif
	}
}
void simZ(std::vector<double> &q,std::vector<double> &lnull, std::vector<double> &lresp,std::vector<double>& z,std::vector<double> &p,std::vector<bool> &filter,int P, int k){
	int i;
	double lq = ::log(q[0]);
	double mlq = ::log(1.0-q[0]);

#ifdef NDEBUG
	printf("Prelim setup in simZ okay\n");
	printf("lq=%f mlq=%f\n",lq,mlq);
#endif

	for(i=0;i < lnull.size(); i++){
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
		printf("z1: %f z2: %f\n",z[i],z[i+P]);
		printf("Loop %d in simZ okay\n",i);
#endif
	}
}
inline double simQ(std::vector<double> &z, int P){
	std::vector<double> ab(0,2);
	double q;
	for(int j=0;j<2;j++){
		for(int i=0;i<P;i++){
			ab[j]=ab[j]+z[j*P+i];
		}
	}
	double q = ::Rf_rbeta(ab[0]+1,ab[1]+1);
	return q;
}
