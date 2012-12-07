#include <Rcpp.h>
//#include <armadillo>
#include <RcppArmadillo.h>
#include <assert.h>
#include <omp.h>
#include <Rmath.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <R_ext/Utils.h>
#include "MIMOSA.h"
//#define DEBUG2
//#define NDEBUG
/*
 * 16 parameters
 */
RcppExport SEXP fitMCMC(SEXP _stim, SEXP _unstim, SEXP _alphas, SEXP _alphau, SEXP _q, SEXP _z,SEXP _iter, SEXP _burn, SEXP _thin, SEXP _tune,SEXP _outfile, SEXP _filter, SEXP _UPPER, SEXP _LOWER,SEXP _FILTER, SEXP _FAST,SEXP _EXPRATE){
	BEGIN_RCPP
	ntune=0;
	//TODO add argument to pass the complete list of unstimulated samples.
	using namespace Rcpp;
	using namespace arma;
	using namespace std;
	bool fixed = false;

	Rcpp::RNGScope globalscope;
	//			printf("%f %f %f %f\n",exp(normconstIBeta(15,3000,30,3000)),(normconstMC(30,3000,15,3000)),(exp(normconstIBeta(3000,30,3000,15))),(normconstMC(3000,15,3000,30)));
	//			exit(0);
	/*
	 * Copy R variables to standard vectors
	 */
	Rcpp::LogicalVector rfilter(_filter);
	std::vector<bool> filter(rfilter.length(),0);
	filter.resize(rfilter.length());
	copy(rfilter.begin(),rfilter.end(),filter.begin());

	FILTER = Rcpp::as<bool> (_FILTER);
	FAST = Rcpp::as<bool> (_FAST);
	double UPPER = Rcpp::as<double>(_UPPER);
	double LOWER = Rcpp::as<double>(_LOWER);
	bool REJECT = false;
	EXPRATE = Rcpp::as<double>(_EXPRATE);
	//printf("exponential prior rate = %f\n",EXPRATE);

	std::vector< double > alphas = Rcpp::as<std::vector<double> >(_alphas);
	std::vector< double > alphau = Rcpp::as<std::vector<double> 	>(_alphau);
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

	Rprintf("Creating %s\n",outfile.data());
  
  
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
	assert(alphas.size()==stim.ncol());
	assert(stim.nrow()==unstim.nrow());
	assert(stim.ncol()==unstim.ncol());
	assert(alphas.size()==alphau.size());
	/*
	 * Dimensions of the problem
	 */
	const int k = stim.ncol();
	const int P = stim.nrow();

	/*
	 * Output variables
	 */
	std::vector<double> ps(P,0),pu(P,0);
	std::vector <double> llnull(P,0);
	std::vector <double> llresp(P,0);
	std::vector <double> llnullNew(P,0);
	std::vector <double> llrespNew(P,0);


	std::vector <double> stdsum_data_alphas(P*k,0);
	std::vector <double> stdsum_data_alphau(P*k,0);

	arma::vec asi, aui;
	/*
	 * Armadillo objects will store the simulations from the tuning phase
	 * so we can easily compute variances and so forth.
	 */
	arma::mat Ms(int(TUNING),alphas.size());
#ifdef DEBUG2
	Rprintf("init Ms\n");
#endif
	arma::mat Mu(int(TUNING),alphau.size());
#ifdef DEBUG2
	Rprintf("init Mu\n");
#endif

	/*
	 * Precompute and preallocate a few things we'll need
	 */
	std::vector <double> sum_stim_unstim(P*k,0);
	std::transform(stdstim.begin(),stdstim.end(),stdunstim.begin(),sum_stim_unstim.begin(),plus<double>());


	std::vector<double> newalphas(alphas.size(),0);
	newalphas.resize(alphas.size());
	std::vector<double> newalphau(alphau.size(),0);
	newalphau.resize(alphau.size());


	NumericVector sigmas(alphas.size(),10.0);
	NumericVector sigmau(alphau.size(),10.0);

	arma::vec rateS(alphas.size());
#ifdef DEBUG2
	Rprintf("init rateS\n");
#endif
	arma::vec rateU(alphau.size());
#ifdef DEBUG2
	Rprintf("init rateU\n");
#endif
	rateS.fill(RATE);
	rateU.fill(RATE);
	NumericVector accepts(alphas.size(),0.0);
	NumericVector acceptu(alphas.size(),0.0);
	std::vector<double> cll(P,0.0);
	std::vector<double> p(P,0.0);
	std::vector<double> cz(P,0.0);


	double prior=0,newprior=0;
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



	for(int i=0;i<alphas.size();i++){
		headers<<"alphas."<<i<<"\t";
	}
	for(int i=0;i<alphau.size();i++){
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
		//::Rprintf("%d ",iteration);
		for(j=0;j<k;j++){
			/*
			 * prior for the current alphas_j
			 */

//			if(FILTER&&k==2&&!FAST){
//				prior=dgeom(alphas[j],EXPRATE);
//			}else{
				prior=::Rf_dexp(alphas[j],1.0/EXPRATE,true);
//			}
			std::copy(alphas.begin(),alphas.end(),newalphas.begin());

			//current null marginal log likelihood
			loglikenull(sum_stim_unstim,alphau,llnull,stdsum_data_alphau,P,k);

			//If alternative is greater compute the ratio of normalizing constants for the alternative marginal log likelihood
			if(FILTER&&k==2&&!FAST){
				normalizingConstant(stdstim,stdunstim,alphas,alphau,llresp,P,k);
			}else{
				//otherwise the two sided marginal log likelihood
				loglikeresp(stdstim,alphas,stdunstim,alphau,llresp,stdsum_data_alphas,stdsum_data_alphau,P,k);
			}

			//compute z1*lnull+z2*lresp+prior
			completeLL(z,llnull,llresp,cll,filter,P,k);
			oldll=std::accumulate(cll.begin(),cll.end(),0.0)+prior;
			assert(isfinite(oldll)==1);


			//simulate alphas_j
			//try integral values for alpha_s beta_s
			if(FILTER&&k==2&&!FAST){
				newalphas[j]=ceil(alphaProposal(alphas,sigmas[j]*rateS[j],j));
			}else{
				newalphas[j]=alphaProposal(alphas,sigmas[j]*rateS[j],j);
			}
			if(newalphas[j]>0){
//				if(FILTER&&k==2&&!FAST){
//					newprior=dgeom(newalphas[j],EXPRATE);
//				}else{
					newprior=::Rf_dexp(newalphas[j],1.0/EXPRATE,true);
//				}

				//don't need to recompute the null marginal log likelihood since it doesn't depend on alphas_j
				std::copy(llnull.begin(),llnull.end(),llnullNew.begin());
				if(FILTER&&k==2&&!FAST){
					//compute one sided marginal log likelihood
					normalizingConstant(stdstim,stdunstim,newalphas,alphau,llrespNew,P,k);
				}else{
					//two sided
					loglikeresp(stdstim,newalphas,stdunstim,alphau,llrespNew,stdsum_data_alphas,stdsum_data_alphau,P,k);
				}
				//compute z1*lnull+z2*lresp+prior
				completeLL(z,llnullNew,llrespNew,cll,filter,P,k);
				newll=std::accumulate(cll.begin(),cll.end(),0.0)+newprior;
			}else{
				REJECT=true;
				newll=nan(0);
			}


			//Rprintf("newll - oldll = %f isfinite()=%d\n",newll-oldll,isfinite(newll-oldll));

			if(!REJECT&&newalphas[j]>0&&(::log(Rf_runif(0.0,1.0)) <= (newll-oldll) )&&(!ISNAN(newll-oldll)&&(!ISNAN(newll))&&isfinite(::exp(newll-oldll)))){
#ifdef NDEBUG
				Rprintf("ACCEPTED alphas_%d %f prob: %f newll %f oldll %f\n",j,newalphas[j],::exp(newll-oldll), newll, oldll);
#endif
				accepts[j]=accepts[j]+1;//increment acceptance count
				alphas[j]=newalphas[j]; //new alphas_j is accepted
				oldll=newll-newprior; //save the new complete data log likelihood (minus the prior for alphas_j) so we don't recompute it for the next step
				std::copy(llnullNew.begin(),llnullNew.end(),llnull.begin()); //ditto for the null and alternative marginal log likelihood
				std::copy(llrespNew.begin(),llrespNew.end(),llresp.begin());
			}else{
				oldll=oldll-prior; //reject so just subtract the alpha-specific prior
#ifdef NDEBUG
				Rprintf("REJECTED alphas_%d %f, deltall: %f newll %f oldll %f\n",j,newalphas[j],(newll-oldll),newll, oldll);
				REJECT=false;
#endif
			}
			assert(isfinite(oldll)==1);



			/*
			 * Simulate alphau_j
			 */
			//prior for the current alphau_j
//			if(FILTER&&k==2&&!FAST){
//				prior=dgeom(alphau[j],EXPRATE);
//			}else{
				prior=Rf_dexp(alphau[j],1.0/EXPRATE,true);
//			}
			oldll=oldll+prior;

			//recompute to see if we still have problems with alpha u beta u estimates
//			loglikenull(stdsum_stim_unstim,stdalphau,stdllnullRes,stdsum_data_alphau,P,k); //new null marginal likelihood.
//			if(FILTER&&k==2&&!FAST){
//				normalizingConstant(stdstim,stdunstim,stdalphas,stdalphau,stdllrespRes,P,k); //new responder marginal LL - one sided
//			}else{
//				//two sided
//				loglikeresp(stdstim,stdalphas,stdunstim,stdalphau,stdllrespRes,stdsum_data_alpha,stdsum_data_alphau,P,k);
//			}
//			completeLL(z,stdllnullRes,stdllrespRes,cll,filter,P,k);
//			oldll=std::accumulate(cll.begin(),cll.end(),0.0)+prior;



			assert(isfinite(oldll)==1);
			//copy the alphau vector to the proposal vector.
			std::copy(alphau.begin(),alphau.end(),newalphau.begin());//copy the current alpha vector to the new alpha vector prior to drawing a sample

			//simulate alphau)j
			if(FILTER&&k==2&&!FAST){
				newalphau[j]=ceil(alphaProposal(alphau,sigmau[j]*rateU[j],j));
			}else{
				newalphau[j]=alphaProposal(alphau,sigmau[j]*rateU[j],j);
			}
			if(newalphau[j]>0){

//				if(FILTER&&k==2&&!FAST){
//					newprior=dgeom(newalphau[j],EXPRATE);
//				}else{
					newprior=Rf_dexp(newalphau[j],1.0/EXPRATE,true);
//				}

				//compute z1*lnull+z2*lresp+prior
				loglikenull(sum_stim_unstim,newalphau,llnullNew,stdsum_data_alphau,P,k); //new null marginal likelihood.
				if(FILTER&&k==2&&!FAST){
					normalizingConstant(stdstim,stdunstim,alphas,newalphau,llrespNew,P,k); //new responder marginal LL - one sided
				}else{
					//two sided
					loglikeresp(stdstim,alphas,stdunstim,newalphau,llrespNew,stdsum_data_alphas,stdsum_data_alphau,P,k);
				}
				//compute z1*lnull+z2*lresp+prior
				completeLL(z,llnullNew,llrespNew,cll,filter,P,k);
				newll=std::accumulate(cll.begin(),cll.end(),0.0)+newprior;
			}else{
				REJECT=true;
				newll=nan(0);
			}

			//	printf("newll - oldll = %f isfinite()=%d\n",newll-oldll,isfinite(newll-oldll));
			if(!REJECT&&newalphau[j]>0&&(log(Rf_runif(0.0,1.0)) <= (newll-oldll) )&&(!ISNAN(newll-oldll)&&(!ISNAN(newll))&&isfinite(::exp(newll-oldll)))){
				//increment acceptance count for alphauj
#ifdef NDEBUG
				Rprintf("ACCEPTED alphau_%d %f, prob ratio %f newll %f oldll %f\n",j,newalphau[j],::exp(newll-oldll),newll, oldll);
#endif
				acceptu[j]=acceptu[j]+1;
				alphau[j]=newalphau[j];//new alphau_j is accepted
				oldll=newll-newprior; //complete data log likelihood (minus the prior)
				//marginal null and alternative log likelihoods for the accepted parameter are saved so we don't have to recompute them
				std::copy(llnullNew.begin(),llnullNew.end(),llnull.begin());
				std::copy(llrespNew.begin(),llrespNew.end(),llresp.begin());
			}else{
				oldll=oldll-prior;
				REJECT=false;
#ifdef NDEBUG
				Rprintf("REJECTED alphau_%d %f, deltall: %f newll %f oldll %f\n",j,newalphau[j],(newll-oldll),newll, oldll);
#endif
			}
			//simulate q (w)
			q=simQ(z,P,k);
			//simulate z
			simZ(q,llnull,llresp,z,p,filter,P,k); //overwrites the current z. A running average is stored in cz

		}

		/*
		 * If we haven't fixed the step sizes yet..
		 */
		if(!fixed){
			/*
			 * Tuning phase
			 */
			//Fill in the Ms and Mu matrices
			asi = conv_to<arma::vec>::from(alphas);
#ifdef DEBUG2
			Rprintf("conv alphas to asi\n");
#endif
			aui = conv_to<arma::vec>::from(alphau);
#ifdef DEBUG2
			Rprintf("conv alphau to aui\n");
#endif
//print some dimensions
#ifdef DEBUG2
			Rprintf("Ms size is %d rows by %d columns and asi size is %d rows by %d columns\n",Ms.n_rows,Ms.n_cols,asi.n_rows,asi.n_cols);
#endif
			Ms.row(iteration%(int)TUNING) = trans(asi);
#ifdef DEBUG2
			Rprintf("assign asi to Ms.row\n");
#endif
			Mu.row(iteration%(int)TUNING) = trans(aui);
#ifdef DEBUG2
			Rprintf("assign aui to Mu.row\n");
#endif
			//Tuning
			if(((iteration+1) % (int)TUNING)==0){
				ntune=ntune+1;
				//Compute the covariances
				arma::mat covarianceS = arma::cov(Ms);
#ifdef DEBUG2
				Rprintf("compute cov Ms\n");
#endif
				arma::mat covarianceU = arma::cov(Mu);
#ifdef DEBUG2
				Rprintf("compute cov Mu\n");
#endif
				arma::vec dS = arma::sqrt(covarianceS.diag());
#ifdef DEBUG2
				Rprintf("sqrt diag S\n");
#endif
				arma::vec dU = arma::sqrt(covarianceU.diag());
#ifdef DEBUG2
				Rprintf("sqrt diag U\n");
#endif
				//Weight and assign
				//Tweak the acceptance rates
				accepts=accepts/TUNING;
				acceptu=acceptu/TUNING;

//				printf("sigmas: ");
//				for(j=0;j<sigmas.length();j++){
//					printf("%f %f ",sigmas[j],sigmau[j]);
//				}
//				printf("\n Acceptance rates:");
//				for(j = 0; j < accepts.length(); j++){
//					printf("%f %f ",accepts[j],acceptu[j]);
//				}
//				printf("\n");
				Rprintf("Tuning: %d\t",ntune);
				if((ntune<=12&&(Rcpp::any(accepts > UPPER).is_true() || Rcpp::any(acceptu > UPPER).is_true() || Rcpp::any(accepts < LOWER).is_true() || Rcpp::any(acceptu < LOWER).is_true()))){
					for(j=0;j<accepts.length();j++){
						//stimulated
						if((accepts[j]/DEFAULT_RATE) > 1.1 || (accepts[j]/DEFAULT_RATE )<0.8){
							rateS[j]=accepts[j]/DEFAULT_RATE;
						}
						//unstimulated
						if((acceptu[j]/DEFAULT_RATE) > 1.1 || (acceptu[j]/DEFAULT_RATE) < 0.8){
							rateU[j]=acceptu[j]/DEFAULT_RATE;
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
						if(ISNAN(sigmas[j])||sigmas[j]==0){
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
						}if(ISNAN(sigmau[j])||sigmau[j]==0){
							R_CheckUserInterrupt();
							sigmau[j]=1;
							iteration=0;
							acceptu[j]=0;
						}
					}

				}else{
					Rprintf("Fixed step size\n");
					Rprintf("sigmas: ");
					for(int j=0;j<sigmas.length();j++){
						Rprintf("%f %f ",sigmas[j],sigmau[j]);
					}
					Rprintf("\n Acceptance rates:");
					for(int j = 0; j < accepts.length(); j++){
						Rprintf("%f %f ",accepts[j],acceptu[j]);
					}
					Rprintf("\n");
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
			Rprintf("--- Done %i iterations ---\n",(int)iteration+1);
			//			printf("\n Acceptance rates:");
			//			for(int j = 0; j < accepts.length(); j++){
			//				printf("%f %f",accepts[j],acceptu[j]);
			//			}
			//			printf("\n");
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
			//when the model is two dimensional..
			if(k==2){
				//and the proportions sampled from each model
				sampleP(sum_stim_unstim,stdstim,stdunstim,alphas,alphau,z,ps,pu,P,k);
				for(int obs=0;obs<P;obs++){
					fprintf(fileP,"%f\t", z[obs]);
				}
				for(int obs=0;obs<P-1;obs++){
					fprintf(fileP,"%f\t%f\t",ps[obs],pu[obs]);
				}
				fprintf(fileP,"%f\t%f\n",ps[P-1],pu[P-1]);
			}
			for(int obs=0;obs<alphas.size();obs++){
				fprintf(file,"%f\t", alphas[obs]);
			}
			for(int obs=0;obs<alphau.size();obs++){
				fprintf(file,"%f\t", alphau[obs]);
			}
			fprintf(file,"%f\n",q);
		}
#ifdef NDEBUG
		Rprintf("alphas: %f %f alphau: %f %f\n",alphas[0],alphas[1],alphau[0],alphau[1]);
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
	ntune=0;
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
double alphaDiscreteProposal(const std::vector<double> &alpha, double d, int i){
	double na;
	d=round(abs(d));
	assert(d>=1);
	na = ::Rf_runif(-d,d);
	na=alpha[i]+na;
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
		//Rprintf("s1 %f, s0 %f, u1 %f u0 %f\n",s[1],s[0],u[1],u[0]);
		//Rprintf("alphas %f, betas %f, alphau %f betau %f\n",alphas[1],alphas[0],alphau[1],alphau[0]);

		numerator=normconstIBeta((double)s[0],(double)s[1],(double)u[0],(double)u[1]);
		//hyperparameters only
		denominator=normconstIBeta((double)alphas[0],(double)alphas[1],(double)alphau[0],(double)alphau[1]);
		//					nummc = normconstMC((double)s[1],(double)s[0],(double)u[1],(double)u[0]);
		//					denommc = normconstMC((double)alphas[1],(double)alphas[0],(double)alphau[1],(double)alphau[0]);
		//	Rprintf("%f/%f\n",nummc,denommc);
		//	Rprintf("%f  %f  %f\n",log(nummc),log(denommc),log(nummc)-log(denommc));


		//Rprintf("numerator: %f, denominator %f\n",exp(numerator), exp(denominator));
		CC=numerator-denominator;
		//printf("alphas %f betas %f alphau %f betau %f\n",alphas[1], alphas[0], alphau[1],alphau[0]);

		//		C=log(nummc)-log(denommc);
		//printf("%d.\t%f\t%f\n",i,nummc,denommc);
		//	printf("C: %f CC: %f  diff: %f\n",C,CC, C-CC);
		//double K=lgamma((double)s[1])+lgamma((double)s[0])-lgamma(s[1]+s[0])+lgamma((double)u[1])+lgamma((double)u[0])-lgamma(u[1]+u[0])-lgamma((double)alphas[1])-lgamma((double)alphas[0])+lgamma(alphas[1]+alphas[0])-lgamma((double)alphau[1])-lgamma((double)alphau[0])+lgamma(alphau[1]+alphau[0]);
		double K = ::Rf_lbeta((double) s[1],(double) s[0])+::Rf_lbeta((double) u[1], (double) u[0])-::Rf_lbeta((double) alphas[1],(double)alphas[0])-::Rf_lbeta((double)alphau[1],(double)alphau[0]);
		if(ISNAN(CC)){
			//			printf("s0: %f s1: %f u0: %f u1: %f as0: %f as1: %f au0: %f au1: %f \n",s[0],s[1],u[0],u[1],alphas[0],alphas[1],alphau[0],alphau[1]);
			//			printf("numerator: %f  denominator %f  C: %f\n",numerator,denominator,numerator-denominator);
			//			printf("log(1-exp(C))=%f\n",log(1-exp(numerator-denominator)));
			nummc = (normconstMC((double)s[1],(double)s[0],(double)u[1],(double)u[0]));
			denommc = (normconstMC((double)alphas[1],(double)alphas[0],(double)alphau[1],(double)alphau[0]));
			C=log(nummc)-log(denommc);
			CC=C;
		}
		//		printf("mc: %f acpprox: %f difference %f\n",C+K,K+CC,C-CC);
		//Rprintf("C:%f\n",exp(CC));
		//Rprintf("ll:%f\n\n",K+CC);

		llresp[i]=K+CC;
		//printf("%d.\t%f\n",i,llresp[i]);
	}
}

double normconstMC(double as, double bs, double au, double bu){
	double res;
	NumericVector r = rbeta(1000,au,bu);
	r=pbeta(r,as,bs,false,false);
	res = std::accumulate(r.begin(),r.end(),0.0)/r.length();
	return res;
}

//estimate_logZus_int<-function(alpha.u,beta.u,alpha.s,beta.s)
//{
//  if(round(beta.u)!=beta.u) # Test if it's an integer
//  {
//    # Compute for integer values around the true value
//    betaU<-beta.u
//    beta.u<-floor(betaU)
//    j<-alpha.u:(alpha.u+beta.u-1)
//    K<- -lbeta(alpha.s,beta.s)+lfactorial(alpha.u+beta.u-1)+lbeta(j+alpha.s,alpha.u+beta.u-1-j+beta.s)-lfactorial(j)-lfactorial(alpha.u+beta.u-1-j)
//    I1<-sum(exp(K))
//
//    beta.u<-ceiling(betaU)
//    j<-alpha.u:(alpha.u+beta.u-1)
//    K<- -lbeta(alpha.s,beta.s)+lfactorial(alpha.u+beta.u-1)+lbeta(j+alpha.s,alpha.u+beta.u-1-j+beta.s)-lfactorial(j)-lfactorial(alpha.u+beta.u-1-j)
//    I2<-sum(exp(K))
//
//    # Interpolation
//    I<-(ceiling(betaU)-betaU)*I1+(betaU-floor(betaU))*I2
//  }
//  else
//  {
//  j<-alpha.u:(alpha.u+beta.u-1)
//  K<- -lbeta(alpha.s,beta.s)+lfactorial(alpha.u+beta.u-1)+lbeta(j+alpha.s,alpha.u+beta.u-1-j+beta.s)-lfactorial(j)-lfactorial(alpha.u+beta.u-1-j)
//  I<-sum(exp(K))
//  }
//  return(I)
//}


double normconstIBeta(double au, double bu, double as, double bs){
	double alphau = round(au);
	double betaU = round(bu);
	double alphas = as;
	double betas = bs;
	double mx;
	std::vector<double> upper((int)betaU,0.0);
	double sum=0;
	double K=-Rf_lbeta(alphas,betas)+lgamma(alphau+betaU);
	//printf("summation length: %d\n",betaU);
	for(int j=alphau;j<=alphau+betaU-1;j++){
		upper[j-alphau]=Rf_lbeta(j+alphas,alphau+betaU-1-j+betas)-lgamma(j+1)-lgamma(alphau+betaU-j)+K;
	}
	std::vector<double>::iterator where=std::max_element(upper.begin(),upper.end());
	mx = *where;
	for(int i=0;i<upper.size();i++){
		sum=sum+::exp(upper[i]-mx);
	}
	sum=log(sum)+mx;
	return sum;
}
//double normconstIBeta(double as, double bs, double au, double bu){
//	double alphas = as;
//	double betas = bs;
//	double alphau = (double)ceil(au);
//	double betau = (double)ceil(bu);
//	//moved the test for negative proposal outside to the main loop
//	//	if(alphas<=0||betas<=0||alphau<=0||betau<=0){
//	//		return 0.0/0.0; //if any parameters are negative return nan.
//	//	}
//	double sum=0,mx=0;
//	double upper = (double) (alphau+betau);
//	std::vector<double> res((int)betau,0.0);
//	//	printf("as=%f bs=%f au=%f bu=%f\n",as,bs,au,bu);
//	//	printf("INTS: as=%f bs=%f au=%f bu=%f\n",alphas,betas,alphau,betau);
//	//	printf("upper: %d, alphau: %d size: %d\n",upper, alphau,res.size());
//	double K = -::Rf_lbeta(alphas,betas)+Rf_lgammafn(alphau+betau)-Rf_lgammafn(alphas+betas+alphau+betau-1);
//#ifdef FOO
//	printf("upper=%f\n",upper);
//	printf("K=%f\n",K);
//	printf("-::Rf_lbeta(alphas,betas)=%f\n",-::Rf_lbeta(alphas,betas));
//	printf("Rf_lgammafn(alphau+betau)=%f\n",Rf_lgammafn(alphau+betau));
//	printf("-Rf_lgammafn(alphas+betas+alphau-1)=%f\n",-Rf_lgammafn(alphas+betas+alphau-1));
//#endif
//	for(int j = (int)alphau;j<=((int)upper-1);j++){
//#ifdef FOO
//		printf("j=%d\n",j);
//		printf("Rf_lgammafn(alphas+j)=%f\n",Rf_lgammafn(alphas+j));
//		printf("Rf_lgammafn(alphau+betau+betas-1)=%f\n",Rf_lgammafn(alphau+betau+betas-j-1));
//		printf("-Rf_lgammafn(j+1)=%f\n",-Rf_lgammafn(j+1));
//		printf("-Rf_lgammafn(alphau+betau-j)=%f\n",-Rf_lgammafn(alphau+betau-j));
//		printf("sm = %f\n",K+Rf_lgammafn(alphas+j)+Rf_lgammafn(alphau+betau+betas-j-1)-Rf_lgammafn(j+1)-Rf_lgammafn(alphau+betau-j));
//#endif
//		res[j-(int)alphau]=K+Rf_lgammafn(alphas+j)+Rf_lgammafn(alphau+betau+betas-j-1)-Rf_lgammafn(j+1)-Rf_lgammafn(alphau+betau-j);
//	}
//
//	//todo normalized sum of exponentials
//	std::vector<double>::iterator where=std::max_element(res.begin(),res.end());
//	mx = *where;
//	for(int i=0;i<res.size();i++){
//		sum=sum+::exp(res[i]-mx);
//	}
//	sum=log(sum)+mx;
//	//printf("%f \n",sum);
//	return(sum);
//}


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

double dgeom(int k,double p){
	assert(k>=1);
	assert(p>=0&&p<=1);
	double olp=log(1-p);
	double lp=log(p);
	double res=olp*(k-1)+lp;
	return res;
}
