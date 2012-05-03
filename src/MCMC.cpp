#include <Rcpp.h>
#include <armadillo>
#include <RcppArmadillo.h>
#include <assert.h>
#include <gsl/gsl_cblas.h>
#include <Rmath.h>
#include <iostream>
#include <fstream>
#include "MIMOSA.h"
#undef NDEBUG

#define RATE 2.4
#define DEFAULT_RATE  0.4
#define UPPER 1.2
#define LOWER 0.8
/*
 * 10 parameters
 */
RcppExport SEXP fitMCMC(SEXP _stim, SEXP _unstim, SEXP _alphas, SEXP _alphau, SEXP _q, SEXP _z,SEXP _iter, SEXP _burn, SEXP _thin, SEXP _tune,SEXP _outfile){
	BEGIN_RCPP
	using namespace Rcpp;
	using namespace arma;
	using namespace std;

	bool fixed = false;
	Rcpp::RNGScope globalscope;
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
	std::string outfile = Rcpp::as<std::string>(_outfile);

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
	NumericVector work(alphas.length(),0.0);
	NumericVector cll(stim.nrow(),0.0);
	NumericVector p(z.nrow(),0.0);
	NumericVector cz(z.nrow(),0.0);
	std::copy((z.column(0)).begin(),(z.column(0)).end(),cz.begin());


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
			prior=::Rf_dexp(alphas[j],1000,true);
			std::copy(alphas.begin(),alphas.end(),nextalphavec.begin());
			nextalphavec[j]=alphaProposal(&alphas,sigmas[j]*rateS[j],j);
			priornext=::Rf_dexp(nextalphavec[j],1000,true);

			loglikenull(&sum_stim_unstim,&alphau,&llnullRes,&sum_data_alphau,&work);
			loglikeresp(&stim,&alphas,&unstim,&alphau,&llrespRes,&sum_data_alpha,&sum_data_alphau,&work);

			//compute z1*lnull+z2*lresp+prior
			completeLL(&z,&llnullRes,&llrespRes,&cll);
			oldll=std::accumulate(cll.begin(),cll.end(),0.0)+prior;

			loglikenull(&sum_stim_unstim,&alphau,&llnullRes,&sum_data_alphau,&work);
			loglikeresp(&stim,&nextalphavec,&unstim,&alphau,&llrespRes,&sum_data_alpha,&sum_data_alphau,&work);

			//compute z1*lnull+z2*lresp+prior
			completeLL(&z,&llnullRes,&llrespRes,&cll);
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
			prior=Rf_dexp(alphau[j],1000,true);
			std::copy(alphau.begin(),alphau.end(),nextalphavec.begin());//copy the current alpha vector to the new alpha vector prior to drawing a sample
			nextalphavec[j]=alphaProposal(&alphau,sigmau[j]*rateU[j],j);
			priornext=Rf_dexp(nextalphavec[j],1000,true);

			/*
			 * Compute log likelihood
			 */
			loglikenull(&sum_stim_unstim,&alphau,&llnullRes,&sum_data_alphau,&work);
			loglikeresp(&stim,&alphas,&unstim,&alphau,&llrespRes,&sum_data_alpha,&sum_data_alphau,&work);
			//compute z1*lnull+z2*lresp+prior
			completeLL(&z,&llnullRes,&llrespRes,&cll);
			oldll=std::accumulate(cll.begin(),cll.end(),0.0)+prior;
			loglikenull(&sum_stim_unstim,&nextalphavec,&llnullRes,&sum_data_alphau,&work);
			loglikeresp(&stim,&alphas,&unstim,&nextalphavec,&llrespRes,&sum_data_alpha,&sum_data_alphau,&work);

			//compute z1*lnull+z2*lresp+prior
			completeLL(&z,&llnullRes,&llrespRes,&cll);
			newll=std::accumulate(cll.begin(),cll.end(),0.0)+priornext;

			if(all(nextalphavec>0).is_true()&&(log(Rf_runif(0.0,1.0)) <= (newll-oldll) )&&(!ISNAN(newll-oldll))){
				//increment acceptance count for alphauj
#ifdef NDEBUG
				printf("ACCEPTED alphau_%d %f, prob %f newll %f oldll %f\n",j,nextalphavec[j],::exp(newll-oldll),newll, oldll);
#endif
				acceptu[j]=acceptu[j]+1;
				alphau[j]=nextalphavec[j];
			}else{
#ifdef NDEBUG
				printf("REJECTED alphau_%d %f, prob: %f newll %f oldll %f\n",j,nextalphavec[j],::exp(newll-oldll),newll, oldll);
#endif
			}

			//simulate q
			q[0]=simQ(&z);

			//simulate z
			loglikenull(&sum_stim_unstim,&alphau,&llnullRes,&sum_data_alphau,&work);
			loglikeresp(&stim,&alphas,&unstim,&alphau,&llrespRes,&sum_data_alpha,&sum_data_alphau,&work);
			simZ(q,&llnullRes,&llrespRes,&z,&p); //overwrites the current z
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
				if((Rcpp::any(accepts/DEFAULT_RATE > UPPER).is_true() || Rcpp::any(acceptu/DEFAULT_RATE > UPPER).is_true() || Rcpp::any(accepts/DEFAULT_RATE < LOWER).is_true() || Rcpp::any(acceptu/DEFAULT_RATE < LOWER).is_true())){
					for(int j=0;j<accepts.length();j++){
						//stimulated
						if((accepts[j] / DEFAULT_RATE) > UPPER || (accepts[j] / DEFAULT_RATE) < LOWER){
							rateS[j]=accepts[j]/DEFAULT_RATE;
						}
						//unstimulated
						if((acceptu[j] / DEFAULT_RATE) > UPPER || (acceptu[j] / DEFAULT_RATE) < LOWER){
							rateU[j]=acceptu[j]/DEFAULT_RATE;
						}
					}
					for(int j=0;j<sigmas.length();j++){
						sigmas[j] = sigmas[j]*0.5+0.5*dS(j);
						sigmau[j] = sigmau[j]*0.5+0.5*dU(j);
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
			printf("--- Done %i iterations ---\n",(int)iteration+1);
		}

		if(iteration>=BURNIN&&iteration%(int)THINNING==0&&fixed){
			realitcounter++;
			for(int j=0;j<cz.length();j++){
				double f = realitcounter/(realitcounter+1.0);
				double foo = z(j,0)/realitcounter;
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
			fprintf(file,"\n");
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
	return Rcpp::List::create(Rcpp::Named("q") = q,
			Rcpp::Named("z") = z,
			Rcpp::Named("alpha.u") = alphau,
			Rcpp::Named("alpha.s") = alphas,
			Rcpp::Named("cz") = cz);
	END_RCPP
}

/*
 * Null component log likelihood
 * data is the data, alpha is the parameters, output is the result, sum_dat_alphau is data+alpha
 */
void loglikenull(const NumericMatrix* data,const NumericVector*  alpha,NumericVector* output, NumericMatrix* sum_dat_alphau,NumericVector* work){
	int i,j;
	double da,a;
	a=lkbeta(alpha,work);
	for(i=0;i<data->nrow();i++){
		for(j=0;j<data->ncol();j++){
			(*sum_dat_alphau)(i,j)=(*data)(i,j)+(*alpha)[j];
		}
		da=lkbeta(sum_dat_alphau->row(i),work);
		(*output)[i]=da-a;
	}
}
/*
 * Responder component log-likelihood
 */
void loglikeresp(const NumericMatrix*  stim, const NumericVector*  alphas,const  NumericVector*  unstim, const NumericVector*  alphau,NumericVector* output, NumericMatrix* sum_dat_alphas,NumericMatrix* sum_dat_alphau, NumericVector* work){
	int i,j;
	double da,db,a,b;

	a=lkbeta(alphas,work);
	b=lkbeta(alphau,work);

	for(i=0;i<stim->nrow();i++){
		for(j=0;j<stim->ncol();j++){
			(*sum_dat_alphas)(i,j)=(*stim)(i,j)+(*alphas)[j];
			(*sum_dat_alphau)(i,j)=(*unstim)(i,j)+(*alphau)[j];
		}
		da=lkbeta(sum_dat_alphas->row(i),work);
		db=lkbeta(sum_dat_alphau->row(i),work);
		(*output)(i)=da+db-a-b;
	}
}

/*
 * K-dimensional beta function
 */
double lkbeta(const NumericVector*  alpha,NumericVector *work){
	double sum_alpha = std::accumulate(alpha->begin(),alpha->end(),0.0);
	double log_gamma_sum_alpha=lgamma(sum_alpha);
	std::transform(alpha->begin(),alpha->end(),work->begin(),op_lgamma);
	double sum_log_gamma_alpha = std::accumulate(work->begin(),work->end(),0.0);
	return sum_log_gamma_alpha - log_gamma_sum_alpha;
}

/*
 * K-dimensional Beta function
 */
double lkbeta(const NumericMatrix::Row  alpha,NumericVector *work){
	double sum_alpha = std::accumulate(alpha.begin(),alpha.end(),0.0);
	double log_gamma_sum_alpha=lgamma(sum_alpha);
	std::transform(alpha.begin(),alpha.end(),work->begin(),op_lgamma);
	double sum_log_gamma_alpha = std::accumulate(work->begin(),work->end(),0.0);
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
double alphaProposal(const NumericVector* alpha, double sigma, int i){
	double na;
	na = ::Rf_rnorm((*alpha)[i],sigma);
	return na;
}

void completeLL(NumericMatrix* z,NumericVector* lnull, NumericVector* lresp,NumericVector* cll){
	int i;
	for(i=0;i< z->nrow();i++){
		(*cll)[i] = (*z)(i,0)*(*lnull)(i)+(*z)(i,1)*(*lresp)(i);
#ifdef NDEBUG
		printf("%f = %f*%f + %f*%f\n",(*cll)[i],(*z)(i,0),(*lnull)(i),(*z)(i,1),(*lresp)(i));
#endif
	}
}
void simZ(NumericVector q,NumericVector *lnull, NumericVector *lresp,NumericMatrix* z,NumericVector *p){
	int i;
	double lq = ::log(q[0]);
	double mlq = ::log(1.0-q[0]);

#ifdef NDEBUG
	printf("Prelim setup in simZ okay\n");
	printf("lq=%f mlq=%f\n",lq,mlq);
#endif

	for(i=0;i < lnull->length(); i++){
		(*lnull)[i]=(*lnull)[i]+lq;
		(*lresp)[i]=(*lresp)[i]+mlq;
		double mx=std::max((*lnull)[i],(*lresp)[i]);
		(*p)[i] = ::exp((*lnull)[i]-::log(::exp((*lnull)[i]-mx)+::exp((*lresp)[i]-mx))-mx);
		(*z)(i,0) = ::Rf_rbinom(1.0,(*p)[i]);
		(*z)(i,1) = 1.0-(*z)(i,0);
#ifdef NDEBUG
		printf("z1: %f z2: %f\n",(*z)(i,0),(*z)(i,1));
		printf("Loop %d in simZ okay\n",i);
#endif
	}
}
inline double simQ(NumericMatrix* z){
	double alpha=std::accumulate((z->column(0)).begin(),(z->column(0)).end(),0.0);
	double beta=std::accumulate((z->column(1)).begin(),(z->column(1)).end(),0.0);
	double q = ::Rf_rbeta(alpha+0.0001,beta+0.0001);
	return q;
}
