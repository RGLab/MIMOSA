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

RcppExport SEXP fitMCMCMultiStim(SEXP _stim, SEXP _unstim, SEXP _alphas, SEXP _alphau, SEXP _q, SEXP _z,SEXP _iter, SEXP _burn, SEXP _thin, SEXP _tune,SEXP _outfile, SEXP _filter, SEXP _UPPER, SEXP _LOWER,SEXP FILTER_, SEXP FAST_,SEXP EXPRATE_, SEXP CLASSES_,SEXP P_, SEXP k_, SEXP c_){
	BEGIN_RCPP
	int ntune=0;
	double ERATE=1000;
	bool FILTER,FAST;
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

	int k = Rcpp::as<int> (k_);
	int P = Rcpp::as<int> (P_);
	int c = Rcpp::as<int> (c_);

	FILTER = Rcpp::as<bool> (FILTER_);
	FAST = Rcpp::as<bool> (FAST_);
	double UPPER = Rcpp::as<double>(_UPPER);
	double LOWER = Rcpp::as<double>(_LOWER);
	bool REJECT = false;
	const double NITERS = Rcpp::as <double> (_iter);
	const double BURNIN = Rcpp::as <double> (_burn);
	const double THINNING = Rcpp::as <double> (_thin);
	const double TUNING = Rcpp::as <double> (_tune);

	mcmcparams parameters(k,P,c,NITERS);
	//set the prior for alphas.. exponential in this case
	parameters.setPrior(::Rf_dexp,EXPRATE_);
	//could also be geometric
	//parameters.setPrior(dgeom);

	parameters.initialize(_stim, _unstim, _alphas, _alphau, CLASSES_, _filter, _q, _z);
	exit(0);


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
	double realitcounter=0;


	/*
	 * Armadillo objects will store the simulations from the tuning phase
	 * so we can easily compute variances and so forth.
	 */
	arma::vec asi, aui;
	arma::mat Ms(int(TUNING),parameters.alphas.size());
	arma::mat Mu(int(TUNING),parameters.alphau.size());


	NumericMatrix sigmas(parameters.c,parameters.k,10.0);
	NumericMatrix sigmau(parameters.c,parameters.k,10.0);

	arma::mat rateS(parameters.c,parameters.k);
	arma::mat rateU(parameters.c,parameters.k);

	rateS.fill(RATE);
	rateU.fill(RATE);

	/*
	 * Construct the header string
	 *
	 */
	std::stringstream headers(stringstream::in|stringstream::out);
	std::stringstream headersP(stringstream::in|stringstream::out);

	//write the z's and p's to the second file
	for(int i=0;i<parameters.P;i++){
		headersP<<"z."<<i<<"\t";
	}
	for(int i=0;i<parameters.P-1;i++){
		headersP<<"ps."<<i<<"\t"<<"pu."<<i<<"\t";
	}
	headersP<<"ps."<<parameters.P-1<<"\t"<<"pu."<<parameters.P-1<<std::endl;



	for(int i=0;i<parameters.alphas.size();i++){
		headers<<"alphas."<<i<<"\t";
	}
	for(int i=0;i<parameters.alphau.size();i++){
		headers<<"alphau."<<i<<"\t";
	}
	headers<<"q"<<std::endl;
	//counters
	int iteration=0,j=0;

	/*
	 * Run the MCMC algorithm
	 */

	while(!parameters.DONE){
		parameters.proposeAlphaS(sigmas,rateS);
		parameters.update_sum_as();
		parameters.updaterespLL();
		parameters.updatenulLL();

		parameters.next();
	}

	for(iteration = 0; iteration < NITERS; iteration++){
			for(j=0;j<k;j++){

				/*
				 * prior for the current alphas_j
				 */

				//				prior=dgeom(alphas[j],EXPRATE);
				//calculate the priors
				parameters.calcPriorS(cl,j,true);
				//simulate propose new alphas
				parameters.proposeAlphaS(cl,j,sigmas,rateS);

				//compute the current and new marginal log likelihood and complete log likelihood
				parameters.loglikenull(cl,j);
				parameters.loglikeresp(cl,j);

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
					newprior=::Rf_dexp(newalphas[j],1.0/ERATE,true);
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
				prior=Rf_dexp(alphau[j],1.0/ERATE,true);
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
					newprior=Rf_dexp(newalphau[j],1.0/ERATE,true);
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
				aui = conv_to<arma::vec>::from(alphau);

				Ms.row(iteration%(int)TUNING) = asi;
				Mu.row(iteration%(int)TUNING) = aui;
				//Tuning
				if(((iteration+1) % (int)TUNING)==0){
					ntune=ntune+1;
					//Compute the covariances
					arma::mat covarianceS = arma::cov(Ms);
					arma::mat covarianceU = arma::cov(Mu);
					arma::vec dS = arma::sqrt(covarianceS.diag());
					arma::vec dU = arma::sqrt(covarianceU.diag());
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
						printf("Fixed step size\n");
						printf("sigmas: ");
						for(int j=0;j<sigmas.length();j++){
							printf("%f %f ",sigmas[j],sigmau[j]);
						}
						printf("\n Acceptance rates:");
						for(int j = 0; j < accepts.length(); j++){
							printf("%f %f ",accepts[j],acceptu[j]);
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
	fflush(file);fflush(fileP);
	fclose(file);fclose(fileP);
	ntune=0;
	return Rcpp::List::create(
			Rcpp::Named("z") = cz,
			Rcpp::Named("stepsizeS") = sigmas,
			Rcpp::Named("stepsizeU") = sigmau);
	END_RCPP
}
