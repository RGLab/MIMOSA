BetaMix <-
		function(d=NULL,maxiter=20,shrink=1,alternative.model="greater",mciter=50,fixedNULL=FALSE,ics=NULL,priorXi=1,inits=NULL,scl=10,K=200){
	#priorXi is a uniform Beta(Xi,Xi) prior on the mixing proportions. By default it is set to 2.
	alternative.model<-match.arg(alternative.model,c("greater","not equal"))
	#result<-replicate(1,{
	LL<--.Machine$double.xmax
	#try many initializations and start from the largest log likelihood
	if(is.null(inits)){		
		inits<-initBetaMix(d,fixedNULL=fixedNULL,ics=ics,alternative=alternative.model,priorXi=priorXi,scl=scl,K=K,mciter=mciter)
		llmax<-CompleteDataLLRcpp(d=d,alpha0=inits$alpha0,alphaS=inits$alphaS,beta0=inits$beta0,betaS=inits$betaS,z=inits$z,w=inits$w,alternative=alternative.model,mciter=mciter)		
	}
	lltraj<-LL
	iter<-1
	#Given Z, estiamate alpha0,beta0,alphaS,betaS
	#optimize over alpha0,beta0, alphaS, betaS for the complete data (conditional on z's)
	repeat{
		fail<-FALSE
		last<-inits;
		n0.est<-(c(inits$alpha0,inits$beta0))
		nS.est<-(c(inits$alphaS,inits$betaS))
		#Update alpha0,beta0 on alternate iterations, and only if they are not fixed.
		#Alternatively try to update both parameters simultaneously using the mean / sample-size parameterization.. this seems to work.
		if(fixedNULL){
			pars<-try(optim(par=(c(nS.est[1]/(nS.est[1]+nS.est[2]), nS.est[2]+nS.est[1])),function(p=par,data=d,Z=inits$z,W=inits$w,ALT=alternative.model,MC=mciter,a0=n0.est[1],b0=n0.est[2])f0m(p=p,d=data,z=Z,w=W,alternative=ALT,mciter=MC,alpha0=a0,beta0=b0),control=list(parscale=c(scl,1)),method="L-BFGS-B",lower=c(1e-6,10),upper=c(0.9999,K/inits$muS)),silent=TRUE)			
		}
		else{
			pars<-try(optim(par=(c(nS.est[1]/(nS.est[1]+nS.est[2]), nS.est[2]+nS.est[1],n0.est[1]/(n0.est[1]+n0.est[2]),n0.est[2]+n0.est[1])),function(p=par,data=d,Z=inits$z,W=inits$w,ALT=alternative.model,MC=mciter)f0(p=p,d=data,z=Z,w=W,alternative=ALT,mciter=MC),method="L-BFGS-B",control=list(parscale=c(scl,1,scl,1)),lower=c(1e-6,10,1e-6,10),upper=c(0.9999,K/inits$muS,0.9999,K/inits$mu0)),silent=TRUE)			
		}
		if(inherits(pars,"try-error")|inherits(try(pars$convergence,silent=TRUE)!=0,"try-error")){
			message("failed to converge estimating alpha0,beta0, alphaS, betaS in BetaMix: returning last best iteration, ",iter)
			iter<-maxiter+1
		}else{
			inits$alphaS<-pars$par[1]*pars$par[2]
			inits$betaS<-(1-pars$par[1])*pars$par[2]
			if(!fixedNULL){
				inits$alpha0<-pars$par[3]*pars$par[4]
				inits$beta0<-(1-pars$par[3])*pars$par[4]
			}
		}
		
		#Update the zs
		if(alternative.model=="greater"){
			m1<-log(inits$w[1])+.Call("MarginalNULL",ns=d[,"ns"],Ns=d[,"Ns"],nu=d[,"nu"],Nu=d[,"Nu"],alpha0=rep(inits$alpha0,nrow(d)),beta0=rep(inits$beta0,nrow(d)),alphaS=rep(inits$alphaS,nrow(d)),betaS=rep(inits$betaS,nrow(d)),inits$w,log=TRUE,package="MIMOSA")
			m2<-log(inits$w[2])+.Call("MarginalGT",ns=d[,"ns"],Ns=d[,"Ns"],nu=d[,"nu"],Nu=d[,"Nu"],alpha0=rep(inits$alpha0,nrow(d)),beta0=rep(inits$beta0,nrow(d)),alphaS=rep(inits$alphaS,nrow(d)),betaS=rep(inits$betaS,nrow(d)),inits$w,log=TRUE,mciter=mciter,package="MIMOSA")
			z<-exp(m1-(log(rowSums(exp(cbind(m1, m2) - pmax(m1, m2))))+pmax(m1,m2)))
		} else{
			m1<-log(inits$w[1])+.Call("MarginalNULL",ns=d[,"ns"],Ns=d[,"Ns"],nu=d[,"nu"],Nu=d[,"Nu"],alpha0=rep(inits$alpha0,nrow(d)),beta0=rep(inits$beta0,nrow(d)),alphaS=rep(inits$alphaS,nrow(d)),betaS=rep(inits$betaS,nrow(d)),inits$w,log=TRUE,package="MIMOSA")
			m2<-log(inits$w[2])+.Call("MarginalNE",ns=d[,"ns"],Ns=d[,"Ns"],nu=d[,"nu"],Nu=d[,"Nu"],alpha0=rep(inits$alpha0,nrow(d)),beta0=rep(inits$beta0,nrow(d)),alphaS=rep(inits$alphaS,nrow(d)),betaS=rep(inits$betaS,nrow(d)),inits$w,log=TRUE,package="MIMOSA")
			z<-exp(m1-(log(rowSums(exp(cbind(m1, m2) - pmax(m1, m2))))+pmax(m1,m2)))
		}
		inits$z<-cbind(z,1-z)		
		
		#Update w
		inits$w<-(colSums(inits$z)+priorXi-1)/sum(colSums(inits$z)+priorXi-1)
		
		#Compute the complete data LL
		ll<-CompleteDataLLRcpp(d=d,alpha0=inits$alpha0,alphaS=inits$alphaS,beta0=inits$beta0,betaS=inits$betaS,z=inits$z,w=inits$w,alternative=alternative.model,mciter=mciter)
		lltraj<-c(lltraj,ll)
		cat(iter);
		
		if(iter>=maxiter){
			iter<-iter+1
			break
			#also test for hyperparameter convergence..
		}else if(((ll-LL)>0)&((abs(ll-LL)/abs(LL))<(.Machine$double.eps)^(1/2))){ 
			break
		}else {
			iter<-iter+1
			LL<-ll
		}
	}
	
	if(!fail){
		result<-BetaMixResult(alternative.model=alternative.model,control=attr(d,"control"),stimulation=attr(d,"stimulation"),cytokine=attr(d,"cytokine"),ll=ll,iter=iter,traj=lltraj,z=inits$z,w=inits$w,alpha0=inits$alpha0,beta0=inits$beta0,alphaS=inits$alphaS,betaS=inits$betaS,data=as.data.frame(d))
	}else{
		return(new("BetaMixResult",data=as.data.frame(d),z=matrix(rep(NA_real_,nrow(d)*2),ncol=2),fdr=rep(NA_real_,nrow(d))))
	}
	return(result);
}

