setGeneric("estimateProportions",function(bmr,method,...){
			standardGeneric("estimateProportions")
		})
setMethod("estimateProportions",signature=c("BetaMixResult"),function(bmr,method="posterior",N=10000,burn=1000,volcano=FALSE){
			method<-match.arg(method,c("posterior","ML"))
			if(method%in%"posterior"){
				r<-t(sapply(1:nrow(bmr@data),function(i){colMeans(gibbsPsPu(bmr@data,inits=list(w=bmr@w,z=bmr@z,alpha0=bmr@alpha0,beta0=bmr@beta0,alphaS=bmr@alphaS,betaS=bmr@betaS),alt=bmr@alternative.model,which=i,N=N,burn=burn,z=NULL,volcano=volcano))}))
			}else{
				r<-data.frame(Pu=bmr@data[,"nu"]/(bmr@data[,"nu"]+bmr@data[,"Nu"]),Ps=bmr@data[,"ns"]/(bmr@data[,"ns"]+bmr@data[,"Ns"]))
			}
			return(r);
		})

gibbsPsPu <-
		function(curdat=NULL,inits=NULL,which=NULL,N=10000,burn=5000,alt,z=NULL,ecdf.approx=FALSE,volcano=FALSE,...){
	threshold<-list(...)$threshold
	if(!is.null(threshold)){
		if(!is.numeric(threshold)){
			stop("whoops.. threshold must be numeric between zero and one.")
		}
	}
	Ns<-curdat[which,"Ns"];ns<-curdat[which,"ns"];Nu<-curdat[which,"Nu"];nu<-curdat[which,"nu"]
	alpha0<-inits$alpha0;beta0<-inits$beta0;alphaS<-inits$alphaS;betaS<-inits$betaS;
	ps<-double(N)
	pu<-double(N)
	PU<-double(N)
	PS<-double(N)
	pu<-rbeta(N,nu+alpha0,Nu+beta0)
	PU<-rbeta(N,nu+ns+alpha0,Nu+Ns+beta0)
	PS<-rbeta(N,nu+ns+alpha0,Nu+Ns+beta0)
	for(i in 1:N){
		#ps[i]<-inits$w[2]*qbeta(runif(1,pbeta(pu[i],alphaS+ns-1,betaS+Ns-1),1),alphaS+ns-1,betaS+Ns-1)+inits$w[1]*pu[i]
		#TODO estimates of ps and pu given model 2 should be different when alternative model changes	
		if(alt=="greater"){
			ps[i]<-qbeta(runif(1,pbeta(pu[i],alphaS+ns,betaS+Ns,log.p=F),1),alphaS+ns,betaS+Ns,log.p=F)
		}else{
			ps[i]<-rbeta(1,alphaS+ns,betaS+Ns)
		}
	}
	#Sometimes qbeta and pbeta don't give a decent answer when very near the tails. Approximate using the ecdf. We'll say we need at least 2000 samples
	if(ecdf.approx){
		if(N<2000){
			stop("Please sample at least 2000 points to use the empirical cdf approximation.")
		}
		e<-ecdf(ps)
		for(i in 1:N){
			ps[i]<-qbeta(runif(1,e(pu[i]),1),alphaS+ns,betaS+Ns)
		}
		
	}
	#pu<-inits$w[1]*PU+inits$w[2]*pu
	#ps<-PS*inits$w[1]+ps*inits$w[2]
	#Given the data, compute which is the most likely generating cluster, and return simulations from that.
	#If z is not given 
	if(volcano){
		return(cbind(pu=inits$z[which,1]*PU+inits$z[which,2]*pu,ps=inits$z[which,1]*PS+inits$z[which,2]*ps))
	}
	if(is.null(z)){
		if(alt=="greater"){
			m1<-exp(.Call("MarginalNULL",ns=ns,Ns=Ns,nu=nu,Nu=Nu,alpha0=rep(inits$alpha0,length(ns)),beta0=rep(inits$beta0,length(ns)),alphaS=rep(inits$alphaS,length(ns)),betaS=rep(inits$betaS,length(ns)),inits$w,log=TRUE,package="MIMOSA"))
			m2<-exp(.Call("MarginalGT",ns=ns,Ns=Ns,nu=nu,Nu=Nu,alpha0=rep(inits$alpha0,length(ns)),beta0=rep(inits$beta0,length(ns)),alphaS=rep(inits$alphaS,length(ns)),betaS=rep(inits$betaS,length(ns)),inits$w,log=TRUE,mciter=50,package="MIMOSA"))
			z<-m1/(m1+m2)
			z<-cbind(z,1-z)
		}else{
			m1<-exp(.Call("MarginalNULL",ns=ns,Ns=Ns,nu=nu,Nu=Nu,alpha0=rep(inits$alpha0,length(ns)),beta0=rep(inits$beta0,length(ns)),alphaS=rep(inits$alphaS,length(ns)),betaS=rep(inits$betaS,length(ns)),inits$w,log=TRUE,package="MIMOSA"))
			m2<-exp(.Call("MarginalNE",ns=ns,Ns=Ns,nu=nu,Nu=Nu,alpha0=rep(inits$alpha0,length(ns)),beta0=rep(inits$beta0,length(ns)),alphaS=rep(inits$alphaS,length(ns)),betaS=rep(inits$betaS,length(ns)),inits$w,log=TRUE,package="MIMOSA"))
			z<-m1/(m1+m2)
			z<-cbind(z,1-z)
		}
		if(is.null(threshold)){
			if(max.col(z)==1){
				return(cbind(pu=PU,ps=PS)[(burn+1):N,])
			}else{
				return(cbind(pu=pu,ps=ps)[(burn+1):N,])
			}
		}else{
			z.fdr<-fdr(inits$z)
			if((z.fdr<threshold)[which]){
				return(cbind(pu=pu,ps=ps)[(burn+1):N,])
			}else{
				return(cbind(pu=PU,ps=PS)[(burn+1):N,])
			}
		}
	}else{
		#if z is given, mix the PS,PU nd ps, pu samples in the appropriate proportions
		pu[z==1]<-PU[z==1]
		ps[z==1]<-PS[z==1]
		return(cbind(pu=pu,ps=ps)[(burn+1):N,])
	}
}


##Samples ps and pu from their conditional distributions
#gibbsPsPu <-
#function(curdat=NULL,inits=NULL,which=NULL,N=10000,burn=5000,alt){
#  Ns<-curdat[which,"Ns"];ns<-curdat[which,"ns"];Nu<-curdat[which,"Nu"];nu<-curdat[which,"nu"]
#  alpha0<-inits$alpha0;beta0<-inits$beta0;alphaS<-inits$alphaS;betaS<-inits$betaS;
#  ps<-double(N)
#  pu<-double(N)
#  PU<-double(N)
#  PS<-double(N)
#  pu<-rbeta(N,nu+alpha0,Nu+beta0)
#  PU<-rbeta(N,nu+ns+alpha0,Nu+Ns+beta0)
#  PS<-rbeta(N,nu+ns+alpha0,Nu+Ns+beta0)
#  for(i in 1:N){
#      #ps[i]<-inits$w[2]*qbeta(runif(1,pbeta(pu[i],alphaS+ns-1,betaS+Ns-1),1),alphaS+ns-1,betaS+Ns-1)+inits$w[1]*pu[i]
#	  #TODO estimates of ps and pu given model 2 should be different when alternative model changes	
#	if(alt=="greater"){
#	  ps[i]<-qbeta(runif(1,pbeta(pu[i],alphaS+ns,betaS+Ns),1),alphaS+ns,betaS+Ns)
#  	}else{
#		ps[i]<-rbeta(1,alphaS+ns-1,betaS+Ns)
#	}
#  }
#  #pu<-inits$w[1]*PU+inits$w[2]*pu
#  #ps<-PS*inits$w[1]+ps*inits$w[2]
#  #Given the data, compute which is the most likely generating cluster, and return simulations from that.
#  if(alt=="greater"){
#		m1<-exp(.Call("MarginalNULL",ns=ns,Ns=Ns,nu=nu,Nu=Nu,alpha0=rep(inits$alpha0,length(ns)),beta0=rep(inits$beta0,length(ns)),alphaS=rep(inits$alphaS,length(ns)),betaS=rep(inits$betaS,length(ns)),inits$w,log=TRUE,package="flowModels"))
#		m2<-exp(.Call("MarginalGT",ns=ns,Ns=Ns,nu=nu,Nu=Nu,alpha0=rep(inits$alpha0,length(ns)),beta0=rep(inits$beta0,length(ns)),alphaS=rep(inits$alphaS,length(ns)),betaS=rep(inits$betaS,length(ns)),inits$w,log=TRUE,mciter=50,package="flowModels"))
#		z<-m1/(m1+m2)
#		z<-cbind(z,1-z)
#  }else{
#		m1<-exp(.Call("MarginalNULL",ns=ns,Ns=Ns,nu=nu,Nu=Nu,alpha0=rep(inits$alpha0,length(ns)),beta0=rep(inits$beta0,length(ns)),alphaS=rep(inits$alphaS,length(ns)),betaS=rep(inits$betaS,length(ns)),inits$w,log=TRUE,package="flowModels"))
#		m2<-exp(.Call("MarginalNE",ns=ns,Ns=Ns,nu=nu,Nu=Nu,alpha0=rep(inits$alpha0,length(ns)),beta0=rep(inits$beta0,length(ns)),alphaS=rep(inits$alphaS,length(ns)),betaS=rep(inits$betaS,length(ns)),inits$w,log=TRUE,package="flowModels"))
#		z<-m1/(m1+m2)
#		z<-cbind(z,1-z)
#  }
#  
#  if(max.col(z)==1){
#	  cbind(pu=PU,ps=PS)[(burn+1):N,]
#  }else{
#  	cbind(pu=pu,ps=ps)[(burn+1):N,]
#  }
#}

#Given z, returns the q-values for each observation. The q-value for observation j is the false discovery rate for observations with q<q_j

#' Compute the fdr (q-value) from posterior probabilities
#' 
#' Given the z's from a MIMOSA model, calculates the q-values for each observation.
#' @rdname fdr
#' @param z matrix of posterior probabilties, or a \code{MIMOSAResult}, or \code{MIMOSAResultList}
#' @return a vector of q-values or a list of vectors of q-values.
#' @export fdr
fdr<-function(z){
  UseMethod("fdr")
}

#'@rdname fdr
#'@method fdr matrix
#'@S3method fdr matrix
fdr.matrix<-function(z){
	fdr<-rep(0,nrow(z)); o<-order(z[,2],decreasing=T); 
	fdr[o]<-(cumsum(z[o,1])/1:nrow(z))
	return(fdr)
}

#'@rdname fdr
#'@method fdr MIMOSAResult
#'@S3method fdr MIMOSAResult
fdr.MIMOSAResult <- function(z){
  fdr(z@z)
}

#'@rdname fdr
#'@method fdr MIMOSAResultList
#'@S3method fdr MIMOSAResultList
fdr.MIMOSAResultList <- function(z){
  r<-lapply(z,function(y)fdr(y))
  names(r)<-names(z)
  return(data.frame(r))
}
