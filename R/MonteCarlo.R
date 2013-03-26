# Implements an MCMC sampler for the MIMOSA
# Author: finak
###############################################################################


alphaProp1<-function(alpha,sigma,i){
	alpha[i]<-alpha[i]+rnorm(1,sd=sigma)
	alpha
}


simAlpha.s<-function(alpha.s,alpha.u,z,S,llnull,llresp,i,rate){
	a<-llnull(c(alpha.s,alpha.u))
	b<-llresp(c(alpha.s,alpha.u))
	
	old<-sum(a*z[,1]+b*z[,2])+dexp(alpha.s[i],0.0001,log=TRUE)
	alpha.s.prop<-alphaProp1(alpha.s,S,i)	
	if(any(alpha.s.prop<0)){
		return(list(alpha.s,rate))
	}
	
	a<-llnull(c(alpha.s.prop,alpha.u))
	b<-llresp(c(alpha.s.prop,alpha.u))
	
	new<-sum(a*z[,1]+b*z[,2])+dexp(alpha.s.prop[i],0.0001,log=TRUE)
	
	if(log(runif(1))<=(new-old)&is.finite(new-old)){
		rate[i]<-rate[i]+1
		return(list(alpha.s.prop,rate))
	}else{
		return(list(alpha.s,rate))
	}
	
}

simAlpha.u<-function(alpha.s,alpha.u,z,S,llnull,llresp,i,rate){
	a<-llnull(c(alpha.s,alpha.u))
	b<-llresp(c(alpha.s,alpha.u))
	
	old<-sum(a*z[,1]+b*z[,2])+dexp(alpha.u[i],0.0001,log=TRUE)
	alpha.u.prop<-alphaProp1(alpha.u,S,i)
	if(any(alpha.u.prop<0)){
		return(list(alpha.u,rate))
	}
	
	a<-llnull(c(alpha.s,alpha.u.prop))
	b<-llresp(c(alpha.s,alpha.u.prop))
	
	new<-sum(a*z[,1]+b*z[,2])+dexp(alpha.u.prop[i],0.0001,log=TRUE)
	
	if(log(runif(1))<=(new-old)&is.finite(new-old)){
		rate[i]<-rate[i]+1
		return(list(alpha.u.prop,rate))
	}else{
		return(list(alpha.u,rate))
	}
	
}


#TODO: prior on Q so that rbeta is positive
simQ<-function(z){
	rbeta(1,sum(z[,1])+1,sum(z[,2])+1)
}

simZ<-function(q,alpha.s,alpha.u,llnull,llresp){
	a<-llnull(c(alpha.s,alpha.u))+log(q) 
	b<-llresp(c(alpha.s,alpha.u))+log(1-q)
	den<-apply(cbind(a,b),1,function(x)log(sum(exp(x-max(x))))+max(x))
	p<-exp(a-den)
	z<-sapply(p,function(p)sample(c(0,1),1,prob=c(1-p,p),replace=FALSE))
	return(cbind(z,1-z))
}


icsdata2mvicsdata<-function(x){
	if(any(class(x)%in%"icsdata")){
		xx<-list(n.stim=x[,c("Ns","ns")],n.unstim=x[,c("Nu","nu")])
		attr(xx,"pData")<-attributes(x)$pData
		class(xx)<-c(class(xx),"mvicsdata")
		return(xx)
	}else{
		return(x)
	}
}

#'Fit the MIMOSA model via MCMC
#'
#'This is an internal function that fits the MIMOSA model via MCMC. It is called
#'from \code{MIMOSA}
#'
#'@param data a \code{list} with elements names "n.stim" and "n.unstim", the
#'  stimulated and unstimulated counts. Must be at least of dimension 2.
#'@param inits the initialization parameters for the MCMC routine. Can be
#'  initialized from \code{MDMix} with \code{initonly=TRUE}.
#'@param iter the number of Mote Carlo iterations
#'@param burn the number of burn-in iterations
#'@param thin The thinning interval
#'@param tune the number of iterations used for tuning the step size
#'@param outfile the output file name
#'@param alternative either "greater" or "not equal" for fitting the one-sided
#'  or two-sided MIMOSA model, respectively.
#'@param UPPER tuning parameter for the upper bound on the acceptance ratio of
#'  each paramter
#'@param LOWER tuning parmeter for the lower bound on the acceptance ratio of
#'  each paramter
#'@param FAST \code{TRUE,FALSE}. Use the heuristic (FAST=TRUE) for fitting a
#'  one-sided model rather than recomputing the normalization constant via MCMC
#'  for each step.
#'@param EXPRATE the mean of the prior distribution for the model hyperparameters.
#'  @rdname fitMCMC
#'  @name .fitMCMC
#'  @export
.fitMCMC<-function(data,inits=NULL,iter=250000, burn=50000, thin=1,tune=100,outfile=basename(tempfile(tmpdir=".",fileext=".dat")),alternative="greater",UPPER=0.5,LOWER=0.15,FAST=TRUE,EXPRATE=1e-4,pXi=1){
	alternative<-match.arg(alternative,c("greater","not equal"))
	data<-icsdata2mvicsdata(data)
	if(is.null(inits)){
		r<-MDMix(data);
		inits<-list(alpha.s=r@par.stim,alpha.u=r@par.unstim,q=r@w[1],z=round(r@z))
	}

	#If the alternative hypothesis is one-sided, then compute a filter for pu>ps and pass that to the MCMC code
	if(alternative=="greater"){
		ps<-t(do.call(cbind,apply(data$n.stim,1,function(x)(data.frame(prop.table(x))[-1L,,drop=FALSE]))))
		pu<-t(do.call(cbind,apply(data$n.unstim,1,function(x)(data.frame(prop.table(x))[-1L,,drop=FALSE]))))
		
		filter<-sapply(1:nrow(ps),function(i)all(ps[i,]<pu[i,]))
		FILTER=TRUE
	}else{
		filter<-rep(FALSE,nrow(data$n.stim))
		FILTER<-FALSE
	}
	result<-.Call("fitMCMC",as.matrix(data$n.stim),as.matrix(data$n.unstim),as.vector(inits$alpha.s),as.vector(inits$alpha.u),as.vector(inits$q),as.matrix(inits$z),as.vector(iter),as.vector(burn),as.vector(thin),as.numeric(tune),as.character(outfile),as.vector(filter),as.numeric(UPPER),as.numeric(LOWER),FILTER,FAST,as.numeric(EXPRATE),as.numeric(pXi), package="MIMOSA")
	if(inherits(result,"character")){
		return(result)
	}
	result$z<-cbind(result$z,1-result$z)
	result$getmcmc<-function(x=outfile){
		mcmc(read.table(x,sep="\t",header=T));
	}
  #TODO rewrite this to be more robust with large files
  #   result$getP<-function(x=paste(outfile,"P",sep="")){
  #     nc<-length(strsplit(readLines(x,n=1),"\t")[[1]])
  #     which.cols<-nc/3+1:nc
  #     #awk each co  lumn and read it
  #     quantiles<-sapply(which.cols[1:10],function(i){
  #       system(sprintf("awk '{print $%s}' %s > column",i,x))
  #       quantile(read.table("column",header=TRUE)[,1],c(0.025,0.5,0.975))
  #     })
  #     #return the quantiles 
  #   }
	result$getP<-function(x=paste(outfile,"P",sep=""),thin=1){
    if(thin>1){
      nc<-length(strsplit(readLines(x,1),"\t")[[1]])
      thins<-paste("p",paste(rep(";n",thin-1),collapse=""),sep="")
      s<-sprintf("sed -n '%s' %s|cut -f %s-%s",thins,x,(nc/3+1),nc)
      con<-pipe(s)
      d<-do.call(rbind,lapply(strsplit(readLines(con),"\t")[-1L],as.numeric))
      colnames(d)<-strsplit(readLines(x,1),"\t")[[1]][(nc/3+1):nc]
      d<-split(as.list(data.frame(d)),gl(nc/3,2))
      close(con)
    }else{
		  d<-mcmc(read.table(x,sep="\t",header=T));
		  nc<-ncol(d)
		  d<-split(as.list(data.frame(d[,(nc/3+1):nc])),gl(nc/3,2))
    }
		d
	}
	attr(result,"class")<-c(attr(result,"class"),"MDMixResult")
	attr(result,"pData")<-attr(data,"pData")
	result$n.stim<-data$n.stim
	result$n.unstim<-data$n.unstim
	result
}
