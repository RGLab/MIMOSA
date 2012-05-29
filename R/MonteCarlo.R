# TODO: Having issues with newton-raphson and optim for the dirichlet-multinomial model
# Implement a MCMC sampler.
# Author: finak
###############################################################################


#sProposal<-function(s,sigma=10){
#	return(rnorm(1,mean=s,sd=sigma))
#}
#alphaProposal<-function(alpha,sigma=diag(rep(10,length(alpha))),C=2.38/sqrt(length(alpha))){
#	#sapply(alpha,function(x)rnorm(1,mean=x,sd=sigma))
#	mvrnorm(mu=alpha,Sigma=sigma*C^2)
#}
alphaProp1<-function(alpha,sigma,i){
	alpha[i]<-alpha[i]+rnorm(1,sd=sigma)
	alpha
}

#lddirichlet<-function(x,alpha){
#	logD<-lkbeta(alpha)
#	s<-sum(log(x)*(alpha-1))
#	s-logD
#}
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
	
	new<-sum(a*z[,1]+b*z[,2])+dexp(alpha.s.prop[i],0.0001,log=TRUE)
	
	if(log(runif(1))<=(new-old)&is.finite(new-old)){
		rate[i]<-rate[i]+1
		return(list(alpha.u.prop,rate))
	}else{
		return(list(alpha.u,rate))
	}
	
}

#wProposal<-function(w,c=0.2){
#	W<-t(scale(runif((length(w)-1),-c,c),scale=FALSE))
#}

#TODO: prior on Q so that rbeta is positive

simQ<-function(z){
	rbeta(1,sum(z[,1])+0.0001,sum(z[,2])+0.0001)
}

simZ<-function(q,alpha.s,alpha.u,llnull,llresp){
	a<-llnull(c(alpha.s,alpha.u))+log(q) 
	b<-llresp(c(alpha.s,alpha.u))+log(1-q)
	den<-apply(cbind(a,b),1,function(x)log(sum(exp(x-max(x))))+max(x))
	p<-exp(a-den)
	z<-sapply(p,function(p)sample(c(0,1),1,prob=c(1-p,p),replace=FALSE))
	return(cbind(z,1-z))
}

#simWs<-function(w.s,w.u,z,S.stim,S.unstim,w.stepU,llnull,llresp){
#	a<-llnull(c(w.s*S.stim,w.u*S.unstim))
#	b<-llresp(c(w.s*S.stim,w.u*S.unstim))
#	
#	old<-sum(a*z[,1]+b*z[,2]+lddirichlet(matrix(w.s,nrow=1),rep(1,length(w.s))))
#	w.s.prop<-wProposal(w.s,c=w.stepU)
#	if(any(w.s.prop<0|w.s.prop>1)){
#		return(w.s)
#	}
#	
#	a<-llnull(c(w.s.prop*S.stim,w.u*S.unstim))
#	b<-llresp(c(w.s.prop*S.stim,w.u*S.unstim))
#	
#	new<-sum(a*z[,1]+b*z[,2]+lddirichlet(matrix(w.s.prop,nrow=1),rep(1,length(w.s.prop))))
#	
#	if(log(runif(1))<=(new-old)&is.finite(new-old)){
#		return(w.s.prop)
#	}else{
#		return(w.s)
#	}
#}
#
#simWu<-function(w.s,w.u,z,S.stim,S.unstim,w.stepS,llnull,llresp){
#	a<-llnull(c(w.s*S.stim,w.u*S.unstim))
#	b<-llresp(c(w.s*S.stim,w.u*S.unstim))
#	
#	old<-sum(a*z[,1]+b*z[,2]+lddirichlet(matrix(w.u,nrow=1),rep(1,length(w.u))))
#	w.u.prop<-wProposal(w.u,c=w.stepS)
#	if(any(w.u.prop<0|w.u.prop>1)){
#		return(w.u)
#	}
#	
#	a<-llnull(c(w.s*S.stim,w.u.prop*S.unstim))
#	b<-llresp(c(w.s*S.stim,w.u.prop*S.unstim))
#	
#	new<-sum(a*z[,1]+b*z[,2]+lddirichlet(matrix(w.u.prop,nrow=1),rep(1,length(w.u.prop))))
#	
#	#wierd, should be new/old, but that doesn't work out right
#	if(log(runif(1))<=(new-old)&is.finite(new-old)){
#		return(w.u.prop)
#	}else{
#		return(w.u)
#	}
#}
#
#simS.stim<-function(w.s,w.u,z,S.stim,S.unstim,sigma,llnull,llresp){
#	a<-llnull(c(w.s*S.stim,w.u*S.unstim))
#	b<-llresp(c(w.s*S.stim,w.u*S.unstim))
#	
#	old<-sum(a*z[,1]+b*z[,2]+dexp(S.stim,rate=1e-3,log=TRUE))
#	S.stim.prop<-sProposal(S.stim,sigma=sigma)
#	
#	a<-llnull(c(w.s*S.stim.prop,w.u*S.unstim))
#	b<-llresp(c(w.s*S.stim.prop,w.u*S.unstim))
#	
#	new<-sum(a*z[,1]+b*z[,2]+dexp(S.stim.prop,rate=1e-3,log=TRUE))
#	
#	if(log(runif(1))<=(new-old)&is.finite(new-old)){
#		return(S.stim.prop)
#	}else{
#		return(S.stim)
#	}
#	
#}
#simS.unstim<-function(w.s,w.u,z,S.stim,S.unstim,sigma,llnull,llresp){
#	a<-llnull(c(w.s*S.stim,w.u*S.unstim))
#	b<-llresp(c(w.s*S.stim,w.u*S.unstim))
#	
#	old<-sum(a*z[,1]+b*z[,2]+dexp(S.unstim,rate=1e-3,log=TRUE))
#	S.unstim.prop<-sProposal(S.unstim,sigma=sigma)
#	
#	a<-llnull(c(w.s*S.stim,w.u*S.unstim.prop))
#	b<-llresp(c(w.s*S.stim,w.u*S.unstim.prop))
#	
#	new<-sum(a*z[,1]+b*z[,2]+dexp(S.unstim.prop,rate=1e-3,log=TRUE))
#	
#	if(log(runif(1))<=(new-old)&is.finite(new-old)){
#		return(S.unstim.prop)
#	}else{
#		return(S.unstim)
#	}
#}

fitMCMC<-function(data=NULL,inits=NULL,iter=5000,burn=2000,thin=1){
	if(burn>=iter){
		burn<-round(iter/2)
	}
	if(is.null(data)){
		stop("Must provide data")
	}
	if(is.null(inits)){
		stop("Must provide initializing values")
	}
	llnull<-compiler::cmpfun(MIMOSA:::makeLogLikeNULLComponent(data$n.stim,data$n.unstim))
	llresp<-compiler::cmpfun(MIMOSA:::makeLogLikeRespComponent(data$n.stim,data$n.unstim))
	
	alpha.s=inits$alpha.s
	alpha.u=inits$alpha.u
	q=inits$q
	z=inits$z
	Ss<-Su<-rep(10,length(alpha.s))
	rate.s<-rate.u<-rep(0,length(alpha.s))
	f<-matrix(0,nrow=iter,ncol=length(alpha.s)+length(alpha.u)+1)
	cn<-c(paste("alpha.s",1:length(alpha.s),sep=""),paste("alpha.u",1:length(alpha.u),sep=""),"q")
	colnames(f)<-cn
	tune=200
	fixed<-FALSE
	i<-1
	repeat{
		for(j in seq_along(alpha.s)){
			res<-MIMOSA:::simAlpha.s(alpha.s,alpha.u,z,Ss[j],llnull,llresp,j,rate.s)
			alpha.s<-res[[1]]
			rate.s[j]<-res[[2]][j]
			res<-MIMOSA:::simAlpha.u(alpha.s,alpha.u,z,Su[j],llnull,llresp,j,rate.u)
			alpha.u<-res[[1]]
			rate.u[j]<-res[[2]][j]
		}
		if(i%%tune==0&!fixed){
			if(any(rate.u/tune>0.55|rate.u/tune<0.35|rate.s/tune>0.55|rate.s/tune<0.35)){
				cat(rate.s/tune," ",rate.u/tune,"\n")
				w.u<-which(rate.u/tune>0.55|rate.u/tune<0.35)
				w.s<-which(rate.s/tune>0.55|rate.s/tune<0.35)
				for(j in w.s){
					Ss[j]<-Ss[j]*0.5+0.5*sd(f[(i-tune+1):(i-1),j])
					Ss[j]<-Ss[j]*((rate.s[j]/tune)/0.44)^2
				}
				for(j in w.u){
					Su[j]<-Su[j]*0.5+0.5*sd(f[(i-tune+1):(i-1),length(rate.u)+j])
					Su[j]<-Su[j]*((rate.u[j]/tune)/0.44)^2
				}
				Ss[Ss==0]<-1
				Su[Su==0]<-1
				rate.u<-rate.s<-rep(0,length(rate.s))
			}else {
				cat("Fixing: ",rate.s/tune," ",rate.u/tune,"\n")
				fixed<-TRUE
				assign("i",1)
			}
			assign("i",1)
		}
		z<-MIMOSA:::simZ(q=q,alpha.s=alpha.s,alpha.u=alpha.u,llnull,llresp)
		if(i==burn&fixed){
			sz=z
		}else if(i>burn&fixed){
			sz=(sz+z/(i-burn))*((i-burn)/(i-burn+1))	
		}
		q<-MIMOSA:::simQ(z=z)
		f[i,]<-c(alpha.s,alpha.u,q)
		i<-i+1
		if(i==iter&fixed){
			break
		}
	}
	write.table(f,file="mcmc.dat")
	return(list(z=sz,mcmc=mcmc(f[seq(burn,iter,by=thin),])))
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

.fitMCMC<-function(data,inits,iter, burn, thin,tune=100,outfile="mcmc.dat",alternative="greater",FAST=FALSE,LOWER=0.15,UPPER=0.5){
	alternative<-match.arg(alternative,c("greater","not equal"))
	data<-icsdata2mvicsdata(data)
	#If the alternative hypothesis is one-sided, then compute a filter for pu>ps and pass that to the MCMC code
	if(alternative=="greater"){
		ps<-t(do.call(cbind,apply(data$n.stim,1,function(x)(data.frame(prop.table(x))[-1L,,drop=FALSE]))))
		pu<-t(do.call(cbind,apply(data$n.unstim,1,function(x)(data.frame(prop.table(x))[-1L,,drop=FALSE]))))		
		filter<-sapply(1:nrow(ps),function(i)all(ps[i,]<pu[i,]))
		if(!FAST){
			FILTER<-TRUE
		}else{
			FILTER<-FALSE
		}
	}else{
		filter<-rep(FALSE,nrow(data$n.stim))
		FILTER<-FALSE;
		FAST<-FALSE
	}
	result<-.Call("fitMCMC",as.matrix(data$n.stim),as.matrix(data$n.unstim),as.vector(inits$alpha.s),as.vector(inits$alpha.u),as.vector(inits$q),as.matrix(inits$z),as.vector(iter),as.vector(burn),as.vector(thin),as.numeric(tune),as.character(outfile),as.vector(filter),as.logical(FAST),as.logical(FILTER),as.numeric(LOWER),as.numeric(UPPER),package="MIMOSA")
	if(inherits(result,"character")){
		return(result)
	}
	result$z<-cbind(result$z,1-result$z)
	result$getmcmc<-function(x=outfile){
		mcmc(read.table(x,sep="\t",header=T));
	}
	attr(result,"class")<-c(attr(result,"class"),"MDMixResult")
	attr(result,"pData")<-attr(data,"pData")
	result
}
#test<-function(iter,burn,w.stepS=0.01,w.stepU=0.001,sigmaS=20,sigmaU=20){
#	set.seed(5)
#	alpha.s=c(1000-120,40,40,40,40)
#	alpha.u=c(1000,10,10,10,10)
#	dat<-MIMOSA:::simMD(alpha.s=c(1000-120,40,40,40,40),alpha.u=c(1000,10,10,10,10),w=0.25)
#	
#	
#	llnull<-cmpfun(MIMOSA:::makeLogLikeNULLComponent(dat$n.stim,dat$n.unstim))
#	llresp<-cmpfun(MIMOSA:::makeLogLikeRespComponent(dat$n.stim,dat$n.unstim))
#	q<-runif(1)
#	z<-MIMOSA:::simZ(q,w.s,w.u,S.stim,S.unstim,llnull,llresp)
#	f<-matrix(0,nrow=iter,ncol=length(w.s)+length(w.u)+1)
#	cn<-c(paste("alpha.s",1:length(w.s),sep=""),paste("alpha.u",1:length(w.u),sep=""),"q")
#	colnames(f)<-cn
#	for(i in 1:iter){
#		#every 100 iterations do some tuning until we reach a good acceptance ratio
#		if(i%%tune==0){
#			#browser()
#			rate.s<-rate.s/tune
#			rate.u<-rate.u/tune
#			if(!(rate.s<=0.26&rate.s>=0.22&rate.u>=0.22&rate.u<=0.26)&!fixed){
#				cat("rate.s=",rate.s," rate.u=",rate.u,"\n");
#				Cs=Cs*abs(qnorm(0.24/2)/qnorm(rate.s/2))
#				Cu=Cu*abs(qnorm(0.24/2)/qnorm(rate.u/2))
#				S.s<-S.s*0.75+0.25*(cov(f[(i-tune+1):i,1:length(alpha.s)]))
#				S.u<-S.u*0.75+0.25*(cov(f[(i-tune+1):i,(length(alpha.s)+1):(length(alpha.s)+length(alpha.u))]))
#				rate.s<-rate.u<-0
#			} else if(!fixed) {
#				cat("rate.s=",rate.s," rate.u=",rate.u,"\n");
#				cat("Fixing the step size and sampling ", iter, " samples with a burn-in of ", burn,"\n")
#				fixed<-TRUE
#				i<-1
#			}else{
#				cat(".")
#			}
#		}
#		ret<-MIMOSA:::simAlpha.s(alpha.s,alpha.u,z,S.s,C=Cs,llnull,llresp,rate.s)
#		alpha.s<-ret[[1]]
#		rate.s<-ret[[2]]
#		ret<-MIMOSA:::simAlpha.u(alpha.s,alpha.u,z,S.u,C=Cu,llnull,llresp,rate.u)
#		alpha.u<-ret[[1]]
#		rate.u<-ret[[2]]
#		z<-MIMOSA:::simZ(q=q,w.s=alpha.s*sum(alpha.s),w.u=alpha.u*sum(alpha.u),S.stim=sum(alpha.s),S.unstim=sum(alpha.u),llnull,llresp)
#		if(i==burn){
#				sz=z
#		}else if(i>burn){
#			sz=(sz+z/(i-1))*((i-1)/i)	
#		}
#		q<-MIMOSA:::simQ(z=z)
#		f[i,]<-c(alpha.s,alpha.u,q)
#	}
#	write.table(f,file="mcmc.dat")
#	return(list(z=sz,mcmc=mcmc(f[burn:iter,])))
#}