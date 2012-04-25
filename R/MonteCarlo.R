# TODO: Having issues with newton-raphson and optim for the dirichlet-multinomial model
# Implement a Gibbs sampler.
# Author: finak
###############################################################################


sProposal<-function(s,sigma=10){
	return(rnorm(1,mean=s,sd=sigma))
}
alphaProposal<-function(alpha){
	sapply(alpha,function(x)rnorm(1,mean=x,sd=10))
}

lddirichlet<-function(x,alpha){
	logD<-lkbeta(alpha)
	s<-sum(log(x)*(alpha-1))
	s-logD
}
simAlpha.s<-function(alpha.s,alpha.u,z){
	a<-llnull(c(alpha.s,alpha.u))
	b<-llresp(c(alpha.s,alpha.u))
	
	old<-sum(a*z[,1]+b*z[,2]+sum(dexp(alpha.s,rate=0.001,log=TRUE)))
	alpha.s.prop<-alphaProposal(alpha.s)
	if(any(alpha.s.prop<0)){
		return(alpha.s)
	}
	
	a<-llnull(c(alpha.s.prop,alpha.u))
	b<-llresp(c(alpha.s.prop,alpha.u))
	
	new<-sum(a*z[,1]+b*z[,2]+sum(dexp(alpha.s.prop,rate=0.001,log=TRUE)))
	
	if(log(runif(1))<=(new-old)&is.finite(new-old)){
		return(alpha.s.prop)
	}else{
		return(alpha.s)
	}
	
}
simAlpha.u<-function(alpha.s,alpha.u,z){
	a<-llnull(c(alpha.s,alpha.u))
	b<-llresp(c(alpha.s,alpha.u))
	
	old<-sum(a*z[,1]+b*z[,2]+sum(dexp(alpha.u,rate=0.001,log=TRUE)))
	alpha.u.prop<-alphaProposal(alpha.u)
	if(any(alpha.u.prop<0)){
		return(alpha.u)
	}
	
	a<-llnull(c(alpha.s,alpha.u.prop))
	b<-llresp(c(alpha.s,alpha.u.prop))
	
	new<-sum(a*z[,1]+b*z[,2]+sum(dexp(alpha.u.prop,rate=0.001,log=TRUE)))
	
	if(log(runif(1))<=(new-old)&is.finite(new-old)){
		return(alpha.u.prop)
	}else{
		return(alpha.u)
	}
	
}

wProposal<-function(w,c=0.2){
	w+scale(runif(length(w),-c,c),scale=FALSE)
}

#TODO: prior on Q so that rbeta is positive

simQ<-function(z){
	rbeta(1,sum(z[,1])+0.01,sum(z[,2])+0.01)
}

simZ<-function(q,w.s,w.u,S.stim,S.unstim){
	a<-llnull(c(w.s*S.stim,w.u*S.unstim))+log(q) 
	b<-llresp(c(w.s*S.stim,w.u*S.unstim))+log(1-q)
	den<-apply(cbind(a,b),1,function(x)log(sum(exp(x-max(x))))+max(x))
	p<-exp(a-den)
	z<-sapply(p,function(p)sample(c(0,1),1,prob=c(1-p,p),replace=FALSE))
	return(cbind(z,1-z))
}

simWs<-function(w.s,w.u,z,S.stim,S.unstim){
	a<-llnull(c(w.s*S.stim,w.u*S.unstim))
	b<-llresp(c(w.s*S.stim,w.u*S.unstim))
	
	old<-sum(a*z[,1]+b*z[,2]+lddirichlet(matrix(w.s,nrow=1),rep(1,length(w.s))))
	w.s.prop<-wProposal(w.s,c=w.step)
	if(any(w.s.prop<0|w.s.prop>1)){
		return(w.s)
	}
	
	a<-llnull(c(w.s.prop*S.stim,w.u*S.unstim))
	b<-llresp(c(w.s.prop*S.stim,w.u*S.unstim))
	
	new<-sum(a*z[,1]+b*z[,2]+lddirichlet(matrix(w.s.prop,nrow=1),rep(1,length(w.s.prop))))
	
	if(log(runif(1))<=(new-old)&is.finite(new-old)){
		return(w.s.prop)
	}else{
		return(w.s)
	}
}

simWu<-function(w.s,w.u,z,S.stim,S.unstim){
	a<-llnull(c(w.s*S.stim,w.u*S.unstim))
	b<-llresp(c(w.s*S.stim,w.u*S.unstim))
	
	old<-sum(a*z[,1]+b*z[,2]+lddirichlet(matrix(w.u,nrow=1),rep(1,length(w.u))))
	w.u.prop<-wProposal(w.u,c=w.step)
	if(any(w.u.prop<0|w.u.prop>1)){
		return(w.u)
	}
	
	a<-llnull(c(w.s*S.stim,w.u.prop*S.unstim))
	b<-llresp(c(w.s*S.stim,w.u.prop*S.unstim))
	
	new<-sum(a*z[,1]+b*z[,2]+lddirichlet(matrix(w.u.prop,nrow=1),rep(1,length(w.u.prop))))
	
	#wierd, should be new/old, but that doesn't work out right
	if(log(runif(1))<=(new-old)&is.finite(new-old)){
		return(w.u.prop)
	}else{
		return(w.u)
	}
}

simS.stim<-function(w.s,w.u,z,S.stim,S.unstim){
	a<-llnull(c(w.s*S.stim,w.u*S.unstim))
	b<-llresp(c(w.s*S.stim,w.u*S.unstim))
	
	old<-sum(a*z[,1]+b*z[,2]+dexp(S.stim,rate=0.001,log=TRUE))
	S.stim.prop<-sProposal(S.stim,sigma=sigma)
	
	a<-llnull(c(w.s*S.stim.prop,w.u*S.unstim))
	b<-llresp(c(w.s*S.stim.prop,w.u*S.unstim))
	
	new<-sum(a*z[,1]+b*z[,2]+dexp(S.stim.prop,rate=0.001,log=TRUE))
	
	if(log(runif(1))<=(new-old)&is.finite(new-old)){
		return(S.stim.prop)
	}else{
		return(S.stim)
	}
	
}
simS.unstim<-function(w.s,w.u,z,S.stim,S.unstim){
	a<-llnull(c(w.s*S.stim,w.u*S.unstim))
	b<-llresp(c(w.s*S.stim,w.u*S.unstim))
	
	old<-sum(a*z[,1]+b*z[,2]+dexp(S.unstim,rate=0.001,log=TRUE))
	S.unstim.prop<-sProposal(S.unstim,sigma=sigma)
	
	a<-llnull(c(w.s*S.stim,w.u*S.unstim.prop))
	b<-llresp(c(w.s*S.stim,w.u*S.unstim.prop))
	
	new<-sum(a*z[,1]+b*z[,2]+dexp(S.unstim.prop,rate=0.001,log=TRUE))
	
	if(log(runif(1))<=(new-old)&is.finite(new-old)){
		return(S.unstim.prop)
	}else{
		return(S.unstim)
	}
}

test<-function(x){
	set.seed(5)
	dat<-MIMOSA:::simMD(alpha.s=c(1000-120,40,40,40,40),alpha.u=c(1000,10,10,10,10),w=0.25)
	
	
	llnull<-MIMOSA:::makeLogLikeNULLComponent(dat$n.stim,dat$n.unstim)
	llresp<-MIMOSA:::makeLogLikeRespComponent(dat$n.stim,dat$n.unstim)
	
	S.stim<-1000
	S.unstim<-1000
	w.u<-rowMeans(apply(dat$n.unstim,1,prop.table))
	w.s<-rowMeans(apply(dat$n.stim,1,prop.table))
	q<-runif(1)
	z<-MIMOSA:::simZ(q,w.s,w.u,S.stim,S.unstim)
	sigma<-10
	w.step<-0.01
	iter<-1000000
	f<-file("mcmc.dat",open="wt",blocking=FALSE)
	writeLines(paste(c(paste("ws",1:length(w.s),sep="",collapse="\t"),paste("wu",1:length(w.u),sep="",collapse="\t"),"S.stim","S.unstim",paste("z",1:nrow(dat[[1]]),sep="",collapse="\t"),"q","\n"),collapse="\t"),f)
	writeLines(text=as.character(c(w.s,w.u,S.unstim,S.stim,z[,1],q,"\n")),sep="\t",f)
	#coda needs the order reversed. Bizarre
	for(i in (iter-1):1){
		w.s<-MIMOSA:::simWs(w.s=w.s,w.u=w.u,z=z,S.stim=S.stim,S.unstim=S.unstim)
		w.u<-MIMOSA:::simWu(w.s=w.s,w.u=w.u,z=z,S.stim=S.stim,S.unstim=S.unstim)
		S.unstim<-MIMOSA:::simS.unstim(w.s=w.s,w.u=w.u,z=z,S.stim=S.stim,S.unstim=S.unstim)
		S.stim<-MIMOSA:::simS.stim(w.s=w.s,w.u=w.u,z=z,S.stim=S.stim,S.unstim=S.unstim)
		z<-MIMOSA:::simZ(q=q,w.s=w.s,w.u=w.u,S.stim=S.stim,S.unstim=S.unstim)
		q<-MIMOSA:::simQ(z=z)
		writeLines(text=as.character(c(w.s,w.u,S.unstim,S.stim,z[,1],q,"\n")),sep="\t",f)
	}
	close(f)
	saved
}