# TODO: Implements the multinomial-dirichelt model for two cytokines
# 
# Author: Greg Finak
###############################################################################



#Observed data
#n.stim - vector of counts from the stimulated sample
#n.unstim - vector of counts from the unstimulated sample
#alpha.unstim
#alpha.stim

#compute the log of the factorial term
#' 
#' @param x observed data (counts)
#' @returnType 
#' @return 
#' @author Greg Finak
#' @export
lfctrl<-function(x){
	lfactorial(sum(x))-sum(lfactorial(x))
}

#compute the log of the k-dimensional beta function
#' compute the log of the k-dimensional beta function
#' @param alpha k-vector of beta function parameters
#' @returnType numeric
#' @return 
#' @author Greg Finak
#' @export
lkbeta<-function(alpha){
	sum(lgamma(alpha))-lgamma(sum(alpha))
}



makeLogLikeDirichlet<-function(data)
{
	P<-nrow(data)
	loglike<-function(x)
	{
		ll<-P*lgamma(sum(x))-P*sum(lgamma(x))+sum((x-1)*log(data))
		return(ll)
	}
}
makeGradDirichlet<-function(data){
	P<-nrow(data)
	data<-matrix(data,ncol=4)
	grad<-function(x)
	{
		g<-P*digamma(sum(x))-P*(digamma(x))+P*colMeans(log(data))
		return(g)
	}
}
makeHessDirichlet<-function(data){
	P<-nrow(data)
	data<-matrix(data,ncol=4)
	hess<-function(x)
	{
		O<-P*trigamma(sum(x))
		H<-matrix(O,nrow=length(x),ncol=length(x))
		H<-H-P*diag(trigamma(x))
		return(H)
	}
}




makeLogLikeNULLComponent<-function(data.stim,data.unstim){
	N.s<-rowSums(data.stim)
	N.u<-rowSums(data.unstim)
	data<-data.stim+data.unstim
	loglike<-function(x){
		(lgamma(sum(x))-lgamma(rowSums(t(x+t(data))))+rowSums(t(t(lgamma(t(x+t(data))))-lgamma(x)))+lgamma(N.s+1)+lgamma(N.u+1) - rowSums(lgamma(data.stim+1)+lgamma(data.unstim+1)))
	}
	return(loglike)
}
makeLogLikeRespComponent<-function(data.stim,data.unstim){
	N.s<-rowSums(data.stim)
	N.u<-rowSums(data.unstim)
	stim.ind<-1:4
	unstim.ind<-5:8
	data<-cbind(data.stim,data.unstim)
	loglike<-function(x){
		rowSums(t(t(lgamma(t(x[stim.ind]+t(data.stim)))+lgamma(t(x[unstim.ind]+t(data.unstim)))-lgamma(data.stim+1)-lgamma(data.unstim+1))-lgamma(x[stim.ind])-lgamma(x[unstim.ind])))+
				lgamma(sum(x[stim.ind]))+lgamma(sum(x[unstim.ind]))-lgamma(rowSums(t(x[unstim.ind]+t(data.unstim))))-lgamma(rowSums(t(x[stim.ind]+t(data.stim))))+lgamma(N.s+1)+lgamma(N.u+1)
	}
	return(loglike)
}



makeGradientNULLComponent<-function(data.stim,data.unstim,z=NULL){
	data<-data.stim+data.unstim
	if(is.null(z)){
		z<-matrix(1,nrow=nrow(data.unstim),ncol=2)
	}
	grad<-function(x){
		gr<-(digamma(sum(x))-digamma(colSums((t(data)+x)))+t(digamma(t(x+t(data))))-digamma(x))
		gr%*%z[,1]
	}
}
makeGradientHalfComponent<-function(data,z=NULL){
	if(is.null(z)){
		z<-matrix(1,nrow=nrow(data),ncol=1)
	}
	grad<-function(x){
		gr<-(digamma(sum(x))-digamma(colSums((t(data)+x)))+t(digamma(t(x+t(data))))-digamma(x))
		gr%*%z
	}
}

makeGradientRespComponent<-function(data.stim,data.unstim,z=NULL){
	stim.ind<-1:4
	unstim.ind<-5:8
	if(is.null(z)){
		z<-matrix(1,nrow=nrow(data.stim),ncol=2)
	}
	gs<-makeGradientHalfComponent(data.stim,z=z[,2]);gu<-makeGradientHalfComponent(data.unstim,z=z[,1]);
	grad<-function(x){
		c(gs(x[stim.ind]),gu(x[unstim.ind]))
	}
	return(grad)
}

makeHessianNULLComponent<-function(data.stim,data.unstim,z){
	data<-data.stim+data.unstim
	if(is.null(z)){
		z<-matrix(1,nrow=nrow(data.unstim),ncol=2)
	}
	hess<-function(x){
		H<-matrix(sum(z[,1]*(trigamma(sum(x))-trigamma(colSums(t(data)+x)))),ncol=length(x),nrow=length(x))
		D<-(trigamma(x+t(data))-trigamma(x))%*%z[,1]
		diag(H)<-diag(H)+D
		return(H)
	}
	return(hess)
}

makeHessianHalfComponent<-function(data,z=NULL){
	if(is.null(z)){
		z<-matrix(1,nrow=nrow(data),ncol=1)
	}
	hess<-function(x){
		H<-matrix(sum(z*(trigamma(sum(x))-trigamma(colSums(t(data)+x)))),ncol=length(x),nrow=length(x))
		D<-(trigamma(x+t(data))-trigamma(x))%*%z
		diag(H)<-diag(H)+D
		return(H)
	}
}
makeHessianRespComponent<-function(data.stim,data.unstim,z){
	stim.ind<-1:4
	unstim.ind<-5:8
	if(is.null(z)){
		z<-matrix(1,nrow=data.stim,ncol=2)
	}
	hs<-makeHessianHalfComponent(data.stim,z[,2]);hu<-makeHessianHalfComponent(data.unstim,z[,1]);
	hess<-function(x){
		H<-matrix(0,nrow=length(c(stim.ind,unstim.ind)),ncol=length(c(stim.ind,unstim.ind)))
		H[stim.ind,stim.ind]<-hs(x[stim.ind]);
		H[unstim.ind,unstim.ind]<-hu(x[unstim.ind]);
		return(H)
	}
	return(hess)
}


#a,b are the non-responder and responder hyperparameters
testComponentOptimization <- function(a, b)
{
	a<-c(100000,40,20,100)
	b<-c(100000,900,100,100)
	pu<-rdirichlet(100,a)
	ps<-rdirichlet(50,a)
	ps<-rbind(ps,rdirichlet(50,b))
	nu<-t(sapply(seq_along(1:nrow(pu)),function(i)rmultinom(1,runif(1,1e5,1.5e5),pu[i,])))
	ns<-t(sapply(seq_along(1:nrow(ps)),function(i)rmultinom(1,runif(1,1e5,1.5e5),ps[i,])))
	
	llnull<-makeLogLikeNULLComponent(ns,nu)
	llresp<-makeLogLikeRespComponent(ns,nu)
	gnull<-makeGradientNULLComponent(ns,nu)
	gresp<-makeGradientRespComponent(ns,nu)
	hessnull<-makeHessianNULLComponent(ns,nu)
	hessresp<-makeHessianRespComponent(ns,nu)
	
	guess<-c(colMeans(ns/rowSums(ns)),colMeans(nu/rowSums(nu)))
	iter<-1
	ll<-rep(0,100000)
	repeat{
		#new<-guess-ginv(hessresp(as.vector(guess)))%*%gresp(as.vector(guess))
		t<-try(solve(hessresp(as.vector(guess)))%*%gresp(as.vector(guess))+c(rep(0,4),solve(hessnull(as.vector(guess[5:8])))%*%gnull(as.vector(guess[5:8]))))
		inherits(t,"try-error")
		t<-ginv(hessresp(as.vector(guess)))%*%gresp(as.vector(guess))+c(rep(0,4),ginv(hessnull(as.vector(guess[5:8]))%*%gnull(as.vector(guess[5:8]))))
		new<-guess-t
		ll[iter]<--sum(llresp(as.vector(new)))
		new<-abs(new)
		if(sum(abs(new-as.matrix(guess)))<sqrt(.Machine$double.eps)|iter>999){
			break
		}
		guess<-new
		iter<-iter+1
	}
}





#EM algorithm fitting the multinomial dirichlet mixture
#' 
#' @param data The observed data
#' @param modelmatrix a model matrix specifying which components should be computed
#' @returnType 
#' @return 
#' @author Greg Finak
#' @export
MDMix<-function(data=NULL,modelmatrix=NULL){
	inits<-initMDMix(data=data,modelmatrix=modelmatrix)
	
}

#initialization of z's and parameters for the multinomial dirichlet mixture
#' 
#' @param data The observed data
#' @param modelmatrix a model matrix specifying which components should be computed
#' @returnType 
#' @return 
#' @author Greg Finak
#' @export
initMDMix<-function(data=NULL,modelmatrix=NULL){
	unstim<-data$n.unstim
	stim<-data$n.stim
	
	#fisher's exact test of all the marginals
	mm<-do.call(cbind,lapply(2:4,function(i)apply(cbind(rowSums(unstim[,-i]),unstim[,i],rowSums(stim[,-i]),stim[,i]),1,function(x)fisher.test(matrix(x,2),alternative="greater")$p.value)))<0.01
	
	#observations with no significant marginals belong to the null component.
	mm<-apply(mm,1,function(x)all(!x))
	
	#The rest are from the responder component
	#construct the z-matrix
	z<-matrix(0,length(mm),2)
	z[mm,1]<-1
	z[!mm,2]<-1
	
	#estimate hyperparamters.
	pu.u<-unstim/rowSums(unstim)
	pu.s<-stim[which(z[,1]==1),]/rowSums(stim[which(z[,1]==1),])
	pu<-(rbind(pu.u,pu.s))
	ps.s<-stim[z[,2]==1,]/rowSums(stim[z[,2]==1,])
	alpha.u<-colMeans(pu)
	alpha.s<-colMeans(ps.s)
	
	
	guess<-c(alpha.s,alpha.u)
	w<-colSums(z)/sum(z)

	#EM
	LL<-NULL
	repeat{
		if(is.null(LL)){
			last<--.Machine$double.xmax
		}
		#update parameters
		llnull<-makeLogLikeNULLComponent(stim,unstim)
		llresp<-makeLogLikeRespComponent(stim,unstim)
		gnull<-makeGradientNULLComponent(stim,unstim,z)
		gresp<-makeGradientRespComponent(stim,unstim,z)
		hessnull<-makeHessianNULLComponent(stim,unstim,z)
		hessresp<-makeHessianRespComponent(stim,unstim,z)
		
		iter<-1
		ll<-rep(0,10000)
		repeat{
			#t<-try(solve(hessresp(as.vector(guess)))%*%gresp(as.vector(guess))+c(rep(0,4),solve(hessnull(as.vector(guess[5:8])))%*%gnull(as.vector(guess[5:8]))),silent=TRUE)
			#if(inherits(t,"try-error"))
				t<-ginv(hessresp(as.vector(guess)))%*%gresp(as.vector(guess))+c(rep(0,4),ginv(hessnull(as.vector(guess[5:8]))%*%gnull(as.vector(guess[5:8]))))
			new<-guess-t
			ll[iter]<--sum(llnull(new[5:8])*z[,1]+llresp(as.vector(new))*z[,2])
			new<-abs(new)
			if(sum(abs(new-as.matrix(guess)))<.Machine$double.eps^(1/3)|iter>9999){
				break
			}
			guess<-new
			iter<-iter+1
		}
		#compute z's and w's
		
		z1<-(exp(llnull(new[5:8]))*w[1])/(exp(llnull(new[5:8]))*w[1]+exp(llresp(new))*w[2])
		z2<-1-z1
		z<-cbind(z1,z2)
		w<-colSums(z)/sum(z)
		cll<--sum(llnull(new[5:8])*z[,1]+llresp(new)*z[,2])
		LL<-c(LL,cll)
		browser()
		if(abs(last-cll)<1e-3){
			break;
		}
		last<-cll
	}
	return(list(llnull=llnull,llresp=llresp,gresp=gresp,hresp=hessresp,gnull=gnull,hnull=hessnull,z=z,LL=LL,par=new))
}

estAlpha<-function(est,dat=NULL) {
	NLL = -sum(log(ddirichlet(dat,est)))
	return(NLL)	
}

#' extracts bifunctional cytokine data from and ICS object given the two marginals (A, B) and A||B for stimualted and unstimulated. Used for the multinomial dirichlet model. ORDER OF CYTOKINES MATTERS!
#' @param ics 
#' @param cytokineA 
#' @param cytokineB 
#' @param or 
#' @param stim 
#' @param control 
#' @param subset 
#' @param shrink 
#' @param scl 
#' @returnType 
#' @return 
#' @author Greg Finak
#' @export
extractDataMultinomDir<-function(ics=NULL,cytokineA=NULL,cytokineB=NULL,or=NULL,stim=NULL,control=NULL,subset=NULL,shrink=1,scl=1){
	if(any(c(is.null(ics),is.null(cytokineA),is.null(stim),is.null(cytokineB),is.null(or),is.null(control),is.null(subset)))){
		stop("All arguments must be non--null");
	}
	A<-extractData(ics,control=control,stim=stim,subset=c(subset,cytokineA))
	B<-extractData(ics,control=control,stim=stim,subset=c(subset,cytokineB))
	OR<-extractData(ics,control=control,stim=stim,subset=c(subset,or))
	AND<-A+B-OR
	AND[,"Ns"]<-A[,"Ns"]+A[,"ns"]-AND[,"ns"]
	AND[,"Nu"]<-A[,"Nu"]+A[,"nu"]-AND[,"nu"]
	ds<-AND[,"ns"]
	bs<-A[,"ns"]-AND[,"ns"]
	cs<-B[,"ns"]-ds
	as<-A[,"Ns"]-cs
	
	du<-AND[,"nu"]
	bu<-A[,"nu"]-AND[,"nu"]
	cu<-B[,"nu"]-du
	au<-A[,"Nu"]-cu
	
	n.unstim<-cbind("u1"=au,"u2"=bu,"u3"=cu,"u4"=du)
	n.stim<-cbind("s1"=as,"s2"=bs,"s3"=cs,"s4"=ds)
	r<-(list(n.stim,n.unstim))
	names(r)<-c("n.stim","n.unstim")
	attr(r,"cytokines")<-c(cytokineA,cytokineB)
	attr(r,"subset")<-subset
	attr(r,"stim")<-stim
	class(r)<-c("list","MDlist")
	return(r)
}
setOldClass("MDlist")

setMethod(show,"MDlist",function(x){
	cat("Cytokines ",attr(x,"cytokines"),"\n")
	cat("Stimulation ",attr(x,"stim"),"\n")
	cat("Subset ", attr(x,"subset"),"\n")
	cat("Number of obs: ",nrow(x[[1]]),"\n")
})

print.MDlist<-function(x){
			cat("Cytokines ",attr(x,"cytokines"),"\n")
			cat("Stimulation ",attr(x,"stim"),"\n")
			cat("Subset ", attr(x,"subset"),"\n")
			cat(nrow(x[[1]])," observations","\n")
}
