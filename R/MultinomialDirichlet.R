# TODO: Implements the multinomial-dirichelt model for two cytokines
# 
# Author: Greg Finak
###############################################################################



#Observed data
#n.stim - vector of counts from the stimulated sample
#n.unstim - vector of counts from the unstimulated sample
#alpha.unstim
#alpha.stim

##compute the log of the k-dimensional beta function
##' compute the log of the k-dimensional beta function
##' @param alpha k-vector of beta function parameters
##' @returnType numeric
##' @return 
##' @author Greg Finak
##' @export
lkbeta<-function(alpha){
	sum(lgamma(alpha))-lgamma(sum(alpha))
}


makeLogLikeNULLComponent<-function(data.stim,data.unstim){
	data<-data.stim+data.unstim
	ll<-function(x){
		a<-x[(length(x)/2+1):length(x)]
		apply(data,1,function(y)lkbeta(y+a))-lkbeta(a)
		#-rowSums(lfactorial(data.stim))-rowSums(lfactorial(data.unstim))+lfactorial(rowSums(data.stim))+lfactorial(rowSums(data.unstim))
	}
}
makeLogLikeRespComponent<-function(data.stim,data.unstim){
	ll<-function(x){
		a<-x[(length(x)/2+1):length(x)]
		b<-x[1:(length(x)/2)]
		apply(data.stim,1,function(y)lkbeta(y+b))+apply(data.unstim,1,function(y)lkbeta(y+a))-lkbeta(b)-lkbeta(a)
		#-rowSums(lfactorial(data.stim))-rowSums(lfactorial(data.unstim))+lfactorial(rowSums(data.stim))+lfactorial(rowSums(data.unstim))
	}
}

#makeLogLikeNULLComponent<-function(data.stim,data.unstim){
#	N.s<-rowSums(data.stim)
#	N.u<-rowSums(data.unstim)
#	data<-data.stim+data.unstim
#	loglike<-function(x){
#		x<-x[(ncol(data.stim)+1):(ncol(data.stim)*2)]
#		lgamma(sum(x))-lgamma(colSums(t(data)+x))+colSums(lgamma(t(data)+x)-lgamma(x))+lgamma(N.s+1)+lgamma(N.u+1)-colSums(t(lgamma(data.stim+1)+lgamma(data.unstim+1)))
#	}
#	return(loglike)
#}
#makeLogLikeRespComponent<-function(data.stim,data.unstim){
#	N.s<-rowSums(data.stim)
#	N.u<-rowSums(data.unstim)
#	stim.ind<-1:ncol(data.stim)
#	unstim.ind<-(ncol(data.stim)+1):((ncol(data.stim)*2))
#	loglike<-function(x){
#		colSums(lgamma(x[unstim.ind]+t(data.unstim))+lgamma(x[stim.ind]+t(data.stim))-lgamma(x[stim.ind])-lgamma(x[unstim.ind])-lgamma(t(data.stim)+1)-lgamma(t(data.unstim)+1))+lgamma(sum(x[stim.ind]))+lgamma(sum(x[unstim.ind]))-lgamma(colSums(t(data.unstim)+x[unstim.ind]))-lgamma(colSums(t(data.stim)+x[stim.ind]))+lgamma(N.s+1)+lgamma(N.u+1)
#	}
#	return(loglike)
#}



makeGradientNULLComponent<-function(data.stim,data.unstim,z=NULL){
	data<-data.stim+data.unstim
	if(is.null(z)){
		z<-matrix(1,nrow=nrow(data.unstim),ncol=2)
	}
	grad<-function(x){
		x<-x[(ncol(data.stim)+1):(ncol(data.stim)*2)]
		gr<-(t(digamma(sum(x))-digamma(colSums((t(data)+x)))+digamma(t(x+t(data))))-digamma(x))
		gr<-gr%*%z[,1]
		return(c(rep(0,ncol(data.stim)),gr))
	}
}
makeGradientHalfComponent<-function(data,z=NULL){
	if(is.null(z)){
		z<-matrix(1,nrow=nrow(data),ncol=1)
	}
	grad<-function(x){
		gr<-(t(digamma(sum(x))-digamma(colSums((t(data)+x)))+(digamma(t(x+t(data)))))-digamma(x))
		gr%*%as.vector(z)
	}
}

makeGradientRespComponent<-function(data.stim,data.unstim,z=NULL){
	stim.ind<-1:ncol(data.stim)
	unstim.ind<-(ncol(data.stim)+1):(ncol(data.stim)*2)
	if(is.null(z)){
		z<-matrix(1,nrow=nrow(data.stim),ncol=2)
	}
	gs<-makeGradientHalfComponent(data.stim,z=z[,2]);gu<-makeGradientHalfComponent(data.unstim,z=z[,2]);
	grad<-function(x){
		c(gs(x[stim.ind]),gu(x[unstim.ind]))
	}
	return(grad)
}

makeHessianNULLComponent<-function(data.stim,data.unstim,z=NULL){
	data<-data.stim+data.unstim
	if(is.null(z)){
		z<-matrix(1,nrow=nrow(data.unstim),ncol=2)
	}
	hess<-function(x){
		x<-x[(ncol(data.stim)+1):(ncol(data.stim)*2)]
#		H<-matrix(sum(z[,1]*(trigamma(sum(x))-trigamma(colSums(t(data)+x)))),ncol=length(x),nrow=length(x))
#		D<-(trigamma(x+t(data))-trigamma(x))%*%z[,1]
#		diag(H)<-diag(H)+D
#		HH<-matrix(0,ncol(data)*2,ncol(data)*2)
#		HH[(ncol(data.stim)+1):(ncol(data.stim)*2),(ncol(data.stim)+1):(ncol(data.stim)*2)]<-H
#		return(HH)
		H<- (trigamma(sum(x)) - trigamma(colSums(t(data) +x)))
		D<-(trigamma(x + t(data)) - trigamma(x))
		HH<-matrix(0,ncol(data)*2,ncol(data)*2)
		H<-sapply(1:length(H),function(i){D[,i]<-D[,i]+H[i];M<-matrix(0,ncol(data)*2,ncol(data)*2);inds<-(ncol(data.stim)+1):(ncol(data.stim)*2);M[inds,inds]<-H[i];diag(M)<-c(rep(0,ncol(data)),D[,i]);list(M*z[i,1])})
		for(i in seq_along(H)){
			HH<-HH+H[[i]]
		}
		return(HH)
	}
	return(hess)
}

makeHessianHalfComponent<-function(data,z=NULL){
	if(is.null(z)){
		z<-matrix(1,nrow=nrow(data),ncol=1)
	}
	hess<-function(x){
#		H<-matrix(sum(z*(trigamma(sum(x))-trigamma(colSums(t(data)+x)))),ncol=length(x),nrow=length(x))
#		D<-(trigamma(x+t(data))-trigamma(x))%*%as.vector(z)
#		diag(H)<-diag(H)+D
#		return(H)
		H<- (trigamma(sum(x)) - trigamma(colSums(t(data) +x)))
		D<-(trigamma(x + t(data)) - trigamma(x))
		HH<-matrix(0,ncol(data),ncol(data))
		H<-sapply(1:length(H),function(i){D[,i]<-D[,i]+H[i];M<-matrix(H[i],ncol(data),ncol(data));diag(M)<-D[,i];list(M*z[i])})
		for(i in seq_along(H)){
			HH<-HH+H[[i]]
		}
		return(HH)
	}
}
makeHessianRespComponent<-function(data.stim,data.unstim,z=NULL){
	stim.ind<-1:ncol(data.stim)
	unstim.ind<-(ncol(data.stim)+1):(ncol(data.stim)*2)
	if(is.null(z)){
		z<-matrix(1,nrow=nrow(data.stim),ncol=2)
	}
	hs<-makeHessianHalfComponent(data.stim,z[,2]);hu<-makeHessianHalfComponent(data.unstim,z[,2]);
	hess<-function(x){
		H<-matrix(0,nrow=length(c(stim.ind,unstim.ind)),ncol=length(c(stim.ind,unstim.ind)))
		H[stim.ind,stim.ind]<-hs(x[stim.ind]);
		H[unstim.ind,unstim.ind]<-hu(x[unstim.ind]);
		return(H)
	}
	return(hess)
}

simMD<-function(alpha.s=c(100,50,10,10),alpha.u=c(100,10,10,10),N=100,w=0.5,n=2){
	pu<-rdirichlet(1,alpha.u)
	pu<-matrix(pu,nrow=N,ncol=length(pu),byrow=T)
	ps<-matrix(pu[1,],nrow=N*(1-w),ncol=length(alpha.s),byrow=T)
	ps<-rbind(ps,rdirichlet(round(N*w),alpha.s))
	NU<-runif(N,1.5*10^n,1.5*10^n)
	NS<-runif(N,1.5*10^n,1.5*10^n)
	nu<-t(sapply(seq_along(1:N),function(i)rmultinom(1,NU[i],pu[i,])))
	ns<-t(sapply(seq_along(1:N),function(i)rmultinom(1,NS[i],ps[i,])))
	data<-list(n.stim=ns,n.unstim=nu)
	return(data)
}



# EM fitting
#' 
#' @param data The observed data
#' @param modelmatrix a model matrix specifying which components should be computed
#' @returnType 
#' @return 
#' @author Greg Finak
#' @export
MDMix<-function(data=NULL,modelmatrix=NULL){
	unstim<-data$n.unstim
	stim<-data$n.stim
	
	#fisher's exact test of all the marginals
	if(ncol(unstim)==4){
		mm<-do.call(cbind,lapply(2:4,function(i)apply(cbind(rowSums(unstim[,-i]),unstim[,i],rowSums(stim[,-i]),stim[,i]),1,function(x)fisher.test(matrix(x,2),alternative="two.sided")$p.value)))<(0.01)/3
		mm<-apply(mm,1,function(x)all(!x))
	}else{ #two-D case
		mm<-sapply(1:nrow(unstim),function(i)fisher.test(rbind(unstim[i,],stim[i,]),alternative="two.sided")$p.value<(0.01))
		mm<-!mm
	}
	#observations with no significant marginals belong to the null component.
	
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
	if(any(is.nan(alpha.s)))
		alpha.s[is.nan(alpha.s)]<-1
	if(any(is.nan(alpha.u)))
		alpha.s[is.nan(alpha.u)]<-1
	
	
	guess<-c(alpha.s,alpha.u)
	w<-colSums(z)/sum(z)
	#EM
	LL<-NULL
	repeat{
		if(is.null(LL)){
			last<-.Machine$double.xmax
		}
		#update parameters
		llnull<-makeLogLikeNULLComponent(stim,unstim)
		llresp<-makeLogLikeRespComponent(stim,unstim)
		gnull<-makeGradientNULLComponent(stim,unstim,z)
		gresp<-makeGradientRespComponent(stim,unstim,z)
		hessnull<-makeHessianNULLComponent(stim,unstim,z)
		hessresp<-makeHessianRespComponent(stim,unstim,z)
		
#		optFun<-function(x,z){
#			ll<- -sum(llnull(x)*z[,1]+llresp(x)*z[,2])
#			#grad<-gnull(x)+gresp(x)
#			#hess<-hessnull(x)+hessresp(x)
#			#list(value=ll,gradient=grad,hessian=hess)
#		}
		
		iter<-2
		ll<-rep(0,1000)
		ll[1]<-.Machine$double.xmax
		lastguess<-guess;
		repeat{
			t<-try(solve(hessresp(guess)+hessnull(guess),gnull(guess)+gresp(guess)),silent=TRUE)
			if(inherits(t,"try-error"))
				t<-ginv(hessresp(guess)+hessnull(guess))%*%(gnull(guess)+gresp(guess)) #uses SVD
			new<-guess-t
#			ll[iter]<- -sum(llnull(new)*z[,1]+llresp(new)*z[,2])
			if((all(abs(new-guess)/abs(guess)<1e-4))|(iter>999)){
				guess<-new
				break
			}
			guess<-new
			iter<-iter+1
		}
#		t<-trust(optFun,parinit=guess,rinit=1,rmax=10^6,minimize=TRUE,z=z)
#		guess<-new<-t$argument

		#compute z's and w's
	
		den<-apply(cbind(log(w[1])+llnull(new), log(w[2])+llresp(new)), 1, function(x)log(sum(exp(x-max(x))))+max(x))
		z2<-exp((llresp(new)+log(w[2]))-(den))
		z<-cbind(1-z2,z2)	
		w<-colSums(z)/sum(z)
		cll<- -sum(llnull(new)*z[,1]+llresp(new)*z[,2])
		if((abs((last-cll)/last)<1e-8)&cll<last){
			break;
		}else if(cll>last){
			new<-lastguess
			break;
		}
		LL<-c(LL,cll)
		cat(cll,"\n")
		last<-cll
	}
	gnull<-makeGradientNULLComponent(stim,unstim,z)
	gresp<-makeGradientRespComponent(stim,unstim,z)
	hessnull<-makeHessianNULLComponent(stim,unstim,z)
	hessresp<-makeHessianRespComponent(stim,unstim,z)

	return(list(llnull=llnull,llresp=llresp,gresp=gresp,hresp=hessresp,gnull=gnull,w=w,hnull=hessnull,z=z,LL=LL,par=new))
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

setMethod(show,"MDlist",function(object){
	cat("Cytokines ",attr(object,"cytokines"),"\n")
	cat("Stimulation ",attr(object,"stim"),"\n")
	cat("Subset ", attr(object,"subset"),"\n")
	cat("Number of obs: ",nrow(object[[1]]),"\n")
})

print.MDlist<-function(x){
			cat("Cytokines ",attr(x,"cytokines"),"\n")
			cat("Stimulation ",attr(x,"stim"),"\n")
			cat("Subset ", attr(x,"subset"),"\n")
			cat(nrow(x[[1]])," observations","\n")
}
