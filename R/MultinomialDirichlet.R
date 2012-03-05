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

#Compute the marginal-log-likelihood for the null distribution (vector of length P)
#' 
#' @param alpha.unstim vector of hyperparameters for the unstimulated sample
#' @param n.stim vector of counts from the stimulated sample
#' @param n.unstim vector of counts from the unstimulated sample
#' @returnType 
#' @return 
#' @author Greg Finak
#' @export
MDnull<-function(alpha.unstim,n.stim,n.unstim){
	lkbeta(alpha.unstim+n.stim+n.unstim)-lkbeta(alpha.unstim)+lfctrl(n.stim)+lfctrl(n.unstim)
}

#Compute the marginal log-likelihood for the alternative distribution (vector of length P)
#' 
#' @param alpha.unstim vector of hyperparameters for the unstimulated sample
#' @param alpha.stim vector of hyperparameters for the stimulated sample
#' @param n.stim vector of counts from the stimulated sample
#' @param n.unstim vector of counts from the unstimulated sample
#' @returnType 
#' @return 
#' @author Greg Finak
#' @export
MDalternative<-function(alpha.unstim,alpha.stim,n.stim,n.unstim){
	lkbeta(alpha.unstim+n.unstim)-lkbeta(alpha.unstim)-lkbeta(alpha.stim)+lkbeta(alpha.stim+n.stim)+lfctrl(n.stim)+lfctrl(n.unstim)
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
	alpha.u<-colSums(data$n.unstim)/sum(data$n.unstim)
#	1.	Nelson, W. JSTOR: The American Statistician, Vol. 26, No. 3 (Jun., 1972), pp. 22-27. The American statistician (1972).
	p<-(data$n.unstim/rowSums(data$n.unstim))
	q<-(data$n.stim/rowSums(data$n.stim))
	ddirichlet(q,data$n.stim)/(ddirichlet(q,data$n.stim)+ddirichlet(q,data$n.unstim))
	rho<-p/q
	rowSums(data$n.stim)
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
