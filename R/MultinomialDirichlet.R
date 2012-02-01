# TODO: Implements the multinomial-dirichelt model for two cytokines
# 
# Author: finak
###############################################################################



#Observed data
#n.stim - vector of counts from the stimulated sample
#n.unstim - vector of counts from the unstimulated sample
#alpha.unstim
#alpha.stim

#compute the log of the factorial term
lfctrl<-function(x){
	lfactorial(sum(x))-sum(lfactorial(x))
}

#compute the log of the k-dimensional beta function
lkbeta<-function(alpha){
	sum(lgamma(alpha))-lgamma(sum(alpha))
}

#Compute the marginal-log-likelihood for the null distribution (vector of length P)
MDnull<-function(alpha.unstim,n.stim,n.unstim){
	lkbeta(alpha.unstim+n.stim+n.unstim)-lkbeta(alpha.unstim)+lfctrl(n.stim)+lfctrl(n.unstim)
}

#Compute the marginal log-likelihood for the alternative distribution (vector of length P)
MDalternative<-function(alpha.unstim,alpha.stim,n.stim,n.unstim){
	lkbeta(alpha.unstim+n.unstim)-lkbeta(alpha.unstim)-lkbeta(alpha.stim)+lkbeta(alpha.stim+n.stim)+lfcrl(n.stim)+lfctrl(n.unstim)
}


#EM algorithm for the multinomial dirichlet mixture
MDMix<-function(data=NULL,modelmatrix=NULL){
	inits<-initMDMix(data=data,modelmatrix=modelmatrix)
	
}

#initialization of z's and parameters for the multinomial dirichlet mixture
initMDMix<-function(data=NULL,modelmatrix=NULL){
	
}