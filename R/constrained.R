# TODO: Constrained, reparameterized dirichlet multinomial model
# 
# Author: finak
###############################################################################


#The constrained log likelihood for the non responder component
#pars is the vector of parameters (n, lambda, w1, w2, ...)
#'@importFrom pracma hessian grad
constrainedNRLL<-function(data.stim,data.unstim,z=NULL){
	data<-t(data.stim+data.unstim)
	if(is.null(z)){
		z<-matrix(1,nrow=ncol(data),ncol=2)
	}
	function(pars){
		#TODO make general for multiple components (i.e. more than two columns in z).
		s<-pars[1]
		lambda<-pars[2]
		w<-pars[-c(1:2)]
		sw<-s*w
		sws<-sum(sw)
		data<-sw+data
		sum((sum(lgamma(sw))-rowSums(lgamma(t(data)))+lgamma(colSums(data))-lgamma(sws)+lambda*(sum(w)-1))*z[,1])
	}
}




#Gradient for constrained log likelihood for the non responder component, parameterized as n*s

gradCostrainedNRLL<-function(data.stim,data.unstim,z=NULL){
	data<-t(data.stim+data.unstim)
	if(is.null(z)){
		z<-matrix(1,nrow=ncol(data),ncol=2)
	}
	function(pars){
		#TODO make general for multiple components (i.e. more than two columns in z).
		#calculate things once
		n<-pars[1]
		lambda<-pars[2]
		s<-pars[-c(1:2)]
		nc<-ncol(data)
		ns<-n*s
		dns<-digamma(ns)
		ss<-sum(s)
		sns<-sum(ns)
		dsns<-digamma(sns)
		data<-ns+data
		csd<-colSums(data)
		digcsd<-digamma(csd)
		digdata<-digamma(data)
		
		gn<-colSums(s*digamma(ns)-s*digamma(data))+ss*digamma(csd)-ss*digamma(sns)
		gs<-((n*dns-n*digdata)+n*digcsd-n*dsns+lambda)
		gl<-rep(ss-1,nc)
		gr<-rbind(gn,gl,gs)
		return(gr%*%z[,1])
	}	
}

#Hessian of the constrained log likelihood for the non responder component
hessConstrainedNRLL<-function(data.stim,data.unstim,z=NULL){
	data<-t(data.stim+data.unstim)
	if(is.null(z)){
		z<-matrix(1,nrow=ncol(data),ncol=2)
	}
	function(pars){
		#TODO make general for multiple components (i.e. more than two columns in z).
		#Calculate everything once
		n<-pars[1]
		lambda<-pars[2]
		s<-pars[-c(1:2)]
		ns<-n*s
		s2<-s*s
		sns<-sum(ns)
		n2<-n^2
		ss<-sum(s)
		data<-ns+data
		csd<-colSums(data)
		td<-trigamma(data)
		tdata<-trigamma(data)
		tns<-trigamma(ns)
		tcsd<-trigamma(csd)
		tsns<-trigamma(sns)
		ddata<-digamma(data)
		dsns<-digamma(sns)
		dcsd<-digamma(csd)
		dns<-digamma(ns)
		
		#Partial deriviatives
		d2n<- (colSums(s2*tns-s2*tdata)+ss^2*tcsd-ss^2*tsns)%*%z[,1]
		d2l<-0
		d2s<-n2*((tns-t(t(td)-tcsd+tsns))%*%z[,1])
		D<-diag(c(d2n,d2l,d2s))
		dndl<-0
		dnds<-t(t(dns+ns*tns - ddata - ns*tdata) + dcsd +n*ss*tcsd - dsns-n*ss*tsns)%*%z[,1]
		dsds<-(n2*tcsd-n2*tsns)%*%z[,1]
		dsdl<-rep(1,ncol(data))%*%z[,1]
		lt<-c(dndl,dnds,rep(dsdl,length(s)),rep(dsds,sum(1:(length(s)-1))))
		
		#Construct hessian
		D[lower.tri(D)]<-lt
		D<-t(D)
		D[lower.tri(D)]<-lt
		return(D)
	}
}

#log likelihood for constrained model responder component
constrainedRLL<-function(data.stim,data.unstim,z=NULL){
	dats<-t(data.stim)
	datu<-t(data.unstim)
	if(is.null(z)){
		z<-matrix(1,nrow=ncol(datu),ncol=2)
	}
	function(pars){
		s1<-pars[1]
		s2<-pars[2]
		lambda1<-pars[3]
		lambda2<-pars[4]
		w<-pars[-c(1:4)]
		l<-length(w)
		ws<-w[(1:(l/2))] #stim parameters
		wu<-w[(((l/2)+1):l)] #unstim parameters
		dats<-dats+ws#stim hyperparameters
		datu<-datu+wu   #unstim hyperparameters
		
		swu<-s1*ws
		sws<-s2*wu
		
		gsws<-lgamma(sws)
		gswu<-lgamma(swu)
		
		gdats<-lgamma(dats)
		gdatu<-lgamma(datu)
		
		csdats<-colSums(dats)
		csdatu<-colSums(datu)
		
		gcsdats<-lgamma(csdats)
		gcsdatu<-lgamma(csdatu)
		
		ssws<-sum(sws)
		sswu<-sum(swu)
		
		gssws<-lgamma(ssws)
		gsswu<-lgamma(sswu)
		
		ll<- (colSums(gsws+gswu-gdats-gdatu)+gcsdats+gcsdatu-gssws-gsswu+lambda1*(sum(wu)-1)+lambda2*(sum(ws)-1))%*%z[,2]
		return(ll)
	}
}

#gradient for the log likelihood of the constrained responder component
gradConstrainedRLL<-function(data.stim,data.unstim,z=NULL){
	
}

#hessian for the log likelihood of the constrained responder component
hessConstrainedRLL<-function(data.stim,data.unstim,z=NULL){
	
}

#Construct the gradient, hessian and likelihood functions
#numerical gradient and hessian
test<-function(x){
    set.seed(5)
	gr<-gradCostrainedNRLL(dat[[1]],dat[[2]])
	ll<-constrainedNRLL(dat[[1]],dat[[2]])
	hess<-hessConstrainedNRLL(dat[[1]],dat[[2]])
	#Closed form. Should be equal to above, but a bit faster.
	time2<-system.time({
				dat<-MIMOSA:::simMD(alpha.s=c(100,10,10,10),alpha.u=c(100,10,10,10),w=0)
				repeat{
					new<-guess-solve(hess(guess))%*%gr(guess)
					if(norm(new-guess)<0.01){
						guess.closedform<-guess
						names(guess.closedform)<-c("N","lambda","w1","w2","w3","w4")
						break
					}
					else
						guess<-new
				}
			})
	
time1<-system.time({
			guess<-c(mean(rowSums(dat[[1]])),1,prop.table(colSums(dat[[1]])))
			repeat{
				new<-guess-solve(hessian(ll,guess))%*%grad(ll,guess)
				if(norm(new-guess)<0.01){
					guess.numerical<-guess
					names(guess.numerical)<-c("N","lambda","w1","w2","w3","w4")
					break
					
				}
				else
					guess<-new
			}
		})

cat("Speed increase: ",time1[3]/time2[3], " fold.\n")
cat("Estimates from numerical derivatives:\n \t",guess.numerical[-2L],"\n")
cat("Estimates from closed form derivatives:\n\t ",guess.closedform[-2L],"\n")
optres<-try(optim(fn=ll,par=c(mean(rowSums(dat[[1]])),1,prop.table(colSums(dat[[1]]))),method="L-BFGS-B",gr=gr,lower=c(0,-Inf,1e-10,1e-10,1e-10,1e-10),upper=c(Inf,Inf,1,1,1,1)),silent=TRUE)
if(inherits(optres,"try-error")){
	message("optim failed")
	return(list(optim=optres,Newton.Estimates=guess.numerical))
}
else return(list(optim=optres,Newton.Estimates=guess.numerical))
}


#responder component test
#llr<-constrainedRLL(dat[[1]],dat[[2]])
#guess<-c(150,150,1,1,prop.table(colMeans(dat[[1]])),prop.table(colMeans(dat[[2]]))) 
#llr(guess)
#guess<-guess-0.1*solve(hessian(LL,guess))%*%grad(LL,guess)

