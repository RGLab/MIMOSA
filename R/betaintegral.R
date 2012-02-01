#deprecated the numerical integration
betaintegral <-
function(alphaS,betaS,alpha0,beta0,Nu,nu,Ns,ns,p=1000){
 l<-length(Nu)
 #set the upper limit of integration if the density is concentrated aroud 0.. (no need to go all the way to 1)
   ms<-(alphaS+ns)/(alphaS+betaS+ns+Ns)
   m0<-(alpha0+nu)/(alpha0+beta0+nu+Nu)
   vs<-sqrt(((alphaS+ns)*(betaS+Ns))/((alphaS+ns+betaS+Ns)^2*(alphaS+ns+betaS+ns+1)))
   v0<-sqrt(((alpha0+nu)*(beta0+Nu))/((alpha0+nu+beta0+Nu)^2*(alpha0+nu+beta0+nu+1)))
   mx<-max(ms+7*vs,m0+7*v0)
   pu<-seq(0,mx,l=p)
   dpu<-diff(pu)[1];
  res<-double(l)
 res<-.C("betaintegral_c",as.double(alphaS),as.double(betaS),as.double(alpha0),as.double(beta0),as.integer(Nu),as.integer(nu),as.integer(Ns),as.integer(ns),as.double(pu),as.double(dpu),as.double(res),as.integer(l),as.integer(p));
  return(res[[11]]);
}


#New code does Monte-Carlo integration
betaintegral_R<-function(alphaS,betaS,alpha0,beta0,Nu,nu,Ns,ns,l=1000){
	res<-lbeta(ns+alphaS, betaS+Ns)+lbeta(nu+alpha0,beta0+Nu)+sapply(1:length(Ns),function(i){
  	log(mean(pbeta(rbeta(l,alpha0+nu[i],beta0+Nu[i]),alphaS+ns[i],betaS+Ns[i],lower.tail=FALSE,log=F)));
	})
	return(res)
}

integrand<-function(alpha0,beta0,alphaS,betaS,d){
	function(s){dbeta(s,alpha0+d[,"nu"],beta0+d[,"Nu"])*pbeta(s,alphaS+d[,"ns"],betaS+d[,"Ns"])}
}

#explicit
betaintegral_R2<-function(alphaS,betaS,alpha0,beta0,Nu,nu,Ns,ns,pu,dpu){
	res=lbeta(ns+alphaS, betaS+Ns)+lbeta(nu+alpha0,beta0+Nu);
	for(i in 1:length(ns)){
		acc=0;
		acc=sum(dbeta(pu,alpha0+ nu[i],beta0+ Nu[i])*(1-pbeta(pu,alphaS+ns[i],betaS+ Ns[i]))*dpu)
		res[i]=res[i]+log(acc)
	}
	res
}