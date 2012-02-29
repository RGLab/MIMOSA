#densityPlot.R
#
# Created on: Dec 6, 2011
#     Author: finak

#densityplot of the posterior distributions of ps and pu
#function(curdat=NULL,inits=NULL,which=NULL,N=10000,burn=5000,alt=NULL){

densityplot.BetaMixResult<-function(x,data,full=F,N=5000,...){
	if(class(x)!="BetaMixResult"){
		stop("Class of x must be BetaMixResult")
	}
	if(!is.numeric(data)){
		stop("data must be an integer vector specifying which samples to plot (by index, 1..n")
	}
	burn<-50
	inits<-list(z=x@z,w=x@w,alphaS=x@alphaS,alpha0=x@alpha0,betaS=x@betaS,beta0=x@beta0)
	if(full){
		z<-sample(c(1,2),N+burn,x@w,replace=T)
	}else{
		z<-NULL
	}
	p<-lapply(as.list(data),function(y)gibbsPsPu(curdat=x@data,alt=x@alternative.model,inits=inits,which=y,burn=burn,z=z,N=N,...))
	p<-do.call(rbind,p)
	p<-data.frame(p,subject=gl(length(data),N-burn,labels=rownames(x@data)[data]))
	p<-melt(p,measured=c("pu","ps"),id="subject")
	p<-rename(p,list(subject="subject",variable="estimate",value="p"))
	densityplot(~p|subject,p,auto.key=TRUE,groups=estimate,pch='.',main=expression(paste("Posterior Distributions of    ", p[s] ,"  and  ", p[u] ,"  | subject")),n=1000,lwd=2,scales=list(relation="free"),...)
}
