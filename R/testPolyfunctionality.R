#testPolyfunctionality.R
#
# Created on: Dec 6, 2011
#     Author: finak

testPolyfunctionality<-function(ics=NULL,cytokineA=NULL,cytokineB=NULL,or=NULL,stim=NULL,control=NULL,subset=NULL,shrink=1,scl=1){
	if(any(c(is.null(ics),is.null(cytokineA),is.null(stim),is.null(cytokineB),is.null(or),is.null(control),is.null(subset)))){
		stop("All arguments must be non--null");
	}
	A<-extractData(ics,control=control,stim=stim,subset=c(subset,cytokineA))
	B<-extractData(ics,control=control,stim=stim,subset=c(subset,cytokineB))
	OR<-extractData(ics,control=control,stim=stim,subset=c(subset,or))
	AND<-A+B-OR
	AND[,"Ns"]<-A[,"Ns"]+A[,"ns"]-AND[,"ns"]
	AND[,"Nu"]<-A[,"Nu"]+A[,"nu"]-AND[,"nu"]
	attr(AND,"stimulation")<-attr(A,"stimulation")
	attr(AND,"control")<-attr(A,"control")
	attr(AND,"cytokine")<-union(attr(A,"cytokine"),attr(B,"cytokine"))
	#resA<-BetaMix(A,shrink=shrink,scl=9,K=5000)
	#resB<-BetaMix(B,shrink=shrink,scl=9,K=5000)
	ns<-cbind(A[,1]-AND[,1],B[,1]-AND[,1],AND[,1])
	nu<-cbind(A[,2]-AND[,2],B[,2]-AND[,2],AND[,2])
	ns<-cbind(rowSums(A[,c(1,3)])-OR[,1],ns)
	nu<-cbind(rowSums(A[,c(2,4)])-OR[,2],nu)
	colnames(ns)<-c("--","+-","-+","++")
	colnames(nu)<-c("--","+-","-+","++")
	dat<-list(n.stim=ns,n.unstim=nu)
	poly<-MIMOSA:::.fitMCMC(dat,inits=MIMOSA:::MDMix(dat,initonly=TRUE),100000,50000,thin=2)
	return(resAND)
}