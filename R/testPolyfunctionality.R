#testPolyfunctionality.R
#
# Created on: Dec 6, 2011
#     Author: finak

testPolyfunctionality<-function(ics=NULL,cytokineA=NULL,cytokineB=NULL,or=NULL,stim=NULL,control=NULL,subset=NULL,shrink=1,scl=1){
	if(any(c(is.null(ics),is.null(cytokineA),is.null(stim),is.null(cytokineB),is.null(or),is.null(control),is.null(subset)))){
		stop("All arguments must be non--null");
	}
	A<-flowContrasts(ics,control=control,stim=stim,subset=c(subset,cytokineA))
	B<-flowContrasts(ics,control=control,stim=stim,subset=c(subset,cytokineB))
	OR<-flowContrasts(ics,control=control,stim=stim,subset=c(subset,or))
	AND<-A+B-OR
	AND[,"Ns"]<-A[,"Ns"]+A[,"ns"]-AND[,"ns"]
	AND[,"Nu"]<-A[,"Nu"]+A[,"nu"]-AND[,"nu"]
	attr(AND,"stimulation")<-attr(A,"stimulation")
	attr(AND,"control")<-attr(A,"control")
	attr(AND,"cytokine")<-union(attr(A,"cytokine"),attr(B,"cytokine"))
	#resA<-BetaMix(A,shrink=shrink,scl=9,K=5000)
	#resB<-BetaMix(B,shrink=shrink,scl=9,K=5000)
	resAND<-BetaMix(AND,shrink=shrink,scl=scl,K=5000)
	return(resAND)
}