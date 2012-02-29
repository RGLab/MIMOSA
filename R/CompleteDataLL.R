
CompleteDataLLRcpp<-function(d,alpha0,beta0,alphaS,betaS,z,w,alternative="greater",mciter=50){
	if(length(alpha0)==1){
		alpha0<-rep(alpha0,nrow(d));
	}
	if(length(beta0)==1){
		beta0<-rep(beta0,nrow(d));
	}
	if(length(alphaS)==1){
		alphaS<-rep(alphaS,nrow(d));
	}
	if(length(betaS)==1){
		betaS<-rep(betaS,nrow(d));
	}
	result<-.Call("CompleteDataLLRcpp",d[,"ns"],d[,"Ns"],d[,"nu"],d[,"Nu"],alpha0,beta0,alphaS,betaS,z,w,alternative,mciter=mciter,package="MIMOSA");
	return(result);
}
CompleteDataLL<-function(d,alpha0,beta0,alphaS,betaS,z,w,alternative="greater",mciter=50){
	if(alternative=="greater"){
		result<-MIMOSA:::MarginalGT(d[,"Ns"],d[,"ns"],d[,"Nu"],d[,"nu"],alpha0,beta0,alphaS,betaS)*z[,2]
	}else{
		result<-MIMOSA:::MarginalNE(d[,"Ns"],d[,"ns"],d[,"Nu"],d[,"nu"],alpha0,beta0,alphaS,betaS)*z[,2]
	}
	result<-result+MIMOSA:::MarginalNULL(d[,"Ns"],d[,"ns"],d[,"Nu"],d[,"nu"],alpha0,beta0)*z[,1]
	
	return(sum(result));
}
