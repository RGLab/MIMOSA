#volcanoPlot.R
#
# Created on: Nov 29, 2011
#     Author: finak

setGeneric("volcanoPlot",function(bmr,...){
			standardGeneric("volcanoPlot");
})

setMethod("volcanoPlot",signature("BetaMixResult"),function(bmr,threshold=0.01,N=500,...){
  ml<-estimateProportions(bmr,method="ML")
  post<-estimateProportions(bmr,method="posterior",volcano=TRUE,N=N,...)
  bmr.nlfdr<--log10(bmr@fdr)
  
  d<-bmr@data[,c("ns","nu","Ns","Nu")]
  if(bmr@alternative.model=="greater"){
 	 fisher.p<-apply(d,1,function(x){fisher.test(as.table(matrix(x,ncol=2)),alternative="greater")$p})
 }else{
	 fisher.p<-apply(d,1,function(x){fisher.test(as.table(matrix(x,ncol=2)),alternative="two.sided")$p})
 }
 
  fisher.nlfdr<--log10(p.adjust(fisher.p,"fdr"))
  
  yup<-max(c(fisher.nlfdr,bmr.nlfdr))
  if(is.infinite(yup)){
	yup<-.Machine$double.max.exp
  }
  fisher.nlfdr[is.infinite(fisher.nlfdr)]<-.Machine$double.max.exp
  bmr.nlfdr[is.infinite(bmr.nlfdr)]<-.Machine$double.max.exp
  
  ylow<-0
  
  xup<-max(c(apply(ml,1,diff),apply(post,1,diff)))
  xlow<-min(c(0,apply(ml,1,diff),apply(post,1,diff)))
  
  bmr.colors<-(bmr.nlfdr > -log10(threshold))+1
  
  fisher.colors<-(fisher.nlfdr > -log10(threshold))+1
  fisher.colors[fisher.colors==2]<-3
  
  print(plot(apply(ml,1,diff),fisher.nlfdr,xlim=c(xlow,xup),ylim=c(ylow,yup),pch=2,col=fisher.colors,main=strwrap(paste("Difference Between ",bmr@stimulation," and ",bmr@control," Proportions",paste(bmr@cytokine,collapse="/")),...),xlab="Ps-Pu",ylab="-log10(q-value)"))
  points(post[,"ps"]-post[,"pu"],bmr.nlfdr,pch=3,col=bmr.colors)
  
  abline(v=0,lty=2)
  legend(x="bottomright",pch=c(2,3),c("ML Estimate","Posterior Estimate"))
  invisible(cbind(betabin=length(which(bmr@fdr<threshold)),fisher=length(which(p.adjust(fisher.p,"fdr")<threshold))))
})

#TODO fix the volcano plot when fdr=0. -log10(0) = Inf and is not plotted.
#TODO also adjust the visible range so that all the data is visible, not just the ML or posterior estimates