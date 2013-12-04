# contourplot.R Created on: Dec 2, 2011 Author: finak
contourPlot <- function(bmr, which, nsamps = 1000, ecdf.approx = TRUE, 
                        ...) {
  z <- sample(c(1, 2), nsamps + 1000, prob = bmr@w, replace = T)
  inits <- list(alpha0 = bmr@alpha0, beta0 = bmr@beta0, alphaS = bmr@alphaS, 
                betaS = bmr@betaS, w = bmr@w, z = bmr@z)
  s <- gibbsPsPu(curdat = bmr@data, inits = inits, alt = bmr@alternative.model, 
                 which = which, N = nsamps + 1000, burn = 1000, z = z, ecdf.approx = ecdf.approx)
  S <- t(apply(s, 1, function(x) {
    cbind(ns = rbinom(1, bmr@data[which, "ns"] + bmr@data[which, 
                                                          "Ns"], x[2]), nu = rbinom(1, bmr@data[which, "nu"] + bmr@data[which, 
                                                                                                                        "Nu"], x[1]))
  }))
  colnames(S) <- c("ns", "nu")
  m <- MASS::kde2d(S[, 1], S[, 2], 30)
  m$z <- sqrt(m$z)
  contour(m, xlab = "ns", ylab = "nu", ...)
  points(S[, 1], S[, 2], pch = ".", cex = 2)
  points(x = bmr@data[which, "ns"], y = bmr@data[which, "nu"], pch = 20, 
         col = "red")
}


## Old deprecated code TODO FIX THIS!!!
## contourplot<-function(d,k,mciter=50,n=20){
## p<-gibbsPsPu(curdat=d@data,inits=list(alpha0=d@alpha0,alphaS=d@alphaS,beta0=d@beta0,betaS=d@betaS,z=d@z,w=d@w),which=k,alt='greater')
## cm<-colMeans(p)
## lims<-rbind(cm-apply(p,2,sd)*5,cm+apply(p,2,sd)*5) k<-k; #upper
## limit to lower limit.
## x1<-seq(max(lims[1,'ps'],0),min(lims[2,'ps'],1))
## y1<-seq(max(lims[1,'pu'],0),min(lims[2,'pu'],1))
## xup<-d@data[k,'ns']*2 yup<-d@data[k,'nu']*2 #If counts are
## really small, go from zero to 10.  if(xup<10){ xup<-10 }
## if(yup<10){ yup<-10 } x1<-seq(0,xup,l=100) y1<-seq(0,yup,l=100)
## gr<-expand.grid(x1,y1)
## G<-data.frame(gr,d@data[k,'Ns']+d@data[k,'ns']-gr[,1],d@data[k,'Nu']+d@data[k,'nu']-gr[,2])
## colnames(G)<-c('ns','nu','Ns','Nu') #this should be tot adjusted
## by alpha0 and beta0 depending on the generating model, which
## depends on z if(d@z[k,1]>0.5){ #ps=pu model
## tot.x1<-d@data[k,'ns']+d@data[k,'Ns']+d@alpha0+d@beta0+d@data[k,'nu']+d@data[k,'Nu']
## tot.y1<-d@data[k,'nu']+d@data[k,'Nu']+d@alpha0+d@beta0+d@data[k,'ns']+d@data[k,'Ns']
## #integer counts of ns and nu
## x1.i<-unique(round(x1*tot.x1-d@alpha0))
## y1.i<-unique(round(y1*tot.y1-d@alpha0)) #unique indices #sort
## out which model this came from and compute ns and nu
## x1.u<-x1[match(unique(round(x1*tot.x1-d@alpha0)),round(x1*tot.x1-d@alpha0))]
## y1.u<-y1[match(unique(round(y1*tot.y1-d@alpha0)),round(y1*tot.y1-d@alpha0))]
## gr<-expand.grid(x1.i,y1.i) colnames(gr)<-c('ns','nu')
## G<-t(apply(gr,1,function(x)cbind(ns=x[1],nu=x[2],Ns=tot.x1-x[1]-d@alpha0,Nu=tot.y1-x[2]-d@alpha0)))
## }else{ #ps>pu model
## tot.x1<-d@data[k,'ns']+d@data[k,'Ns']+d@alphaS+d@betaS
## tot.y1<-d@data[k,'nu']+d@data[k,'Nu']+d@alpha0+d@beta0 #integer
## counts of ns and nu x1.i<-unique(round(x1*tot.x1-d@alphaS))
## y1.i<-unique(round(y1*tot.y1-d@alpha0)) #unique indices #sort
## out which model this came from and compute ns and nu
## x1.u<-x1[match(unique(round(x1*tot.x1-d@alphaS)),round(x1*tot.x1-d@alphaS))]
## y1.u<-y1[match(unique(round(y1*tot.y1-d@alpha0)),round(y1*tot.y1-d@alpha0))]
## gr<-expand.grid(x1.i,y1.i) colnames(gr)<-c('ns','nu')
## G<-t(apply(gr,1,function(x)cbind(ns=x[1],nu=x[2],Ns=tot.x1-x[1]-d@alphaS,Nu=tot.y1-x[2]-d@alpha0)))
## } colnames(G)<-c('ns','nu','Ns','Nu')
## x1<-x1/(d@data[k,'ns']+d@data[k,'Ns'])
## y1<-y1/(d@data[k,'nu']+d@data[k,'Nu']) #comment out R based code
## #M2<-sqrt(matrix((exp(MarginalNULL(G[,'Ns'],G[,'ns'],G[,'Nu'],G[,'nu'],d@alpha0,d@beta0))*d@w[1]+exp(MarginalGT(G[,'Ns'],G[,'ns'],G[,'Nu'],G[,'nu'],d@alpha0,d@beta0,d@alphaS,d@betaS))*d@w[2]),nrow=length(x1.u)))
## m1<-.Call('MarginalNULL',ns=G[,'ns'],Ns=G[,'Ns'],nu=G[,'nu'],Nu=G[,'Nu'],alpha0=rep(d@alpha0,nrow(G)),beta0=rep(d@beta0,nrow(G)),alphaS=rep(d@alphaS,nrow(G)),betaS=rep(d@betaS,nrow(G)),d@w,log=TRUE,package='flowModels');
## m2<-.Call('MarginalGT',ns=G[,'ns'],Ns=G[,'Ns'],nu=G[,'nu'],Nu=G[,'Nu'],alpha0=rep(d@alpha0,nrow(G)),beta0=rep(d@beta0,nrow(G)),alphaS=rep(d@alphaS,nrow(G)),betaS=rep(d@betaS,nrow(G)),d@w,log=TRUE,mciter=mciter,package='flowModels');
## M<-sqrt(exp(m1)+exp(m2)) M<-matrix(M,nrow=length(x1))
## #par(mfrow=c(1,2))
## contour(x=x1,y=y1,z=M,n=n,xlab=expression(n[s]),ylab=expression(n[u]))
## points(x=d@data[k,'ns'],y=d@data[k,'nu'],pch=20,col='red') } 
