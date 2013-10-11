load("~/Downloads/testMIMOSA.rda")
E<-ConstructMIMOSAExpressionSet(data,measure.columns=c("NSUB","CYTNUM"),.variables=.(SIMPLEID),reference=NULL)
result<-MIMOSA(NSUB+CYTNUM~RefTreat+SIMPLEID,E,ref=RefTreat%in%"Reference",subset=RefTreat%in%"Treatment")
df1<-data.frame(pData(result2[[1]]),z=result[[1]]@z[,2],dp=prop.table(as.matrix(result2[[1]]@result@data$n.stim),1)[,2]-prop.table(as.matrix(result2[[1]]@result@data$n.unstim),1)[,2])

result2<-MIMOSA(NSUB+CYTNUM~RefTreat+SIMPLEID,E,ref=RefTreat%in%"Reference",subset=RefTreat%in%"Treatment",method="EM")
data.frame(pData(result2[[1]]),result2[[1]]@z[,2],dp=prop.table(as.matrix(result2[[1]]@result@data$n.stim),1)[,2]-prop.table(as.matrix(result2[[1]]@result@data$n.unstim),1)[,2])

result3<-MIMOSA(NSUB+CYTNUM~RefTreat+SIMPLEID,E,ref=RefTreat%in%"Reference",subset=RefTreat%in%"Treatment",run.parallel=TRUE)
df1<-data.frame(pData(result2[[1]]),result3[[1]]@z[,2],dp=prop.table(as.matrix(result2[[1]]@result@data$n.stim),1)[,2]-prop.table(as.matrix(result2[[1]]@result@data$n.unstim),1)[,2])

data.frame(pData(result2[[1]]),em=result[[1]]@z[,2],mcmc=result2[[1]]@z[,2],dpr=result[[1]]@z[,2]-result2[[1]]@z[,2],dp=prop.table(as.matrix(result2[[1]]@result@data$n.stim),1)[,2]-prop.table(as.matrix(result2[[1]]@result@data$n.unstim),1)[,2])

mpc<-read.csv("~/Downloads/MacvsPC.csv")
mpc$NSUB<-as.numeric(as.character(mpc$NSUB))
mpc$NSUB_REF<-as.numeric(as.character(mpc$NSUB_REF))
mpc$CYTNUM<-as.numeric(as.character(mpc$CYTNUM))
mpc$CYTNUM_REF<-as.numeric(as.character(mpc$CYTNUM_REF))
raw<-with(subset(mpc,computer=="PC"),data.frame(pr=prop.table(cbind(NSUB,CYTNUM),1)[,2]-prop.table(cbind(NSUB_REF,CYTNUM_REF),1)[,2],effect=as.numeric(as.character(effect)),SIMPLEID=SIMPLEID))

library(data.table)
raw<-data.table(raw)
df1<-data.table(df1)
setkey(df1,SIMPLEID)
setkey(raw,SIMPLEID)
merge(raw,df1)
