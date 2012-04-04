#Classes.R
#
# Created on: Nov 30, 2011
#     Author: finak

#classes
#ICS class holds the data
setClass("ICS",representation=list(pos="integer",neg="integer",fname="factor",fcsfile="character",parent="factor",antigen="factor",ID="factor",rest="data.frame",.Data="array"),
		validity=function(object){
			test<-all.equal(length(object@pos),length(object@neg),length(object@fname),length(object@fcsfile),length(object@parent),length(object@antigen),length(object@ID))
			if(test){
				return(TRUE)
			}else{
				return("All object slots are not the same length")
			}
		}
)

#Holds the output of BetaMix
setClass("BetaMixResult",representation=list(alternative.model="character",cytokine="character",control="character",stimulation="character",ll="numeric",traj="numeric",iter="numeric",z="matrix",w="numeric",alpha0="numeric",beta0="numeric",alphaS="numeric",betaS="numeric",fdr="numeric",data="data.frame"))

setClass("MDMixResult",representation=list(llnull="function",llresp="function",gresp="function",w="numeric",z="matrix",hresp="function",gnull="function",ll="numeric",hnull="function",par.unstim="numeric",par.stim="numeric",data="list"))

setClassUnion("MixResult",c("BetaMixResult","MDMixResult"))
#Constructors
ICS<-function(x=data.frame(pos=NA_integer_,neg=NA_integer_,fname=NA_character_,fcsfile=NA_character_,parent=NA_character_,antigen=NA_character_,ID=NA_character_),extra.identifiers=NULL){
	rest<-setdiff(colnames(x),c("pos","neg","fname","fcsfile","parent","antigen","ID"))
	value<-new("ICS",pos=x$pos,neg=x$neg,fname=x$fname,fcsfile=x$fcsfile,parent=x$parent,antigen=x$antigen,ID=x$ID,rest=x[,rest],.Data=list())
	value@.Data<-prepCounts(value,extra.identifiers=extra.identifiers)
	value
}

BetaMixResult<-function(alternative.model=NA_character_,cytokine=NA_character_,control=NA_character_,stimulation=NA_character_,ll,traj,iter,z,w,alpha0,beta0,alphaS,betaS,data){
	new("BetaMixResult",alternative.model=alternative.model,control=control,stimulation=stimulation,cytokine=cytokine,ll=ll,traj=traj,iter=iter,z=z,w=w,alpha0=alpha0,beta0=beta0,alphaS=alphaS,betaS=betaS,fdr=fdr(z),data=data)	
}

