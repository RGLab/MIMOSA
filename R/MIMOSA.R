#'MIMOSA: Mixture Models for Single Cell Assays
#'
#'MIMOSA implements mxitures of Dirichlet-multinomial or Beta-binomial models 
#'for paired count data from single--cell assays that typically arise 
#'in immunological studies. It can be used for ICS (Intracellular Cytokine
#'Staining) assays to detect vaccine responders, for example, or to detect
#'changes in proportions of cells expressing a gene, such as in Fluidigm Biomark
#'Single--cell gene expression.
#'
#'@docType package
#'@useDynLib MIMOSA
#'@rdname MIMOSA-package
#'@import Formula Biobase
#'@importFrom MASS ginv
#'@importClassesFrom methods array character data.frame factor integer matrix numeric
#'@name MIMOSA-package
#'@references Finak, Greg; McDavid, Andrew; Chattopadhyay, Pratip; Dominguez, Maria; De Rosa, Steve; Roederer, Mario; Gottardo, Raphael
#'  Mixture Models for Single Cell Assays with Applications to Vaccine Studies
#'  eprint arXiv:1208.5809 \url{http://arxiv.org/abs/1208.5809}
NULL


#' Stimulated and unstimulated T-cell counts for an ICS assay
#' 
#' A data set containing T-cell counts for various stimulations and cytokines in an ICS assay.
#' 
#' \itemize{
#' \item pos. The positive cell counts
#' \item neg. The negative cell counts
#' \item fname. The feature name (cytokine) measured
#' \item parent. The parent T-cell population
#' \item antigen. The antigen stimulation for this sample
#' \item ID. The subject ID
#' }
#'
#'@docType data
#'@keywords datasets
#'@format A data frame with 3960 rows
#'@name ICS
NULL

#'Fit a MIMOSA Model
#'
#'This method fits a MIMOSA model to count data stored in an ExpressionSet
#'object.
#'
#'@details The ExpressionSet should be fully annotated with featureData and
#'  phenoData. For ICS data, for example, features would be positive and
#'  negative counts for different cytokine producing cell subsets (i.e.
#'  IFNg_pos, IFNg_neg) The formula lhs should contain features and the rhs
#'  should contain phenotypic variable. See the vignette for an example.
#'  
#'@param formula describing the features on the lhs and the phenodata on the
#'  rhs, supporting extended formula interface with conditioning.
#'@param data an \code{ExpressionSet} object with features on rows and samples
#'  (labelled with phenoData) on columns.
#'@return an object of type \code{MIMOSAResult}
#'@aliases MIMOSA,formula,ExpressionSet-method
#'@export
setGeneric("MIMOSA",def=function(formula,data,...){
  standardGeneric("MIMOSA")
})


setMethod("MIMOSA",c("formula","ExpressionSet"),definition=function(formula,data,ref,method="mcmc",subset,...){
  if(!exists("ref")){
    stop("ref must contain expressions that define subsets to be compared")
  }
  method<-match.arg(method,c("mcmc","EM"))
  if(!inherits(formula,"Formula")){
    formula<-Formula(formula)
  }   
  
  mf.ref<-model.frame(formula,data)[eval(substitute(ref),pData(data)),,drop=FALSE]
  mf.test<-model.frame(formula,data)[eval(substitute(subset),pData(data)),,drop=FALSE]
  
  if(length(formula)[2]>1){
    spl.ref<-model.part(formula,mf.ref,rhs=2)
    spl.test<-model.part(formula,mf.test,rhs=2)
  }
  interact<-function(x,...){
    if(length(x)==1){
      return(factor(x[[1]]))
    }
    if(length(x)==2){
      interaction(x[[1]],x[[2]],drop=TRUE,...)
    }else{
      interaction(x[[1]],Recall(x[[-1L]],drop=TRUE),...)
    }
  }
  if(length(formula)[2]>1){
    ref<-split(model.part(formula,mf.ref,lhs=1),interact(spl.ref))
    test<-split(model.part(formula,mf.test,lhs=1),interact(spl.test))
    pd<-split(model.part(formula,mf.test,rhs=1:2),interact(spl.test))
    result<-vector("list",length(test))
  }else{
    ref<-list(model.part(formula,mf.ref,lhs=1))
    test<-list(model.part(formula,mf.test,lhs=1))
    pd<-list(model.part(formula,mf.test,rhs=1))
    result<-vector("list",1)
  }
  for(i in 1:length(test)){
    #recycle the reference
    j<-data.frame(1:length(test),1:length(ref))[i,2]
    fitme<-list(n.stim=test[[i]],n.unstim=ref[[j]])
    if(method%in%"mcmc"){
      res<-.fitMCMC(fitme,inits=MDMix(fitme,initonly=TRUE),...)
      res$params<-res$params<-apply(res$getmcmc(),2,function(x)quantile(x,c(0.025,0.5,0.975),na.rm=TRUE))
      if(ncol(fitme[[1]])==2){
        res$p<-lapply(res$getP(),function(x)do.call(rbind,lapply(x,function(x)quantile(x,c(0.025,0.5,0.975)))))
      }
      attr(res,"pData")<-new("AnnotatedDataFrame",pd[[i]])
      res<-MCMCResult(object=res)
    }
    else if (method%in%"EM"){
      res<-MDMix(fitme)
      res@pd<-new("AnnotatedDataFrame",pd[[i]])
    }
    res<-new("MIMOSAResult",result=res)
    result[[i]]<-res
  }
  return(result)
})
