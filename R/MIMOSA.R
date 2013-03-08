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


setMethod("MIMOSA",c("formula","ExpressionSet"),definition=function(formula,data,ref,method="mcmc",subset,getP=FALSE,p.thin=1,run.parallel=FALSE,...){
  if(!exists("ref")){
    stop("ref must contain expressions that define subsets to be compared")
  }
  method<-match.arg(method,c("mcmc","EM"))
  if(!inherits(formula,"Formula")){
    formula<-Formula(formula)
  }   
  
  mf.ref<-model.frame(formula,data,na.action=NULL)
  mf.test<-model.frame(formula,data,na.action=NULL)
  
  if(!missing(ref)){
    mf.ref<-mf.ref[eval(substitute(ref),pData(data)),,drop=FALSE]
  }
  if(!missing(subset)){
    mf.test<-mf.test[eval(substitute(subset),pData(data)),,drop=FALSE]
  }
  
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
  if(!run.parallel){
    for(i in 1:length(test)){
      #recycle the reference
      j<-data.frame(1:length(test),1:length(ref))[i,2]
      fitme<-list(n.stim=test[[i]],n.unstim=ref[[j]])
      #remove NAs
      nas<-unique(rbind(which(is.na(fitme[[1]]),T),which(is.na(fitme[[2]]),T))[,1])
      if(length(nas)>0){
        fitme[[1]]<-fitme[[1]][-nas,]
        fitme[[2]]<-fitme[[2]][-nas,]
        pd[[i]]<-pd[[i]][-nas,]
      }
      if(method%in%"mcmc"){
        res<-.fitMCMC(fitme,inits=MDMix(fitme,initonly=TRUE),...)
        res$params<-res$params<-apply(res$getmcmc(),2,function(x)quantile(x,c(0.025,0.5,0.975),na.rm=TRUE))
        if(ncol(fitme[[1]])==2&getP){
          res$p<-lapply(res$getP(thin=p.thin),function(x)do.call(rbind,lapply(x,function(x)quantile(x,c(0.025,0.5,0.975),na.rm=TRUE))))
        }else{
          res$p<-list()
        }
        attr(res,"pData")<-new("AnnotatedDataFrame",pd[[i]])
        outfile<-get("outfile",environment(res$getmcmc))
        unlink(outfile)
        unlink(paste(outfile,"P",sep="")) 
        res<-MCMCResult(object=res)
      }
      else if (method%in%"EM"){
        res<-MDMix(fitme)
        res@pd<-new("AnnotatedDataFrame",pd[[i]])
      }
      res<-new("MIMOSAResult",result=res)
      result[[i]]<-res
    }
  }else if(run.parallel&any(grepl("multicore",loadedNamespaces()))){
    result<-mclapply(1:length(test),function(i){
      #recycle the reference
      j<-data.frame(1:length(test),1:length(ref))[i,2]
      fitme<-list(n.stim=test[[i]],n.unstim=ref[[j]])
      #remove NAs
      nas<-unique(rbind(which(is.na(fitme[[1]]),T),which(is.na(fitme[[2]]),T))[,1])
      if(length(nas)>0){
        fitme[[1]]<-fitme[[1]][-nas,]
        fitme[[2]]<-fitme[[2]][-nas,]
        pd[[i]]<-pd[[i]][-nas,]
      }
      if(method%in%"mcmc"){
        res<-.fitMCMC(fitme,inits=MDMix(fitme,initonly=TRUE),...)
        res$params<-res$params<-apply(res$getmcmc(),2,function(x)quantile(x,c(0.025,0.5,0.975),na.rm=TRUE))
        if(ncol(fitme[[1]])==2&getP){
          res$p<-lapply(res$getP(thin=p.thin),function(x)do.call(rbind,lapply(x,function(x)quantile(x,c(0.025,0.5,0.975),na.rm=TRUE))))
        }else{
          res$p<-list()
        }
        attr(res,"pData")<-new("AnnotatedDataFrame",pd[[i]])
        outfile<-get("outfile",environment(res$getmcmc))
        unlink(outfile)
        unlink(paste(outfile,"P",sep=""))
        res<-MCMCResult(object=res)
      }
      else if (method%in%"EM"){
        res<-MDMix(fitme)
        res@pd<-new("AnnotatedDataFrame",pd[[i]])
      }
      res<-new("MIMOSAResult",result=res)
      res
    })
  }else{
    stop("Can't run parallel MIMOSA. Must load multicore package.")
  }
  return(result)
})


#'pData extract the phenoData table from a MIMOSA result
#'@details Extracts the phenoData data.frame from a MIMOSAResult object
#'
#'@param object is the MIMOSAResult returned from a call to MIMOSA
#'@return an object of type \code{data.frame}
#'@importMethodsFrom Biobase pData
#'@aliases pData,MIMOSAResult-methods
#'@rdname pData-methods
setMethod("pData","MIMOSAResult",function(object){
  pData(object@result)
})

#'@rdname pData-methods
#'@aliases pData,MDMixResult-methods
setMethod("pData","MDMixResult",function(object){
  pData(object@pd)
})

#'@rdname pData-methods
#'@aliases pData,MCMCResult-methods
setMethod("pData","MCMCResult",function(object){
  pData(object@phenoData)
})


#'roc computes an ROC curve
#'
#'@param p is the probability of a positive result
#'@param truth is a logical with the true positive and negative results
roc<-function(p,truth){
  s<-seq(0,1,l=1000)
  table<-t(sapply(s,function(th){
    prop.table(table(test=factor(p<=th,levels=c("FALSE","TRUE")),truth=truth),margin=2)["TRUE",]
  }))
  colnames(table)<-c("FPR","TPR")
  table
}

#'fdrComparison calculates the observed vs expected false discovery rate
#'@param fdr is the expected false discovery rate
#'@param truth is a logical with the true positive and negative results
#'@export
fdrComparison<-function(fdr,truth){
  truth<-truth[order(fdr)]
  true.fdr<-cumsum(1-truth)/1:length(truth)
  fdr<-sort(fdr)
  cbind(fdr,true.fdr)
}


#' Construct an ExpressionSet for MIMOSA
#' 
#' Starting from a reshaped data frame in the correct format, construct
#' an ExpressionSet object that can be used with MIMOSA.
#' 
#' @details The featureCols will be used to construct feature names, and these
#'  columns will be dropped from the exprs matrix. The column names are assumed
#'  to have names that contain "_" characters separating phenotypic
#'  characteristics. These would be generated automatically
#'  if the data frame was constrcuted with "reshape". They
#'  are used to construct the phenoData for the expression set
#' @param df a data.frame that is in the correct form
#' @param featureCols the indices of the columns that identify features.
#' @export
MIMOSAExpressionSet<-function(df,featureCols){
  #feature<-df[,featureCols,drop=FALSE] 
  
  featuredata<-attributes(df)$rdimnames[[1]]
  pdata<-attributes(df)$rdimnames[[2]]
  #datanames<-colnames(melt(df)) #slow
  #datanames<-colnames(df)[-featureCols]
  #browser()
  #pnames<-datanames[(which(datanames%in%"value")+1):length(datanames)]
  
  df<-df[,-featureCols,drop=FALSE] #adata
  
  #pdata<-data.frame(do.call(rbind,strsplit(colnames(df),"_")))
  
  #fnames<-apply(as.data.frame(feature),1,function(x)paste(x,collapse="_"))
  
  #colnames(pdata)<-pnames
  #rownames(feature)<-fnames
  rownames(df)<-rownames(featuredata)
  #rownames(pdata)<-colnames(df)
  
  E<-ExpressionSet(as.matrix(as.data.frame(df)),featureData=AnnotatedDataFrame(as.data.frame(featuredata)),phenoData=AnnotatedDataFrame(pdata))
  return(E)
}


#'Replicate the reference observations across all treatment groups using ddply
#'
#'@details Should be passed to .fun argument of ddply.
#'
#'@usage ddply(data, .(VAR1,VAR2,VAR3),.fun=setReference,ref=VAR4%in%"Reference",cols=c("VAR5","VAR6"))
#'@param dat the piece we will work on
#'@param ref an expression evaluating to a logical vector that identifies the reference class within the piece.
#'@param cols A character vector of the names of the columns that hold our measurements
#@param annotations a characte vector of additional annotations to return for the data
#'@export
setReference<-function(dat,ref=NULL,cols=NULL,annotations=NULL,default.formula=component~...){
  REFERENCE<-try(eval(ref,dat))
  if(inherits(REFERENCE,"try-error")){
    REFERENCE<-eval(substitute(ref),dat)
  }
  if(inherits(REFERENCE,"try-error")){
    stop("Cannot evaluate ref argument")
  }
  #if((nrow(dat)!=2)|nlevels(factor(REFERENCE))==1){
#	return(NULL)
 # }
  if(sum(REFERENCE)!=1){
	return(NULL)
  }
  MEASUREMENTS<-do.call(cbind,with(dat,mget(cols,envir=as.environment(-1L))))
  #MEASUREMENTS<-model.frame(as.formula(paste("~",paste(cols,collapse="+"))),dat)
  NEWNAMES<-paste(cols,"REF",sep="_")
  REF.MEAS<-MEASUREMENTS[REFERENCE,,drop=FALSE]
  TREAT.MEAS<-MEASUREMENTS[!REFERENCE,,drop=FALSE]
  if(nrow(TREAT.MEAS)==0){
	return(NULL)
  }
  REF.MEAS<-matrix(REF.MEAS,nrow=dim(TREAT.MEAS)[1],ncol=dim(TREAT.MEAS)[2],byrow=TRUE)
  colnames(REF.MEAS)<-NEWNAMES
  retme<-cbind(TREAT.MEAS,REF.MEAS)
  if(!is.null(annotations)){
    ANNOTATIONS<-do.call(data.frame,with(dat,mget(annotations,envir=as.environment(-1L))))[!REFERENCE,,drop=FALSE]
    #ANNOTATIONS<-model.frame(as.formula(paste("~",paste(annotations,collapse="+"))),dat)[!REFERENCE,,drop=FALSE]
    retme<-cbind(ANNOTATIONS,retme)
  }
  retme
}

##'A wrapper for constructing an Expression Set for MIMOSA
##'
##'Calls a series of other functions that will reshape and refactor the data frame into the right format for use by MIMOSA
##'We provide some default arguments as examples. Currently slow, and very much prototype code.
##'@param \code{thisdata} is the input data frame
##'@param \code{reference} is an \code{expression} that evaluates to a \code{logical} vector which specifices the observations in the data frame that are to be used for the negative control or reference set
##'@param \code{measure.columns} is a \code{chracter} vector that specifies which columns hold the observed counts
##'@param \code{other.annotations} is a \code{character} vector that specifies which additional columns in the data frame should be included in the returned data. By default we take everything, but you could specify only relevant phenotypic information.
##'@param \code{default.cast.formula} is a \code{formula} that tells reshape how to recast the data frame so that rows corresponde to different measured components and columns correspond to samples. By default \code{component~...} will put the components as the rows (i.e. positive and negative cell counts) and all measured phenotypic information on the columns.
##'@param \code{.variable} is a dotted list that specifies the variable names (columns of the data frame) by which to group the data when organzing stimulated and unstimulated observations. i.e. PTID x ANTIGEN x TCELLSUBSET x TESTDT, or something else for your own data.
##'@param \code{featureCols} is a \code{numeric} vector that specifies the indices of the columns to be used to name the features. If the casting formula is \code{component~...} then there is only one feature column (and it is the first one), so \code{featureCols = 1}, by default.
##'@export
ConstructMIMOSAExpressionSet<-function(thisdata,reference=STAGE%in%"CTRL"&PROTEIN%in%"Media+cells",measure.columns=c("Neg","Pos"),other.annotations=setdiff(colnames(thisdata),measure.columns),default.cast.formula=component~...,.variables=.(PTID,TESTDT,ASSAYID,PLATEID),featureCols=1,ref.append.replace="NEG"){
  #if reference is null, then we already have the data in a form we need
  if(!is.null(substitute(reference))){
  #Set the Reference Class to be the negative control Media+cells for each ptid/date/assayid/plate
    thisdata<-ddply(thisdata,.variables=.variables,setReference,ref=substitute(reference),cols=measure.columns,annotations=other.annotations,default.formula=default.cast.formula)
  }else{
    colnames(thisdata)<-gsub(ref.append.replace,"_REF",colnames(thisdata))
  }
  #transform the data so that features are on the rows and samples are on the columns
  #build a function to reshape the returned data and assign it to the calling environment
  MIMOSAReshape<-function(mydata=NULL,default.formula=NULL,cols=measure.columns){
    if(!grepl("RefTreat",paste(deparse(default.formula),collapse=""))){
	default.formula<-as.formula(paste(paste(deparse(default.formula),collapse=""),"RefTreat",sep="+"))
    }
	NEWNAMES<-paste(cols,"REF",sep="_")
    mydata<-melt(mydata,measure.var=c(cols,NEWNAMES))
    mydata$RefTreat<-factor(grepl("_REF$",mydata$variable),labels=c("Treatment","Reference"))
    mydata<-rename(mydata,c("variable"="component","value"="count"))
    mydata$component<-factor(gsub("_REF$","",mydata$component))
    mydata<-cast(melt(mydata,measure="count"),default.formula)
    mydata
  }
  thisdata<-MIMOSAReshape(mydata=thisdata,default.formula=default.cast.formula,cols=measure.columns)
  MIMOSAExpressionSet(thisdata,featureCols=featureCols)
}
#'Sumarize the replicates in an elispot epitope mapping data set from SCHARP
#'
#' Thus function will summarize the replicate observations for the negative controls and stimulations in an epitope mapping data set from SCHARP. The user provides information on grouping structure, cells per well, replicate observation columns and so forth. The function is meant to be passed to ddply
#'@param x is the piece of data currently being worked on. It is a unique combination of the grouping variables passed to ddply
#'@param CONTROL is an expression that evaluates to a logical vector specifying which rows of the data frame are control observations.
#'@param WELLMULTIPLIER is a numeric argument that specifies the factor by which to multiply the denominator (in the \code{cells.per.well} argument). i.e., if there are 100K cells per well, but the cells.per.well column has 100, the WELLMULTIPLIER would be 1000
#'@param replicates is a character vector identifying the column names that contain the replicated observations. These will be summed.
#'@param cells.per.well is a character vector that identifies the column containing the number of cells per well.
#'@param SAMPLE is an expression evaluating to a logical vector that identifies the rows of the data frame which are antigen or peptide stimualtions rather than control samples.
#'@export
match.elispot.antigens<-function(x,CONTROL=STAGE%in%"CTRL"&PROTEIN%in%"Media+cells",WELLMULTIPLIER=1000,replicates=c("REP1","REP2","REP3","REP4","REP5","REP6"),cells.per.well="CELLWELL",SAMPLE=!STAGE%in%"CTRL"){
  control.sub<-eval(substitute(CONTROL),x)
  ctrl<-subset(x,control.sub)
  CPOS<-sum(ctrl[,replicates],na.rm=TRUE)
  l<-length(na.omit(t(ctrl[,replicates,drop=FALSE])))
  CNEG<-WELLMULTIPLIER*l*ctrl[,cells.per.well]
  
  samples<-eval(substitute(SAMPLE),x)
  samps<-subset(x,samples)
  
  PNEG<-NULL
  PPOS<-NULL
  for(i in 1:nrow(samps)){
    ppos<-na.omit(t(samps[i,replicates,drop=FALSE]))
    l<-length(ppos)
    pneg<-samps[i,cells.per.well]*l*WELLMULTIPLIER
    ppos<-sum(ppos)
    PNEG<-c(PNEG,pneg)
    PPOS<-c(PPOS,ppos)
  }
  ret<-rbind(data.frame(samps,Neg=PNEG,Pos=PPOS),data.frame(ctrl,Neg=CNEG,Pos=CPOS))
}
