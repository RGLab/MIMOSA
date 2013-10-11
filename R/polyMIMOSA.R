#' Multivariate MIMOSA data structure to store single-cell data.
#' 
#' @name SingleCellData
#' @docType class
#' @param sc is a list of single-cell expression values, similar to a list of flowFrames with measured channels on the columns and cells on the rows.
#' These are generally the `positive` cells identified by some gating method.
#' @param counts is a vector of total cell counts for each sample
#' @param meta is a data.frame of metadata for the samples
#' @export
SingleCellData<-function(sc,counts,meta){
  if(!all(length(sc)==c(length(counts),dim(meta)[1]))){
    stop("Dimensions of arguments don't match.")
  }
  obj<-vector('list',3)
  obj$obs<-sc
  obj$totalcounts<-counts
  obj$meta<-meta
  class(obj)<-"SingleCellData"
  obj
}

#'Print a SingleCellData object
#'
#'@name print.SingleCellData
#'@export 
print.SingleCellData<-function(x, ...){
  cat(sprintf("A SingleCellData object with %s samples\n",length(obs(x))));
}

#'Accessor for observations
#'
#'@name obs
#'@export
obs<-function(x){
  UseMethod("obs")
}

#'@S3method obs SingleCellData
obs.SingleCellData<-function(x){
  x$obs;
}

#'Accessor for total counts
#'
#'@name counts
#'@export
counts<-function(x){
  UseMethod("counts")
}

#'@S3method counts SingleCellData
counts.SingleCellData<-function(x){
  x$totalcounts;
}

#'Accessor for metadata
#'
#'@name metadata
#'@export
metadata<-function(x){
  UseMethod("metadata")
}

#'@S3method metadata SingleCellData
metadata.SingleCellData<-function(x){
  x$meta;
}

#'Length of a SingleCellData object
#'
#'@name length.SingleCellData
#'@export 
length.SingleCellData<-function(x){
  length(obs(x));
}

setOldClass("SingleCellData")

#'Subset via bracket operator
#'@name [.SingleCellData
#'@export
`[.SingleCellData`<-function(x,i,...){
  if(class(i)=="factor")
    i<-as.character(i)
  if(!missing(i)&(class(i)=="character")){
    obs(x)<-obs(x)[i]
    counts(x)<-counts(x)[i]
    metadata(x)<-metadata(x)[i,]
  }
  x
}

#'Replacement method for observations
#'@export
setGeneric("obs<-",function(object,value){
  standardGeneric("obs<-")
})

setMethod("obs<-",c("SingleCellData"),function(object,value){
  object$obs<-value
  object
})

#'Replacement method for metadata
#'@export
setGeneric("metadata<-",function(object,value){
  standardGeneric("metadata<-")
})

setMethod("metadata<-",c("SingleCellData"),function(object,value){
  object$meta<-value
  object
})

#'Replacement method for counts
#'@export
setGeneric("counts<-",function(object,value){
  standardGeneric("counts<-")
})

setMethod("counts<-",c("SingleCellData"),function(object,value){
  object$totalcounts<-value
  object
})

setMethod("sampleNames",c("SingleCellData"),function(object){
    names(obs(object))
})
