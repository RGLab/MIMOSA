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
#'@useDynLib MIMOSA, .registration=TRUE
#'@rdname MIMOSA-package
#'@import Formula
#'@importFrom testthat test_dir expect_that context 
#'@importFrom MASS ginv
#'@importClassesFrom methods array character data.frame factor integer matrix numeric
#'@importFrom reshape rename
#'@importFrom plyr ddply
#'@importFrom methods new
#'@importFrom Rcpp evalCpp
#'@name MIMOSA-package
#'@seealso \code{\link{MIMOSA}}, \code{\link{ConstructMIMOSAExpressionSet}}
#'@references Greg Finak, Andrew McDavid, Pratip Chattopadhyay, Maria Dominguez, Stephen C De Rosa, Mario Roederer, Raphael Gottardo
#'  Mixture Models for Single Cell Assays with Applications to Vaccine Studies
#'  Biostatistics, 2013, \url{http://biostatistics.oxfordjournals.org/content/early/2013/07/24/biostatistics.kxt024.abstract}
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
#'@param ... additional arguments
#'@return an object of type \code{MIMOSAResult}
#'@aliases MIMOSA,formula,ExpressionSet-method
#'@importFrom data.table key 
#'@importFrom coda mcmc
#'@examples
#' data(ICS)
#' E<-ConstructMIMOSAExpressionSet(ICS,
#'   reference=ANTIGEN%in%'negctrl',measure.columns=c('CYTNUM','NSUB'),
#'   other.annotations=c('CYTOKINE','TCELLSUBSET','ANTIGEN','UID'),
#'   default.cast.formula=component~UID+ANTIGEN+CYTOKINE+TCELLSUBSET,
#'   .variables=.(TCELLSUBSET,CYTOKINE,UID),
#'   featureCols=1,ref.append.replace='_REF')
#'
#' result<-MIMOSA(NSUB+CYTNUM~UID+TCELLSUBSET+CYTOKINE|ANTIGEN,
#'     data=E, method='EM',
#'     subset=RefTreat%in%'Treatment'&ANTIGEN%in%'ENV',
#'     ref=ANTIGEN%in%'ENV'&RefTreat%in%'Reference')
#'@seealso \code{\link{MIMOSA-package}} \code{\link{ConstructMIMOSAExpressionSet}} \code{\link{MIMOSAResult}}
#'@export
setGeneric("MIMOSA", def = function(formula, data, ...) {
    standardGeneric("MIMOSA")
})


setMethod("MIMOSA", c("formula", "ExpressionSet"), definition = function(formula,
    data, ref = RefTreat %in% "Reference", RT = TRUE, method = "mcmc", subset = RefTreat %in%
        "Treatment", getP = FALSE, p.thin = 1, run.parallel = FALSE, cleanup = TRUE,
    ...) {
    if (!exists("ref")) {
        stop("ref must contain expressions that define subsets to be compared")
    }
    if (!inherits(formula, "Formula")) {
        formula <- Formula(formula)
    }
    # Does the formula contain RefTreat?
    if (!any(attr(terms(formula), "term.labels") %in% "RefTreat") & RT) {
        warning("Formula does not contain the RefTreat variable. It will be added automatically. Set RT=FALSE to disable this.")
        formula <- Formula(formula(gsub("\\|", "+ RefTreat |", paste0(deparse(formula),
            collapse = ""))))
    }

    method <- match.arg(method, c("mcmc", "EM"))
   # if (Sys.info()["sysname"] == "Darwin") {
    #    run.parallel <- FALSE
    #    warning("Parallel is disabled on the Mac due to a bug that is yet to be tracked down.\n Your analysis will run in serial.\n Grab a coffee and wait..\n")
    #}
    mf.ref <- model.frame(formula, data, na.action = NULL)
    mf.test <- model.frame(formula, data, na.action = NULL)

    if (!missing(ref)) {
        mf.ref <- mf.ref[eval(substitute(ref), pData(data)), , drop = FALSE]
    }
    if (!missing(subset)) {
        mf.test <- mf.test[eval(substitute(subset), pData(data)), , drop = FALSE]
    }

    if (RT) {
        if (nlevels(factor(mf.ref$RefTreat)) != 1) {
            mf.ref <- mf.ref[mf.ref$RefTreat %in% "Reference", ]
        }
        if (nlevels(factor(mf.test$RefTreat)) != 1) {
            mf.test <- mf.test[mf.test$RefTreat %in% "Treatment", ]
        }
    }

    if (length(formula)[2] > 1) {
        spl.ref <- model.part(formula, mf.ref, rhs = 2)
        spl.test <- model.part(formula, mf.test, rhs = 2)
    }
    interact <- function(x, ...) {
        if (length(x) == 1) {
            return(factor(x[[1]]))
        }
        if (length(x) == 2) {
            interaction(x[[1]], x[[2]], drop = TRUE, ...)
        } else {
            interaction(x[[1]], Recall(x[-1L]), ...)
        }
    }
    if (length(formula)[2] > 1) {
        ref <- split(model.part(formula, mf.ref, lhs = 1), factor(interact(spl.ref)))
        test <- split(model.part(formula, mf.test, lhs = 1), factor(interact(spl.test)))
        pd <- split(model.part(formula, mf.test, rhs = 1:2), factor(interact(spl.test)))
        result <- vector("list", length(test))
        # Here should filter stuff that has empty categories.
        filter_empty <- do.call(c, lapply(ref, function(x) dim(x)[1] != 0))
        ref <- ref[filter_empty]
        test <- test[filter_empty]
        pd <- pd[filter_empty]
    } else {
        ref <- list(model.part(formula, mf.ref, lhs = 1))
        test <- list(model.part(formula, mf.test, lhs = 1))
        pd <- list(model.part(formula, mf.test, rhs = 1))
        result <- vector("list", 1)
        filter_empty <- do.call(c, lapply(ref, function(x) dim(x)[1] != 0))
        ref <- ref[filter_empty]
        test <- test[filter_empty]
        pd <- pd[filter_empty]
    }
    if (!run.parallel) {
        for (i in 1:length(test)) {
            # recycle the reference
            j <- data.frame(1:length(test), 1:length(ref))[i, 2]
            fitme <- list(n.stim = test[[i]], n.unstim = ref[[j]])
            # remove NAs
            nas <- unique(rbind(which(is.na(fitme[[1]]), TRUE), which(is.na(fitme[[2]]),
                TRUE))[, 1])
            # also remove zeros
            if (length(nas) > 0) {
                fitme[[1]] <- fitme[[1]][-nas, ]
                fitme[[2]] <- fitme[[2]][-nas, ]
                pd[[i]] <- pd[[i]][-nas, ]
            }
            # don't fit empty data! These are now filtered out earlier..
            if (nrow(pd[[i]]) > 0) {
                if (method %in% "mcmc") {
                  res <- .fitMCMC(fitme, inits = MDMix(fitme, initonly = TRUE), ...)
                  res$params <- res$params <- apply(res$getmcmc(), 2, function(x) quantile(x,
                    c(0.025, 0.5, 0.975), na.rm = TRUE))
                  if (ncol(fitme[[1]]) == 2 & getP) {
                    res$p <- lapply(res$getP(thin = p.thin), function(x) do.call(rbind,
                      lapply(x, function(x) quantile(x, c(0.025, 0.5, 0.975), na.rm = TRUE))))
                  } else {
                    res$p <- list()
                  }
                  if(ncol(fitme[[1]]) == 2){
                      res$IndMat<-res$getRespInd()
                  }
                  attr(res, "pData") <- new("AnnotatedDataFrame", pd[[i]])
                  colnames(res$IndMat)<-rownames(pd[[i]])
                  outfile <- get("outfile", environment(res$getmcmc))
                  if (cleanup) {
                    unlink(outfile)
                    unlink(paste(outfile, "P", sep = ""))
                  }
                  res <- MCMCResult(object = res)
                } else if (method %in% "EM") {
                  res <- try(MDMix(fitme))
                  if (inherits(res, "try-error")) {
                    message("Failed to fit subset ", i, " with EM. Trying mcmc")
                    res <- .fitMCMC(fitme, inits = MDMix(fitme, initonly = TRUE),
                      ...)
                    res$params <- res$params <- apply(res$getmcmc(), 2, function(x) quantile(x,
                      c(0.025, 0.5, 0.975), na.rm = TRUE))
                    if (ncol(fitme[[1]]) == 2 & getP) {
                      res$p <- lapply(res$getP(thin = p.thin), function(x) do.call(rbind,
                        lapply(x, function(x) quantile(x, c(0.025, 0.5, 0.975), na.rm = TRUE))))
                    } else {
                      res$p <- list()
                    }
                    attr(res, "pData") <- new("AnnotatedDataFrame", pd[[i]])
                    outfile <- get("outfile", environment(res$getmcmc))
                    if (cleanup) {
                      unlink(outfile)
                      unlink(paste(outfile, "P", sep = ""))
                    }
                    res <- MCMCResult(object = res)
                    if(!all.equal(res@IndMat,data.frame())){
                        colnames(res@IndMat)<-rownames(pData(res))
                    }
                  } else {
                    res@pd <- new("AnnotatedDataFrame", pd[[i]])
                  }
                }
                res <- new("MIMOSAResult", result = res)
                result[[i]] <- res
            } else {
                result[[i]] <- NA
            }
        }
    } else if (run.parallel & any(grepl("parallel", loadedNamespaces()))) {
        result <- mclapply(1:length(test), function(i) {
            # recycle the reference
            j <- data.frame(1:length(test), 1:length(ref))[i, 2]
            fitme <- list(n.stim = test[[i]], n.unstim = ref[[j]])
            # remove NAs
            nas <- unique(rbind(which(is.na(fitme[[1]]), TRUE), which(is.na(fitme[[2]]),
                TRUE))[, 1])
            if (length(nas) > 0) {
                fitme[[1]] <- fitme[[1]][-nas, ]
                fitme[[2]] <- fitme[[2]][-nas, ]
                pd[[i]] <- pd[[i]][-nas, ]
            }
            if (nrow(pd[[i]]) > 0) {
                if (method %in% "mcmc") {
                  res <- .fitMCMC(fitme, inits = MDMix(fitme, initonly = TRUE), ...)
                  res$params <- res$params <- apply(res$getmcmc(), 2, function(x) quantile(x,
                    c(0.025, 0.5, 0.975), na.rm = TRUE))
                  if (ncol(fitme[[1]]) == 2 & getP) {
                    res$p <- lapply(res$getP(thin = p.thin), function(x) do.call(rbind,
                      lapply(x, function(x) quantile(x, c(0.025, 0.5, 0.975), na.rm = TRUE))))
                  } else {
                    res$p <- list()
                  }
                  attr(res, "pData") <- new("AnnotatedDataFrame", pd[[i]])
                  outfile <- get("outfile", environment(res$getmcmc))
                  unlink(outfile)
                  unlink(paste(outfile, "P", sep = ""))
                  res <- MCMCResult(object = res)
                } else if (method %in% "EM") {
                  res <- try(MDMix(fitme))
                  res <- try(MDMix(fitme))
                  if (inherits(res, "try-error")) {
                    message("Failed to fit subset ", i, " with EM. Trying mcmc")
                    res <- .fitMCMC(fitme, inits = MDMix(fitme, initonly = TRUE),
                      ...)
                    res$params <- res$params <- apply(res$getmcmc(), 2, function(x) quantile(x,
                      c(0.025, 0.5, 0.975), na.rm = TRUE))
                    if (ncol(fitme[[1]]) == 2 & getP) {
                      res$p <- lapply(res$getP(thin = p.thin), function(x) do.call(rbind,
                        lapply(x, function(x) quantile(x, c(0.025, 0.5, 0.975), na.rm = TRUE))))
                    } else {
                      res$p <- list()
                    }
                    attr(res, "pData") <- new("AnnotatedDataFrame", pd[[i]])
                    outfile <- get("outfile", environment(res$getmcmc))
                    unlink(outfile)
                    unlink(paste(outfile, "P", sep = ""))
                    res <- MCMCResult(object = res)
                  } else {
                    res@pd <- new("AnnotatedDataFrame", pd[[i]])
                  }
                }
                res <- new("MIMOSAResult", result = res)
                return(res)
            } else {
                return(NA)
            }
        })
    } else {
        stop("Can't run parallel MIMOSA. Must load multicore package.")
    }
    class(result) <- c("MIMOSAResultList", "list")
    if (length(result) > 0) {
        splitted <- strsplit(paste0(deparse(formula[[3]]), collapse = ""), "|", fixed = TRUE)[[1]]
        if (length(splitted) > 1) {
            depf <- strsplit(gsub(" ", "", splitted[[2]]), "+", fixed = TRUE)[[1]]
        } else {
            depf <- NULL
        }
        n_vars <- length(depf)
        # NOTE this is probably wrong..
        if (n_vars > 0) {
          #browser()
            #names(result) <- levels(factor(interaction(as.data.frame(lapply(pData(result)[,
            #    depf, with = FALSE], factor)))))
            names(result) <- names(test)
        } else {
            names(result) <- "result"
        }
    }
    return(result)
})

#'@importFrom Biobase `pData<-`
#'@title pData
#'@docType methods
#'@rdname pData-methods
#'@aliases `pData<-`,MCMCResult, ANY
#'@aliases `pData<-`,MIMOSAResult, ANY
#'@param value the phenoData table to be assigned to the object.
#'@export
setMethod("pData<-","MCMCResult",function(object,value){
    pData(object@result@phenoData)<-value
    object
})

#'@export
#'@rdname pData-methods
setMethod("pData<-","MIMOSAResult",function(object,value){
    pData(object@result@phenoData)<-value
    object
})

#'pData extract the phenoData table from a MIMOSA result
#'@details Extracts the phenoData data.frame from a MIMOSAResult object
#'
#'@param object is the MIMOSAResult returned from a call to MIMOSA
#'@return an object of type \code{data.frame}
#'@importMethodsFrom Biobase pData
#'@aliases pData,MIMOSAResult-method
#'@method pData MIMOSAResult
#'@rdname pData-methods
#'@export
setMethod("pData", "MIMOSAResult", function(object) {
    pData(object@result)
})


#'@method pData MDMixResult
#'@aliases pData,MDMixResult-method
#'@rdname pData-methods
#'@export
setMethod("pData", "MDMixResult", function(object) {
    pData(object@pd)
})

#'@aliases pData,MCMCResult-method
#'@method pData MCMCResult
#'@rdname pData-methods
#'@export
setMethod("pData", "MCMCResult", function(object) {
    pData(object@phenoData)
})


# 'roc computes an ROC curve ' '@param p is the probability of a positive result
# '@param truth is a logical with the true positive and negative results
roc <- function(p, truth) {
    s <- seq(0, 1, l = 1000)
    table <- t(sapply(s, function(th) {
        prop.table(table(test = factor(p <= th, levels = c("FALSE", "TRUE")), truth = truth),
            margin = 2)["TRUE", ]
    }))
    colnames(table) <- c("FPR", "TPR")
    table
}

# 'fdrComparison calculates the observed vs expected false discovery rate '@param
# fdr is the expected false discovery rate '@param truth is a logical with the
# true positive and negative results
fdrComparison <- function(fdr, truth) {
    truth <- truth[order(fdr)]
    true.fdr <- cumsum(1 - truth)/1:length(truth)
    fdr <- sort(fdr)
    cbind(fdr, true.fdr)
}


#' Construct an ExpressionSet for MIMOSA
#'
#' Starting from a reshaped data frame in the correct format, construct
#' an ExpressionSet object that can be used with MIMOSA.
#'
#' @details The featureCols will be used to construct feature names, and these
#'  columns will be dropped from the exprs matrix. The column names are assumed
#'  to have names that contain '_' characters separating phenotypic
#'  characteristics. These would be generated automatically
#'  if the data frame was constrcuted with 'reshape'. They
#'  are used to construct the phenoData for the expression set
#' @param df a data.frame that is in the correct form
#' @param featureCols the indices of the columns that identify features.
#' @examples
#' E<-ConstructMIMOSAExpressionSet(ICS,
##'   reference=ANTIGEN%in%'negctrl',measure.columns=c('CYTNUM','NSUB'),
##'   other.annotations=c('CYTOKINE','TCELLSUBSET','ANTIGEN','UID'),
##'   default.cast.formula=component~UID+ANTIGEN+CYTOKINE+TCELLSUBSET,
##'   .variables=.(TCELLSUBSET,CYTOKINE,UID),
##'   featureCols=1,ref.append.replace='_REF')
#' @export
MIMOSAExpressionSet <- function(df, featureCols) {

    featuredata <- attributes(df)$rdimnames[[1]]
    pdata <- attributes(df)$rdimnames[[2]]
    df <- df[, -featureCols, drop = FALSE]  #adata

    rownames(df) <- rownames(featuredata)

    E <- ExpressionSet(as.matrix(as.data.frame(df)), featureData = AnnotatedDataFrame(as.data.frame(featuredata)),
        phenoData = AnnotatedDataFrame(pdata))
    return(E)
}


# 'Replicate the reference observations across all treatment groups using ddply '
# '@details Should be passed to .fun argument of ddply.  ' '@param dat the piece
# we will work on '@param ref an expression evaluating to a logical vector that
# identifies the reference class within the piece.  '@param cols A character
# vector of the names of the columns that hold our measurements '@param
# annotations A character vector of additional annotation columns '@param
# default.formula a default formula to be used for casting the data.
setReference <- function(dat, ref = NULL, cols = NULL, annotations = NULL, default.formula = component ~
    ...) {
    REFERENCE <- try(eval(ref, dat))
    if (inherits(REFERENCE, "try-error")) {
        REFERENCE <- eval(substitute(ref), dat)
    }
    if (inherits(REFERENCE, "try-error")) {
        stop("Cannot evaluate ref argument")
    }
    # if((nrow(dat)!=2)|nlevels(factor(REFERENCE))==1){ return(NULL) }
    if (sum(REFERENCE) != 1) {
        return(NULL)
    }
    MEASUREMENTS <- do.call(cbind, with(dat, mget(cols, envir = as.environment(-1L))))
    # MEASUREMENTS<-model.frame(as.formula(paste('~',paste(cols,collapse='+'))),dat)
    NEWNAMES <- paste(cols, "REF", sep = "_")
    REF.MEAS <- MEASUREMENTS[REFERENCE, , drop = FALSE]
    TREAT.MEAS <- MEASUREMENTS[!REFERENCE, , drop = FALSE]
    if (nrow(TREAT.MEAS) == 0) {
        return(NULL)
    }
    REF.MEAS <- matrix(REF.MEAS, nrow = dim(TREAT.MEAS)[1], ncol = dim(TREAT.MEAS)[2],
        byrow = TRUE)
    colnames(REF.MEAS) <- NEWNAMES
    retme <- cbind(TREAT.MEAS, REF.MEAS)
    if (!is.null(annotations)) {
        ANNOTATIONS <- do.call(data.frame, with(dat, mget(annotations, envir = as.environment(-1L))))[!REFERENCE,
            , drop = FALSE]
        # ANNOTATIONS<-model.frame(as.formula(paste('~',paste(annotations,collapse='+'))),dat)[!REFERENCE,,drop=FALSE]
        retme <- cbind(ANNOTATIONS, retme)
    }
    retme
}

##'A wrapper for constructing an Expression Set for MIMOSA
##'
##'Calls a series of other functions that will reshape and refactor the data frame into the right format for use by MIMOSA
##'Standardized for use with internal SCHARP data sets.
##'We provide some default arguments as examples. Currently slow, and very much prototype code.
##'@param thisdata is the input data frame
##'@param reference is an \code{expression} that evaluates to a \code{logical} vector which specifices the observations in the data frame that are to be used for the negative control or reference set
##'@param measure.columns is a \code{chracter} vector that specifies which columns hold the observed counts
##'@param other.annotations is a \code{character} vector that specifies which additional columns in the data frame should be included in the returned data. By default we take everything, but you could specify only relevant phenotypic information.
##'@param default.cast.formula is a \code{formula} that tells reshape how to recast the data frame so that rows corresponde to different measured components and columns correspond to samples. By default \code{component~...} will put the components as the rows (i.e. positive and negative cell counts) and all measured phenotypic information on the columns.
##'@param .variables is a dotted list that specifies the variable names (columns of the data frame) by which to group the data when organzing stimulated and unstimulated observations. i.e. PTID x ANTIGEN x TCELLSUBSET x TESTDT, or something else for your own data.
##'@param featureCols is a \code{numeric} vector that specifies the indices of the columns to be used to name the features. If the casting formula is \code{component~...} then there is only one feature column (and it is the first one), so \code{featureCols = 1}, by default.
##'@param ref.append.replace the terminating character string in the column names of the negative controls. It will be replaces with _REF for 'reference'
##'@examples
##' data(ICS)
##' E<-ConstructMIMOSAExpressionSet(ICS,
##'   reference=ANTIGEN%in%'negctrl',measure.columns=c('CYTNUM','NSUB'),
##'   other.annotations=c('CYTOKINE','TCELLSUBSET','ANTIGEN','UID'),
##'   default.cast.formula=component~UID+ANTIGEN+CYTOKINE+TCELLSUBSET,
##'   .variables=.(TCELLSUBSET,CYTOKINE,UID),
##'   featureCols=1,ref.append.replace='_REF')
##'
##'@export
ConstructMIMOSAExpressionSet <- function(thisdata, reference = quote(STAGE %in% "CTRL" &
    PROTEIN %in% "Media+cells"), measure.columns = c("Neg", "Pos"), other.annotations = setdiff(colnames(thisdata),
    measure.columns), default.cast.formula = component ~ ..., .variables = quote(.(PTID,
    TESTDT, ASSAYID, PLATEID)), featureCols = 1, ref.append.replace = "_NEG") {
    # if reference is null, then we already have the data in a form we need
    if (!is.null(substitute(reference))) {
        # Set the Reference Class to be the negative control Media+cells for each
        # ptid/date/assayid/plate
        thisdata <- ddply(thisdata, .variables = .variables, setReference, ref = substitute(reference),
            cols = measure.columns, annotations = other.annotations, default.formula = default.cast.formula)
    } else {
        colnames(thisdata) <- gsub(ref.append.replace, "_REF", colnames(thisdata))
    }
    if (nrow(thisdata) == 0)
        stop("All your data was filtered when reshaping due to non-unique pairs of samples. Perhaps you need to aggregate negative controls")

    # transform the data so that features are on the rows and samples are on the
    # columns build a function to reshape the returned data and assign it to the
    # calling environment
    MIMOSAReshape <- function(mydata = NULL, default.formula = NULL, cols = measure.columns) {
        if (!grepl("RefTreat", paste(deparse(default.formula), collapse = ""))) {
            default.formula <- as.formula(paste(paste(deparse(default.formula), collapse = ""),
                "RefTreat", sep = "+"))
        }
        NEWNAMES <- paste(cols, "REF", sep = "_")
        mydata <- melt(mydata, measure.var = c(cols, NEWNAMES))
        mydata$RefTreat <- factor(grepl("_REF$", mydata$variable), labels = c("Treatment",
            "Reference"))
        mydata <- reshape::rename(mydata, c(variable = "component", value = "count"))
        mydata$component <- factor(gsub("_REF$", "", mydata$component))
        mydata <- cast(melt(mydata, measure = "count"), default.formula)
        mydata
    }
    thisdata <- MIMOSAReshape(mydata = thisdata, default.formula = default.cast.formula,
        cols = measure.columns)
    MIMOSAExpressionSet(thisdata, featureCols = featureCols)
}


# 'Sumarize the replicates in an elispot epitope mapping data set from SCHARP ' '
# Thus function will summarize the replicate observations for the negative
# controls and stimulations in an epitope mapping data set from SCHARP. The user
# provides information on grouping structure, cells per well, replicate
# observation columns and so forth. The function is meant to be passed to ddply
# '@param x is the piece of data currently being worked on. It is a unique
# combination of the grouping variables passed to ddply '@param CONTROL is an
# expression that evaluates to a logical vector specifying which rows of the data
# frame are control observations.  '@param WELLMULTIPLIER is a numeric argument
# that specifies the factor by which to multiply the denominator (in the
# \code{cells.per.well} argument). i.e., if there are 100K cells per well, but
# the cells.per.well column has 100, the WELLMULTIPLIER would be 1000 '@param
# replicates is a character vector identifying the column names that contain the
# replicated observations. These will be summed.  '@param cells.per.well is a
# character vector that identifies the column containing the number of cells per
# well.  '@param SAMPLE is an expression evaluating to a logical vector that
# identifies the rows of the data frame which are antigen or peptide stimualtions
# rather than control samples.
match.elispot.antigens <- function(x, CONTROL = quote(STAGE %in% "CTRL" & PROTEIN %in%
    "Media+cells"), WELLMULTIPLIER = 1000, replicates = c("REP1", "REP2", "REP3",
    "REP4", "REP5", "REP6"), cells.per.well = "CELLWELL", SAMPLE = quote(!STAGE %in%
    "CTRL")) {
    control.sub <- eval(eval(substitute(CONTROL)), x)
    ctrl <- subset(x, control.sub)
    CPOS <- sum(ctrl[, replicates], na.rm = TRUE)
    l <- length(na.omit(t(ctrl[, replicates, drop = FALSE])))
    CNEG <- WELLMULTIPLIER * l * ctrl[, cells.per.well, drop = TRUE]

    samples <- eval(eval(substitute(SAMPLE)), x)
    samps <- subset(x, samples)

    PNEG <- NULL
    PPOS <- NULL
    for (i in 1:nrow(samps)) {
        ppos <- na.omit(t(samps[i, replicates, drop = FALSE]))
        l <- length(ppos)
        pneg <- samps[i, cells.per.well] * l * WELLMULTIPLIER
        ppos <- sum(ppos)
        PNEG <- c(PNEG, pneg)
        PPOS <- c(PPOS, ppos)
    }
    if (any(c(nrow(samps), nrow(ctrl)) %in% 0))
        return(NULL) else ret <- rbind(data.frame(samps, Neg = PNEG, Pos = PPOS), data.frame(ctrl,
        Neg = CNEG, Pos = CPOS))
    ret
}

setOldClass("MIMOSAResultList")

#'Print a MIMOSAResultList
#'
#'Print a summary of the list of results returned by a call to \code{MIMOSA}
#'@rdname print
#'@method print MIMOSAResultList
#'@param x a \code{MIMOSAResultList}
#'@param ... additional arguments passed down
print.MIMOSAResultList <- function(x, ...) {
    cat(sprintf("A MIMOSAResultList with %s models for\n", length(x)))
    cat(sprintf("%s ", names(x)))
}


#'Show a MIMOSAResultList
#'
#'Show a summary of a MIMOSAResultList.
#'@rdname show
#'@title show
#'@param object \code{MIMOSAResultList}
#'@name show
#'@aliases show,MIMOSAResult-method
setMethod("show", "MIMOSAResult", function(object) {
    cat("A MIMOSA Model with ")
    cat(sprintf("%s observations.\n", nrow(object@z)))
    cat(sprintf("Response rate (w) of %.4g percent.\n", 100 * object@w[2]))
})


#'@rdname pData-methods
#'@importFrom data.table rbindlist
#'@method pData MIMOSAResultList
pData.MIMOSAResultList <- function(object) {
    do.call(rbind,lapply(object, pData))
}


#'@rdname pData-methods
#'@method pData MIMOSAResultList
#'@aliases pData,MIMOSAResultList-method
setMethod("pData", "MIMOSAResultList", pData.MIMOSAResultList)

#'Extract the posterior probabilities of response from a MIMOSA model
#'
#'@rdname MIMOSA-accessors
#'@param x output from a MIMOSA model
#'@return a \code{matrix} of posterior probabilities
#'@examples
#'  data(ICS)
##' E<-ConstructMIMOSAExpressionSet(ICS,
##'   reference=ANTIGEN%in%'negctrl',measure.columns=c('CYTNUM','NSUB'),
##'   other.annotations=c('CYTOKINE','TCELLSUBSET','ANTIGEN','UID'),
##'   default.cast.formula=component~UID+ANTIGEN+CYTOKINE+TCELLSUBSET,
##'   .variables=.(TCELLSUBSET,CYTOKINE,UID),
##'   featureCols=1,ref.append.replace='_REF')
##'
##' result<-MIMOSA(NSUB+CYTNUM~UID+TCELLSUBSET+CYTOKINE|ANTIGEN,
##'     data=E, method='EM',
##'     subset=RefTreat%in%'Treatment'&ANTIGEN%in%'ENV',
##'     ref=ANTIGEN%in%'ENV'&RefTreat%in%'Reference')
##'     getZ(result)
#'@export
getZ <- function(x) {
    UseMethod("getZ")
}

#'@rdname MIMOSA-accessors
#'@method getZ MIMOSAResultList
#'@export
getZ.MIMOSAResultList <- function(x) {
    as.matrix(do.call(rbind, lapply(x, getZ)))
}

#'@rdname MIMOSA-accessors
#'@method getZ MIMOSAResult
#'@export
getZ.MIMOSAResult <- function(x) {
    z <- x@z
    colnames(z) <- c("Pr.Nonresponse", "Pr.response")
    z
}

#'Extract the component weights from a MIMOSA model
#'
#'@rdname MIMOSA-accessors
#'@return a \code{vector} of component weights
#'@examples
#'data(ICS)
##' E<-ConstructMIMOSAExpressionSet(ICS,
##'   reference=ANTIGEN%in%'negctrl',measure.columns=c('CYTNUM','NSUB'),
##'   other.annotations=c('CYTOKINE','TCELLSUBSET','ANTIGEN','UID'),
##'   default.cast.formula=component~UID+ANTIGEN+CYTOKINE+TCELLSUBSET,
##'   .variables=.(TCELLSUBSET,CYTOKINE,UID),
##'   featureCols=1,ref.append.replace='_REF')
##'
##' result<-MIMOSA(NSUB+CYTNUM~UID+TCELLSUBSET+CYTOKINE|ANTIGEN,
##'     data=E, method='EM',
##'     subset=RefTreat%in%'Treatment'&ANTIGEN%in%'ENV',
##'     ref=ANTIGEN%in%'ENV'&RefTreat%in%'Reference')
##' getW(result)
#'@export
getW <- function(x) {
    UseMethod("getW")
}

#'@rdname MIMOSA-accessors
#'@method getW MIMOSAResultList
#'@export
getW.MIMOSAResultList <- function(x) {
    w <- data.frame(lapply(x, getW))
    colnames(w) <- names(x)
    w
}

#'@rdname MIMOSA-accessors
#'@method getW MIMOSAResult
#'@export
getW.MIMOSAResult <- function(x) {
    w <- x@w
    names(w) <- c("w.nonresp", "w.resp")
    w
}

countsTable <- function(object, proportion = FALSE) {
}
#'Extract the table of counts from a MIMOSA model
#'
#'@param object a \code{MIMOSAResult}
#'@param proportion \code{logical} return the counts or the proportions
#'@return a \code{data.frame} of counts to which the model was fit.
#'@rdname countsTable
#'@docType methods
#'@return a \code{data.frame} of counts for the stimulated and unstimulated samples
#'@examples
##' data(ICS)
##' E<-ConstructMIMOSAExpressionSet(ICS,
##'   reference=ANTIGEN%in%'negctrl',measure.columns=c('CYTNUM','NSUB'),
##'   other.annotations=c('CYTOKINE','TCELLSUBSET','ANTIGEN','UID'),
##'   default.cast.formula=component~UID+ANTIGEN+CYTOKINE+TCELLSUBSET,
##'   .variables=.(TCELLSUBSET,CYTOKINE,UID),
##'   featureCols=1,ref.append.replace='_REF')
##'
##' result<-MIMOSA(NSUB+CYTNUM~UID+TCELLSUBSET+CYTOKINE|ANTIGEN,
##'     data=E, method='EM',
##'     subset=RefTreat%in%'Treatment'&ANTIGEN%in%'ENV',
##'     ref=ANTIGEN%in%'ENV'&RefTreat%in%'Reference')
##' head(countsTable(result))
##' head(countsTable(result,proportion=TRUE))
#'@export
setGeneric("countsTable")

#'@rdname countsTable
#'@method countsTable MIMOSAResult
#'@aliases countsTable,MIMOSAResult-method
#'@export
setMethod("countsTable", "MIMOSAResult", definition = function(object, proportion = FALSE) {
    countsTable(object@result, proportion = proportion)
})

#'@rdname countsTable
#'@method countsTable MCMCResult
#'@aliases countsTable,MCMCResult-method
#'@export
setMethod("countsTable", "MCMCResult", definition = function(object, proportion = FALSE) {
    if (proportion) {
        m <- (cbind(prop.table(as.matrix(object@n.stim), 1), prop.table(as.matrix(object@n.unstim),
            1)))
    } else {
        m <- (as.matrix(cbind(object@n.stim, object@n.unstim)))
    }
    nc <- ncol(m)
    colnames(m) <- paste0(colnames(m), rep(c("", "_REF"), each = nc/2))
    return(m)
})

#'@rdname countsTable
#'@method countsTable MDMixResult
#'@aliases countsTable,MDMixResult-method
#'@export
setMethod("countsTable", "MDMixResult", definition = function(object, proportion = FALSE) {
    if (proportion == TRUE) {
        m <- (do.call(cbind, lapply(object@data, function(x) prop.table(as.matrix(x),
            1))))
    } else {
        m <- (as.matrix(do.call(cbind, object@data)))
    }
    nc <- ncol(m)
    colnames(m) <- paste0(colnames(m), rep(c("", "_REF"), each = nc/2))
    return(m)
})

#'@rdname countsTable
#'@method countsTable MIMOSAResultList
#'@export
countsTable.MIMOSAResultList <- function(object, proportion = FALSE) {
    as.matrix(do.call(rbind, lapply(object, function(x) countsTable(x, proportion = proportion))))
}

#'@rdname countsTable
#'@method countsTable MIMOSAResultList
#'@aliases countsTable,MIMOSAResultList-method
#'@export
setMethod("countsTable", "MIMOSAResultList", countsTable.MIMOSAResultList)

#'Volcano plot for a MIMOSA model
#'
#'Plots effect size vs posterior probablilty of response from a MIMOSAResultList, faceting by the conditioning variables.
#'
#'@rdname volcanoPlot
#'@param x A \code{MIMOSAResultList}
#'@param effect_expression an \code{expression} that defines the effect size. Usually a function of the stimulated and unstimulated proportions from \code{countsTable(x,proportion=TRUE)}
#'@param facet_var an \code{expression} defining the faceting in ggplot parlance. i.e. \code{~ faceting + variables}
#'@param threshold a \code{numeric} value between [0,1] for coloring significant observations (based on the q-value)
#'@examples
#'data(ICS)
##' E<-ConstructMIMOSAExpressionSet(ICS,
##'   reference=ANTIGEN%in%'negctrl',measure.columns=c('CYTNUM','NSUB'),
##'   other.annotations=c('CYTOKINE','TCELLSUBSET','ANTIGEN','UID'),
##'   default.cast.formula=component~UID+ANTIGEN+CYTOKINE+TCELLSUBSET,
##'   .variables=.(TCELLSUBSET,CYTOKINE,UID),
##'   featureCols=1,ref.append.replace='_REF')
##'
##' result<-MIMOSA(NSUB+CYTNUM~UID+TCELLSUBSET+CYTOKINE|ANTIGEN,
##'     data=E, method='EM',
##'     subset=RefTreat%in%'Treatment'&ANTIGEN%in%'ENV',
##'     ref=ANTIGEN%in%'ENV'&RefTreat%in%'Reference')
##' volcanoPlot(result,CYTNUM-CYTNUM_REF)
#'@importFrom ggplot2 ggplot geom_point theme_bw aes_string scale_y_continuous
#'@seealso \code{\link{countsTable}}
#'@export
volcanoPlot <- function(x, effect_expression = NA, facet_var = NA, threshold = 0.01) {
    UseMethod("volcanoPlot")
}

#'@method volcanoPlot MIMOSAResultList
#'@importFrom data.table data.table
#'@export
volcanoPlot.MIMOSAResultList <- function(x, effect_expression = NA, facet_var = NA,
    threshold = 0.01) {
    err <- FALSE
    effect_expression <- deparse(substitute(effect_expression))
    if (effect_expression %in% "NA")
        err <- TRUE
    if (err) {
        stop("Must provide an expression for the effect size (i.e. CYTNUM-CYTNUM_REF)")
    }

    q <- -log10(unlist(fdr(x), use.names = FALSE))
    pspu <- countsTable(x, proportion = TRUE)
    p.stim <- getZ(x)
    pd <- pData(x)
    df <- data.table(q, pspu, p.stim, pd, signif = q > -log10(threshold))
    if (!is.na(effect_expression)) {
        p <- ggplot(df) + aes_string(x = effect_expression, y = "Pr.response", col = "signif.fdr") +
            geom_point() + theme_bw() + scale_y_continuous("Probability of Response")
    }
    if (is.formula(facet_var)) {
        p <- p + facet_grid(facet_var)
    }
    p
}

