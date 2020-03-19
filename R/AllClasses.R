# Classes.R Created on: Nov 30, 2011 Author: finak

# ' The output of fitting Beta-Binomial EM implementation \code{BetaMix}.  ' '
# BetaMix will return an object of this class.  '
setClass("BetaMixResult", representation = list(alternative.model = "character", 
    cytokine = "character", control = "character", stimulation = "character", ll = "numeric", 
    traj = "numeric", iter = "numeric", z = "matrix", w = "numeric", alpha0 = "numeric", 
    beta0 = "numeric", alphaS = "numeric", betaS = "numeric", fdr = "numeric", data = "data.frame", 
    pd = "AnnotatedDataFrame"))
# ' The output of fitting Dirichlet-Multinomial EM implementation \code{MDMix}.
# ' ' MDMix will return an object of this class.  ' ' @rdname MDMixResult '
# @docType class
setClass("MDMixResult", representation = list(llnull = "function", llresp = "function", 
    gresp = "function", w = "numeric", z = "matrix", hresp = "function", gnull = "function", 
    ll = "numeric", hnull = "function", par.unstim = "numeric", par.stim = "numeric", 
    data = "list", pd = "AnnotatedDataFrame"))



# 'Constructor for the BetaMixResult class ' 'BetaMix uses this constructor to
# create a BetaMixResult object ' '@rdname BetaMixResult '@docType function
# '@aliases BetaMixResult-class '@aliases BetaMixResult '@param alternative.model
# character, either 'greater' or 'not equal' '@param cytokine 'character' '@param
# control 'character' '@param stimulation 'character' '@param ll log-likelihood
# '@param traj ll trajectory during fitting '@param iter number of iterations
# '@param z the posterior probabilties '@param w the responder and non-responder
# fractions '@param alpha0 model parameters '@param alphaS model parameters
# '@param beta0 model parameters '@param betaS model parameters '@param data the
# dat '@param pd the phenoData
BetaMixResult <- function(alternative.model = NA_character_, cytokine = NA_character_, 
    control = NA_character_, stimulation = NA_character_, ll, traj, iter, z, w, alpha0, 
    beta0, alphaS, betaS, data, pd = new("AnnotatedDataFrame")) {
    new("BetaMixResult", alternative.model = alternative.model, control = control, 
        stimulation = stimulation, cytokine = cytokine, ll = ll, traj = traj, iter = iter, 
        z = z, w = w, alpha0 = alpha0, beta0 = beta0, alphaS = alphaS, betaS = betaS, 
        fdr = fdr(z), data = data, pd = pd)
}

# ' The output of fitting the Dirichlet--Multinomial MIMOSA model via MCMC ' '
# Fitting the model using the \code{MIMOSA} function returns an object of this
# class.  ' ' @docType class ' @rdname MCMCResult
setClass("MCMCResult", representation = list(z = "matrix", n.stim = "data.frame", 
    n.unstim = "data.frame", params = "matrix", p = "list", phenoData = "AnnotatedDataFrame",IndMat="data.frame"))

# ' Constructor for a \code{MCMCResult} object ' ' This constructor is used by
# the \code{MIMOSA} model fitting method.  ' ' @rdname MCMCResult ' @param
# object of class \code{MDMixResult},\code{MCMCResult},or \code{BetaMixResult}
MCMCResult <- function(object = NULL) {
    if (class(object) == "MDMixResult") {
        if (is.null(names(object))) {
            stop("Wrong input class")
        } else {
            return(new("MCMCResult", z = object$z, n.stim = object$n.stim, n.unstim = object$n.unstim, 
                params = object$params, p = object$p, phenoData = attr(object, "pData"),IndMat=if(!is.null(object$IndMat)){object$IndMat}else{data.frame()}))
        }
    }
}

# ' Represents the output of any implementation of MIMOSA ' ' Virtual class that
# combines the results of an \code{MDMixResult}, \code{MCMCResult} and
# \code{BetaMixResult} object.  ' ' @name MixResult
setClassUnion("MixResult", c("BetaMixResult", "MDMixResult", "MCMCResult"))

#' Stores the result of a MIMOSA fitted model
#' 
#' \code{MIMOSA} returns an object of \code{MIMOSAResult} irrespective of which
#' method / implementation is used to fit the data.
#' 
#' @docType class
setClass("MIMOSAResult", representation = list(result = "MixResult", z = "matrix", 
    w = "numeric", IndMat = "data.frame"))

setMethod("initialize", "MIMOSAResult", function(.Object, ...) {
    .Object <- callNextMethod(.Object, ...)
    .Object@z = .Object@result@z
    .Object@w = colMeans(.Object@z)
    if(class(.Object@result)=="MCMCResult"){
        .Object@IndMat = .Object@result@IndMat
    }else{
        .Object@IndMat=data.frame()
    }
    return(.Object)
})



 
