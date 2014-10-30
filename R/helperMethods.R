

estimateProportions2 <- function(x, method = "mode") {
    UseMethod("estimateProportions2")
}

#'@importFrom modeest mlv
estimateProportions2.MDMixResult <- function(x, method = "mode") {
    match.arg(method, c("mode", "mean"))
    posterior <- x$getP()
    if (method == "mean") {
        posterior.diff <- do.call(rbind, lapply(posterior, function(x) mean(apply(do.call(cbind,
            x)[, 2:1], 1, diff))))
        posterior.lods <- do.call(rbind, lapply(posterior, function(x) {
            d <- (apply(do.call(cbind, x)[, 2:1], 1, function(x) diff(log(x))))
            mean(d[is.finite(d)])
        }))
        ml <- cbind(pu = prop.table(as.matrix(x$n.unstim), margin = 1)[, 2], ps = prop.table(as.matrix(x$n.stim),
            margin = 1)[, 2])
        ml.differences <- apply(ml, 1, diff)
        ml.logodds <- apply(ml, 1, function(x) diff(log(x)))
        return(list(posterior.logodds = posterior.lods, posterior.differences = posterior.diff,
            ml.differences = ml.differences, ml.logodds = ml.logodds))
    }
    if (method == "mode") {
        posterior.diff <- do.call(rbind, lapply(posterior, function(x) mlv(apply(do.call(cbind,
            x)[, 2:1], 1, diff), method = "mfv")[[1]]))
        posterior.lods <- do.call(rbind, lapply(posterior, function(x) mlv(apply(do.call(cbind,
            x)[, 2:1], 1, function(x) diff(log(x))), method = "mfv")[[1]]))
        ml <- cbind(pu = prop.table(as.matrix(x$n.unstim), margin = 1)[, 2], ps = prop.table(as.matrix(x$n.stim),
            margin = 1)[, 2])
        ml.differences <- apply(ml, 1, diff)
        ml.logodds <- apply(ml, 1, function(x) diff(log(x)))
        return(list(posterior.logodds = posterior.lods, posterior.differences = posterior.diff,
            ml.differences = ml.differences, ml.logodds = ml.logodds))
    }
}


##' Arcsinh transform for ggplot2
##'
##' Arcsinh transform for use with coord_trans in ggplot2
##' @title asinh_trans
##' @param c \code{numeric} cofactor for asinh trasnform. Default 1.
##' @return transform
##' @import scales
##' @export
##' @author Greg Finak
asinh_trans <- function(c){
    trans_new(name = "asinh", transform = function(x) asinh(x * c), inverse = function(x) sinh(x)/c)
}

##' Combine two or more MIMOSAResultList objects
##'
##' Combines two or more MIMOSAResultList objects. The method is light on error checking so the results should be from the same
##' MIMOSAExpressionSet object.
##' @title Combine MIMOSAResultList objects
##' @param x \code{MIMOSAResultList}
##' @param y \code{MIMOSAResultList}
##' @param ... additional \code{MIMOSAResultList} objects
##' @return a \code{MIMOSAResultList}
##' @export
##' @name combine.MIMOSA
##' @author Greg Finak
combine.MIMOSA <- function(x,y,...){
    if(length(list(...)) > 0L){
        combine.MIMOSA(x,do.call(combine.MIMOSA,list(y,...)))
    }
    else{
        if(class(x)[1]!=class(y)[1]){
            stop("objects must be the same class but are '", class(x), "', '", class(y), "'")
        }
        # actual combine code
        x<-c(x,y)
        class(x) <- c("MIMOSAResultList","list")
        return(x)
    }
}


##' Boxplots of MIMOSA
##'
##' Generate boxplots for MIMOSA positivity calls.
##' @title boxplotMIMOSAResultList
##' @param data \code{MIMOSAResultList}
##' @param title \code{character} Title of the plot.
##' @param x_axis_category \code{name} the column of the phenoData frame for the x-axis of the boxplots.
##' @param cofactor \code{integer} cofactor of the arcsinhTransform for the y axis.
##' @param line \code{logical} whether or not to connect points from the same subject
##' @param threshold \code{numeric} the FDR threshold (q-value) at which to classify responders as a separate category.
##' @return \code{ggplot} object.
##' @importFrom data.table setnames
##' @export
##' @author Greg Finak
boxplotMIMOSAResultList <- function (data, title = "A Boxplot", x_axis_category = NULL,
    cofactor = 5000,line=TRUE,threshold=0.005)
{
    x_axis_category <- deparse(substitute(x_axis_category))
    if (x_axis_category == "x_axis_category") {
        stop("must provide an 'x_axis_category' variable for boxplots")
    }
    d <- (cbind(ldply(data, pData), fdr(data), countsTable(data,
        proportion = TRUE)))
    x_axis_category <- get(x_axis_category, d)
    setnames(d, colnames(d)[(ncol(d) - 3):ncol(d)], c("ParentProportion",
        "Proportion", "ParentProportion_REF", "Proportion_REF"))
    plt <- ggplot(environment = environment(), d, aes(x = x_axis_category,
        y = Proportion - Proportion_REF)) + geom_boxplot(aes(fill = fdr <
        threshold), outlier.colour = NA, position = "identity") +
        coord_trans(y = asinh_trans(cofactor)) + theme_bw() +
        facet_wrap(~.id) + ggtitle(Kmisc::wrap(title, 40)) +
        geom_jitter(aes(color = fdr < threshold), position = position_jitter(width = 0.01,
            height = 0))  +
        scale_fill_brewer(palette = "Pastel1") + scale_color_brewer(palette = "Set1")
    if(line){
        plt <- plt + geom_line(aes(group = PTID), lty = 2)
    }
    plt
}

.wrap <- function (x, width = 8, ...)
{
    return(paste(strwrap(x, width, ...), collapse = "\n"))
}

proportions.icsdata <- function(object) {
    r <- t(apply(object, 1, function(x) cbind(prop.table(x[c("Ns", "ns")])[2], prop.table(x[c("Nu",
        "nu")])[2])))
    colnames(r) <- c("ps", "pu")
    return(r)
}
proportions <- function(object) {
    UseMethod("proportions")
}

getpData <- function(x) {
    UseMethod("getpData")
}
setpData <- function(x, y) {
    UseMethod("setpData")
}

getpData.icsdata <- function(x) {
    pData(attr(x, "pData"))
}
setpData.icsdata <- function(x, y) {
    pData(attr(x, "pData")) <- y
    x
}

estimate_logZus <- function(alpha.u, beta.u, alpha.s, beta.s, B, lower.tail = TRUE) {
    if (lower.tail) {
        x <- rbeta(B, alpha.s, beta.s)
        I <- mean(pbeta(x, alpha.u, beta.u))
    } else if (!lower.tail) {
        x <- rbeta(B, alpha.u, beta.u)
        I <- mean(pbeta(x, alpha.s, beta.s, lower.tail = FALSE))
    }
    return(log(I))
}

refactorPData <- function(x) {
    if (class(x) == "BetaMixResult") {
        cat <- sapply(pData(x), is.factor)
        pData(x)[cat] <- lapply(pData(x)[cat], factor)
        x
    } else {
        stop("class must be BetaMixResult")
    }
}

setGeneric("Data", function(object) standardGeneric("Data"))
setMethod("Data", "BetaMixResult", function(object) object@data[, c("Ns", "ns", "Nu",
    "nu")])

huberFilter <- function(object, sd = 2) {
    if (any(class(object) == "icsdata")) {
        stim <- object[, c("Ns", "ns")]
        unstim <- object[, c("Nu", "nu")]
        difference <- apply(stim, 1, function(x) prop.table(x)[2]) - apply(unstim,
            1, function(x) prop.table(x)[2])
        h <- huber(difference)
        obj <- object[difference >= h$mu - 2 * h$s, ]
        a <- attributes(object)
        a$row.names <- attr(a, "row.names")
        pd <- attr(object, "pData")[difference >= h$mu - 2 * h$s, ]
        attr(obj, "pData") <- pd
        obj
    } else if (class(object) == "BetaMixResult") {
        huberFilter(object@data, sd = sd)
    } else {
        message("Other data classes not yet supported.. soon")
    }
}


setMethod("pData", "BetaMixResult", function(object) {
    pData(object@pd)
})

setMethod("pData<-", c("BetaMixResult", "data.frame"), function(object, value) {
    pData(object@pd) <- value
    object
})
# The mean,sample-size parameterization of the beta-binomial
f0 <- function(p, z, d, w, alternative = alternative, mciter = mciter) {
    -CompleteDataLLRcpp(d = d, alpha0 = p[3] * p[4], beta0 = (1 - p[3]) * p[4], alphaS = p[1] *
        p[2], betaS = (1 - p[1]) * p[2], z = z, w = w, alternative = alternative,
        mciter = mciter)
}
f0m <- function(p, z, d, w, alternative = alternative, mciter = mciter, alpha0, beta0) {
    -CompleteDataLLRcpp(d = d, alpha0 = alpha0, beta0 = beta0, alphaS = p[1] * p[2],
        betaS = (1 - p[1]) * p[2], z = z, w = w, alternative = alternative, mciter = mciter)
}

# estimate_logZus(alpha.u=10,beta.u=100,alpha.s=1,beta.s=100,B=1000000,lower.tail=FALSE)
# estimate_logZus(alpha.u=10,beta.u=100,alpha.s=1,beta.s=100,B=1000000,lower.tail=TRUE)

f1 <- function(x, alpha.u, beta.u, alpha.s, beta.s) {
    dbeta(x, alpha.u, beta.u) * pbeta(x, alpha.s, beta.s, lower.tail = FALSE)
}
f2 <- function(x, alpha.u, beta.u, alpha.s, beta.s) {
    dbeta(x, alpha.s, beta.s) * pbeta(x, alpha.u, beta.u, lower.tail = TRUE)
}

f1.log <- function(x, alpha.u, beta.u, alpha.s, beta.s) {
    exp(dbeta(x, alpha.u, beta.u, log = TRUE) + pbeta(x, alpha.s, beta.s, log.p = TRUE,
        lower.tail = FALSE))
}
f2.log <- function(x, alpha.u, beta.u, alpha.s, beta.s) {
    exp(dbeta(x, alpha.s, beta.s, log = TRUE) + pbeta(x, alpha.u, beta.u, log.p = TRUE,
        lower.tail = TRUE))
}

# integrate(f1,0,1,alpha.u=10,beta.u=100,alpha.s=1,beta.s=100)
# integrate(f2,0,1,alpha.u=10,beta.u=100,alpha.s=1,beta.s=100)


# calcSens<-function(S,delta,shrink){
# data.frame(do.call(rbind,apply(S,2,function(x){do.call(rbind,lapply(x,function(x){x<-t(x)[2:1,2:1];rs<-rowSums(x);cs<-colSums(x);data.frame(ppv=x[1,1]/rs[1],npv=x[2,2]/rs[2],sensitivity=x[1,1]/cs[1],specificity=x[2,2]/cs[2])}))})),test=rep(gl(2,1,labels=c('Fisher','BetaBinomial')),length(delta)),Threshold=rep(delta,each=2),shrinkage=rep(shrink,length(delta)*2))
# } summarize<-function(r=NULL,tpr=0.1){ thresholds=seq(0,1,by=0.001) d<-r@data;
# N<-nrow(d)
# fisher.p<-apply(d,1,function(x)fisher.test(matrix(x,2,byrow=T),alternative='greater')$p)
# #Table for ROC curve
# sensitivity<-rbind(sapply(thresholds,function(x)list(table(truth=factor(c(rep(FALSE,floor(N*(1-tpr))),rep(TRUE,N-floor(N*(1-tpr)))),levels=c('FALSE','TRUE')),fisher=factor(fisher.p<x,levels=c('FALSE','TRUE'))))),sapply(thresholds,function(x)list(table(truth=factor(c(rep(FALSE,floor(N*(1-tpr))),rep(TRUE,N-floor(N*(1-tpr)))),levels=c('FALSE','TRUE')),betabin=factor(r@z[,1]<x,levels=c('FALSE','TRUE'))))))
# sensitivity<-calcSens(sensitivity,thresholds,1) #False discovery rates
# fdr<-rbind(sapply(thresholds,function(x)list(table(truth=factor(c(rep(FALSE,floor(N*(1-tpr))),rep(TRUE,N-floor(N*(1-tpr)))),levels=c('FALSE','TRUE')),fisher=factor(p.adjust(fisher.p,'fdr')<x,levels=c('FALSE','TRUE'))))),sapply(thresholds,function(x)list(table(truth=factor(c(rep(FALSE,floor(N*(1-tpr))),rep(TRUE,N-floor(N*(1-tpr)))),levels=c('FALSE','TRUE')),betabin=factor(r@fdr<x,levels=c('FALSE','TRUE'))))))
# fdr<-calcSens(fdr,thresholds,1) #caculate the AUC from the Z's
# truth<-factor(c(rep(FALSE,N*(1-tpr)),rep(TRUE,N*(tpr))),levels=c('FALSE','TRUE'))
# sm.betabin<-0 sm.fisher<-0 for(i in which(truth=='TRUE')){ for(j in
# which(truth!='TRUE')){ sm.betabin<-sm.betabin+as.numeric(r@z[i,2]>r@z[j,2])
# ##less than since these are p-values
# sm.fisher<-sm.fisher+as.numeric(fisher.p[i]<fisher.p[j]) } }
# AUC.betabin<-sm.betabin/(length(which(truth=='TRUE'))*length(which(truth!='TRUE')))
# AUC.fisher<-sm.fisher/(length(which(truth=='TRUE'))*length(which(truth!='TRUE')))
# return(list(r,sensitivity,fdr,AUC.betabin,AUC.fisher)) }

simulate2 <- function(obs = 1000, As, Bs, A0, B0, NS, N0, w2, alternative = "greater",
    truncated = FALSE) {
    # ps and pu for null component
    match.arg(alternative, c("greater", "not equal"))

    null <- round((1 - w2) * obs)
    stim <- obs - round((1 - w2) * obs)
    p <- matrix(NA, ncol = 2, nrow = null + stim)
    i <- 1
    if (!truncated) {
        while (i <= null) {
            ps <- pu <- rbeta(1, A0, B0)
            p[i, ] <- c(ps, pu)
            i <- i + 1
        }
        while (i <= null + stim) {
            ps <- rbeta(1, As, Bs)
            pu <- rbeta(1, A0, B0)
            if (alternative == "greater" & ps > pu) {
                p[i, ] <- c(ps, pu)
                i <- i + 1
            }
            if (alternative == "not equal" & ps != pu) {
                p[i, ] <- c(ps, pu)
                i <- i + 1
            }
        }
    } else {
        # Simulate from truncated normal...
        mp0 <- A0/(A0 + B0)
        mps <- As/(As + Bs)
        v0 <- sqrt((A0 * B0)/((A0 + B0)^2 * (A0 + B0 + 1)))
        vs <- sqrt((As * Bs)/((As + Bs)^2 * (As + Bs + 1)))
        while (i <= null) {
            ps <- pu <- rnorm(1, mp0, v0)
            if (ps >= 0 & ps <= 1) {
                p[i, ] <- c(ps, pu)
                i <- i + 1
            }
        }
        while (i <= null + stim) {
            pu <- rnorm(1, mp0, v0)
            ps <- rnorm(1, mps, vs)
            if (ps > pu & ps >= 0 & pu >= 0 & pu <= 1 & ps <= 1) {
                p[i, ] <- c(ps, pu)
                i <- i + 1
            }
        }
    }
    colnames(p) <- c("ps", "pu")
    p <- data.frame(p)

    ### Now simulate the counts given the p's
    d <- data.frame(ns = rbinom(obs, NS, prob = p[, "ps"]), nu = rbinom(obs, N0,
        prob = p[, "pu"]))
    d <- data.frame(d, Ns = NS - d[, "ns"], Nu = N0 - d[, "nu"])
    attr(d, "control") <- "control"
    attr(d, "stimulation") <- "simulated data"
    attr(d, "cytokine") <- "simulated data"
    return(d)
}



plotPriors <- function(fit, a0 = NULL, b0 = NULL, as = NULL, bs = NULL, l = 1000,
    ...) {
    null.fit <- qbeta(c(0.001, 0.999), fit@alpha0, fit@beta0)
    stim.fit <- qbeta(c(0.001, 0.999), fit@alphaS, fit@betaS)
    if (is.null(as) | is.null(bs) | is.null(a0) | is.null(b0)) {
        lower <- min(c(null.fit, stim.fit))
        upper <- max(c(null.fit, stim.fit))
        skip = TRUE
    } else {
        true.null <- qbeta(c(0.001, 0.999), a0, b0)
        true.stim <- qbeta(c(0.001, 0.999), as, bs)
        lower <- min(c(null.fit, stim.fit, true.null, true.stim))
        upper <- max(c(null.fit, stim.fit, true.null, true.stim))
        skip = FALSE
    }

    if (fit@alternative.model == "greater") {
        s <- seq(lower, upper, l = 1000)
        null <- dbeta(s, fit@alpha0, fit@beta0) * pbeta(1 - s, fit@betaS, fit@alphaS)
        stim <- dbeta(s, fit@alphaS, fit@betaS) * pbeta(s, fit@alpha0, fit@beta0)
        plot(s, null, col = "green", lty = 1, type = "l")
        lines(s, stim, col = "red", lty = 1)
    } else {
        s <- seq(lower, upper, l = 1000)
        null <- dbeta(s, fit@alpha0, fit@beta0)
        stim <- dbeta(s, fit@alphaS, fit@betaS)
        plot(s, null, col = "green", lty = 1, type = "l")
        lines(s, stim, col = "red", lty = 1)
    }
    if (!skip) {
        s <- seq(lower, upper, l = 1000)
        if (fit@alternative.model == "greater") {
            null <- dbeta(s, a0, b0) * pbeta(1 - s, bs, as)
            stim <- dbeta(s, as, bs) * pbeta(s, a0, b0)
        } else {
            null <- dbeta(s, a0, b0)
            stim <- dbeta(s, as, bs)
        }
        lines(s, null, col = "green", lty = 2)
        lines(s, stim, col = "red", lty = 2)
        legend("topright", c("estimated null", "estimated stimulated", "true null",
            "true stimulated"), lty = c(1, 1, 2, 2), col = c("green", "red", "green",
            "red"))
    } else {
        legend("topright", c("estimated null", "estimated stimulated"), lty = c(1,
            1), col = c("green", "red"))
    }
    title(main = "Prior Distributions")

}

# Accessor and convenience methods setMethod('show','ICS',function(object){
# cat(length(levels(object@ID)),' Samples\n')
# print(head(as(object,'data.frame'))) cat('Antigens:
# ',levels(object@antigen),'\n') cat('Populations:
# ',levels(object@parent),'\n') cat('Cytokines: ',levels(object@fname),'\n')
# #Walk the nested list in .Data and summarize the levels
# head<-mode(object@.Data) this<-object@.Data levels<-NULL nlevels<-0;
# while(head=='list'){ nlevels<-nlevels+1 levels<-c(levels,list(names(this)))
# this<-this[[1]] head<-mode(this) } #Omit reporting the last level corresponding
# to patients cat(nlevels-1,' Levels of nesting\n') for(i in 1:(nlevels-1)){
# cat(rep(' ',i),i,'. ',levels[[i]],'\n') } }) #convert ICS to data.frame
# setAs('ICS','data.frame',def=function(from){
# df<-data.frame(pos=from@pos,neg=from@neg,fname=from@fname,fcsfile=from@fcsfile,parent=from@parent,antigen=from@antigen,ID=from@ID,from@rest)
# if(ncol(from@rest)>0){ colnames(df)[8:(7+ncol(from@rest))]<-colnames(from@rest)
# } return(df) })

# prepCounts populates one of the slots in the ICS structure.
# setGeneric('prepCounts',function(ics,...){ standardGeneric('prepCounts') })

fa <- function(u, M) {
    u * M
}
fb <- function(u, M) {
    M * (1 - u)
}
fu <- function(a, b) {
    a/(a + b)
}
fm <- function(a, b) {
    a + b
}

# setMethod('prepCounts',signature='ICS',function(ics,extra.identifiers=NULL){
# message('Hold on.. melting and casting.') x<-as(ics,'data.frame')
# if(!exists('extra.identifiers')) extra.identifiers<-NULL
# x<-melt(x,id=c('ID','antigen','parent','fcsfile','fname',extra.identifiers),measure=c('pos','neg'))
# form<- paste(paste('antigen~variable|', paste(paste(extra.identifiers,sep = '
# ', collapse = '+'), ifelse(is.null(extra.identifiers),'', '+'))),'parent +
# fname + ID') x<-cast(x,as.formula(form),fun.aggregate=sum) return(x) }) #show
# method for BetaMixResult
setMethod("show", "BetaMixResult", function(object) {
    cat("BetaMixResult for contrast: ", object@stimulation, "-", object@control,
        "\n")
    cat("Number of Samples: ", nrow(object@data), "\n")
    cat("Log Likelihood: ", object@ll, "\n")
    cat("cytokine", object@cytokine, "\n")
    cat("w: Pr(ps=pu)=", object@w[1], ifelse(object@alternative.model == "greater",
        " Pr(ps>pu)=", " Pr(ps!=pu)="), object@w[2], "\n")
})

setMethod("summary", "BetaMixResult", function(object, ...) {
    show(object)
    threshold <- list(...)$threshold
    if (is.null(threshold)) {
        threshold <- 0.01
    }
    cat("Positivity: ", table(factor(object@fdr <= threshold, levels = c("FALSE",
        "TRUE")))[2], " of ", length(object@fdr), " at ", threshold * 100, "% FDR\n")
})
# plot method for BetaMixResult, will plot the log likelihood trajectory of the
# model fitting process
plot.BetaMixResult <- function(x, y, ...) {
    plot(x@traj[-1], main = "Complete Data Log Likelihood Trajectory", xlab = "Iteration",
        ylab = "Complete Data Log Likelihood", type = "o")
}




