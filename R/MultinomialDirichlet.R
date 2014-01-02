# TODO: Implements the multinomial-dirichelt model for two
# cytokines Author: Greg Finak



# Observed data n.stim - vector of counts from the stimulated
# sample n.unstim - vector of counts from the unstimulated sample
# alpha.unstim alpha.stim


lkbeta <- function(alpha) {
  sum(lgamma(alpha)) - lgamma(sum(alpha))
}


makeLogLikeNULLComponent <- function(data.stim, data.unstim) {
  data <- data.stim + data.unstim
  ll <- function(x) {
    a <- x[(length(x)/2 + 1):length(x)]
    apply(data, 1, function(y) lkbeta(y + a)) - lkbeta(a)  #-
    # rowSums(lfactorial(data.stim))-rowSums(lfactorial(data.unstim))+lfactorial(rowSums(data.stim))+lfactorial(rowSums(data.unstim))
  }
}
makeLogLikeRespComponent <- function(data.stim, data.unstim) {
  ll <- function(x) {
    a <- x[(length(x)/2 + 1):length(x)]
    b <- x[1:(length(x)/2)]
    apply(data.stim, 1, function(y) lkbeta(y + b)) + apply(data.unstim, 
                                                           1, function(y) lkbeta(y + a)) - lkbeta(b) - lkbeta(a)  #-
    # rowSums(lfactorial(data.stim))-rowSums(lfactorial(data.unstim))+lfactorial(rowSums(data.stim))+lfactorial(rowSums(data.unstim))
  }
}



makeGradientNULLComponent <- function(data.stim, data.unstim, z = NULL) {
  data <- data.stim + data.unstim
  if (is.null(z)) {
    z <- matrix(1, nrow = nrow(data.unstim), ncol = 2)
  }
  grad <- function(x) {
    x <- x[(ncol(data.stim) + 1):(ncol(data.stim) * 2)]
    gr <- (t(digamma(sum(x)) - digamma(colSums((t(data) + x))) + 
               digamma(t(x + t(data)))) - digamma(x))
    gr <- gr %*% z[, 1]
    return(c(rep(0, ncol(data.stim)), gr))
  }
}
makeGradientHalfComponent <- function(data, z = NULL) {
  if (is.null(z)) {
    z <- matrix(1, nrow = nrow(data), ncol = 1)
  }
  grad <- function(x) {
    gr <- (t(digamma(sum(x)) - digamma(colSums((t(data) + x))) + 
               (digamma(t(x + t(data))))) - digamma(x))
    gr %*% as.vector(z)
  }
}

makeGradientRespComponent <- function(data.stim, data.unstim, z = NULL) {
  stim.ind <- 1:ncol(data.stim)
  unstim.ind <- (ncol(data.stim) + 1):(ncol(data.stim) * 2)
  if (is.null(z)) {
    z <- matrix(1, nrow = nrow(data.stim), ncol = 2)
  }
  gs <- makeGradientHalfComponent(data.stim, z = z[, 2])
  gu <- makeGradientHalfComponent(data.unstim, z = z[, 2])
  grad <- function(x) {
    c(gs(x[stim.ind]), gu(x[unstim.ind]))
  }
  return(grad)
}

makeHessianNULLComponent <- function(data.stim, data.unstim, z = NULL) {
  data <- data.stim + data.unstim
  if (is.null(z)) {
    z <- matrix(1, nrow = nrow(data.unstim), ncol = 2)
  }
  hess <- function(x) {
    x <- x[(ncol(data.stim) + 1):(ncol(data.stim) * 2)]
    H <- (trigamma(sum(x)) - trigamma(colSums(t(data) + x)))
    D <- (trigamma(x + t(data)) - trigamma(x))
    HH <- matrix(0, ncol(data) * 2, ncol(data) * 2)
    H <- sapply(1:length(H), function(i) {
      D[, i] <- D[, i] + H[i]
      M <- matrix(0, ncol(data) * 2, ncol(data) * 2)
      inds <- (ncol(data.stim) + 1):(ncol(data.stim) * 2)
      M[inds, inds] <- H[i]
      diag(M) <- c(rep(0, ncol(data)), D[, i])
      list(M * z[i, 1])
    })
    for (i in seq_along(H)) {
      HH <- HH + H[[i]]
    }
    return(HH)
  }
  return(hess)
}

makeHessianHalfComponent <- function(data, z = NULL) {
  if (is.null(z)) {
    z <- matrix(1, nrow = nrow(data), ncol = 1)
  }
  hess <- function(x) {
    
    H <- (trigamma(sum(x)) - trigamma(colSums(t(data) + x)))
    D <- (trigamma(x + t(data)) - trigamma(x))
    HH <- matrix(0, ncol(data), ncol(data))
    H <- sapply(1:length(H), function(i) {
      D[, i] <- D[, i] + H[i]
      M <- matrix(H[i], ncol(data), ncol(data))
      diag(M) <- D[, i]
      list(M * z[i])
    })
    for (i in seq_along(H)) {
      HH <- HH + H[[i]]
    }
    return(HH)
  }
}

makeHessianRespComponent <- function(data.stim, data.unstim, z = NULL) {
  stim.ind <- 1:ncol(data.stim)
  unstim.ind <- (ncol(data.stim) + 1):(ncol(data.stim) * 2)
  if (is.null(z)) {
    z <- matrix(1, nrow = nrow(data.stim), ncol = 2)
  }
  hs <- makeHessianHalfComponent(data.stim, z[, 2])
  hu <- makeHessianHalfComponent(data.unstim, z[, 2])
  hess <- function(x) {
    H <- matrix(0, nrow = length(c(stim.ind, unstim.ind)), ncol = length(c(stim.ind, 
                                                                           unstim.ind)))
    H[stim.ind, stim.ind] <- hs(x[stim.ind])
    H[unstim.ind, unstim.ind] <- hu(x[unstim.ind])
    return(H)
  }
  return(hess)
}

#'@importFrom MCMCpack rdirichlet
simMD <- function(alpha.s = c(100, 50, 10, 10), alpha.u = c(100, 10, 
                                                            10, 10), N = 100, w = 0.5, n = 2, alternative = "greater") {
  nnull <- round((1 - w) * N)
  nresp <- N - round((1 - w) * N)
  pu <- rdirichlet(nnull, alpha.u)
  if (nnull > 0) {
    ps <- pu
  } else {
    # empty matrix when all responders
    ps <- matrix(ncol = length(alpha.s), nrow = 0)
  }
  i <- nnull + 1
  ps <- rbind(ps, matrix(0, nrow = nresp, ncol = length(alpha.s)))
  pu <- rbind(pu, matrix(0, nrow = nresp, ncol = length(alpha.u)))
  if (alternative == "greater") {
    while (i <= nnull + nresp) {
      p.s <- rdirichlet(1, alpha.s)
      p.u <- rdirichlet(1, alpha.u)
      if (any(p.s[-1L] > p.u[-1L])) {
        ps[i, ] <- p.s
        pu[i, ] <- p.u
        i <- i + 1
      }
    }
  } else {
    p.s <- rdirichlet(nresp, alpha.s)
    p.u <- rdirichlet(nresp, alpha.u)
    i <- (nnull + 1):(nnull + nresp)
    ps[i, ] <- p.s
    pu[i, ] <- p.u
  }
  NU <- runif(N, 1.5 * 10^n, 1.5 * 10^n)
  NS <- runif(N, 1.5 * 10^n, 1.5 * 10^n)
  nu <- t(sapply(seq_along(1:N), function(i) rmultinom(1, NU[i], 
                                                       pu[i, ])))
  ns <- t(sapply(seq_along(1:N), function(i) rmultinom(1, NS[i], 
                                                       ps[i, ])))
  data <- list(n.stim = ns, n.unstim = nu)
  return(data)
}



#' EM fitting of the Multinomial Dirichlet MIMOSA model.
#' 
#' This function fits the multinomial dirichelt MIMOSA model using EM. It can also be used to initialize the model parameters for the MCMC model.
#' @param data The observed data
#' @param modelmatrix a model matrix specifying which components should be computed
#' @param alternative either 'greater' or 'not equal' to fit the one-sided or two-sided model.
#' @param initonly \code{TRUE} or \code{FALSE} to return just the initialization parameters.
#' @return An object of class \code{MDMixResult}
#' @author Greg Finak
#' TODO filtering of pu>ps needs to be corrected here.
MDMix <- function(data = NULL, modelmatrix = NULL, alternative = "greater", 
                  initonly = FALSE) {
  data <- icsdata2mvicsdata(data)
  unstim <- data$n.unstim
  stim <- data$n.stim
  match.arg(alternative, c("greater", "not equal"))
  if (alternative == "greater") {
    alternative <- "greater"
    filt <- apply(stim, 1, function(x) prop.table(x)[2]) < apply(unstim, 
                                                                 1, function(x) prop.table(x)[2])
  } else if (alternative == "not equal") {
    alternative <- "two.sided"
    filt <- rep(FALSE, nrow(stim))
  }
  # fisher's exact test of all the marginals
  if (ncol(unstim) > 2) {
    # TODO double check that Fisher's gives the right result here for
    # one sided test
    mm <- do.call(cbind, lapply(2:ncol(unstim), function(i) apply(cbind(rowSums(unstim[, 
                                                                                       -i]), unstim[, i], rowSums(stim[, -i]), stim[, i]), 1, 
                                                                  function(x) fisher.test(matrix(x, 2), alternative = alternative)$p.value)))
    mm <- matrix(p.adjust(mm, "fdr") < 0.05, ncol = ncol(unstim) - 
                   1)
    # returns all non-significant
    mm <- apply(mm, 1, function(x) all(!x))
  } else {
    # two-D case
    mm <- sapply(1:nrow(unstim), function(i) fisher.test(matrix(unlist(c(unstim[i, 
                                                                                ], stim[i, ])), ncol = 2, byrow = TRUE), alternative = alternative)$p.value)
    mm <- p.adjust(mm, "fdr") < 0.05
    mm <- !mm
  }
  # observations with no significant marginals belong to the null
  # component.
  
  # The rest are from the responder component construct the z-matrix
  z <- matrix(0, length(mm), 2)
  z[mm, 1] <- 1
  z[!mm, 2] <- 1
  
  # estimate hyperparamters.
  pu <- prop.table(colMeans(unstim))
  ps <- prop.table(colMeans(stim[which(z[, 2] == 1), , drop = FALSE]))
  alpha.u <- round(colMeans(unstim))
  alpha.s <- round(colMeans(stim[which(z[, 2] == 1), , drop = FALSE]))
  
  
  if (any(is.nan(alpha.s))) 
    alpha.s[is.nan(alpha.s)] <- 1
  if (any(is.nan(alpha.u))) 
    alpha.u[is.nan(alpha.u)] <- 1
  
  
  alpha.s[alpha.s == 0] <- 1
  alpha.u[alpha.u == 0] <- 1
  
  guess <- c(ps, pu)
  
  w <- colSums(z)/sum(z)
  if (initonly) {
    return(list(q = w[1], z = z, alpha.s = alpha.s, alpha.u = alpha.u))
  }
  # EM
  LL <- NULL
  repeat {
    if (is.null(LL)) {
      last <- .Machine$double.xmax
    }
    # update parameters
    llnull <- makeLogLikeNULLComponent(stim, unstim)
    llresp <- makeLogLikeRespComponent(stim, unstim)
    gnull <- makeGradientNULLComponent(stim, unstim, z)
    gresp <- makeGradientRespComponent(stim, unstim, z)
    hessnull <- makeHessianNULLComponent(stim, unstim, z)
    hessresp <- makeHessianRespComponent(stim, unstim, z)
    
    
    
    iter <- 2
    ll <- rep(0, 1000)
    # make negative
    ll[1] <- .Machine$double.xmax
    lastguess <- guess
    repeat {
      t <- try(solve(hessresp(guess) + hessnull(guess), gnull(guess) + 
                       gresp(guess)), silent = TRUE)
      if (inherits(t, "try-error")) 
        t <- try(ginv(hessresp(guess) + hessnull(guess)) %*% 
                   (gnull(guess) + gresp(guess)),silent=TRUE)  #uses SVD
      if (inherits(t,"try-error")){
          H<-hessresp(guess)+hessnull(guess)
          G<-gresp(guess)+gnull(guess)
          H<-H+G%*%t(G)
          t<-try(solve(H)%*%G,silent=TRUE)
      }
      if (inherits(t, "try-error")) {
        stop("Hessian not positive definite. Trying MCMC")
      }
      new <- guess - t
      ll[iter] <- -sum(llnull(new) * z[, 1] + llresp(new) * 
                         z[, 2])
      step <- 0.5
      # put a limit on the number of iterations here
      piter <- 1
      while (((ll[iter - 1] - ll[iter]) < -1e-04) & piter < 
               10) {
        step <- step/2
        new <- 0.5 * (guess - t)
        ll[iter] <- -sum(llnull(new) * z[, 1] + llresp(new) * 
                           z[, 2])
        piter <- piter + 1
      }
      if ((all(abs(new - guess)/abs(guess) < 1e-04)) | (iter > 
                                                          999)) {
        guess <- new
        if (any(guess < 0)) {
          warning("Parameter estimates became negative!")
        }
        break
      }
      guess <- new
      iter <- iter + 1
    }
    
    
    den <- apply(cbind(log(w[1]) + llnull(new), log(w[2]) + llresp(new)), 
                 1, function(x) log(sum(exp(x - max(x)))) + max(x))
    z2 <- exp((llresp(new) + log(w[2])) - (den))
    
    if (any(filt) & alternative == "greater") {
      ## Fix z's for pu>ps when alternative is greater
      z2[filt] <- 0
    }
    z <- cbind(1 - z2, z2)
    w <- colSums(z)/sum(z)
    cll <- -sum(sapply((llnull(new) + log(w[1])) * z[, 1], function(x) ifelse(is.nan(x), 
                                                                              0, x)) + sapply((llresp(new) + log(w[2])) * z[, 2], function(x) ifelse(is.nan(x), 
                                                                                                                                                     0, x)))
    if ((abs((last - cll)) < 0.01) & cll >= last) {
      break
    }
    
    LL <- c(LL, cll)
    # cat(cll,'\n')
    last <- cll
  }
  gnull <- makeGradientNULLComponent(stim, unstim, z)
  gresp <- makeGradientRespComponent(stim, unstim, z)
  hessnull <- makeHessianNULLComponent(stim, unstim, z)
  hessresp <- makeHessianRespComponent(stim, unstim, z)
  if (any(new < 0)) {
    warning("Failed to converge: negative parameter estimates")
    m <- "Failed to converge: negative parameter estimates"
    class(m) <- "try-error"
    return(m)
  }
  return(new("MDMixResult", llnull = llnull, llresp = llresp, gresp = gresp, 
             hresp = hessresp, gnull = gnull, w = w, hnull = hessnull, 
             z = z, ll = LL, par.stim = new[1:(length(new)/2)], par.unstim = new[(length(new)/2 + 
                                                                                    1):length(new)], data = data))
}

# #' extracts bifunctional cytokine data from and ICS object given
# the two marginals (A, B) and A||B for stimualted and
# unstimulated. Used for the multinomial dirichlet model. ORDER OF
# CYTOKINES MATTERS!  #' @param ics #' @param cytokineA #' @param
# cytokineB #' @param or #' @param stim #' @param control #'
# @param subset #' @param shrink #' @param scl #' @returnType #'
# @return #' @author Greg Finak #' @export
# extractDataMultinomDir<-function(ics=NULL,cytokineA=NULL,cytokineB=NULL,or=NULL,stim=NULL,control=NULL,subset=NULL,shrink=1,scl=1){
# if(any(c(is.null(ics),is.null(cytokineA),is.null(stim),is.null(cytokineB),is.null(or),is.null(control),is.null(subset)))){
# stop('All arguments must be non--null'); }
# A<-extractData(ics,control=control,stim=stim,subset=c(subset,cytokineA))
# B<-extractData(ics,control=control,stim=stim,subset=c(subset,cytokineB))
# OR<-extractData(ics,control=control,stim=stim,subset=c(subset,or))
# AND<-A+B-OR AND[,'Ns']<-A[,'Ns']+A[,'ns']-AND[,'ns']
# AND[,'Nu']<-A[,'Nu']+A[,'nu']-AND[,'nu'] ds<-AND[,'ns']
# bs<-A[,'ns']-AND[,'ns'] cs<-B[,'ns']-ds as<-A[,'Ns']-cs
# du<-AND[,'nu'] bu<-A[,'nu']-AND[,'nu'] cu<-B[,'nu']-du
# au<-A[,'Nu']-cu n.unstim<-cbind('u1'=au,'u2'=bu,'u3'=cu,'u4'=du)
# n.stim<-cbind('s1'=as,'s2'=bs,'s3'=cs,'s4'=ds) ##rownames
# rownames(n.unstim)<-rownames(A) rownames(n.stim)<-rownames(A)
# r<-(list(n.stim,n.unstim)) names(r)<-c('n.stim','n.unstim')
# attr(r,'cytokines')<-c(cytokineA,cytokineB)
# attr(r,'subset')<-subset attr(r,'stim')<-stim
# class(r)<-c('list','MDlist') return(r) }


# setOldClass('MDlist') setMethod(show,'MDlist',function(object){
# cat('Cytokines ',attr(object,'cytokines'),'\n')
# cat('Stimulation ',attr(object,'stim'),'\n') cat('Subset ',
# attr(object,'subset'),'\n') cat('Number of obs:
# ',nrow(object[[1]]),'\n') }) print.MDlist<-function(x,...){
# cat('Cytokines ',attr(x,'cytokines'),'\n') cat('Stimulation
# ',attr(x,'stim'),'\n') cat('Subset ', attr(x,'subset'),'\n')
# cat(nrow(x[[1]]),' observations','\n') } 
