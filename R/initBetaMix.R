
initBetaMix <- function(data = NULL, fixedNULL = FALSE, ics = NULL, 
                        alternative = "greater", priorXi = 1, scl = 10, K = 200, mciter = 200) {
  K <- 200
  if (priorXi < 1) {
    stop("The prior for the mixing proportions must be greater than or equal to 1.")
  }
  if (fixedNULL & !is.null(ics)) {
    # Estimate alpha0,beta0 from the full data set.
    d <- ics@.Data
    for (i in attr(data, "cytokine")) {
      d <- d[[i]]
    }
    d <- subset(do.call(rbind, d), quote(antigen == attr(data, 
                                                         "control")))
    
    # Moment estimators
    mu0 <- sum(d[, "pos"])/sum(d[, "pos"] + d[, "neg"])
    V <- var((d[, "pos"])/sum(d[, "pos"] + d[, "neg"]))
    alpha0 <- mu0 * (mu0 * (1 - mu0)/V - 1)
    beta0 <- (1 - mu0) * (mu0 * (1 - mu0)/V - 1)
    if (is.na(alpha0) | is.na(beta0)) {
      alpha0 <- beta0 <- mu0 * 1000
      beta0 <- 1000 - alpha0
    }
    pars <- try(optim(par = (c(alpha0, beta0)), function(p) {
      -sum(lchoose(d[, "pos"] + d[, "neg"], d[, "pos"]) - lbeta((p[1]), 
                                                                (p[2])) + lbeta(d[, "pos"] + (p[1]), d[, "neg"] + 
                                                                                  (p[2])))
    }, method = "L-BFGS-B", lower = c(0.1, 0.1), upper = c(Inf, 
                                                           Inf)))
    if (inherits(pars, "try-error") | try(pars$convergence, silent = TRUE) != 
          0) {
      # Moment estimators if the above fails
      mu0 <- sum(d[, "pos"])/sum(d[, "pos"] + d[, "neg"])
      V <- var((d[, "pos"])/sum(d["pos"] + d[, "neg"]))
      alpha0 <- mu0 * (mu0 * (1 - mu0)/V - 1)
      beta0 <- (1 - mu0) * (mu0 * (1 - mu0)/V - 1)
      if (is.na(alpha0) | is.na(beta0)) 
        alpha0 <- beta0 <- mu0 * 1000
      beta0 <- 1000 - alpha0
    } else {
      alpha0 <- (pars$par[1])
      beta0 <- (pars$par[2])
    }
  } else if (fixedNULL & is.null(ics)) {
    d <- as.data.frame(data[, c("ns", "nu", "Ns", "Nu")])
    V <- var(d[, "nu"]/(d[, "nu"] + d[, "Nu"]))
    mu0 <- mean(d[, "nu"]/(d[, "nu"] + d[, "Nu"]))
    alpha0 <- mu0 * (mu0 * (1 - mu0)/V - 1)
    beta0 <- (1 - mu0) * (mu0 * (1 - mu0)/V - 1)
    if (is.na(alpha0) | is.na(beta0)) {
      alpha0 <- beta0 <- mu0 * 1000
      beta0 <- 1000 - alpha0
    }
    pars <- try(optim(par = c(alpha0, beta0), function(p) {
      -sum(lchoose(d[, "nu"] + d[, "Nu"], d[, "nu"]) - lbeta(p[1], 
                                                             p[2]) + lbeta(d[, "nu"] + p[1], d[, "Nu"] + p[2]))
    }, method = "L-BFGS-B", lower = c(0.1, 0.1), upper = c(Inf, 
                                                           Inf)))
    if (inherits(pars, "try-error") | try(pars$convergence, silent = TRUE) != 
          0) {
      V <- var(d[, "nu"]/(d[, "nu"] + d[, "Nu"]))
      mu0 <- mean(d[, "nu"]/(d[, "nu"] + d[, "Nu"]))
      alpha0 <- mu0 * (mu0 * (1 - mu0)/V - 1)
      beta0 <- (1 - mu0) * (mu0 * (1 - mu0)/V - 1)
      if (is.na(alpha0) | is.na(beta0)) {
        alpha0 <- beta0 <- mu0 * 1000
        beta0 <- 1000 - alpha0
      }
    } else {
      alpha0 <- pars$par[1]
      beta0 <- pars$par[2]
    }
  }
  
  
  if (!fixedNULL) {
    d <- as.data.frame(data[, c("ns", "nu", "Ns", "Nu")])
    V <- var(d[, "nu"]/(d[, "nu"] + d[, "Nu"]))
    mu0 <- mean(d[, "nu"]/(d[, "nu"] + d[, "Nu"]))
    alpha0 <- mu0 * (mu0 * (1 - mu0)/V - 1)
    beta0 <- (1 - mu0) * (mu0 * (1 - mu0)/V - 1)
    if (is.na(alpha0) | is.na(beta0)) {
      alpha0 <- beta0 <- mu0 * 1000
      beta0 <- 1000 - alpha0
    }
    pars <- try(optim(par = (c(alpha0, beta0)), function(p) {
      -sum(lchoose(d[, "nu"] + d[, "Nu"], d[, "nu"]) - lbeta((p[1]), 
                                                             (p[2])) + lbeta(d[, "nu"] + (p[1]), d[, "Nu"] + (p[2])))
    }, method = "L-BFGS-B", lower = c(0.1, 0.1), upper = c(Inf, 
                                                           Inf)), silent = TRUE)
    if (inherits(pars, "try-error") | try(pars$convergence, silent = TRUE) != 
          0) {
      pars <- try(optim(par = log(c(alpha0, beta0)), function(p) {
        -sum(lchoose(d[, "nu"] + d[, "Nu"], d[, "nu"]) - lbeta(p[1], 
                                                               p[2]) + lbeta(d[, "nu"] + exp(p[1]), d[, "Nu"] + 
                                                                               exp(p[2])))
      }, method = "BFGS"), silent = TRUE)
      if (inherits(pars, "try-error") | try(pars$convergence, 
                                            silent = TRUE) != 0) {
        V <- var(d[, "nu"]/(d[, "nu"] + d[, "Nu"]))
        mu0 <- mean(d[, "nu"]/(d[, "nu"] + d[, "Nu"]))
        alpha0 <- mu0 * (mu0 * (1 - mu0)/V - 1)
        beta0 <- (1 - mu0) * (mu0 * (1 - mu0)/V - 1)
        if (is.na(alpha0) | is.na(beta0)) {
          alpha0 <- beta0 <- mu0 * 1000
          beta0 <- 1000 - alpha0
        }
      } else {
        alpha0 <- exp(pars$par[1])
        beta0 <- exp(pars$par[2])
      }
    } else {
      alpha0 <- (pars$par[1])
      beta0 <- (pars$par[2])
    }
  }
  
  # filter out anything that has ns and nu == 0 now estimate
  # responding, stimulated samples
  if (alternative == "greater") 
    alt <- "greater"
  if (alternative == "not equal") 
    alt <- "two.sided"
  d <- as.data.frame(data[, c("ns", "nu", "Ns", "Nu")])
  fisher.p <- apply(d, 1, function(x) fisher.test(matrix(unlist(x), 
                                                         2, byrow = T), alternative = alt)$p)
  fisher.p.w <- p.adjust(fisher.p, "fdr") < 0.05
  ord <- order(fisher.p, decreasing = F)
  l <- length(which(fisher.p.w[ord]))
  
  ## Odd situation when all (or many) observations are significant.
  ## Taking only the top 10 for the null is fine, but initialization
  ## fails when computing the complete LL for the data since the
  ## observations assigned to the null are so inconsistent with the
  ## model that optim fails.  The workaround, for now, is to use ALL
  ## observations when ALL are significant.
  if (l < length(fisher.p)) 
    l <- min(l, 10)
  
  if (l > 1) 
  {
    wh <- ord[1:l]
    # We'll estimate the w's from the response rate
    w <- c(1 - (sum(fisher.p.w) + priorXi - 1)/(length(fisher.p) + 
                                                  priorXi - 1), (sum(fisher.p.w) + priorXi - 1)/(length(fisher.p) + 
                                                                                                   priorXi - 1))
    
    # Set the z's to the responders / non-responders called by
    # fisher's test
    z <- matrix(0, nrow = length(fisher.p), ncol = 2)
    z[wh, 2] <- 1
    z[setdiff(1:length(fisher.p), wh), 1] <- 1  #90% in the null component\t\t
    
    
    if (alternative == "greater") {
      V <- var(d[wh, "ns"]/(d[wh, "ns"] + d[wh, "Ns"]))
      muS <- mean(d[wh, "ns"]/(d[wh, "ns"] + d[wh, "Ns"]))
      alphaS <- muS * (muS * (1 - muS)/V - 1)
      betaS <- (1 - muS) * (muS * (1 - muS)/V - 1)
      if (is.na(alphaS) | is.na(betaS)) {
        alphaS <- alpha0 * 2
        betaS <- beta0 + alpha0 - alphaS  #double the mean, keep 'sample size' fixed
      }
      # f0 is the mean,sample size parameterization, both components f0m
      # is for the stimulated component only f0 is defined in
      # helperFunctions.R
      if (fixedNULL) {
        # upper bounds k/mu ensure sufficient sample size
        pars <- try(optim(par = (c(alphaS/(alphaS + betaS), 
                                   betaS + alphaS)), function(p = par, data = d, 
                                                              Z = z, W = w, ALT = alternative, MC = mciter, 
                                                              a0 = alpha0, b0 = beta0) f0m(p = p, d = data, 
                                                                                           z = Z, w = W, alternative = ALT, mciter = MC, 
                                                                                           alpha0 = a0, beta0 = b0), method = "L-BFGS-B", 
                          lower = c(1e-06, 10), control = list(parscale = c(scl, 
                                                                            1)), upper = c(0.9999, K/muS)), silent = TRUE)
      } else {
        pars <- try(optim(par = (c(alphaS/(alphaS + betaS), 
                                   betaS + alphaS, alpha0/(alpha0 + beta0), beta0 + 
                                     alpha0)), function(p = par, data = d, Z = z, 
                                                        W = w, ALT = alternative, MC = mciter) f0(p = p, 
                                                                                                  d = data, z = Z, w = W, alternative = ALT, mciter = MC), 
                          method = "L-BFGS-B", lower = c(1e-06, 10, 1e-06, 
                                                         10), control = list(parscale = c(scl, 1, scl, 
                                                                                          1)), upper = c(0.9999, K/mu0, 0.9999, K/muS)), 
                    silent = TRUE)
      }
      if (inherits(pars, "try-error") | inherits(try(pars$convergence, 
                                                     silent = TRUE) != 0, "try-error")) {
        stop("failed to converge estimating initial alphaS,betaS in initBetaMix")
      }
      alphaS <- pars$par[1] * pars$par[2]
      betaS <- (1 - pars$par[1]) * pars$par[2]
      if (!fixedNULL) {
        alpha0 <- pars$par[3] * pars$par[4]
        beta0 <- (1 - pars$par[3]) * pars$par[4]
      }
    } else if (alternative == "not equal") {
      V <- var(d[wh, "ns"]/(d[wh, "ns"] + d[wh, "Ns"]))
      muS <- mean(d[wh, "ns"]/(d[wh, "ns"] + d[wh, "Ns"]))
      alphaS <- muS * (muS * (1 - muS)/V - 1)
      betaS <- (1 - muS) * (muS * (1 - muS)/V - 1)
      if (is.na(alphaS) | is.na(betaS)) {
        alphaS <- (1 - alpha0/(alpha0 + beta0)) * ((alpha0 + 
                                                      beta0)/2)  #same mean, but half the sample size (increased variance)
        betaS <- ((alpha0 + beta0)/2) - alphaS
      }
      
      if (fixedNULL) {
        pars <- try(optim(par = (c(alphaS/(alphaS + betaS), 
                                   betaS + alphaS)), function(p = par, data = d, 
                                                              Z = z, W = w, ALT = alternative, MC = mciter, 
                                                              a0 = alpha0, b0 = beta0) f0m(p = p, d = data, 
                                                                                           z = Z, w = W, alternative = ALT, mciter = MC, 
                                                                                           alpha0 = a0, beta0 = b0), method = "L-BFGS-B", 
                          control = list(parscale = c(scl, 1)), lower = c(1e-06, 
                                                                          10), upper = c(0.9999, K/muS)), silent = TRUE)
      } else {
        pars <- try(optim(par = (c(alphaS/(alphaS + betaS), 
                                   betaS + alphaS, alpha0/(alpha0 + beta0), beta0 + 
                                     alpha0)), function(p = par, data = d, Z = z, 
                                                        W = w, ALT = alternative, MC = mciter) f0(p = p, 
                                                                                                  d = data, z = Z, w = W, alternative = ALT, mciter = MC), 
                          control = list(parscale = c(scl, 1, scl, 1)), 
                          method = "L-BFGS-B", lower = c(1e-06, 10, 1e-06, 
                                                         10), upper = c(0.9999, K/muS, 0.9999, K/mu0)), 
                    silent = TRUE)
      }
      if (inherits(pars, "try-error") | inherits(try(pars$convergence != 
                                                       0), "try-error")) {
        stop("failed to converge estimating initial alpha0,beta0 in initBetaMix")
      }
      alphaS <- pars$par[1] * pars$par[2]
      betaS <- (1 - pars$par[1]) * pars$par[2]
      if (!fixedNULL) {
        alpha0 <- pars$par[3] * pars$par[4]
        beta0 <- (1 - pars$par[3]) * pars$par[4]
      }
    }
    
  }  ##Case where no samples are significant (ie fisher's exact test has no significant p-values)
  else if (alternative == "greater") {
    alphaS <- (1 - (alpha0/(alpha0 + beta0))) * 2  #twice the mean
    betaS <- beta0 + alpha0 - alphaS  #same sample size
    # TODO figure out the right thing to do here.. initialize z and w
    z <- matrix(0, nrow = nrow(d), ncol = 2)
    z[, 1] <- 1
    z[, 2] <- 0
    w <- colSums(z)/sum(z)
    muS <- mu0
  } else if (alternative == "not equal") {
    ## non-informative? TODO also figure out the right initialization
    ## here.
    alphaS <- 1.01
    betaS <- 1.01
    muS <- mu0
    # initialize z and w
    z <- matrix(0, nrow = nrow(d), ncol = 2)
    z[, 1] <- 1
    z[, 2] <- 0
    w <- colSums(z)/sum(z)
  }
  if (alphaS < 0 | betaS < 0) {
    stop("cant' initlialize in initBetaMix")
  }
  if (alpha0 < 0 | beta0 < 0) {
    stop("cant' initlialize in initBetaMix")
  }
  return(list(z = z, w = w, alpha0 = alpha0, beta0 = beta0, alphaS = alphaS, 
              betaS = betaS, muS = muS, mu0 = mu0))
}

