# posteriorSample.R Created on: Dec 6, 2011 Author: finak


posteriorSample <- function(bmr, which, nsamps = 1000) {
  
  inits <- list(alpha0 = bmr@alpha0, beta0 = bmr@beta0, alphaS = bmr@alphaS, 
                betaS = bmr@betaS, w = bmr@w, z = bmr@z)
  s <- gibbsPsPu(curdat = bmr@data, inits = inits, alt = bmr@alternative.model, 
                 which = which, N = nsamps + 1000, burn = 1000)
  S <- t(apply(s, 1, function(x) {
    cbind(ns = rbinom(1, bmr@data[which, "ns"] + bmr@data[which, 
                                                          "Ns"], x[2]), nu = rbinom(1, bmr@data[which, "nu"] + bmr@data[which, 
                                                                                                                        "Nu"], x[1]))
  }))
  colnames(S) <- c("ns", "nu")
  S
} 
