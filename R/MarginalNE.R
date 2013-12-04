MarginalNE <- function(Ns, ns, Nu, nu, alpha0, beta0, alphaS, betaS) {
  res <- lchoose(Ns + ns, ns) + lchoose(Nu + nu, nu) - lbeta(alpha0, 
                                                             beta0) - lbeta(alphaS, betaS) + lbeta(alpha0 + nu, beta0 + 
                                                                                                     Nu) + lbeta(alphaS + ns, betaS + Ns)
  res
}

