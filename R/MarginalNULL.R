MarginalNULL <- function(Ns, ns, Nu, nu, alpha0, beta0) {
  lchoose(Ns + ns, ns) + lchoose(Nu + nu, nu) - lbeta(alpha0, beta0) + 
    lbeta(ns + nu + alpha0, Ns + Nu + beta0)
}

