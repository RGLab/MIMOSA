
MarginalGT <- function(Ns, ns, Nu, nu, alpha0, beta0, alphaS, betaS, mciter = 1000) {
    # s1<-log(mean(1-pbeta(rbeta(mciter,alpha0,beta0),alphaS,betaS)))
    # nc.prior<-log(integrate(f1,0,1,alpha.u=alpha0,beta.u=beta0,alpha.s=alphaS,beta.s=betaS)$value)
    nc.prior <- (estimate_logZus(alpha.u = alpha0, beta.u = beta0, alpha.s = alphaS, 
        beta.s = betaS, B = mciter, lower.tail = FALSE))
    nc.post <- (sapply(seq_along(ns), function(i) {
        estimate_logZus(alpha.u = nu[i] + alpha0, beta.u = Nu[i] + beta0, alpha.s = ns[i] + 
            alphaS, beta.s = Ns[i] + betaS, B = mciter, lower.tail = FALSE)
    }))
    fctr <- nc.post - nc.prior
    lchoose(Ns + ns, ns) + lchoose(Nu + nu, nu) - lbeta(alpha0, beta0) - lbeta(alphaS, 
        betaS) + lbeta(ns + alphaS, betaS + Ns) + lbeta(nu + alpha0, beta0 + Nu) + 
        fctr
} 
