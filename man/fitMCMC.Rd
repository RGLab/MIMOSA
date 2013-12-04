\name{.fitMCMC}
\alias{.fitMCMC}
\title{Fit the MIMOSA model via MCMC}
\usage{
  .fitMCMC(data, inits = NULL, iter = 250000, burn = 50000,
    thin = 1, tune = 100,
    outfile = basename(tempfile(tmpdir = ".", fileext = ".dat")),
    alternative = "greater", UPPER = 0.5, LOWER = 0.15,
    FAST = TRUE, EXPRATE = 1e-04, pXi = 1, seed = 10)
}
\arguments{
  \item{data}{a \code{list} with elements names 'n.stim'
  and 'n.unstim', the stimulated and unstimulated counts.
  Must be at least of dimension 2.}

  \item{inits}{the initialization parameters for the MCMC
  routine. Can be initialized from \code{MDMix} with
  \code{initonly=TRUE}.}

  \item{iter}{the number of Mote Carlo iterations}

  \item{burn}{the number of burn-in iterations}

  \item{thin}{The thinning interval}

  \item{tune}{the number of iterations used for tuning the
  step size}

  \item{outfile}{the output file name}

  \item{alternative}{either 'greater' or 'not equal' for
  fitting the one-sided or two-sided MIMOSA model,
  respectively.}

  \item{UPPER}{tuning parameter for the upper bound on the
  acceptance ratio of each paramter}

  \item{LOWER}{tuning parmeter for the lower bound on the
  acceptance ratio of each paramter}

  \item{FAST}{\code{TRUE,FALSE}. Use the heuristic
  (FAST=TRUE) for fitting a one-sided model rather than
  recomputing the normalization constant via MCMC for each
  step.}

  \item{EXPRATE}{the mean of the prior distribution for the
  model hyperparameters.}

  \item{pXi}{is the prior on the w. 1 by default (as in
  beta(1,1) i.e. uniform).}

  \item{seed}{\code{numeric} random seed}
}
\description{
  This is an internal function that fits the MIMOSA model
  via MCMC. It is called from \code{MIMOSA}
}

