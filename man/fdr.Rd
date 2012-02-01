\name{fdr}
\alias{fdr}
\title{
Compute the fdr for each observation after fitting the beta binomial mixture model to ICS flow data for positivity calls.
}
\description{
Computes the false discovery rate for each observation fit to a beta binomial mixture model for ICS flow data.
}
\usage{
fdr(inits)
}
\arguments{
  \item{inits}{
    \code{inits} is a list of fitted model parameters, an element of the list structure returned by \code{BetaMix}. 
}
}
\author{
Greg Finak <gfinak@fhcrc.org>
}
\keyword{ htest }
