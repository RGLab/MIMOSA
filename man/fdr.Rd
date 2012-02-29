\name{fdr}
\alias{fdr}
\title{
Compute the fdr for each observation after fitting the beta binomial mixture model to ICS flow data for positivity calls.
}
\description{
Computes the false discovery rate for each observation fit to a beta binomial mixture model for ICS flow data.
}
\usage{
fdr(z)
}
\arguments{
  \item{z}{
    \code{z} a \code{matrix} (n x 2) of the z's (posterior probabilities) that an observation (row) is responding to stimulation (column 2) or not responding to stimulation (column 1). These are probabilities so each row should sum to 1. 
}
}
\author{
Greg Finak <gfinak@fhcrc.org>
}
\keyword{ htest }
