\name{MDMix}
\alias{MDMix}
\title{EM fitting of the Multinomial Dirichlet MIMOSA model.}
\usage{
MDMix(data = NULL, modelmatrix = NULL, alternative = "greater",
  initonly = FALSE)
}
\arguments{
  \item{data}{The observed data}

  \item{modelmatrix}{a model matrix specifying which
  components should be computed}

  \item{alternative}{either 'greater' or 'not equal' to fit
  the one-sided or two-sided model.}

  \item{initonly}{\code{TRUE} or \code{FALSE} to return
  just the initialization parameters.}
}
\value{
An object of class \code{MDMixResult}
}
\description{
This function fits the multinomial dirichelt MIMOSA model
using EM. It can also be used to initialize the model
parameters for the MCMC model.
}
\author{
Greg Finak TODO filtering of pu>ps needs to be corrected
here.
}

