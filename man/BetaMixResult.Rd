\docType{class}
\name{BetaMixResult}
\alias{BetaMixResult}
\alias{BetaMixResult-class}
\title{The output of fitting Beta-Binomial EM implementation \code{BetaMix}.}
\usage{
  BetaMixResult(alternative.model = NA_character_,
    cytokine = NA_character_, control = NA_character_,
    stimulation = NA_character_, ll, traj, iter, z, w,
    alpha0, beta0, alphaS, betaS, data,
    pd = new("AnnotatedDataFrame"))
}
\arguments{
  \item{alternative.model}{character, either 'greater' or
  'not equal'}

  \item{cytokine}{'character'}

  \item{control}{'character'}

  \item{stimulation}{'character'}

  \item{ll}{log-likelihood}

  \item{traj}{ll trajectory during fitting}

  \item{iter}{number of iterations}

  \item{z}{the posterior probabilties}

  \item{w}{the responder and non-responder fractions}

  \item{alpha0}{model parameters}

  \item{alphaS}{model parameters}

  \item{beta0}{model parameters}

  \item{betaS}{model parameters}

  \item{data}{the dat}

  \item{pd}{the phenoData}
}
\description{
  BetaMix will return an object of this class.

  BetaMix uses this constructor to create a BetaMixResult
  object
}

