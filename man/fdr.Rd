\name{fdr}
\alias{fdr}
\alias{fdr.MIMOSAResult}
\alias{fdr.MIMOSAResultList}
\alias{fdr.matrix}
\title{Compute the fdr (q-value) from posterior probabilities}
\usage{
  fdr(x)

  \method{fdr}{matrix} (z)

  \method{fdr}{MIMOSAResult} (x)

  \method{fdr}{MIMOSAResultList} (x)
}
\arguments{
  \item{z}{matrix of posterior probabilties, or a
  \code{MIMOSAResult}, or \code{MIMOSAResultList}}
}
\value{
  a vector of q-values or a list of vectors of q-values.
}
\description{
  Given the z's from a MIMOSA model, calculates the
  q-values for each observation.
}

