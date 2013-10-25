\name{getZ}
\alias{getW}
\alias{getW.MIMOSAResult}
\alias{getW.MIMOSAResultList}
\alias{getZ}
\alias{getZ.MIMOSAResult}
\alias{getZ.MIMOSAResultList}
\title{Extract the posterior probabilities of response from a MIMOSA model}
\usage{
  getZ(x)

  getZ.MIMOSAResultList(x)

  getZ.MIMOSAResult(x)

  getW(x)

  getW.MIMOSAResultList(x)

  getW.MIMOSAResult(x)
}
\arguments{
  \item{x}{output from a MIMOSA model}
}
\value{
  a \code{matrix} of posterior probabilities

  a \code{vector} of component weights
}
\description{
  Extract the posterior probabilities of response from a
  MIMOSA model

  Extract the component weights from a MIMOSA model
}

