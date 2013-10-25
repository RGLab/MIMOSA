\docType{methods}
\name{countsTable}
\alias{countsTable}
\alias{countsTable,MCMCResult-method}
\alias{countsTable,MDMixResult-method}
\alias{countsTable,MIMOSAResult-method}
\alias{countsTable,MIMOSAResultList-method}
\alias{countsTable.MIMOSAResultList}
\title{Extract the table of counts from a MIMOSA model}
\usage{
  \method{countsTable}{MIMOSAResultList} (object,
    proportion = FALSE)
}
\arguments{
  \item{object}{a \code{MIMOSAResult}}

  \item{proportion}{\code{logical} return the counts or the
  proportions}
}
\value{
  a \code{data.frame} of counts to which the model was fit.

  a \code{data.frame} of counts for the stimulated and
  unstimulated samples
}
\description{
  Extract the table of counts from a MIMOSA model
}

