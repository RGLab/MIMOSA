\name{setReference}
\alias{setReference}
\title{Replicate the reference observations across all treatment groups using ddply}
\usage{
  setReference(dat, ref = NULL, cols = NULL,
    annotations = NULL, default.formula = component ~ ...)
}
\arguments{
  \item{dat}{the piece we will work on}

  \item{ref}{an expression evaluating to a logical vector
  that identifies the reference class within the piece.}

  \item{cols}{A character vector of the names of the
  columns that hold our measurements}

  \item{annotations}{A character vector of additional
  annotation columns}

  \item{default.formula}{a default formula to be used for
  casting the data.}
}
\description{
  Replicate the reference observations across all treatment
  groups using ddply
}
\details{
  Should be passed to .fun argument of ddply.
}

