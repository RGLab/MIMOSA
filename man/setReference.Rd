\name{setReference}
\alias{setReference}
\title{Replicate the reference observations across all treatment groups using ddply}
\usage{
  ddply(data,
    .(VAR1,VAR2,VAR3),.fun=setReference,ref=VAR4%in%"Reference",cols=c("VAR5","VAR6"))
}
\arguments{
  \item{dat}{the piece we will work on}

  \item{ref}{an expression evaluating to a logical vector
  that identifies the reference class within the piece.}

  \item{cols}{A character vector of the names of the
  columns that hold our measurements}

  \item{annotations}{a characte vector of additional
  annotations to return for the data}
}
\description{
  Replicate the reference observations across all treatment
  groups using ddply
}
\details{
  Should be passed to .fun argument of ddply.
}

