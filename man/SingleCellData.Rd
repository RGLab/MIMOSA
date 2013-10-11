\docType{class}
\name{SingleCellData}
\alias{SingleCellData}
\title{Multivariate MIMOSA data structure to store single-cell data.}
\usage{
  SingleCellData(sc, counts, meta)
}
\arguments{
  \item{sc}{is a list of single-cell expression values,
  similar to a list of flowFrames with measured channels on
  the columns and cells on the rows. These are generally
  the `positive` cells identified by some gating method.}

  \item{counts}{is a vector of total cell counts for each
  sample}

  \item{meta}{is a data.frame of metadata for the samples}
}
\description{
  Multivariate MIMOSA data structure to store single-cell
  data.
}

