\name{MIMOSA}
\alias{MIMOSA}
\alias{MIMOSA,formula,ExpressionSet-method}
\title{Fit a MIMOSA Model}
\arguments{
  \item{formula}{describing the features on the lhs and the
  phenodata on the rhs, supporting extended formula
  interface with conditioning.}

  \item{data}{an \code{ExpressionSet} object with features
  on rows and samples (labelled with phenoData) on
  columns.}
}
\value{
  an object of type \code{MIMOSAResult}
}
\description{
  This method fits a MIMOSA model to count data stored in
  an ExpressionSet object.
}
\details{
  The ExpressionSet should be fully annotated with
  featureData and phenoData. For ICS data, for example,
  features would be positive and negative counts for
  different cytokine producing cell subsets (i.e.
  IFNg_pos, IFNg_neg) The formula lhs should contain
  features and the rhs should contain phenotypic variable.
  See the vignette for an example.
}

