\name{match.elispot.antigens}
\alias{match.elispot.antigens}
\title{Sumarize the replicates in an elispot epitope mapping data set from SCHARP}
\usage{
  match.elispot.antigens(x,
    CONTROL = STAGE \%in\% "CTRL" & PROTEIN \%in\% "Media+cells",
    WELLMULTIPLIER = 1000,
    replicates = c("REP1", "REP2", "REP3", "REP4", "REP5", "REP6"),
    cells.per.well = "CELLWELL",
    SAMPLE = !STAGE \%in\% "CTRL")
}
\arguments{
  \item{x}{is the piece of data currently being worked on.
  It is a unique combination of the grouping variables
  passed to ddply}

  \item{CONTROL}{is an expression that evaluates to a
  logical vector specifying which rows of the data frame
  are control observations.}

  \item{WELLMULTIPLIER}{is a numeric argument that
  specifies the factor by which to multiply the denominator
  (in the \code{cells.per.well} argument). i.e., if there
  are 100K cells per well, but the cells.per.well column
  has 100, the WELLMULTIPLIER would be 1000}

  \item{replicates}{is a character vector identifying the
  column names that contain the replicated observations.
  These will be summed.}

  \item{cells.per.well}{is a character vector that
  identifies the column containing the number of cells per
  well.}

  \item{SAMPLE}{is an expression evaluating to a logical
  vector that identifies the rows of the data frame which
  are antigen or peptide stimualtions rather than control
  samples.}
}
\description{
  Thus function will summarize the replicate observations
  for the negative controls and stimulations in an epitope
  mapping data set from SCHARP. The user provides
  information on grouping structure, cells per well,
  replicate observation columns and so forth. The function
  is meant to be passed to ddply
}

