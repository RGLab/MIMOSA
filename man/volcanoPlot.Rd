\name{volcanoPlot}
\alias{volcanoPlot}
\title{Volcano plot for a MIMOSA model}
\usage{
  volcanoPlot(x, effect_expression = NA, facet_var = NA,
    threshold = 0.01)
}
\arguments{
  \item{x}{A \code{MIMOSAResultList}}

  \item{effect_expression}{an \code{expression} that
  defines the effect size. Usually a function of the
  stimulated and unstimulated proportions from
  \code{countsTable(x,proportion=TRUE)}}

  \item{facet_var}{an \code{expression} defining the
  faceting in ggplot parlance. i.e. \code{~ faceting +
  variables}}

  \item{threshold}{a \code{numeric} value between [0,1] for
  coloring significant observations (based on the q-value)}
}
\description{
  Plots effect size vs posterior probablilty of response
  from a MIMOSAResultList, faceting by the conditioning
  variables.
}

