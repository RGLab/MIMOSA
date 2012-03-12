\name{likelihoodSurfacePlot-methods}
\docType{methods}
\alias{likelihoodSurfacePlot-methods}
\alias{likelihoodSurfacePlot}
\alias{likelihoodSurfacePlot,BetaMixResult,integer,integer,integer,integer,numeric-method}
\title{ Method to plot a likelihood surface from a \code{BetaMixResult}}
\description{
 Plots the likelihood surface for one observation from a \code{BetaMixResult}, as well as the observed counts for the treated and untreated samples. The x,y axes are the counts from the treated and untreated samples.
}
\section{Methods}{
\describe{

\item{\code{signature(bmr = "BetaMixResult", grid.x = "integer", grid.y = "integer", whch = "integer", nlevels = "integer", fdr = "numeric")}}{

The user provides the \code{BetaMixResult}, the upper limit of the grid points in the x and y dimensions (\code{grid.x}, \code{grid.y}), representing the stimulated counts, and the index of the sample to plot (\code{whch}). \code{nlevels} specifies the number of contour levels to plot, and \code{fdr} is the false discovery rate threshold below which observations should be colored as significant on the plot.

}
}}
\keyword{methods}
