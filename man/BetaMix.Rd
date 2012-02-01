\name{BetaMix}
\alias{BetaMix}
\title{
Use the EM algorithm to fit the beta--binomial mixture for ICS assays.
}
\description{
Fits the beta--binomial mixture to data from flow cytometry ICS assays, or other two--dimensional contingency table data with multiple observations.
}
\usage{
BetaMix(d = NULL, maxiter = 100, shrink = 1,alternative.model="greater", mciter=50)
}
\arguments{
  \item{d}{
  A \code{data.frame} or \code{matrix} with columns \code{"Ns","ns","Nu","nu"}, the treated (s), untreated (u), positive (n) and negative (N) counts, and rows corresopnding to different observations (individuals). In an ICS setting \code{Ns,ns} would be the negative and positive cell counts for the stimulated sample and one cytokine. 
}
  \item{maxiter}{
The maximum number of EM iterations. (\code{numeric})
}
  \item{shrink}{
A shrinkage factor to be applied to the beta distribution hyperparameters. (\code{numeric}).
}
	\item{alternative.model}{
	A \code{character} vector. Either "greater" or "not equal", specifying the alternative model. Default "greater".	
}

\item{mciter}{
Number of iterations for monte carlo integration
}
}
\value{
A \code{list} of fitted model parameters and additional information.
 \item{inits}{A list of fitted model parameters. See \code{\link{initBetaMix}}}
\item{ll}{The log--likelihood of the fitted model}
\item{iter}{The number of iterations required to converge}
\item{lltraj}{The trajectory of log--likelihood from interation 1 to iteration iter.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Greg Finak <gfinak@fhcrc.org>
}


\seealso{
See also \code{\link{initBetaMix}}, \code{\link{flowModels-package}}
}

\keyword{ models }
