\name{BetaMix}
\alias{BetaMix}
\title{
Use the EM algorithm to fit the beta--binomial mixture for ICS assays.
}
\description{
Fits the beta--binomial mixture to data from flow cytometry ICS assays, or other two--dimensional contingency table data with multiple observations.
}
\usage{
BetaMix(d = NULL, maxiter = 20, shrink = 1,alternative.model="greater", mciter=50, fixedNULL=FALSE, ics=NULL, priorXi=1,inits=NULL, scl=10,K=200)
}
\arguments{
  \item{d}{
  A \code{data.frame} or \code{matrix} with columns \code{"Ns","ns","Nu","nu"}, the treated (s), untreated (u), positive (n) and negative (N) counts, and rows corresopnding to different observations (individuals). In an ICS setting \code{Ns,ns} would be the negative and positive cell counts for the stimulated sample and one cytokine. 
}
  \item{maxiter}{
The maximum number of EM iterationsi for fitting the model. (\code{numeric})
}
  \item{shrink}{
A shrinkage factor multiplying the beta distribution hyperparameters to reduce the variance of the prior distributions. Defaults to 1, no shrinkage. (\code{numeric}).
}
	\item{alternative.model}{
	A \code{character} vector. Either "greater" or "not equal", specifying the alternative model. Default "greater".	
}

\item{mciter}{
Number of iterations for monte carlo integration in the case of \code{alternative.model="greater"}.
}

\item{fixedNULL}{
A \code{logical} indicating whether the null component hyperparameters should be estimated from all the untreated observations and held fixed during model fitting. This can help sometimes if there are few observations for a given condition 
}
\item{ics}{
An object of type \code{ICS}, the object from which \code{d} is extracted, must be provided if \code{fixedNULL=TRUE}. 
}
\item{priorXi}{
A prior on the stimulated and unstimulated component weights in the mixture model. (\code{numeric}), \code{priorXi=1} is equivalent to no prior. 
}
\item{inits}{
Optionally you can specify initialization values for the algorithm. \code{inits} should be a \code{list} with elements \code{alpha0, beta0,alphaS,betaS,z,w}. 
}
\item{scl}{
Tuning parameter for model fitting passed to the \code{optim} function. Fitting is done using the \code{L-BFGS-B} method. \code{scl} should be a \code{numeric} that rescales the parameter estimates to tune the step size of the algorithm.
}
\item{K}{
A tuning parameter for model fitting passed to the \code{optim} function. \code{K} can be thought of as an upper limit on the equivalent sample size of the estimated priors (i.e. \code{alpha0+beta0} or \code{alphaS+betaS}). It controls the upper limit of the search for the hyperparameters used by optim. It is mostly useful to increase this in order to avoid estimates of \code{alpha0,beta0} less than 1. 
}
}
\value{
A \code{list} of fitted model parameters and additional information.
 \item{inits}{A list of fitted model parameters. See \code{\link{initBetaMix}}}
\item{ll}{The log--likelihood of the fitted model}
\item{iter}{The number of iterations required to converge}
\item{lltraj}{The trajectory of log--likelihood from interation 1 to iteration iter.}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Greg Finak <gfinak@fhcrc.org>
}


\seealso{
See also \code{\link{initBetaMix}}, \code{\link{flowModels-package}}
}

\keyword{ models }
