\name{betaintegral}
\alias{betaintegral}
\title{
Calculates the integral for computing the marginal log--likelihood when ps>pu.
}
\description{
The marginal log likelihood of the alternative model (ps>pu) cannot be computed in closed form and requires numerical integration.
}
\usage{
betaintegral(alphaS, betaS, alpha0, beta0, Nu, nu, Ns, ns,p=1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alphaS}{
	See \code{\link{MarginalGT}}
}
  \item{betaS}{
	See \code{\link{MarginalGT}}
}
  \item{alpha0}{
	See \code{\link{MarginalGT}}
}
  \item{beta0}{
	See \code{\link{MarginalGT}}
}
  \item{Nu}{
	See \code{\link{MarginalGT}}
}
  \item{nu}{
	See \code{\link{MarginalGT}}
}
  \item{Ns}{
	See \code{\link{MarginalGT}}
}
  \item{ns}{
	See \code{\link{MarginalGT}}
}
	\item{p}{
	The number of divisions for integration
	}
}
\value{
Returns the value of the integral. \code{numeric}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Greg Finak <gfinak@fhcrc.org>
}
\note{
Should not be called by the user. Used internally.
}

\seealso{
	See \code{\link{MarginalGT}}
}

\keyword{ models }
