\name{MarginalGT}
\alias{MarginalGT}
\title{
Marginal log likelihood for the alternative model (ps>pu)
}
\description{
This function computes the marginal log likelihood for an alternative beta-binomial component of the mixture model (treated proportion is greater than untreated).
}
\usage{
MarginalGT(Ns, ns, Nu, nu, alpha0, beta0, alphaS, betaS,mciter=1000)
}
\arguments{
  \item{Ns}{
The negative counts in the treated condition. A \code{numeric} vector.
}
  \item{ns}{
The positive counts in the treated condition. A \code{numeric} vector.
}
  \item{Nu}{
The negative counts in the untreated condition. A \code{numeric} vector.
}
  \item{nu}{
The positive counts in the untreated condition. A \code{numeric} vector.
}
  \item{alpha0}{
The alpha0 parameter of the beta distribution for the untreated condition. A \code{numeric}.
}
  \item{beta0}{
The beta0 parameter of the beta distribution for the untreated condition. A \code{numeric}
}
  \item{alphaS}{
The alphaS parameter of the beta distribution for the treated condition. A \code{numeric}.
}
  \item{betaS}{
The betaS parameter of the beta distribution for the treated condition. A \code{numeric}
}
\item{mciter}{
An \code{integer} specifying the number of monte-carlo iterations for integration of this marginal likelihood.
}
}
\details{
This function computes the marginal log likelihood of an alternative model component in the beta binomial mixture. It is called internally by the package during model fitting.
}
\value{
A \code{numeric}. The marginal log--likelihood of the data, given the parameters.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Greg Finak <gfinak@fhcrc.org>
}
\note{
This function should not need to be called by the user.
}


%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\keyword{ models }
