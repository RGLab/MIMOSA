\name{MarginalNULL}
\alias{MarginalNULL}
\title{
Marginal log likelihood for the null model (ps==pu)
}
\description{
This function computes the marginal log likelihood for the null beta-binomial component of the mixture model (no difference between conditions).
}
\usage{
MarginalNULL(Ns, ns, Nu, nu, alpha0, beta0)
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
The alpha0 parameter of the beta distribution. A \code{numeric}.
}
  \item{beta0}{
The beta0 parameter of the beta distribution. A \code{numeric}
}
}
\details{
This function computes the marginal log likelihood of the null model component in the beta binomial mixture. It is called internally by the package during model fitting.
}
\value{
A \code{numeric}. The marginal log--likelihood of the data, given the parameters.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Greg Finak <gfinak@fhcrc.org>
}
\note{
This function should not need to be called by the user.
}


\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\keyword{ models }
