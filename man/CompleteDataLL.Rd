\name{CompleteDataLL}
\alias{CompleteDataLL}
\title{
Compute the complete data log likelihood for the beta binomial mixture model for ICS positivity calls
}
\description{
Computes the complete data log likelihood for the beta binomial mixture model for ICS positivity calls.
}
\usage{
CompleteDataLL(d, alpha0, beta0, alphaS, betaS, z, w, alternative="greater",mciter=50)
}
\arguments{
  \item{d}{
A \code{data.frame} or \code{matrix} of observations. Rows are observations. Columns are named \code{"Ns","ns","Nu","nu"}, the negative stimulated, positive stimulated, negative unstimulated, and positive unstimulated counts.
}
  \item{alpha0}{
	The alpha0 parameter for the null model beta distribution.
}
  \item{beta0}{
	The beta0 parameter for the null model beta distribution.
}
  \item{alphaS}{
	The alphaS parameter for the alternative model beta distribution.
}
  \item{betaS}{
	The betaS parameter for the alternative model beta distribution.
}
  \item{z}{
	The probability of assignment of each observation to to each component of the model. A \code{N x 2} matrix, where \code{N} is the nubmer of observations. Column 1 corresponds to the null model. Rows sum to 1.
}
  \item{w}{
	\code{numeric} vector of length 2. The mixing proportions for the model components. Sum to 1.
}
  \item{alternative}{
	\code{character} vector specifying the alternative model. Either of \code{"greater", "not equal"}.
}
\item{mciter}{
An \code{integer} specifying the number of monte-carlo iterations for integration when \code{alternative="greater"}.
}
}
\details{
	Computes the complete data log-likelihood of the beta binomial mixture model via calls to MarginalNULL and MarginalGT.
}
\value{
	A \code{numeric} vector. The complete--data log likelihood of the model, given \code{d}.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Greg Finak <gfinak@fhcrc.org>
}
\note{
Should not be called by the user. Called internally by BetaMix.
}


\seealso{
\code{\link{BetaMix}}
}
\keyword{ models }
