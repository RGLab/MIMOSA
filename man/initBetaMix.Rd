\name{initBetaMix}
\alias{initBetaMix}
\title{
Initialize the parameters of the beta binomial mixture model.
}
\description{
Initializes the parameters of the beta binomial mixture model prior to starting the EM algorithm.
}
\usage{
initBetaMix(data = NULL)
}
\arguments{
  \item{data}{
A numeric \code{matrix} or \code{data.frame} with column names \code{"Ns","ns","Nu","nu"}, the negative (N) and positive (n), treated (s) and untreated (u) counts, and rows corresponding to observations. 
}
}
\details{
	The model parameters are initialized to \code{alpha0=beta0=alphaS=betaS=1}, representing a uniform prior over the \code{ps} and \code{pu} treated and untreated proportions. The \code{z}s are initialized to a random assignment of observations to model components with equal probability of membership in each component.
}
\value{
A list of initial model parameters
\item{z}{The initial assignment of observations to model components}
\item{w}{The initial weight of the model components in the mixture}
\item{alpha0}{The initial alpha0 value}
\item{beta0}{The initial beta0 value}
\item{alphaS}{The initial alphaS value}
\item{betaS}{The initial betaS value}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Greg Finak <gfinak@fhcrc.org
}
\note{
This function is called internally from the EM fitting code. It should not be called by the user.
}


\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\keyword{ models }