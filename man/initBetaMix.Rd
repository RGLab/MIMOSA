\name{initBetaMix}
\alias{initBetaMix}
\title{
Initialize the parameters of the beta binomial mixture model.
}
\description{
Initializes the parameters of the beta binomial mixture model prior to starting the EM algorithm.
}
\usage{
initBetaMix(data = NULL,fixedNULL=FALSE, ics=NULL,alternative="greater",priorXi=1,scl=10,K=200,mciter=200)
}
\arguments{
  \item{data}{
A numeric \code{matrix} or \code{data.frame} with column names \code{"Ns","ns","Nu","nu"}, the negative (N) and positive (n), treated (s) and untreated (u) counts, and rows corresponding to observations. 
}
\item{fixedNULL}{
\code{logical} specifies whether the null component should be estimated from all untreated samples and the hyperparameters fixed during the EM routine. Useful sometimes if there aren't many observations.
}
\item{ics}{
An object of class \code{ICS} from which \code{data} was extracted. Must be provided if \code{fixedNULL=TRUE}.
}
\item{alternative}{
The model. (\code{character}). Either \code{"greater"}, or \code{"not equal"}.
}
\item{priorXi}{
A \code{numeric}. The prior on the component weights. Default of \code{priorXi=1} is no prior.
}
\item{scl}{
Described in \code{\link{BetaMix}}.
}
\item{K}{
Described in \code{\link{BetaMix}}.
}
\item{mciter}{
The number of monte carlo iterations for integration when the model is "greater".
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
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Greg Finak <gfinak@fhcrc.org
}
\note{
This function is called internally from the EM fitting code. It should not be called by the user.
}


\seealso{
\code{\link{BetaMix}}
}

\keyword{ models }
