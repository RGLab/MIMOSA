\name{gibbsPsPu}
\alias{gibbsPsPu}
\title{
Return the posterior estimate of the proportion of positive events for the treated and untreated samples given the data for one individual.
}
\description{
Uses Gibbs sampling from the posterior predictive distribution of the beta binomial mixture model to estimate the proportion of positive events for the treated and untreated samples given a single observation and the fitted model parameters.
}
\usage{
gibbsPsPu(curdat = NULL, inits = NULL, which = NULL, N = 10000,burn=5000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{curdat}{
A \code{numeric} matrix or data.frame of response X treatment counts with column names \code{"Ns","ns","Nu","nu"}.
}
  \item{inits}{
A \code{list} of the fitted model parameters returned from \code{BetaMix}

}
  \item{which}{
	A \code{numeric} index into the rows of \code{curdat}, of length 1, specifying the observation for which you want to estimate the ps and pu proportions. 
}
  \item{N}{
	The number of samples to draw.
}
	\item{burn}{
	The number of burn--in samples.	
	}
}
\details{
Samples from the posterior predictive distribution of the model are taken via gibbs sampling to estimate ps (the proportion of positive events in the treated case) and pu (the proportion of positive events in the untreated case).
}
\value{
A \code{numeric} matrix with two columns named \code{"ps","pu"}, and \code{N-burn} rows. The samples from the posterior predictive distribution of ps and pu.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Greg Finak <gfinak@fhcrc.org>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{BetaMix}}
}

\keyword{ models }
