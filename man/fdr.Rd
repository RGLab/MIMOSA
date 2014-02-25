\name{fdr}
\alias{fdr}
\alias{fdr.MIMOSAResult}
\alias{fdr.MIMOSAResultList}
\alias{fdr.matrix}
\title{Compute the fdr (q-value) from posterior probabilities}
\usage{
fdr(z)

\method{fdr}{matrix}(z)

\method{fdr}{MIMOSAResult}(z)

\method{fdr}{MIMOSAResultList}(z)
}
\arguments{
  \item{z}{matrix of posterior probabilties, or a
  \code{MIMOSAResult}, or \code{MIMOSAResultList}}
}
\value{
a vector of q-values or a list of vectors of q-values.
}
\description{
Given the z's from a MIMOSA model, calculates the q-values
for each observation.
}
\examples{
data(ICS)
E<-ConstructMIMOSAExpressionSet(ICS,
                              reference=ANTIGEN\%in\%'negctrl',measure.columns=c('CYTNUM','NSUB'),
                               other.annotations=c('CYTOKINE','TCELLSUBSET','ANTIGEN','UID'),
                               default.cast.formula=component~UID+ANTIGEN+CYTOKINE+TCELLSUBSET,
                               .variables=.(TCELLSUBSET,CYTOKINE,UID),
                               featureCols=1,ref.append.replace='_REF')
result<-MIMOSA(NSUB+CYTNUM~UID+TCELLSUBSET+CYTOKINE|ANTIGEN,
   data=E, method='EM',
   subset=RefTreat\%in\%'Treatment'&ANTIGEN\%in\%'ENV',
   ref=ANTIGEN\%in\%'ENV'&RefTreat\%in\%'Reference')
qvalues<-fdr(result)
}

