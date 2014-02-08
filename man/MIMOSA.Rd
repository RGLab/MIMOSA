\name{MIMOSA}
\alias{MIMOSA}
\alias{MIMOSA,formula,ExpressionSet-method}
\title{Fit a MIMOSA Model}
\usage{
MIMOSA(formula, data, ...)
}
\arguments{
  \item{formula}{describing the features on the lhs and the
  phenodata on the rhs, supporting extended formula
  interface with conditioning.}

  \item{data}{an \code{ExpressionSet} object with features
  on rows and samples (labelled with phenoData) on
  columns.}

  \item{...}{additional arguments}
}
\value{
an object of type \code{MIMOSAResult}
}
\description{
This method fits a MIMOSA model to count data stored in an
ExpressionSet object.
}
\details{
The ExpressionSet should be fully annotated with
featureData and phenoData. For ICS data, for example,
features would be positive and negative counts for
different cytokine producing cell subsets (i.e.  IFNg_pos,
IFNg_neg) The formula lhs should contain features and the
rhs should contain phenotypic variable. See the vignette
for an example.
}
\examples{
data(ICS)
E<-ConstructMIMOSAExpressionSet(ICS,
  reference=ANTIGEN\%in\%"negctrl",measure.columns=c("CYTNUM","NSUB"),
  other.annotations=c("CYTOKINE","TCELLSUBSET","ANTIGEN","UID"),
  default.cast.formula=component~UID+ANTIGEN+CYTOKINE+TCELLSUBSET,
  .variables=.(TCELLSUBSET,CYTOKINE,UID),
  featureCols=1,ref.append.replace="_REF")

result<-MIMOSA(NSUB+CYTNUM~UID+TCELLSUBSET+CYTOKINE|ANTIGEN,
    data=E, method="EM",
    subset=RefTreat\%in\%"Treatment"&ANTIGEN\%in\%"ENV",
    ref=ANTIGEN\%in\%"ENV"&RefTreat\%in\%"Reference")
}
\seealso{
\code{\link{MIMOSA-package}}
\code{\link{ConstructMIMOSAExpressionSet}}
\code{\link{MIMOSAResult}}
}

