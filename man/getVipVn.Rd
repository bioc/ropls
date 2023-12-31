% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allGenerics.R, R/opls-methods.R
\name{getVipVn}
\alias{getVipVn}
\alias{getVipVn,opls-method}
\title{getVipVn method for (O)PLS(-DA) models}
\usage{
getVipVn(object, orthoL = FALSE)

\S4method{getVipVn}{opls}(object, orthoL = FALSE)
}
\arguments{
\item{object}{An S4 object of class \code{opls}, created by \code{opls}
function.}

\item{orthoL}{Logical: Should the orthogonal VIP be returned (default is
FALSE and the predictive VIP is returned)}
}
\value{
Numeric vector with a length equal to the number of variables and a
number of columns equal to the number of components
}
\description{
(Orthogonal) VIP of the (O)PLS(-DA) model
}
\examples{

data(sacurine)
attach(sacurine)

sacurine.plsda <- opls(dataMatrix,
                       sampleMetadata[, "gender"])

getVipVn(sacurine.plsda)

detach(sacurine)

}
\references{
Galindo-Prieto B., Eriksson L. and Trygg J. (2014). Variable
influence on projection (VIP) for orthogonal projections to latent
structures (OPLS). Journal of Chemometrics 28, 623-632.
}
\author{
Etienne Thevenot, \email{etienne.thevenot@cea.fr}
}
