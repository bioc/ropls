% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allGenerics.R, R/opls-methods.R
\name{getOpls}
\alias{getOpls}
\alias{getOpls,SummarizedExperiment-method}
\alias{getOpls,MultiAssayExperiment-method}
\title{Getting the models from a SummarizedExperiment or a MultiAssayExperiment object}
\usage{
getOpls(object)

\S4method{getOpls}{SummarizedExperiment}(object)

\S4method{getOpls}{MultiAssayExperiment}(object)
}
\arguments{
\item{object}{An S4 object of class \code{SummarizedExperiment} or \code{MultiAssayExperiment},
once processed by the \code{opls} method}
}
\value{
List of opls models contained in the SummarizedExperiment object(s)
}
\description{
The models are extracted as a list
}
\examples{


# Getting the sacurine data set as a SummarizedExperiment
data(sacurine)
sac.se <- sacurine[["se"]]

# Building the PLS-DA model
sac.se <- opls(sac.se, "gender")

# Getting the models
sac_opls.ls <- getOpls(sac.se)
names(sac_opls.ls)

# Plotting the score plot from the PLS-DA model of the gender response
plot(sac_opls.ls[["gender_PLSDA"]], typeVc = "x-score")

}
\author{
Etienne Thevenot, \email{etienne.thevenot@cea.fr}
}
