% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{aminoacids}
\alias{aminoacids}
\title{Amino-Acids Dataset}
\format{
A data frame with the following parameters:
\itemize{
\item AA amino acid
\item PIE lipophilicity constant
of the AA side chain
\item PIF lipophilicity constant of the AA
side chain
\item DGR free energy of transfer of an AA side chain
from protein interior to water
\item SAC water-accessible
surface area of AA's calculated by MOLSV
\item MR molecular
refractivity
\item Lam polarity parameter
\item Vol
molecular volume of AA's calculated by MOLSV
\item DDGTS free
energy of unfolding of the tryptophane synthase a unit of bacteriophage T4
lysosome
}
}
\source{
'aminoacids' dataset.
}
\value{
Data frame (numeric type except the first column, which can be
transformed into row names) with 19 rows and the 9 columns contaning
information about amino acids. For details see the 'Format' section above.
}
\description{
Quantitative structure property relationship (QSPR)
}
\references{
Wold et al. (2001). PLS-regression: a basic tool of
chemometrics. Chemometrics and Intelligent Laboratory Systems. 58:109-130.
}
\keyword{datasets}
