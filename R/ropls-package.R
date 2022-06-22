#' PCA, PLS(-DA) and OPLS(-DA) for multivariate analysis and feature selection
#' of omics data
#'
#' Latent variable modeling with Principal Component Analysis (PCA) and Partial
#' Least Squares (PLS) are powerful methods for visualization, regression,
#' classification, and feature selection of omics data where the number of
#' variables exceeds the number of samples and with multicollinearity among
#' variables. Orthogonal Partial Least Squares (OPLS) enables to separately
#' model the variation correlated (predictive) to the factor of interest and
#' the uncorrelated (orthogonal) variation. While performing similarly to PLS,
#' OPLS facilitates interpretation. Successful applications of these
#' chemometrics techniques include spectroscopic data such as Raman
#' spectroscopy, nuclear magnetic resonance (NMR), mass spectrometry (MS) in
#' metabolomics and proteomics, but also transcriptomics data. In addition to
#' scores, loadings and weights plots, the package provides metrics and
#' graphics to determine the optimal number of components (e.g. with the R2 and
#' Q2 coefficients), check the validity of the model by permutation testing,
#' detect outliers, and perform feature selection (e.g. with Variable
#' Importance in Projection or regression coefficients). Input formats include matrix, data frame, SummarizedExperiment, MultiAssayExperiment, ExpressionSet, MultiDataSet. The package can be
#' accessed via a user interface on the Workflow4Metabolomics online
#' resource built upon the Galaxy environment.
#'
#' @import Biobase methods
#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#' @importClassesFrom MultiDataSet MultiDataSet
#' @importFrom grDevices dev.new dev.off pdf rainbow
#' @importFrom graphics abline axTicks axis barplot boxplot layout legend lines mtext par pie plot points rect text
#' @importFrom stats complete.cases cor cov median qchisq qf qnorm sd var
#' @importFrom utils head object.size read.table str tail write.table
#' @name ropls-package
#' @aliases ropls ropls-package
#' @docType package
#' @author Etienne A. Thevenot (CEA)
#'
#' Maintainer: Etienne A. Thevenot <etienne.thevenot@@cea.fr>
#' @keywords package
NULL



















