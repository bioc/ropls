####    getMset (oplsMultiDataSet)   ####

#' getMset method
#'
#' Extracts the complemented MultiDataSet when opls has been applied to a MultiDataSet
#'
#' @aliases getMset getMset, oplsMultiDataSet-method
#' @param object An S4 object of class \code{oplsMultiDataSet}, created by \code{opls}
#' function applied to a MultiDataSet
#' @param ... Currently not used
#' @return An S4 object of class \code{MultiDataSet}.
#' @examples
#' # Building a MultiDataSet object
#' # Loading the 'NCI60_4arrays' from the 'omicade4' package
#' data("NCI60_4arrays", package = "omicade4")
#' # Selecting two of the four datasets
#' setNamesVc <- c("agilent", "hgu95")
#' # Creating the MultiDataSet instance
#' nciMset <- MultiDataSet::createMultiDataSet()
#' # Adding the two datasets as ExpressionSet instances
#' for (setC in setNamesVc) {
#'   # Getting the data
#'   exprMN <- as.matrix(NCI60_4arrays[[setC]])
#'   pdataDF <- data.frame(row.names = colnames(exprMN),
#'                         cancer = substr(colnames(exprMN), 1, 2),
#'                         stringsAsFactors = FALSE)
#'   fdataDF <- data.frame(row.names = rownames(exprMN),
#'                         name = rownames(exprMN),
#'                         stringsAsFactors = FALSE)
#'   # Building the ExpressionSet
#'   eset <- Biobase::ExpressionSet(assayData = exprMN,
#'                                  phenoData = new("AnnotatedDataFrame",
#'                                                  data = pdataDF),
#'                                  featureData = new("AnnotatedDataFrame",
#'                                                    data = fdataDF),
#'                                  experimentData = new("MIAME",
#'                                                       title = setC))
#'   # Adding to the MultiDataSet
#'   nciMset <- MultiDataSet::add_eset(nciMset, eset, dataset.type = setC,
#'                                     GRanges = NA, warnings = FALSE)
#' }
#' # Summary of the MultiDataSet
#' nciMset
#' # Principal Component Analysis of each data set
#' nciPca <- ropls::opls(nciMset)
#' # Getting the MultiDataSet with additional information in pData (scores)
#' and fData (loadings; also VIP, coefficients in case of PLS modeling) data frames
#' nciMset <- ropls::getMset(nciPca)
#' @rdname getMset
#' @export
setMethod("getMset", "oplsMultiDataSet",
          function(object) {
            Mset <- MultiDataSet::createMultiDataSet()
            for (setI in 1:length(object@oplsLs)) {
              Mset <- MultiDataSet::add_eset(Mset,
                                             ropls::getEset(object@oplsLs[[setI]]),
                                             dataset.type = names(object@oplsLs)[setI],
                                             GRanges = NA,
                                             overwrite = TRUE,
                                             warnings = FALSE)
            }
            return(Mset)
          })