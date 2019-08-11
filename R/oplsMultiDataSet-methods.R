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
#' # Building the MultiDataSet ('brgeData' Bioconductor package)
#' data("brge_gexp", package = "brgedata")
#' brge_gexp <- brge_gexp[1:1000, ]
#' data("brge_prot", package = "brgedata")
#' brgeMset <- MultiDataSet::createMultiDataSet()
#' brgeMset <- MultiDataSet::add_eset(brgeMset, brge_gexp, dataset.type = "expression",
#' dataset.name = "gexp", warnings = FALSE)
#' brgeMset <- MultiDataSet::add_eset(brgeMset, brge_prot, dataset.type = "proteomics",
#' dataset.name = "prot", warnings = FALSE)
#' # Summary of the MultiDataSet
#' brgeMset
#' # Buidling PCA models for each dataset
#' brgePca <- ropls::opls(brgeMset)
#' # Getting the updated MultiDataSet instance (with scores, loadings, VIP, etc.)
#' brgeMset <- ropls::getMset(brgePca)
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