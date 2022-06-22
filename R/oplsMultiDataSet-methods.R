####    getMset (oplsMultiDataSet)   ####

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