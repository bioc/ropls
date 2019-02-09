#' @export
setGeneric("opls",
           function(x, ...) standardGeneric("opls"))

#' @rdname tested
#' @export
setGeneric("tested",
           function(object, ...) standardGeneric("tested"))

#' @rdname getSummaryDF
#' @export
setGeneric("getSummaryDF",
           function(object, ...) {standardGeneric("getSummaryDF")})

#' @rdname getPcaVarVn
#' @export
setGeneric("getPcaVarVn",
           function(object, ...) {standardGeneric("getPcaVarVn")})

#' @rdname getScoreMN
#' @export
setGeneric("getScoreMN",
           function(object, ...) {standardGeneric("getScoreMN")})

#' @rdname getLoadingMN
#' @export
setGeneric("getLoadingMN",
           function(object, ...) {standardGeneric("getLoadingMN")})

#' @rdname getWeightMN
#' @export
setGeneric("getWeightMN",
           function(object, ...) {standardGeneric("getWeightMN")})

#' @rdname getVipVn
#' @export
setGeneric("getVipVn",
           function(object, ...) {standardGeneric("getVipVn")})

#' @rdname getSubsetVi
#' @export
setGeneric("getSubsetVi",
           function(object, ...) {standardGeneric("getSubsetVi")})

#' @rdname getEset
#' @export
setGeneric("getEset",
           function(object, ...) {standardGeneric("getEset")})

#' @rdname checkW4M
#' @export
setGeneric("checkW4M", function(eset, ...) {standardGeneric("checkW4M")})

#' @rdname toW4M
#' @export
setGeneric("toW4M", function(eset, ...) standardGeneric("toW4M"))
