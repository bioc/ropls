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

#' @rdname getMset
#' @export
setGeneric("getMset",
           function(object, ...) {standardGeneric("getMset")})

#' view
#'
#' Numeric and graphical display of a matrix, a dataframe, an ExpressionSet or a SummarizedExperiment
#'
#' @param x object to be viewed
#' @param printL should the numerical summary be printed?
#' @param plotL should the graphical image be displayed?
#' @param mainC character: plot main title
#' @param subC character(1): plot subtitle
#' @param paletteC character: color palette; either 'heat' [default], 'revHeat', 'grey', 'revGrey', 'palette', 'ramp'
#' @param rowAllL logical: should all rownames be displayed or only the first and
#' last ones?
#' @param rowCexN numeric: size of row labels [default: 1]
#' @param rowMarN numeric: row margin [default: 5.1]
#' @param rowLabC character: label for the y (row) axis
#' @param rowTruncI integer: number of character for truncation of rownames (default,
#' 0, means no truncation)
#' @param colAllL logical: should all column names be displayed or only the first and
#' last ones?
#' @param colCexN numeric: size of column labels [default: 1]
#' @param colMarN numeric: column margin [default: 3.1]
#' @param colLabC character: label for the x (column) axis
#' @param colTruncI integer: number of character for truncation of colnames (default,
#' 0, means no truncation)
#' @param drawScaleL logical: should the color scale be drawn? [default: TRUE]
#' @param delimitReplicatesL logical: should lines be added to the image to delimit
#' replicates in row or column names?
#' @param standardizeL Logical: should columns be standardized for display? 
#' (i.e. subtracting the mean and dividing by the standard deviation) [default: FALSE]
#' @param fig.pdfC character: either 'interactive' [default] or the name of the pdf file to save the figure
#' @param ... Currently not used.
#' @return this method has no output
#' @examples
#' library(ropls)
#' # Get the sacurine dataset
#' data(sacurine)
#' # Display the data matrix
#' view(sacurine[["dataMatrix"]])
#' view(sacurine[["dataMatrix"]][, 1:40], mainC = "'Sacurine' dataset", rowAllL = TRUE,
#' colAllL = TRUE, colTruncI = 13, colMarN = 7)
#' view(sacurine[["dataMatrix"]][, 1:40], mainC = "'Sacurine' dataset", paletteC = "ramp")
#' # Display the sample metadata (dataframe)
#' view(sacurine[["sampleMetadata"]])
#' # Display the ExpressionSet
#' view(sacurine[["eset"]])
#' # Display the SummarizedExperiment
#' view(sacurine[["se"]])
#' @rdname view
#' @export
setGeneric("view",
           function(x, ...) {standardGeneric("view")})

#' @rdname checkW4M
#' @export
setGeneric("checkW4M", function(eset, ...) {standardGeneric("checkW4M")})

#' @rdname toW4M
#' @export
setGeneric("toW4M", function(eset, ...) standardGeneric("toW4M"))
