#### view (ExpressionSet) ####

#' view
#'
#' Numeric and graphical display of exprs, pData and fData slots
#' from an ExpressionSet object
#'
#' @param x data frame to be viewed
#' @param printL should the numerical summary be printed?
#' @param plotL should the graphical image be displayed?
#' @param mainC character: plot main title
#' @param paletteC character: color palette; either 'heat' [default], 'revHeat', 'grey', 'revGrey', 'palette', 'ramp'
#' @param rowAllL logical: should all rownames be displayed or only the first and
#' last ones?
#' @param rowCexN numeric: size of row labels [default: 1]
#' @param rowMarN numeric: row margin [default: 6.1]
#' @param rowLabC character: label for the y (row) axis
#' @param rowTruncI integer: number of character for truncation of rownames (default,
#' 0, means no truncation)
#' @param colAllL logical: should all column names be displayed or only the first and
#' last ones?
#' @param colCexN numeric: size of column labels [default: 1]
#' @param colMarN numeric: column margin [default: 6.1]
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
#' data(sacurine)
#' sacSet <- Biobase::ExpressionSet(assayData = t(sacurine[["dataMatrix"]]), 
#'                                  phenoData = new("AnnotatedDataFrame", 
#'                                                  data = sacurine[["sampleMetadata"]]), 
#'                                  featureData = new("AnnotatedDataFrame", 
#'                                                    data = sacurine[["variableMetadata"]]),
#'                                  experimentData = new("MIAME", 
#'                                                       title = "sacurine"))
#' view(sacSet)
#' @rdname view
#' @export
setMethod("view", signature(x = "ExpressionSet"),
          function(x,
                   printL = TRUE,
                   plotL = TRUE,
                   mainC = "",
                   paletteC = c("heat",
                                "revHeat",
                                "grey",
                                "revGrey",
                                "palette",
                                "ramp")[1],
                   rowAllL = FALSE,
                   rowCexN = 1,
                   rowMarN = 5.1,
                   rowLabC = "",
                   rowTruncI = 0,
                   colAllL = FALSE,
                   colCexN = 1,
                   colMarN = 1.1,
                   colLabC = "",
                   colTruncI = 0,
                   drawScaleL = TRUE,
                   delimitReplicatesL = FALSE,
                   standardizeL = FALSE,
                   fig.pdfC = "interactive") {
            
            if (is.na(mainC) || mainC == "")
              mainC <- Biobase::experimentData(x)@title
            
            message("'exprs(x)':")
            ropls::view(Biobase::exprs(x),
                        printL = printL,
                        plotL = plotL,
                        mainC = mainC,
                        subC = "exprs",
                        paletteC = paletteC,
                        rowAllL = rowAllL,
                        rowCexN = rowCexN,
                        rowMarN = rowMarN,
                        rowLabC = rowLabC,
                        rowTruncI = rowTruncI,
                        colAllL = colAllL,
                        colCexN = colCexN,
                        colMarN = colMarN,
                        colLabC = colLabC,
                        colTruncI = colTruncI,
                        drawScaleL = drawScaleL,
                        delimitReplicatesL = delimitReplicatesL,
                        standardizeL = standardizeL,
                        fig.pdfC = fig.pdfC)
            message("'pData(x)':")
            ropls::view(Biobase::pData(x),
                        printL = printL,
                        plotL = plotL,
                        mainC = mainC,
                        subC = "pData",
                        paletteC = paletteC,
                        rowAllL = rowAllL,
                        rowCexN = rowCexN,
                        rowMarN = rowMarN,
                        rowLabC = rowLabC,
                        rowTruncI = rowTruncI,
                        colAllL = colAllL,
                        colCexN = colCexN,
                        colMarN = colMarN,
                        colLabC = colLabC,
                        colTruncI = colTruncI,
                        drawScaleL = drawScaleL,
                        delimitReplicatesL = delimitReplicatesL,
                        standardizeL = standardizeL,
                        fig.pdfC = fig.pdfC)
            message("'fData(x)':")
            ropls::view(Biobase::fData(x),
                        printL = printL,
                        plotL = plotL,
                        mainC = mainC,
                        subC = "fData",
                        paletteC = paletteC,
                        rowAllL = rowAllL,
                        rowCexN = rowCexN,
                        rowMarN = rowMarN,
                        rowLabC = rowLabC,
                        rowTruncI = rowTruncI,
                        colAllL = colAllL,
                        colCexN = colCexN,
                        colMarN = colMarN,
                        colLabC = colLabC,
                        colTruncI = colTruncI,
                        drawScaleL = drawScaleL,
                        delimitReplicatesL = delimitReplicatesL,
                        standardizeL = standardizeL,
                        fig.pdfC = fig.pdfC)
            
          })

#### view (data.frame) ####

#' view
#'
#' Numeric and graphical display of a data frame
#'
#' @param x data frame to be viewed
#' @param printL should the numerical summary be printed?
#' @param plotL should the graphical image be displayed?
#' @param mainC character: plot main title
#' @param subC character: plot subtitle
#' @param paletteC character: color palette; either 'heat' [default], 'revHeat', 'grey', 'revGrey', 'palette', 'ramp'
#' @param rowAllL logical: should all rownames be displayed or only the first and
#' last ones?
#' @param rowCexN numeric: size of row labels [default: 1]
#' @param rowMarN numeric: row margin [default: 6.1]
#' @param rowLabC character: label for the y (row) axis
#' @param rowTruncI integer: number of character for truncation of rownames (default,
#' 0, means no truncation)
#' @param colAllL logical: should all column names be displayed or only the first and
#' last ones?
#' @param colCexN numeric: size of column labels [default: 1]
#' @param colMarN numeric: column margin [default: 6.1]
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
#' data(sacurine)
#' view(sacurine[["sampleMetadata"]])
#' @rdname view
#' @export
setMethod("view", signature(x = "data.frame"),
          function(x,
                   printL = TRUE,
                   plotL = TRUE,
                   mainC = "",
                   subC = "",
                   paletteC = c("heat",
                                "revHeat",
                                "grey",
                                "revGrey",
                                "palette",
                                "ramp")[1],
                   rowAllL = FALSE,
                   rowCexN = 1,
                   rowMarN = 5.1,
                   rowLabC = "",
                   rowTruncI = 0,
                   colAllL = FALSE,
                   colCexN = 1,
                   colMarN = 1.1,
                   colLabC = "",
                   colTruncI = 0,
                   drawScaleL = TRUE,
                   delimitReplicatesL = FALSE,
                   standardizeL = FALSE,
                   fig.pdfC = "interactive") {
            
            if (printL)
              ropls::strF(x)
            
            if (plotL) {
              
              if (cumprod(dim(x))[2] == 0) {
                
                warning("Data frame with no row and/or no column cannot be plotted.",
                        immediate. = TRUE,
                        call. = FALSE)
                
              } else {
                
                class.vc <- sapply(x, data.class)
                class.vuc <- unique(class.vc)
                
                if ("logical" %in% class.vuc) {
                  logical.vi <- which(class.vc == "logical")
                  message(length(logical.vi), " data.frame 'logical' column(s) converted to 'numeric' for plotting.")
                  for (j in logical.vi)
                    x[, j] <- as.numeric(x[, j])
                }
                if ("character" %in% class.vuc) {
                  character.vi <- which(class.vc == "character")
                  message(length(character.vi), " data.frame 'character' column(s) converted to 'numeric' for plotting.")
                  for (j in character.vi) {
                    x.fc <- factor(x[, j])
                    x[, j] <- as.numeric(x.fc)
                  }
                }
                if ("factor" %in% class.vuc) {
                  factor.vi <- which(class.vc == "factor")
                  message(length(factor.vi), " data.frame 'factor' column(s) converted to 'numeric' for plotting.")
                  for (j in factor.vi) {
                    x[, j] <- as.numeric(x[, j])
                  }
                }
                
                if (all(sapply(x, data.class) == "numeric")) {
                  imageF(x = as.matrix(x),
                         mainC = mainC,
                         subC = subC,
                         paletteC = paletteC,
                         rowAllL = rowAllL,
                         rowCexN = rowCexN,
                         rowMarN = rowMarN,
                         rowLabC = rowLabC,
                         rowTruncI = rowTruncI,
                         colAllL = colAllL,
                         colCexN = colCexN,
                         colMarN = colMarN,
                         colLabC = colLabC,
                         colTruncI = colTruncI,
                         drawScaleL = drawScaleL,
                         delimitReplicatesL = delimitReplicatesL,
                         standardizeL = standardizeL,
                         fig.pdfC = fig.pdfC)
                } else {
                  warning("Data frame could not be plotted because some columns could not be converted to 'numeric'.",
                          call. = FALSE)
                }
              }
            }
          })


#### view (matrix) ####

#' view
#'
#' Numeric and graphical display of a matrix
#'
#' @param x matrix to be viewed
#' @param printL should the numerical summary be printed?
#' @param plotL should the graphical image be displayed?
#' @param mainC character: plot main title
#' @param subC character: plot subtitle
#' @param paletteC character: color palette; either 'heat' [default], 'revHeat', 'grey', 'revGrey', 'palette', 'ramp'
#' @param rowAllL logical: should all rownames be displayed or only the first and
#' last ones?
#' @param rowCexN numeric: size of row labels [default: 1]
#' @param rowMarN numeric: row margin [default: 6.1]
#' @param rowLabC character: label for the y (row) axis
#' @param rowTruncI integer: number of character for truncation of rownames (default,
#' 0, means no truncation)
#' @param colAllL logical: should all column names be displayed or only the first and
#' last ones?
#' @param colCexN numeric: size of column labels [default: 1]
#' @param colMarN numeric: column margin [default: 6.1]
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
#' data(sacurine)
#' dataMN <- sacurine[["dataMatrix"]]
#' view(dataMN)
#' view(dataMN[, 1:40], mainC = "'Sacurine' dataset", rowAllL = TRUE,
#' colAllL = TRUE, colTruncI = 13, colMarN = 7)
#' view(dataMN[, 1:40], mainC = "'Sacurine' dataset", paletteC = "ramp")
#' sacSet <- Biobase::ExpressionSet(assayData = t(sacurine[["dataMatrix"]]), 
#'                                  phenoData = new("AnnotatedDataFrame", 
#'                                                  data = sacurine[["sampleMetadata"]]), 
#'                                  featureData = new("AnnotatedDataFrame", 
#'                                                    data = sacurine[["variableMetadata"]]),
#'                                  experimentData = new("MIAME", 
#'                                                       title = "sacurine"))
#' view(sacSet)
#' @rdname view
#' @export
setMethod("view", signature(x = "matrix"),
          function(x,
                   printL = TRUE,
                   plotL = TRUE,
                   mainC = "",
                   subC = "",
                   paletteC = c("heat",
                                "revHeat",
                                "grey",
                                "revGrey",
                                "palette",
                                "ramp")[1],
                   rowAllL = FALSE,
                   rowCexN = 1,
                   rowMarN = 5.1,
                   rowLabC = "",
                   rowTruncI = 0,
                   colAllL = FALSE,
                   colCexN = 1,
                   colMarN = 1.1,
                   colLabC = "",
                   colTruncI = 0,
                   drawScaleL = TRUE,
                   delimitReplicatesL = FALSE,
                   standardizeL = FALSE,
                   fig.pdfC = "interactive") {
            
            if (printL)
              strF(x)
            
            if (plotL) {
              
              if (cumprod(dim(x))[2] == 0) {
                
                warning("Matrix with no row and/or no column cannot be plotted.",
                        immediate. = TRUE,
                        call. = FALSE)
                
              } else {
                
                if (mode(x) == "logical") {
                  warning("Matrix converted from 'logical' to 'numeric' mode for plotting.",
                          immediate. = TRUE,
                          call. = FALSE)
                  mode(x) <- "numeric"
                } else if (mode(x) == "character") {
                  warning("Matrix converted from 'character' to 'numeric' mode for plotting.",
                          immediate. = TRUE,
                          call. = FALSE)
                  x <- apply(x, 2, function(y) {
                    y <- factor(y)
                    levels(y) <- seq_along(levels(y))
                    y
                  })
                }
                imageF(x = x,
                       mainC = mainC,
                       subC = subC,
                       paletteC = paletteC,
                       rowAllL = rowAllL,
                       rowCexN = rowCexN,
                       rowMarN = rowMarN,
                       rowLabC = rowLabC,
                       rowTruncI = rowTruncI,
                       colAllL = colAllL,
                       colCexN = colCexN,
                       colMarN = colMarN,
                       colLabC = colLabC,
                       colTruncI = colTruncI,
                       drawScaleL = drawScaleL,
                       delimitReplicatesL = delimitReplicatesL,
                       standardizeL = standardizeL,
                       fig.pdfC = fig.pdfC)
              }
              
            }
            
            invisible(NA)
            
          })


#' Printed summary of an R object
#'
#' Display of the class, mode, size and first...last values from the object; used
#' inside the 'view' wrapper method
#'
#' @param tableMF Input matrix, dataframe or vector
#' @param borderI Number of border (first and last) rows and columns to display
#' @param bigMarkC Big mark separator for summary results
#' @return This function has no output.
#' @seealso \code{\link{str}}, \code{\link{view}}
#' @examples
#'
#' data(sacurine)
#' strF(sacurine[['dataMatrix']])
#' strF(sacurine[['sampleMetadata']])
#'
#' @rdname view
#' @export strF
strF <- function(tableMF,
                 borderI = 2,
                 bigMarkC = ",") {
  
  if (any(class(tableMF) %in% c("character", "integer", "logical", "numeric", "double"))) {
    classC <- "vector"
  } else
    classC <- class(tableMF)
  
  numericL <- mode(tableMF) %in% c("numeric", "integer", "double")
  
  if (!(classC %in% c("vector", "matrix", "data.frame"))) {
    str(tableMF)
    return(invisible(NULL))
  }
  
  .header(tableMF = tableMF,
          borderI = borderI,
          bigMarkC = bigMarkC,
          classC = classC,
          numericL = numericL)
  
  if (tail(cumprod(dim(tableMF)), 1) > 0)
    .table(tableMF = tableMF,
           borderI = borderI,
           classC = classC,
           numericL = numericL)
  
} ## strF

.header <- function(tableMF,
                    borderI,
                    bigMarkC,
                    classC,
                    numericL) {
  
  switch(classC,
         
         vector = {
           
           headerDF <- data.frame(length = format(length(tableMF), big.mark = bigMarkC),
                                  class = class(tableMF),
                                  mode = mode(tableMF),
                                  typeof = typeof(tableMF),
                                  size = format(object.size(tableMF), units = "Mb"),
                                  stringsAsFactors = FALSE)
           
           if (numericL)
             headerDF <- cbind.data.frame(headerDF,
                                          data.frame(min = formatC(min(tableMF, na.rm = TRUE),
                                                                   digits = 2, format = "g"),
                                                     mean = formatC(mean(tableMF, na.rm = TRUE),
                                                                    digits = 2, format = "g"),
                                                     median = formatC(median(tableMF, na.rm = TRUE),
                                                                      digits = 2, format = "g"),
                                                     max = formatC(max(tableMF, na.rm = TRUE),
                                                                   digits = 2, format = "g")))
           
           
           
         }, ## vector
         
         matrix = {
           
           headerDF <- data.frame(dim = paste(format(nrow(tableMF), big.mark = bigMarkC), format(ncol(tableMF), big.mark = bigMarkC), sep = " x "),
                                  class = class(tableMF),
                                  mode = mode(tableMF),
                                  typeof = typeof(tableMF),
                                  size = format(object.size(tableMF), units = "Mb"),
                                  NAs = length(which(is.na(tableMF))))
           
           if (numericL)
             headerDF <- cbind.data.frame(headerDF,
                                          data.frame(min = formatC(min(tableMF, na.rm = TRUE),
                                                                   digits = 2, format = "g"),
                                                     mean = formatC(mean(tableMF, na.rm = TRUE),
                                                                    digits = 2, format = "g"),
                                                     median = formatC(median(tableMF, na.rm = TRUE),
                                                                      digits = 2, format = "g"),
                                                     max = formatC(max(tableMF, na.rm = TRUE),
                                                                   digits = 2, format = "g")))
           
           
         }, ## matrix
         
         data.frame = {
           
           claVc <- sapply(tableMF, data.class)
           
           if (length(claVc) > 2 * borderI)
             claVc <- c(head(claVc, borderI),
                        "...",
                        tail(claVc, borderI))
           
           if (!is.null(names(claVc))) {
             
             if (length(claVc) > 2 * borderI)
               names(claVc)[borderI + 1] <- "..."
             
             claDF <- as.data.frame(t(claVc))
             
             rownames(claDF) <- ""
             
             print(claDF)
             
           } else {
             
             claVc <- paste(claVc, collapse = " ")
             class(claVc) <- "table"
             
             message(claVc)
             
           }
           
           headerDF <- data.frame(nRow = format(dim(tableMF)[1], big.mark = bigMarkC),
                                  nCol = format(dim(tableMF)[2], big.mark = bigMarkC),
                                  size = format(object.size(tableMF), units = "Mb"),
                                  NAs = length(which(is.na(tableMF))))
         }) ## data.frame
  
  rownames(headerDF) <- ""
  
  print(headerDF)
  
}


.table <- function(tableMF,
                   borderI,
                   classC,
                   numericL) {
  
  if (classC %in% c("data.frame", "matrix")) {
    
    dimAbbVl <- dim(tableMF) > 2 * borderI
    
    if (is.data.frame(tableMF)) {
      if (dimAbbVl[2]) {
        bordColVi <- c(1:borderI,
                       (ncol(tableMF) - borderI + 1):ncol(tableMF))
      } else
        bordColVi <- 1:ncol(tableMF)
      
      for (j in bordColVi)
        if (is.factor(tableMF[, j]))
          tableMF[, j] <- as.character(tableMF[, j])
    }
    
    if (all(dimAbbVl)) {
      
      tableMF <- rbind(cbind(tableMF[1:borderI, 1:borderI, drop = FALSE],
                             ... = rep("...", times = borderI),
                             tableMF[1:borderI, (ncol(tableMF) - borderI + 1):ncol(tableMF), drop = FALSE]),
                       rep("...", 2 * borderI + 1),
                       cbind(tableMF[(nrow(tableMF) - borderI + 1):nrow(tableMF), 1:borderI, drop = FALSE],
                             ... = rep("...", times = borderI),
                             tableMF[(nrow(tableMF) - borderI + 1):nrow(tableMF), (ncol(tableMF) - borderI + 1):ncol(tableMF), drop = FALSE]))
      
      if (classC == "matrix") {
        
        if (!is.null(rownames(tableMF)) && any(duplicated(rownames(tableMF)))) {
          rownames(tableMF) <- make.names(rownames(tableMF),
                                          unique = TRUE)
        } else if (is.null(rownames(tableMF)))
          rownames(tableMF) <- c(1:borderI,
                                 "...",
                                 (nrow(tableMF) - borderI + 1):nrow(tableMF))
        
        if (is.null(colnames(tableMF)))
          colnames(tableMF) <- c(1:borderI,
                                 "...",
                                 (ncol(tableMF) - borderI + 1):ncol(tableMF))
        
        tableMF <- as.data.frame(tableMF)
        
      }
      
      rownames(tableMF)[borderI + 1] <- "..."
      
    } else if (dimAbbVl[1]) {
      
      if (is.data.frame(tableMF) && is.null(colnames(tableMF)))
        colnames(tableMF) <- 1:ncol(tableMF)
      
      tableMF <- rbind(tableMF[1:borderI, , drop = FALSE],
                       rep("...", ncol(tableMF)),
                       tableMF[(nrow(tableMF) - borderI + 1):nrow(tableMF), , drop = FALSE])
      
      if (classC == "matrix") {
        
        if (!is.null(rownames(tableMF)) && any(duplicated(rownames(tableMF)))) {
          rownames(tableMF) <- make.names(rownames(tableMF),
                                          unique = TRUE)
        } else if (is.null(rownames(tableMF)))
          rownames(tableMF) <- c(1:borderI,
                                 "...",
                                 (nrow(tableMF) - borderI + 1):nrow(tableMF))
        
        if (is.null(colnames(tableMF)))
          colnames(tableMF) <- 1:ncol(tableMF)
        
        tableMF <- as.data.frame(tableMF)
        
      }
      
      rownames(tableMF)[borderI + 1] <- "..."
      
    } else if (dimAbbVl[2]) {
      
      tableMF <- cbind(tableMF[, 1:borderI, drop = FALSE],
                       ... = rep("...", times = nrow(tableMF)),
                       tableMF[, (ncol(tableMF) - borderI + 1):ncol(tableMF), drop = FALSE])
      
      if (classC == "matrix") {
        
        if (!is.null(rownames(tableMF)) && any(duplicated(rownames(tableMF))))
          rownames(tableMF) <- make.names(rownames(tableMF),
                                          unique = TRUE)
        
        if (is.null(colnames(tableMF)))
          colnames(tableMF) <- c(1:borderI,
                                 "...",
                                 (ncol(tableMF) - borderI + 1):ncol(tableMF))
        
        tableMF <- as.data.frame(tableMF)
        
      }
      
    }
    
  } ## 'data.frame' and 'matrix'
  
  if (classC == "vector") {
    
    tableMF <- tableMF
    
    if (numericL)
      tableMF <- round(tableMF, 3)
    
    if (length(tableMF) > 2 * borderI) {
      
      tableMF <- c(head(tableMF, borderI),
                   "...",
                   tail(tableMF, borderI))
      
      if (!is.null(names(tableMF)))
        names(tableMF)[borderI + 1] <- "..."
      
    }
    
    if (!is.null(names(tableMF))) {
      
      tableMF <- as.data.frame(t(tableMF))
      
      rownames(tableMF) <- ""
      
    } else {
      
      tableMF <- paste(tableMF, collapse = " ")
      
      message(tableMF)
      
      return(invisible(NULL))
      
    }
    
  } ## vector
  
  print(tableMF)
  
} ## tabF


#' Displaying a numerical matrix
#' 
#' Wrapper of the stats::image function used inside the 'view' method
#'
#' @seealso \code{\link{image}}, \code{\link{view}}
#' @examples
#'
#' data(sacurine)
#' imageF(sacurine[['dataMatrix']])
#' 
#' @rdname view
#' @export
imageF <- function(x,
                   mainC = "",
                   subC = "",
                   paletteC = c("heat",
                                "revHeat",
                                "grey",
                                "revGrey",
                                "palette",
                                "ramp")[1],
                   rowAllL = FALSE,
                   rowCexN = 1,
                   rowMarN = 5.1,
                   rowLabC = "",
                   rowTruncI = 0,
                   colAllL = FALSE,
                   colCexN = 1,
                   colMarN = 1.1,
                   colLabC = "",
                   colTruncI = 0,
                   drawScaleL = TRUE,
                   delimitReplicatesL = FALSE,
                   standardizeL = FALSE,
                   fig.pdfC = "interactive") {
  
  if (class(x) != "matrix" || mode(x) != "numeric")
    stop("'x must be a matrix of number for the 'image' plot type", call. = FALSE)
  
  if (delimitReplicatesL && (length(rownames(x)) * length(colnames(x))) == 0)
    stop("Rownames and colnames are required when the delimitReplicatesL argument is TRUE", call. = FALSE)
  
  if (standardizeL) {
    message("Standardization of the columns for plotting.")
    mainC <- paste0(mainC, " (standardized)")
    x <- apply(x, 2, function(y) {
      y <- y - mean(y, na.rm = TRUE)
      sd.n <- sqrt(var(y, na.rm = TRUE))
      if (sd.n < .Machine$double.eps) {
        return(y)
      } else {
        return(y / sd.n)
      }
    })
  }
  
  dimVc <- c("row", "col")
  dimnamesLs <- lapply(dimVc,
                       function(dimC) {
                         .dimNames(x,
                                   which(dimVc == dimC),
                                   ifelse(dimC == "row",
                                          rowTruncI,
                                          colTruncI))
                       })
  names(dimnamesLs) <- dimVc
  
  imageMN <- t(x[nrow(x):1, , drop = FALSE])
  
  paletteVc <- .palette(imageMN, paletteC)
  
  currentParLs <- par()
  for (parC in c("cin", "cra", "csi", "cxy", "din", "page"))
    currentParLs[[parC]] <- NULL
  
  if (fig.pdfC != "interactive") {
    basenameC <- basename(fig.pdfC)
    basenameSplitVc <- unlist(strsplit(basenameC, ".", fixed = TRUE))
    extC <- tail(basenameSplitVc, 1)
    if (extC != "pdf")
      stop("The name of the pdf file should end with a '.pdf' extension.",
           call. = FALSE)
    pdf(fig.pdfC)
  }
  
  par(bg = "white",
      font = 2,
      lwd  = 2)
  
  if (mainC != "" || subC != "") {
    
    omaVn <- rep(0, times = 4)
    
    if (mainC != "")
      omaVn[3] <- 1.6
    
    if (subC != "")
      omaVn[1] <- 1.6
    
    par(oma = omaVn)
    
  }
  
  ## drawing color scale
  
  if (drawScaleL)
    .drawScale(imageMN = imageMN,
               paletteC = paletteC,
               paletteVc = paletteVc,
               colAllL = colAllL,
               colMarN = colMarN)
  
  
  
  ## draw image
  
  .drawImage(imageMN = imageMN,
             mainC = mainC,
             subC = subC,
             paletteVc = paletteVc,
             rowAllL = rowAllL,
             colAllL = colAllL,
             rowCexN = rowCexN,
             colCexN = colCexN,
             rowMarN = rowMarN,
             colMarN = colMarN,
             rowLabC = rowLabC,
             colLabC = colLabC,
             drawScaleL = drawScaleL,
             delimitReplicatesL = delimitReplicatesL,
             dimnamesLs = dimnamesLs)
  
  par(currentParLs)
  
  if (fig.pdfC != "interactive")
    dev.off()
  
}

.dimNames <- function(matMN,
                      dimI,
                      truncateI = 0) {
  
  if (length(dimnames(matMN)[[dimI]]) == 0) {
    
    namesVc <- rep("", times = dim(matMN)[dimI])
    
    namesCharVsNumL <- TRUE
    
  } else {
    
    namesVc <- dimnames(matMN)[[dimI]]
    
    namesCharVsNumL <- suppressWarnings(any(is.na(as.numeric(namesVc))))
    
    if (namesCharVsNumL) {
      
      namesDuplicateVl <- duplicated(namesVc)
      
      namesVc[namesDuplicateVl] <- ""
      
    }
    
    if (truncateI > 0)
      namesVc <- .truncate(namesVc, truncateI)
    
  }
  
  return(list(namesVc = namesVc,
              namesCharVsNumL = namesCharVsNumL))
  
}

.palette <- function(imaMN, typC = c("heat",
                                     "revHeat",
                                     "grey",
                                     "revGrey",
                                     "palette",
                                     "ramp")[1]) {
  
  switch(typC,
         heat = {return(rev(grDevices::rainbow(ceiling(256 * 1.5))[1:256]))},
         revHeat = {return(grDevices::rainbow(ceiling(256 * 1.5))[1:256])},
         grey = {return(grDevices::grey((0:255) / 256))},
         revGrey = {return(rev(grDevices::grey((0:255) / 256)))},
         palette = {return(seq(from = min(imaMN),
                               to = max(imaMN),
                               by = 1))},
         ramp = {return(grDevices::colorRampPalette(c("blue", "orange", "red"),
                                                    space = "rgb")(256)[1:256])})
  
}


.drawScale <- function(imageMN,
                       paletteC,
                       paletteVc,
                       colAllL,
                       colMarN) {
  
  layout(matrix(c(2, 1),
                nrow = 1,
                ncol = 2,
                byrow = TRUE),
         widths = c(8, 2))
  
  par(mar = c(1.1,
              0.6,
              ifelse(colAllL,
                     yes = colMarN,
                     no = 3.1),
              4.1))
  
  ylimVn <- c(0, 256)
  ybottomVn <- 0:255
  ytopVn <- 1:256
  
  if (paletteC == "palette") {
    
    ylimVn <- c(min(imageMN) - 1, max(imageMN))
    ybottomVn <- seq(from = min(imageMN) - 1, to = max(imageMN) - 1, by = 1)
    ytopVn <- seq(from = min(imageMN), to = max(imageMN))
    
  }
  
  plot(x = 0,
       y = 0,
       font.axis = 2,
       font.lab = 2,
       type = "n",
       xlim = c(0, 1),
       ylim = ylimVn,
       xlab = "",
       ylab = "",
       xaxs = "i",
       yaxs = "i",
       xaxt = "n",
       yaxt = "n")
  
  rect(xleft = 0,
       ybottom = ybottomVn,
       xright = 1,
       ytop = ytopVn,
       col = paletteVc,
       border = NA)
  
  axis(at = .prettyAxis(range(imageMN, na.rm = TRUE), 256)$atVn,
       font = 2,
       font.axis = 2,
       labels = .prettyAxis(range(imageMN, na.rm = TRUE), 256)$labelVn,
       las = 1,
       lwd = 2,
       lwd.ticks = 2,
       side = 4,
       xpd = TRUE)
  
  graphics::arrows(par("usr")[2],
                   par("usr")[4],
                   par("usr")[2],
                   par("usr")[3],
                   code = 0,
                   lwd = 2,
                   xpd = TRUE)
  
  graphics::arrows(par("usr")[1],
                   par("usr")[4],
                   par("usr")[1],
                   par("usr")[3],
                   code = 0,
                   lwd = 2,
                   xpd = TRUE)
  
  graphics::box(lwd = 2)
  
  
}


.prettyAxis <- function(axisValuesVn,
                        opLengthN) {
  
  if (NA %in% axisValuesVn) {
    
    warning("NA in axisValuesVn", call. = FALSE)
    
    axisValuesVn <- as.vector(stats::na.omit(axisValuesVn))
    
  }
  
  if (opLengthN < length(axisValuesVn))
    stop("The length of in vector must be inferior to the length of the length parameter.")
  
  if (length(axisValuesVn) < opLengthN) {
    
    axisValuesVn <- seq(from = min(axisValuesVn), to = max(axisValuesVn), length.out = opLengthN)
    
  }
  
  prettyAxisValues <- pretty(axisValuesVn)
  
  prettyLabelsVn <- prettyAtVn <- c()
  
  for (n in 1:length(prettyAxisValues))
    if (min(axisValuesVn) < prettyAxisValues[n] && prettyAxisValues[n] < max(axisValuesVn)) {
      prettyLabelsVn <- c(prettyLabelsVn, prettyAxisValues[n])
      prettyAtVn <- c(prettyAtVn, which(abs(axisValuesVn - prettyAxisValues[n]) == min(abs(axisValuesVn - prettyAxisValues[n])))[1])
    }
  
  prettyAxisLs <- list(atVn = prettyAtVn,
                       labelVn = prettyLabelsVn)
  
  
  return(prettyAxisLs)
  
}


.drawImage <- function(imageMN,
                       mainC,
                       subC,
                       paletteVc,
                       rowAllL,
                       colAllL,
                       rowCexN,
                       colCexN,
                       rowMarN,
                       colMarN,
                       rowLabC,
                       colLabC,
                       drawScaleL,
                       delimitReplicatesL,
                       dimnamesLs) {
  
  par(mar = c(1.1,
              rowMarN,
              ifelse(colAllL,
                     yes = colMarN,
                     no = 3.1),
              ifelse(drawScaleL, yes = 0.3, no = 1.6)))
  
  graphics::image(x = 1:nrow(imageMN),
                  y = 1:ncol(imageMN),
                  z = imageMN,
                  col = paletteVc,
                  font.axis = 2,
                  font.lab = 2,
                  xaxt = "n",
                  yaxt = "n",
                  xlab = "",
                  ylab = "")
  
  .drawAxis(imageMN = imageMN,
            rowAllL = rowAllL,
            colAllL = colAllL,
            rowCexN = rowCexN,
            colCexN = colCexN,
            dimnamesLs = dimnamesLs)
  
  ## xlab
  
  mtext(font = 2,
        line = colMarN - 1,
        side = 3,
        text = colLabC)
  
  ## ylab
  
  mtext(font = 2,
        line = rowMarN - 1,
        side = 2,
        text = rowLabC)
  
  ## additional lines
  
  if (delimitReplicatesL) {
    for (dimI in 1:2) {
      dimAllL <- ifelse(dimI == 1, rowAllL, colAllL)
      if (dimAllL && dimnamesLs[[dimI]][["namesCharVsNumL"]]) {
        if (dimI == 1)  {
          delimN <- ncol(imageMN) - which(dimnamesLs[[dimI]][["namesVc"]] != "") + 1.5
          abline(h = delimN, lwd = 1)
        } else {
          delimN <- which(dimnamesLs[[dimI]][["namesVc"]] != "") - 0.5
          abline(v = delimN, lwd = 1)
        }
      }
    }
  }
  
  
  ## border
  
  graphics::box(lwd = 2)
  
  ## arrows at the end of the axes
  
  graphics::arrows(par("usr")[1],
                   par("usr")[4],
                   par("usr")[1],
                   par("usr")[3],
                   length = 0.1,
                   lwd = 2,
                   xpd = TRUE)
  
  graphics::arrows(par("usr")[1],
                   par("usr")[4],
                   par("usr")[2],
                   par("usr")[4],
                   length = 0.1,
                   lwd = 2,
                   xpd = TRUE)
  
  ## writing figure title
  
  if (mainC != "")
    mtext(adj = 0.5,
          cex = 1.2,
          font = 2,
          line = 0.1,
          outer = TRUE,
          side = 3,
          text = mainC)
  
  ## writing figure caption
  
  if (subC != "")
    mtext(adj = 1,
          cex = 0.9,
          font = 3,
          line = 0.6,
          outer = TRUE,
          side = 1,
          text = paste0(subC, " "))
  
}


.drawAxis <- function(imageMN,
                      rowAllL,
                      colAllL,
                      rowCexN,
                      colCexN,
                      dimnamesLs) {
  
  for (dimI in 1:2) {
    
    dimAllL <- ifelse(dimI == 1, rowAllL, colAllL)
    dimCexN <- ifelse(dimI == 1, rowCexN, colCexN)
    labelAtVn <- 1:dim(imageMN)[dimI]
    
    dimNamesVc <- dimnamesLs[[dimI]][["namesVc"]]
    dimNamesCharVsNumL <- dimnamesLs[[dimI]][["namesCharVsNumL"]]
    
    if (!dimAllL) {
      
      labelIndiceVi <- c(1, dim(imageMN)[3 - dimI])
      # labelIndiceVi <- c(1, dim(imageMN)[dimI])
      
      labelNameVc <- rep("", times = 2)
      
      if (suppressWarnings(any(is.na(as.numeric(c(dimNamesVc[1],
                                                  tail(dimNamesVc , 1))))))) {
        
        labelNameVc <- c(dimNamesVc[1], tail(dimNamesVc, 1))
        
      } else {
        
        labelNameVc <- round(as.numeric(c(dimNamesVc[1], tail(dimNamesVc, 1))),
                             digits = 1)
        
      }
      
      labelIndiceVc <- paste0(rep("[", 2),
                              labelIndiceVi,
                              rep("] ", 2))
      
      labelVc <- paste0(labelIndiceVc,
                        rep("\n", times = 2),
                        labelNameVc)
      
      atVn <- c(1, dim(imageMN)[3 - dimI])
      
      if (dimI == 1) {
        labelIndiceVc <- rev(labelIndiceVc)
        labelVc <- rev(labelVc)
        axis(side = dimI + 1,
             at = atVn,
             font = 2,
             hadj = 0,
             labels = labelVc,
             las = 2,
             line = 3,
             lty = "blank",
             tick = FALSE)
      } else {
        axis(side = dimI + 1,
             at = atVn[1],
             font = 2,
             hadj = 0,
             labels = labelVc[1],
             line = -0.5,
             tick = FALSE)
        axis(side = dimI + 1,
             at = atVn[2],
             font = 2,
             hadj = 1,
             labels = labelVc[2],
             line = -0.5,
             tick = FALSE)
      }
      
    } else {
      # dimAllL
      
      if (dimNamesCharVsNumL) {
        
        labelVc <- dimNamesVc[dimNamesVc != ""]
        
        labelLowIndiceVn <- which(dimNamesVc != "")
        
        labelSpanVn <- diff(c(labelLowIndiceVn, length(dimNamesVc) + 1))
        
        labelAtVn <- labelLowIndiceVn - rep(1, times = length(labelLowIndiceVn)) + labelSpanVn / 2 + rep(0.5, times = length(labelLowIndiceVn))
        
        par(cex = dimCexN)
        
        if (dimI == 1) {
          labelAtVn <- ncol(imageMN) - rev(labelLowIndiceVn) + rep(1, times = length(labelLowIndiceVn)) - rev(labelSpanVn) / 2 + rep(0.5, times = length(labelLowIndiceVn))
          labelVc <- rev(labelVc)
        }
        
        axis(side = dimI + 1,
             at = labelAtVn,
             font = 2,
             labels = labelVc,
             las = 2,
             line = -0.5,
             tick = FALSE)
        
        
      } else {
        
        dimNamesVn <- as.numeric(dimNamesVc)
        
        prettyVn <- pretty(dimNamesVn)
        
        prettyVn <- prettyVn[min(dimNamesVn) <= prettyVn & prettyVn <= max(dimNamesVn)]
        
        indiceVn <- numeric()
        for (k in 1:length(prettyVn)) {
          
          indiceVn[k] <- which(abs(dimNamesVn - prettyVn[k]) == min(abs(dimNamesVn - prettyVn[k])))[1]
          
        }
        
        if (dimI == 1)
          indiceVn <- nrow(imageMN) - rev(indiceVn)
        
        axis(side = dimI + 1,
             at = indiceVn,
             font = 2,
             labels = as.character(prettyVn))
        
      }
      
      par(cex = 1)
      
    }
    
  }
}

.truncate <- function(stringVc, ncharI = 15) {
  
  truncatedVc <- character(length(stringVc))
  
  for (i in 1:length(stringVc)) {
    
    strC <- stringVc[i]
    nchI <- nchar(strC)
    
    if (nchI <= ncharI) {
      truncatedVc[i] <- strC
    } else {
      truncatedVc[i] <- paste0(substr(strC, 1, ncharI - 1), ".")
    }
  }
  
  truncatedVc
  
}


#' fromW4M (deprecated)
#'
#' Creating a ExpressionSet object from the 3 'dataMatrix.tsv',
#' 'sampleMetadata.tsv' and 'variableMetadata.tsv' tabulated files
#'
#' @param dirC Character: directory containing the 3 .tsv files
#' @param namePatternC Character: optional file name pattern common to all three
#' file names (e.g., when you want to distinguish between two sets of files
#' within the same directory)
#' @param fileTableNamesVc Vector of characters: if your file names do not
#' contain the standard 'dataMatrix', 'sampleMetadata', and 'variableMetadata'
#' patterns (e.g. if you use 'profile', 'observation', and 'feature' instead),
#' please indicate them here
#' @param verboseL Logical: should comments be printed?
#' @return ExpressionSet instance
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#' @examples sacSet <- fromW4M(file.path(path.package("ropls"), "extdata"))
#' @rdname fromW4M
#' @export
fromW4M <- function(dirC,
                    namePatternC = "",
                    fileTableNamesVc = c("dataMatrix",
                                         "sampleMetadata",
                                         "variableMetadata"),
                    verboseL = TRUE) {
  
  tabVc <- c("dataMatrix",
             "sampleMetadata",
             "variableMetadata")
  
  if (!file.exists(dirC))
    stop("Directory '", dirC, "' was not found.",
         call. = FALSE)
  
  filVc <- character(length(tabVc))
  names(filVc) <- tabVc
  
  filAllVc <- list.files(dirC,
                         pattern = "^.*\\.tsv$")
  
  ## restricting to files with pattern
  if (namePatternC != "")
    filAllVc <- grep(namePatternC, filAllVc, value = TRUE)
  
  ## restricting to one file for each table
  for (tabC in tabVc) {
    namC <- fileTableNamesVc[tabVc == tabC]
    filC <- grep(namC, filAllVc, value = TRUE)
    if (length(filC) == 0) {
      stop("No file found for the ", tabC, " with ",
           ifelse(namePatternC != "",
                  paste0("'", namC, "' pattern and "), ""),
           "a name including '", namC, "' in the '", dirC,
           "' directory.", call. = FALSE)
    } else if (length(filC) > 1) {
      stop("Several files found for the ", tabC, " with ",
           ifelse(namePatternC != "", paste0("'", namC, "' pattern and "),
                  ""), "a name including '", namC, "' in the '",
           dirC, "' directory.", call. = FALSE)
    } else {
      filVc[tabC] <- filC
      ## R standards for row and column names in matrices and data frames
      .checkRformatF(dirC, filC, verboseL)
    }
  }
  
  ## Loading data
  
  for(tabC in tabVc) {
    
    tabDF <- read.table(file.path(dirC, filVc[tabC]),
                        check.names = FALSE,
                        header = TRUE,
                        row.names = 1,
                        sep = "\t",
                        stringsAsFactors = FALSE)
    switch(tabC,
           dataMatrix = {
             datMN <- as.matrix(tabDF)
           },
           sampleMetadata = {
             samDF <- tabDF
           },
           variableMetadata = {
             varDF <- tabDF
           })
    
  }
  
  chkL <- .checkW4mFormatF(t(datMN), samDF, varDF)
  
  if(chkL) {
    TRUE
  } else
    "Problem with the sample or variable names in the tables (see above)"
  
  eset <- ExpressionSet(assayData = datMN,
                        phenoData = new("AnnotatedDataFrame",
                                        data = samDF),
                        featureData = new("AnnotatedDataFrame",
                                          data = varDF),
                        experimentData = new("MIAME",
                                             title = namePatternC))
  
  validObject(eset)
  
  return(eset)
  
}

# deprecated
.checkRformatF <- function(dirCa, filCa, vrbLa) {
  
  rowVc <- read.table(file.path(dirCa, filCa),
                      check.names = FALSE,
                      header = TRUE,
                      sep = "\t",
                      stringsAsFactors = FALSE)[, 1]
  
  colVc <- unlist(read.table(file.path(dirCa, filCa),
                             check.names = FALSE,
                             nrows = 1,
                             sep = "\t",
                             stringsAsFactors = FALSE))[-1]
  
  if (any(duplicated(rowVc)))
    stop("The following ",
         ifelse(names(filCa) == 'sampleMetadata', 'sample', 'variable'),
         " name(s) is/are duplicated in the ",
         names(filCa),
         ": '",
         paste(rowVc[duplicated(rowVc)], collapse = "', '"), "'",
         call. = FALSE)
  
  if (any(duplicated(colVc)))
    stop("The following ", ifelse(names(filCa) == 'sampleMetadata', 'variable', 'sample'), " name(s) is/are duplicated in the ",
         names(filCa),
         ": '",
         paste(colVc[duplicated(colVc)], collapse = "', '"), "'",
         call. = FALSE)
  
  rowMakVc <- make.names(rowVc, unique = TRUE)
  
  rowDifVl <- rowVc != rowMakVc
  
  if (any(rowDifVl)) {
    rowDifDF <- data.frame(row = 1:length(rowVc),
                           actual = rowVc,
                           preferred = rowMakVc)
    rowDifDF <- rowDifDF[rowDifVl, , drop = FALSE]
    if (vrbLa) {
      warning("The following ",
              ifelse(names(filCa) == 'sampleMetadata', 'sample', 'variable'),
              " name(s) of the ",
              names(filCa),
              " is/are not in the standard R format, which may result in errors when loading the data:")
      print(rowDifDF)
    }
  }
  
  colMakVc <- make.names(colVc, unique = TRUE)
  
  colDifVl <- colVc != colMakVc
  
  if (any(colDifVl)) {
    colDifDF <- data.frame(col = 1:length(colVc),
                           actual = colVc,
                           preferred = colMakVc)
    colDifDF <- colDifDF[colDifVl, , drop = FALSE]
    if (vrbLa) {
      warning("The following ",
              ifelse(names(filCa) == 'sampleMetadata', 'variable', 'sample'),
              " name(s) of the ",
              names(filCa),
              " is/are not in the standard R format, which may result in errors when loading the data:")
      print(colDifDF)
    }
  }
}

# deprecated
.checkW4mFormatF <- function(datMN, samDF, varDF) {
  
  chkL <- TRUE
  
  if(!identical(rownames(datMN), rownames(samDF))) {
    ## checking sample names
    
    chkL <- FALSE
    
    datSamDifVc <- setdiff(rownames(datMN), rownames(samDF))
    
    if(length(datSamDifVc)) {
      cat("\nThe following samples were found in the dataMatrix column names but not in the sampleMetadata row names:\n", sep="")
      print(cbind.data.frame(col = as.numeric(sapply(datSamDifVc, function(samC) which(rownames(datMN) == samC))),
                             name = datSamDifVc))
    }
    
    samDatDifVc <- setdiff(rownames(samDF), rownames(datMN))
    
    if(length(samDatDifVc)) {
      cat("\n\nThe following samples were found in the sampleMetadata row names but not in the dataMatrix column names:\n", sep="")
      print(cbind.data.frame(row = as.numeric(sapply(samDatDifVc, function(samC) which(rownames(samDF) == samC))),
                             name = samDatDifVc))
    }
    
    if(nrow(datMN) != nrow(samDF)) {
      cat("\n\nThe dataMatrix has ", nrow(datMN), " columns (ie samples) whereas the sampleMetadata has ", nrow(samDF), " rows\n", sep="")
    } else if(identical(gsub("^X", "", rownames(datMN)), rownames(samDF))) {
      cat("\n\nThe dataMatrix column names start with an 'X' but not the sampleMetadata row names\n", sep="")
    } else if(identical(gsub("^X", "", rownames(samDF)), rownames(datMN))) {
      cat("\n\nThe sampleMetadata row names start with an 'X' but not the dataMatrix column names\n", sep="")
    } else if(identical(sort(rownames(datMN)), sort(rownames(samDF)))) {
      cat("\n\nThe dataMatrix column names and the sampleMetadata row names are not in the same order:\n", sep="")
      print(cbind.data.frame(indice = 1:nrow(datMN),
                             dataMatrix_columnnames=rownames(datMN),
                             sampleMetadata_rownames=rownames(samDF))[rownames(datMN) != rownames(samDF), , drop = FALSE])
    } else {
      cat("\n\nThe dataMatrix column names and the sampleMetadata row names are not identical:\n", sep="")
      print(cbind.data.frame(indice = 1:nrow(datMN),
                             dataMatrix_columnnames=rownames(datMN),
                             sampleMetadata_rownames=rownames(samDF))[rownames(datMN) != rownames(samDF), , drop = FALSE])
    }
    
  }
  
  if(!identical(colnames(datMN), rownames(varDF))) {
    ## checking variable names
    
    chkL <- FALSE
    
    datVarDifVc <- setdiff(colnames(datMN), rownames(varDF))
    
    if(length(datVarDifVc)) {
      cat("\nThe following variables were found in the dataMatrix row names but not in the variableMetadata row names:\n", sep="")
      print(cbind.data.frame(row = as.numeric(sapply(datVarDifVc, function(varC) which(colnames(datMN) == varC))),
                             name = datVarDifVc))
      
    }
    
    varDatDifVc <- setdiff(rownames(varDF), colnames(datMN))
    
    if(length(varDatDifVc)) {
      cat("\n\nThe following variables were found in the variableMetadata row names but not in the dataMatrix row names:\n", sep="")
      print(cbind.data.frame(row = as.numeric(sapply(varDatDifVc, function(varC) which(rownames(varDF) == varC))),
                             name = varDatDifVc))
    }
    
    if(ncol(datMN) != nrow(varDF)) {
      cat("\n\nThe dataMatrix has ", nrow(datMN), " rows (ie variables) whereas the variableMetadata has ", nrow(varDF), " rows\n", sep="")
    } else if(identical(sort(colnames(datMN)), sort(rownames(varDF)))) {
      cat("\n\nThe dataMatrix row names and the variableMetadata row names are not in the same order:\n", sep="")
      print(cbind.data.frame(row = 1:ncol(datMN),
                             dataMatrix_rownames=colnames(datMN),
                             variableMetadata_rownames=rownames(varDF))[colnames(datMN) != rownames(varDF), , drop = FALSE])
    } else {
      cat("\n\nThe dataMatrix row names and the variableMetadata row names are not identical:\n", sep="")
      print(cbind.data.frame(row = 1:ncol(datMN),
                             dataMatrix_rownames=colnames(datMN),
                             variableMetadata_rownames=rownames(varDF))[colnames(datMN) != rownames(varDF), , drop = FALSE])
    }
  }
  
  return(chkL)
  
}

