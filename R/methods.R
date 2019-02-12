#' @rdname opls
#' @export
setMethod("opls", signature(x = "ExpressionSet"),
          function(x, y = NULL, ...) {

              if(is.null(y)) {
                  opl <- opls(t(exprs(x)), ...)
              } else {
                  if(!is.character(y)) {
                      stop("'y' must be a character when the 'opls' method is applied to an 'ExpressionSet' instance")
                  } else {
                      if(!(y %in% colnames(pData(x)))) {
                          stop("'y' must be the name of a column of the sampleMetadata slot of the 'ExpressionSet' instance")
                      } else {
                          rspFcVcn <- pData(x)[, y]
                          opl <- ropls::opls(t(exprs(x)), rspFcVcn, ...)
                      }
                  }
              }
              
              modC <- opl@typeC
              sumDF <- getSummaryDF(opl)
              scoreMN <- getScoreMN(opl)
              loadingMN <- getLoadingMN(opl)
              
              vipVn <- coeMN <- orthoScoreMN <- orthoLoadingMN <- orthoVipVn <- NULL
              
              if(grepl("PLS", modC)) {
                
                vipVn <- getVipVn(opl)
                coeMN <- coef(opl)
                
                if(grepl("OPLS", modC)) {
                  orthoScoreMN <- getScoreMN(opl, orthoL = TRUE)
                  orthoLoadingMN <- getLoadingMN(opl, orthoL = TRUE)
                  orthoVipVn <- getVipVn(opl, orthoL = TRUE)
                }
                
              }
              
              rspModC <- gsub("-", "", modC)
              if(rspModC != "PCA")
                rspModC <- paste0(make.names(y), "_", rspModC)

              if(sumDF[, "pre"] + sumDF[, "ort"] < 2) {

                tCompMN <- scoreMN
                pCompMN <- loadingMN

              } else {

                if(sumDF[, "ort"] > 0) {
                 tCompMN <- cbind(scoreMN[, 1], orthoScoreMN[, 1])
                  pCompMN <- cbind(loadingMN[, 1], orthoLoadingMN[, 1])
                  colnames(pCompMN) <- colnames(tCompMN) <- c("h1", "o1")
                } else {
                  tCompMN <- scoreMN[, c(1, 2), drop = FALSE]
                  pCompMN <- loadingMN[, c(1, 2), drop = FALSE]
                }

              }

              ## x-scores and prediction

              colnames(tCompMN) <- paste0(rspModC, "_XSCOR-", colnames(tCompMN))
              tCompDF <- as.data.frame(tCompMN)[sampleNames(x), , drop = FALSE]

              if(modC != "PCA") {

                fitMCN <- fitted(opl)

                fitDF <- as.data.frame(fitMCN, stringsAsFactors = FALSE)[sampleNames(x), , drop = FALSE]
                
                colnames(fitDF) <- paste0(rspModC,
                                           "_fitted")
                
                tCompDF <- cbind.data.frame(tCompDF, fitDF)
                
              }

              pData(x) <- Biobase::combine(pData(x), tCompDF)

              ## x-loadings and VIP

              colnames(pCompMN) <- paste0(rspModC, "_XLOAD-", colnames(pCompMN))
              if(!is.null(vipVn)) {
                pCompMN <- cbind(pCompMN, vipVn)
                colnames(pCompMN)[ncol(pCompMN)] <- paste0(rspModC,
                                                           "_VIP",
                                                           ifelse(!is.null(orthoVipVn),
                                                                  "_pred",
                                                                  ""))
                if(!is.null(orthoVipVn)) {
                  pCompMN <- cbind(pCompMN, orthoVipVn)
                  colnames(pCompMN)[ncol(pCompMN)] <- paste0(rspModC,
                                                             "_VIP_ortho")
                }
              }
              if(!is.null(coeMN)) {
                pCompMN <- cbind(pCompMN, coeMN)
                if(ncol(coeMN) == 1)
                  colnames(pCompMN)[ncol(pCompMN)] <- paste0(rspModC, "_coef")
                else
                  colnames(pCompMN)[(ncol(pCompMN) - ncol(coeMN) + 1):ncol(pCompMN)] <- paste0(rspModC, "_", colnames(coeMN), "-COEFF")
              }
              pCompDF <- as.data.frame(pCompMN)[featureNames(x), , drop = FALSE]
              fData(x) <- Biobase::combine(fData(x), pCompDF)

              opl@eset <- x
              
              opl

          })

#' @rdname opls
#' @export
setMethod("opls", signature(x = "data.frame"),
          function(x, ...) {
              if(!all(sapply(x, data.class) == "numeric")) {
                  stop("'x' data frame must contain columns of 'numeric' vectors only", call. = FALSE)
              } else
                  x <- as.matrix(x)
              opl <- opls(x, ...)
              opl
          })

####    opls    ####


#' PCA, PLS(-DA), and OPLS(-DA)
#'
#' PCA, PLS, and OPLS regression, classification, and cross-validation with the
#' NIPALS algorithm
#'
#' @name opls
#' @rdname opls
#' @aliases opls opls,ExpressionSet-method opls,data.frame-method
#' opls,matrix-method
#' @docType methods
#' @param x Numerical data frame or matrix (observations x variables; NAs are
#' allowed); or ExpressionSet object with non empty exprs, for PCA, and
#' phenoData@data, for (O)PLS(-DA), slots
#' @param y Response to be modelled: Either 1) 'NULL' for PCA (default) or 2) a
#' numerical vector (same length as 'x' row number) for single response (O)PLS,
#' or 3) a numerical matrix (same row number as 'x') for multiple response PLS,
#' 4) a factor (same length as 'x' row number) for (O)PLS-DA, or 5) a character
#' indicating the name of the column of the phenoData@data to be used, when x
#' is an ExpressionSet object. Note that, for convenience, character vectors
#' are also accepted for (O)PLS-DA as well as single column numerical (resp.
#' character) matrix for (O)PLS (respectively (O)PLS-DA). NAs are allowed in
#' numeric responses.
#' @param predI Integer: number of components (predictive componenents in case
#' of PLS and OPLS) to extract; for OPLS, predI is (automatically) set to 1; if
#' set to NA [default], autofit is performed: a maximum of 10 components are
#' extracted until (i) PCA case: the variance is less than the mean variance of
#' all components (note that this rule requires all components to be computed
#' and can be quite time-consuming for large datasets) or (ii) PLS case: either
#' R2Y of the component is < 0.01 (N4 rule) or Q2Y is < 0 (for more than 100
#' observations) or 0.05 otherwise (R1 rule)
#' @param orthoI Integer: number of orthogonal components (for OPLS only); when
#' set to 0 [default], PLS will be performed; otherwise OPLS will be peformed;
#' when set to NA, OPLS is performed and the number of orthogonal components is
#' automatically computed by using the cross-validation (with a maximum of 9
#' orthogonal components).
#' @param algoC Default algorithm is 'svd' for PCA (in case of no missing
#' values in 'x'; 'nipals' otherwise) and 'nipals' for PLS and OPLS; when
#' asking to use 'svd' for PCA on an 'x' matrix containing missing values, NAs
#' are set to half the minimum of non-missing values and a warning is generated
#' @param crossvalI Integer: number of cross-validation segments (default is
#' 7); The number of samples (rows of 'x') must be at least >= crossvalI
#' @param log10L Should the 'x' matrix be log10 transformed? Zeros are set to 1
#' prior to transformation
#' @param permI Integer: number of random permutations of response labels to
#' estimate R2Y and Q2Y significance by permutation testing [default is 20 for
#' single response models (without train/test partition), and 0 otherwise]
#' @param scaleC Character: either no centering nor scaling ('none'),
#' mean-centering only ('center'), mean-centering and pareto scaling
#' ('pareto'), or mean-centering and unit variance scaling ('standard')
#' [default]
#' @param subset Integer vector: indices of the observations to be used for
#' training (in a classification scheme); use NULL [default] for no partition
#' of the dataset; use 'odd' for a partition of the dataset in two equal sizes
#' (with respect to the classes proportions)
#' @param printL Logical: Should informations regarding the data set and the
#' model be printed? [default = TRUE]
#' @param plotL Logical: Should the 'summary' plot be displayed? [default =
#' TRUE]
#' @param .sinkC Character: Name of the file for R output diversion [default =
#' NULL: no diversion]; Diversion of messages is required for the integration
#' into Galaxy
#' @param ... Currently not used.
#' @return An S4 object of class 'opls' containing the following slots:
#' \itemize{
#' \item typeC Character: model type (PCA, PLS, PLS-DA, OPLS, or OPLS-DA)
#' \item descriptionMC Character matrix: Description of the data set (number
#' of samples, variables, etc.)
#' \item modelDF Data frame with the model overview (number of components, R2X, R2X(cum), R2Y, R2Y(cum), Q2, Q2(cum),
#' significance, iterations)
#' \item summaryDF Data frame with the model summary (cumulated R2X, R2Y and Q2); RMSEE is the square root of the mean
#' error between the actual and the predicted responses
#' \item subsetVi Integer vector: Indices of observations in the training data set
#' \item pcaVarVn PCA: Numerical vector of variances of length: predI
#' \item vipVn PLS(-DA): Numerical vector of Variable Importance in Projection; OPLS(-DA): Numerical vector of Variable Importance for Prediction (VIP4,p from Galindo-Prieto et al, 2014)
#' \item orthoVipVn OPLS(-DA): Numerical vector of Variable Importance for Orthogonal Modeling (VIP4,o from Galindo-Prieto et al, 2014)
#' \item xMeanVn Numerical vector: variable means of the 'x' matrix
#' \item xSdVn Numerical vector: variable standard deviations of the 'x' matrix
#' \item yMeanVn (O)PLS: Numerical vector: variable means of the 'y' response (transformed into a dummy matrix in case it is of 'character' mode initially)
#' \item ySdVn (O)PLS: Numerical vector: variable standard deviations of the 'y' response (transformed into a dummy matrix in case it is of 'character' mode initially)
#' \item xZeroVarVi Numerical vector: indices of variables with variance < 2.22e-16 which were excluded from 'x' before building the model
#' \item scoreMN Numerical matrix of x scores (T; dimensions: nrow(x) x predI) X = TP' + E; Y = TC' + F
#' \item loadingMN Numerical matrix of x loadings (P; dimensions: ncol(x) x predI) X = TP' + E
#' \item weightMN (O)PLS: Numerical matrix of x weights (W; same dimensions as loadingMN)
#' \item orthoScoreMN OPLS: Numerical matrix of orthogonal scores (Tortho; dimensions: nrow(x) x number of orthogonal components)
#' \item orthoLoadingMN OPLS: Numerical matrix of orthogonal loadings (Portho; dimensions: ncol(x) x number of orthogonal components)
#' \item orthoWeightMN OPLS: Numerical matrix of orthogonal weights (same dimensions as orthoLoadingMN)
#' \item cMN (O)PLS: Numerical matrix of Y weights (C; dimensions: number of responses or number of classes in case of qualitative response) x number of predictive components; Y = TC' + F
#' \item coMN) (O)PLS: Numerical matrix of Y orthogonal weights; dimensions: number of responses or number of classes in case of qualitative response with more than 2 classes x number of orthogonal components
#' \item uMN (O)PLS: Numerical matrix of Y scores (U; same dimensions as scoreMN); Y = UC' + G
#' \item weightStarMN Numerical matrix of projections (W*; same dimensions as loadingMN); whereas columns of weightMN are derived from successively deflated matrices, columns of weightStarMN relate to the original 'x' matrix: T = XW*; W*=W(P'W)inv
#' \item suppLs List of additional objects to be used internally by the 'print', 'plot', and 'predict' methods
#' }
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#' @references Eriksson et al. (2006). Multi- and Megarvariate Data Analysis.
#' Umetrics Academy.  Rosipal and Kramer (2006). Overview and recent advances
#' in partial least squares Tenenhaus (1990). La regression PLS : theorie et
#' pratique. Technip.  Wehrens (2011). Chemometrics with R. Springer.  Wold et
#' al. (2001). PLS-regression: a basic tool of chemometrics
#' @examples
#'
#' #### PCA
#'
#' data(foods) ## see Eriksson et al. (2001); presence of 3 missing values (NA)
#' head(foods)
#' foodMN <- as.matrix(foods[, colnames(foods) != "Country"])
#' rownames(foodMN) <- foods[, "Country"]
#' head(foodMN)
#' foo.pca <- opls(foodMN)
#'
#' #### PLS with a single response
#'
#' data(cornell) ## see Tenenhaus, 1998
#' head(cornell)
#' cornell.pls <- opls(as.matrix(cornell[, grep("x", colnames(cornell))]),
#'                     cornell[, "y"])
#'
#' ## Complementary graphics
#'
#' plot(cornell.pls, typeVc = c("outlier", "predict-train", "xy-score", "xy-weight"))
#'
#' #### PLS with multiple (quantitative) responses
#'
#' data(lowarp) ## see Eriksson et al. (2001); presence of NAs
#' head(lowarp)
#' lowarp.pls <- opls(as.matrix(lowarp[, c("glas", "crtp", "mica", "amtp")]),
#'                    as.matrix(lowarp[, grepl("^wrp", colnames(lowarp)) |
#'                                       grepl("^st", colnames(lowarp))]))
#'
#' #### PLS-DA
#'
#' data(sacurine)
#' attach(sacurine)
#' sacurine.plsda <- opls(dataMatrix, sampleMetadata[, "gender"])
#'
#' #### OPLS-DA
#'
#' sacurine.oplsda <- opls(dataMatrix, sampleMetadata[, "gender"], predI = 1, orthoI = NA)
#'
#' #### Application to an ExpressionSet
#' 
#' exprs <- t(dataMatrix)
#' sacSet <- Biobase::ExpressionSet(assayData = exprs, 
#'                                  phenoData = new("AnnotatedDataFrame", 
#'                                                  data = sampleMetadata), 
#'                                  featureData = new("AnnotatedDataFrame", 
#'                                                    data = variableMetadata),
#'                                  experimentData = new("MIAME", 
#'                                                       title = "sacurine"))
#'                                                       
#' sacPlsda <- opls(sacSet, "gender")
#' sacSet <- getEset(sacPlsda)
#' head(Biobase::pData(sacSet))
#' head(Biobase::fData(sacSet))
#' 
#' detach(sacurine)
#'
#' @rdname opls
#' @export
setMethod("opls", signature(x = "matrix"),
          function(x,
                   y = NULL,
                   predI = NA,
                   orthoI = 0,
                   
                   algoC = c("default", "nipals", "svd")[1],
                   crossvalI = 7,
                   log10L = FALSE,
                   permI = 20,
                   scaleC = c("none", "center", "pareto", "standard")[4],
                   subset = NULL,
                   
                   printL = TRUE,
                   plotL = TRUE,
                   
                   .sinkC = NULL,
                   ...) {
            
            if(!is.null(.sinkC)) ##  Diversion of messages is required for the integration into Galaxy
              sink(.sinkC, append = TRUE)
            
            ## Checking arguments
            
            ## x -> xMN
            
            if(is.data.frame(x)) {
              if(!all(sapply(x, data.class) == "numeric")) {
                stop("'x' data frame must contain columns of 'numeric' vectors only", call. = FALSE)
              } else
                x <- as.matrix(x)
            } else if(is.matrix(x)) {
              if(mode(x) != "numeric")
                stop("'x' matrix must be of 'numeric' mode", call. = FALSE)
            } else
              stop("'x' must be either a data.frame or a matrix", call. = FALSE)
            
            if(any(apply(x, 2, function(colVn) all(is.na(colVn)))))
              stop("'x' contains columns with 'NA' only", call. = FALSE)
            
            if(any(is.infinite(x)))
              stop("'x' contains infinite values", call. = FALSE)
            
            if(any(is.nan(x)))
              stop("'x' contains NaN values", call. = FALSE)
            
            xMN <- x
            
            ## y -> yMCN
            
            yLevelVc <- NULL
            
            if(!is.null(y)) {
              
              if(is.vector(y)) {
                
                if(!(mode(y) %in% c("character", "numeric")))
                  stop("'y' vector must be of 'character' or 'numeric' type", call. = FALSE)
                
                if(length(y) != nrow(xMN))
                  stop("'y' vector length must be equal to the number of rows of 'x'", call. = FALSE)
                
                yMCN <- matrix(y, ncol = 1)
                
              } else if(is.factor(y)) {
                
                if(length(y) != nrow(xMN))
                  stop("'y' factor length must be equal to the number of rows of 'x'", call. = FALSE)
                
                yLevelVc <- levels(y)
                
                yMCN <- matrix(as.character(y), ncol = 1)
                
              } else if(is.matrix(y)) {
                
                if(!(mode(y) %in% c("character", "numeric")))
                  stop("'y' matrix must be of 'character' or 'numeric' type", call. = FALSE)
                
                if(nrow(y) != nrow(xMN))
                  stop("'x' and 'y' must have the same number of rows", call. = FALSE)
                
                if(ncol(y) > 1 && mode(y) != "numeric")
                  stop("Multiple response 'y' matrices must be of numeric type", call. = FALSE)
                
                yMCN <- y
                
              } else
                stop("'y' must be either a vector, a factor, or a matrix", call. = FALSE)
              
            } else
              yMCN <- NULL
            
            
            ## NA in Y only possible for multi-response regression (i.e., Y is a numeric matrix)
            
            if(!is.null(yMCN) &&
               ncol(yMCN) == 1 &&
               any(is.na(drop(yMCN))))
              stop("In case of single response modeling, 'y' must not contain missing values", call. = FALSE)
            
            if(!is.logical(log10L))
              stop("'log10L' must be a logical", call. = FALSE)
            
            if(permI < 0 || (permI - floor(permI)) > 1e-10)
              stop("'permI' must be an integer", call. = FALSE)
            
            if(permI > 0 && (is.null(yMCN) || ncol(yMCN) > 1)) {
              ## warning("Permutation testing available for single response (O)PLS(-DA) models only", call. = FALSE)
              permI <- 0
            }
            
            if(permI > 0 && !is.null(subset)) {
              permI <- 0
              warning("'permI' set to 0 because train/test partition is selected", call. = FALSE)
            }
            
            if(!(algoC %in% c('default', 'nipals', 'svd')))
              stop("'algoC' must be either 'default', 'nipals', or 'svd'", call. = FALSE)
            
            if(algoC == "default")
              algoC <- ifelse(is.null(yMCN) && !any(is.na(c(xMN))), "svd", "nipals")
            
            if(!is.null(yMCN) && algoC != "nipals")
              stop("'nipals' algorithm must be used for (O)PLS(-DA)", call. = FALSE)
            
            if((is.na(orthoI) || orthoI > 0) && is.null(yMCN))
              stop("'y' response cannot be NULL for OPLS(-DA) modeling", call. = FALSE)
            
            if(!is.null(yMCN)) {
              if(is.na(orthoI) || orthoI > 0) {
                if(ncol(yMCN) > 1) {
                  stop("OPLS regression only available for a single 'y' response", call. = FALSE)
                } else if(mode(yMCN) == "character" && length(unique(drop(yMCN))) > 2)
                  stop("OPLS-DA only available for binary classification (use PLS-DA for multiple classes)", call. = FALSE)
              }
            }
            
            if(is.na(orthoI) || orthoI > 0)
              if(is.na(predI) || predI > 1) {
                predI <- 1
                warning("OPLS: number of predictive components ('predI' argument) set to 1", call. = FALSE)
              }
            
            if(!is.na(predI) && !is.na(orthoI) && ((predI + orthoI) > min(dim(xMN))))
              stop("The sum of 'predI' (", predI, ") and 'orthoI' (", orthoI, ") exceeds the minimum dimension of the 'x' data matrix (", min(dim(xMN)), ")" , call. = FALSE)
            
            if(!(length(scaleC) == 1 && scaleC %in% c('none', 'center', 'pareto', 'standard')))
              stop("'scaleC' must be either 'none', 'center', 'pareto', or 'standard'", call. = FALSE)
            
            if(!is.null(subset) && (is.null(yMCN) || ncol(yMCN) > 1))
              stop("train/test partition with 'subset' only available for (O)PLS(-DA) models of a single 'y' response", call. = FALSE)
            
            if(!is.null(subset) &&
               !(mode(subset) == 'character' && subset == 'odd') &&
               !all(subset %in% 1:nrow(xMN)))
              stop("'subset' must be either set to 'odd' or an integer vector of 'x' row numbers", call. = FALSE)
            
            if(crossvalI > nrow(xMN))
              stop("'crossvalI' must be less than the row number of 'x'", call. = FALSE)
            
            
            ## Constants
            
            epsN <- .Machine[["double.eps"]] ## [1] 2.22e-16
            
            
            ## Character to numeric convertion function (for classification)
            
            if(!is.null(yMCN) && mode(yMCN) == "character") {
              
              if(!is.null(yLevelVc)) {
                claVc <- yLevelVc
              } else
                claVc <- sort(unique(drop(yMCN)))
              
              if(length(claVc) == 2) {
                ## binary response kept as a single vector for OPLS-DA computations
                .char2numF <- function(inpMCN,
                                       c2nL = TRUE) {
                  
                  if(c2nL) {
                    
                    outMCN <- inpMCN == claVc[2]
                    mode(outMCN) <- "numeric"
                    
                  } else {
                    
                    outMCN <- matrix(claVc[as.numeric(inpMCN > 0.5) + 1],
                                     ncol = 1,
                                     dimnames = dimnames(inpMCN))
                    
                  }
                  
                  return(outMCN)
                  
                }
                
              } else
                .char2numF <- function(inpMCN,
                                       c2nL = TRUE) {
                  
                  if(c2nL) {
                    
                    outMCN  <- t(sapply(drop(inpMCN),
                                        function(claC) as.numeric(claVc == claC)))
                    colnames(outMCN) <- claVc
                    
                    
                  } else {
                    
                    outMCN <- t(t(apply(inpMCN, 1,
                                        function(rowVn) claVc[which(rowVn == max(rowVn))[1]])))
                    colnames(outMCN) <- "y1"
                    
                  }
                  
                  return(outMCN)
                  
                }
              
            } else
              .char2numF <- NULL
            
            
            ####   Computations   ####
            
            
            ## rownames and colnames
            
            if(is.null(rownames(xMN)))
              if(!is.null(yMCN) && !is.null(rownames(yMCN))) {
                rownames(xMN) <- rownames(yMCN)
              } else
                rownames(xMN) <- paste0("s", 1:nrow(xMN))
            if(is.null(colnames(xMN)))
              colnames(xMN) <- paste0("x", 1:ncol(xMN))
            
            if(!is.null(yMCN)) {
              if(is.null(rownames(yMCN)))
                rownames(yMCN) <- rownames(xMN)
              if(is.null(colnames(yMCN)))
                colnames(yMCN) <- paste0("y", 1:ncol(yMCN))
            }
            
            ## Log10 transformation
            
            if(log10L)
              xMN <- .log10F(xMN)
            
            
            ## Test indices
            
            obsIniVi <- 1:nrow(xMN)
            
            if(!is.null(subset)) {
              subsetL <- TRUE
              if(length(subset) == 1 && subset == "odd") {
                
                if(mode(yMCN) == "numeric")
                  subsetVi <- seq(1, nrow(xMN), by = 2)
                else {
                  subsetVi <- integer()
                  for(claC in unique(drop(yMCN)))
                    subsetVi <- c(subsetVi,
                                  which(drop(yMCN) == claC)[seq(1, sum(drop(yMCN) == claC), by = 2)])
                  subsetVi <- sort(subsetVi)
                }
              } else
                subsetVi <- subset
              if(crossvalI > length(subsetVi))
                stop("'crossvalI' must be less than the number of samples in the subset", call. = FALSE)
            } else {
              subsetL <- FALSE
              subsetVi <- numeric()
            }
            
            
            ## Filtering out zero variance variables
            
            xVarIndLs <- list()
            xVarIndLs[[1]] <- 1:nrow(xMN)
            
            if(subsetL) {
              xVarIndLs[[1]] <- subsetVi
            } ## else if(!is.null(yMCN) && ncol(yMCN) == 1 && nrow(xMN) >= 2 * crossvalI)
            ##   for(cvkN in 1:crossvalI)
            ##       xVarIndLs <- c(xVarIndLs, list(setdiff(1:nrow(xMN), cvaOutLs[[cvkN]])))
            
            xVarVarLs <- lapply(xVarIndLs,
                                function(xVarVi) {
                                  apply(xMN[xVarVi, , drop = FALSE],
                                        2,
                                        function(colVn) var(colVn, na.rm = TRUE))
                                })
            
            xZeroVarVi <- integer()
            for(k in 1:length(xVarVarLs))
              xZeroVarVi <- union(xZeroVarVi, which(xVarVarLs[[k]] < epsN))
            
            if(length(xZeroVarVi) > 0) {
              names(xZeroVarVi) <- colnames(xMN)[xZeroVarVi]
              xMN <- xMN[, -xZeroVarVi, drop = FALSE]
              warning("The variance of the ",
                      length(xZeroVarVi),
                      " following variables is less than ",
                      signif(epsN, 2),
                      " in the full or partial (cross-validation) dataset: these variables will be removed:\n",
                      paste(names(xZeroVarVi), collapse = ", "),
                      call. = FALSE)
            }
            
            
            ## Core
            
            opl <- .coreF(xMN = xMN,
                          yMCN = yMCN,
                          orthoI = orthoI,
                          predI = predI,
                          scaleC = scaleC,
                          algoC = algoC,
                          crossvalI = crossvalI,
                          subsetL = subsetL,
                          subsetVi = subsetVi,
                          .char2numF = .char2numF)
            
            opl@suppLs[["y"]] <- y
            
            if(is.null(opl@suppLs[["yMCN"]])) {
              opl@typeC <- "PCA"
            } else {
              if(ncol(opl@suppLs[["yMCN"]]) > 1 || mode(opl@suppLs[["yMCN"]]) == "numeric")
                opl@typeC <- "PLS"
              else
                opl@typeC <- "PLS-DA"
            }
            if(opl@summaryDF[, "ort"] > 0)
              opl@typeC <- paste("O", opl@typeC, sep = "")
            
            opl@xZeroVarVi <- xZeroVarVi
            ## opl@suppLs[["yLevelVc"]] <- yLevelVc
            
            
            ## Permutation testing (Szymanska et al, 2012)
            
            if(permI > 0) {
              
              modSumVc <- colnames(opl@summaryDF)
              
              permMN <- matrix(0,
                               nrow = 1 + permI,
                               ncol = length(modSumVc),
                               dimnames = list(NULL, modSumVc))
              
              perSimVn <- numeric(1 + permI)
              perSimVn[1] <- 1
              
              
              permMN[1, ] <- as.matrix(opl@summaryDF)
              
              for(k in 1:permI) {
                
                yVcn <- drop(opl@suppLs[["yMCN"]])
                if(!subsetL) {
                  yPerVcn <- sample(yVcn)
                } else {
                  yPerVcn <- numeric(nrow(xMN))
                  refVi <- opl@subsetVi
                  tesVi <- setdiff(1:nrow(xMN), refVi)
                  yPerVcn[refVi] <- sample(yVcn[refVi])
                  yPerVcn[tesVi] <- yVcn[tesVi]
                }
                yPerMCN <- matrix(yPerVcn, ncol = 1)
                
                perOpl <- .coreF(xMN = xMN,
                                 yMCN = yPerMCN,
                                 orthoI = opl@summaryDF[, "ort"],
                                 predI = opl@summaryDF[, "pre"],
                                 scaleC = scaleC,
                                 algoC = algoC,
                                 crossvalI = crossvalI,
                                 subsetL = subsetL,
                                 subsetVi = opl@subsetVi,
                                 .char2numF = .char2numF)
                
                permMN[1 + k, ] <- as.matrix(perOpl@summaryDF)
                
                perSimVn[1 + k] <- .similarityF(opl@suppLs[["yMCN"]], yPerMCN,
                                                .char2numF = .char2numF,
                                                charL = mode(opl@suppLs[["yMCN"]]) == "character")
                
              }
              
              permMN <- cbind(permMN, sim = perSimVn)
              
              perPvaVn <- c(pR2Y = (1 + length(which(permMN[-1, "R2Y(cum)"] >= permMN[1, "R2Y(cum)"]))) / (nrow(permMN) - 1),
                            pQ2 = (1 + length(which(permMN[-1, "Q2(cum)"] >= permMN[1, "Q2(cum)"]))) / (nrow(permMN) - 1))
              opl@summaryDF[, "pR2Y"] <- perPvaVn["pR2Y"]
              opl@summaryDF[, "pQ2"] <- perPvaVn["pQ2"]
              
              opl@suppLs[["permMN"]] <- permMN
              
            }
            
            ####   Numerical results   ####
            
            opl@descriptionMC <- rbind(samples = ifelse(!subsetL,
                                                        nrow(xMN),
                                                        length(subsetVi)),
                                       X_variables = ncol(xMN),
                                       near_zero_excluded_X_variables = length(opl@xZeroVarVi))
            
            totN <- length(c(xMN))
            nasN <- sum(is.na(c(xMN)))
            
            if(!is.null(opl@suppLs[["yMCN"]])) {
              
              opl@descriptionMC <- rbind(opl@descriptionMC,
                                         Y_variables = ncol(opl@suppLs[["yMCN"]]))
              totN <- totN + length(c(opl@suppLs[["yMCN"]]))
              nasN <- nasN + sum(is.na(c(opl@suppLs[["yMCN"]])))
              
            }
            
            opl@descriptionMC <- rbind(opl@descriptionMC,
                                       missing_values = paste0(nasN, " (", round(nasN / totN * 100), "%)"))
            
            
            ## Raw summary
            
            opl@suppLs[["topLoadI"]] <- 3
            
            if(ncol(xMN) > opl@suppLs[["topLoadI"]]) {
              xVarVn <- apply(xMN, 2, var)
              names(xVarVn) <- 1:length(xVarVn)
              xVarVn <- sort(xVarVn)
              xVarSorVin <- as.numeric(names(xVarVn[seq(1, length(xVarVn), length = opl@suppLs[["topLoadI"]])]))
              opl@suppLs[["xSubIncVarMN"]] <- xMN[, xVarSorVin, drop = FALSE]
            } else
              opl@suppLs[["xSubIncVarMN"]] <- xMN
            
            if(ncol(xMN) <= 100) {
              
              xCorMN <- cor(xMN, use = "pairwise.complete.obs")
              xCorMN[lower.tri(xCorMN, diag = TRUE)] <- 0
              
              if(ncol(xMN) > opl@suppLs[["topLoadI"]]) {
                
                xCorNexDF <- which(abs(xCorMN) >= sort(abs(xCorMN), decreasing = TRUE)[opl@suppLs[["topLoadI"]] + 1],
                                   arr.ind = TRUE)
                
                xCorDisMN <- matrix(0,
                                    nrow = nrow(xCorNexDF),
                                    ncol = nrow(xCorNexDF),
                                    dimnames = list(colnames(xMN)[xCorNexDF[, "row"]],
                                                    colnames(xMN)[xCorNexDF[, "col"]]))
                
                for(k in 1:nrow(xCorDisMN))
                  xCorDisMN[k, k] <- xCorMN[xCorNexDF[k, "row"], xCorNexDF[k, "col"]]
                
              } else
                xCorDisMN <- xCorMN
              
              opl@suppLs[["xCorMN"]] <- xCorDisMN
              
              rm(xCorDisMN)
              
            }
            
            ## Printing
            
            if(printL) {
              show(opl)
              warnings()
            }
            
            ## Plotting
            
            if(plotL)
              plot(opl, typeVc = "summary")
            
            ## Closing connection
            
            if(!is.null(.sinkC)) ## Used in the Galaxy module
              sink()
            
            ## Returning
            
            return(invisible(opl))
            
          })


####    show    ####

#' Show method for 'opls' objects
#'
#' Displays information about the dataset and the model.
#'
#' @aliases show.opls show,opls-method
#' @param object An S4 object of class \code{opls}, created by the \code{opls}
#' function.
#' @return Invisible.
#' @author Philippe Rinaudo and Etienne Thevenot (CEA)
#' @examples
#'
#' data(sacurine)
#' attach(sacurine)
#' sacurine.plsda <- opls(dataMatrix, sampleMetadata[, "gender"])
#'
#' show(sacurine.plsda)
#'
#' detach(sacurine)
#'
#' @rdname show
#' @export
setMethod("show", "opls",
          function(object) {

    cat(object@typeC, "\n", sep = "")

    cat(object@descriptionMC["samples", ],
        " samples x ",
        object@descriptionMC["X_variables", ],
        " variables",
        ifelse(grepl("PLS", object@typeC),
               paste0(" and ", ncol(object@suppLs[["yMCN"]]), " response", ifelse(ncol(object@suppLs[["yMCN"]]) > 1, "s", "")),
               ""), "\n", sep = "")

    cat(object@suppLs[["scaleC"]], " scaling of predictors",
        ifelse(object@typeC == "PCA",
               "",
               paste0(" and ",
                      ifelse(mode(object@suppLs[["yMCN"]]) == "character" && object@suppLs[["scaleC"]] != "standard",
                             "standard scaling of ",
                             ""),
                      "response(s)")), "\n", sep = "")

    if(substr(object@descriptionMC["missing_values", ], 1, 1) != "0")
        cat(object@descriptionMC["missing_values", ], " NAs\n", sep = "")

    if(substr(object@descriptionMC["near_zero_excluded_X_variables", ], 1, 1) != "0")
        cat(object@descriptionMC["near_zero_excluded_X_variables", ],
            " excluded variables (near zero variance)\n", sep = "")

    optDigN <- options()[["digits"]]
    options(digits = 3)
    print(object@summaryDF)
    options(digits = optDigN)

}) ## show


####    print    ####

#' Print method for 'opls' objects
#'
#' Displays information about the dataset and the model.
#'
#' @aliases print.opls print,opls-method
#' @param x An S4 object of class \code{opls}, created by the \code{opls}
#' function.
#' @param ... Currently not used.
#' @return Invisible.
#'
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#' @examples
#'
#' data(sacurine)
#' attach(sacurine)
#' sacurine.plsda <- opls(dataMatrix, sampleMetadata[, "gender"])
#'
#' print(sacurine.plsda)
#'
#' detach(sacurine)
#'
#' @rdname print
#' @export
setMethod("print", "opls",
          function(x, ...) {

    cat("\n1) Data set:\n", sep = "")

    cat(x@descriptionMC["samples", ],
            " samples x ",
            x@descriptionMC["X_variables", ],
            " variables",
            ifelse(grepl("PLS", x@typeC),
                   paste0(" and ", ncol(x@suppLs[["yMCN"]]), " response", ifelse(ncol(x@suppLs[["yMCN"]]) > 1, "s", "")),
                   ""),
        "\n", sep = "")

    cat(x@descriptionMC["missing_values", ], " NAs\n", sep = "")

    cat(x@descriptionMC["near_zero_excluded_X_variables", ], " excluded variables (near zero variance)\n", sep = "")

    cat(x@suppLs[["scaleC"]], " x", ifelse(x@typeC != "PCA", " and y", ""), " scaling\n", sep = "")

    cat("Summary of the ", x@suppLs[["topLoadI"]], " increasing variance spaced raw variables:\n", sep = "")
    print(summary(x@suppLs[["xSubIncVarMN"]]))

    if(!is.null(x@suppLs[["xCorMN"]])) {
        cat("Correlations between the X-variables:\n")
        print(signif(x@suppLs[["xCorMN"]], 2))
        cat("\n", sep = "")
    }

    cat("\n2) Model: ", x@typeC, "\n", sep = "")

    cat("Correlations between variables and first 2 components:\n", sep = "")

    if(x@summaryDF[, "pre"] + x@summaryDF[, "ort"] < 2) {

        warning("A single component model has been selected by cross-validation", call. = FALSE)

        tCompMN <- x@scoreMN
        pCompMN <- x@loadingMN

    } else {

        if(x@summaryDF[, "ort"] > 0) {
            tCompMN <- cbind(x@scoreMN[, 1], x@orthoScoreMN[, 1])
            pCompMN <- cbind(x@loadingMN[, 1], x@orthoLoadingMN[, 1])
            colnames(pCompMN) <- colnames(tCompMN) <- c("h1", "o1")
        } else {
            tCompMN <- x@scoreMN[, 1:2]
            pCompMN <- x@loadingMN[, 1:2]
        }

    }

    cxtCompMN <- cor(x@suppLs[["xModelMN"]], tCompMN,
                     use = "pairwise.complete.obs")

    if(x@suppLs[["topLoadI"]] * 4 < ncol(x@suppLs[["xModelMN"]])) {

        pexVi <- integer(x@suppLs[["topLoadI"]] * ncol(pCompMN) * 2) ## 'ex'treme values

        for(k in 1:ncol(pCompMN)) {

            pkVn <-  pCompMN[, k]

            pexVi[1:(2 * x@suppLs[["topLoadI"]]) + 2 * x@suppLs[["topLoadI"]] * (k - 1)] <- c(order(pkVn)[1:x@suppLs[["topLoadI"]]],
                                                                                rev(order(pkVn, decreasing = TRUE)[1:x@suppLs[["topLoadI"]]]))

        }

    } else
        pexVi <- 1:ncol(x@suppLs[["xModelMN"]])

    pxtCompMN <- cbind(pCompMN,
                       cxtCompMN)

    if(ncol(pCompMN) == 1) {
        colnames(pxtCompMN)[2] <- paste0("cor_", colnames(pxtCompMN)[2])
    } else
        colnames(pxtCompMN)[3:4] <- paste0("cor_", colnames(pxtCompMN)[3:4])

    topLoadMN <- pxtCompMN

    topLoadMN <- topLoadMN[pexVi, , drop = FALSE]

    if(x@suppLs[["topLoadI"]] * 4 < ncol(x@suppLs[["xModelMN"]]) &&
       ncol(pCompMN) > 1) {

        topLoadMN[(2 * x@suppLs[["topLoadI"]] + 1):(4 * x@suppLs[["topLoadI"]]), c(1, 3)] <- NA
        topLoadMN[1:(2 * x@suppLs[["topLoadI"]]), c(2, 4)] <- NA

    }
    print(signif(topLoadMN, 2))

    message("")

    optDigN <- options()[["digits"]]
    options(digits = 3)
    print(x@modelDF)
    options(digits = optDigN)

})


####    plot    ####

#' Plot Method for (O)PLS(-DA)
#'
#' This function plots values based upon a model trained by \code{opls}.
#'
#' @aliases plot.opls plot,opls-method
#' @param x An S4 object of class \code{opls}, created by the \code{opls}
#' function.
#' @param y Currently not used.
#' @param typeVc Character vector: the following plots are available:
#' 'correlation': Variable correlations with the components, 'outlier':
#' Observation diagnostics (score and orthogonal distances), 'overview': Model
#' overview showing R2Ycum and Q2cum (or 'Variance explained' for PCA),
#' 'permutation': Scatterplot of R2Y and Q2Y actual and simulated models after
#' random permutation of response values; 'predict-train' and 'predict-test':
#' Predicted vs Actual Y for reference and test sets (only if Y has a single
#' column), 'summary' [default]: 4-plot summary showing permutation, overview,
#' outlier, and x-score together, 'x-variance': Spread of raw variables
#' corresp. with min, median, and max variances, 'x-loading': X-loadings (the 6
#' of variables most contributing to loadings are colored in red to facilitate
#' interpretation), 'x-score': X-Scores, 'xy-score': XY-Scores, 'xy-weight':
#' XY-Weights
#' @param parAsColFcVn Optional factor character or numeric vector to be
#' converted into colors for the score plot; default is NA [ie colors will be
#' converted from 'y' in case of (O)PLS(-DA) or will be 'black' for PCA]
#' @param parCexN Numeric: amount by which plotting text should be magnified
#' relative to the default
#' @param parCompVi Integer vector of length 2: indices of the two components
#' to be displayed on the score plot (first two components by default)
#' @param parDevNewL Should the graphics be displayed in a new window? 
#' [default]
#' @param parEllipsesL Should the Mahalanobis ellipses be drawn? If 'NA'
#' [default], ellipses are drawn when either a character parAsColVcn is
#' provided (PCA case), or when 'y' is a character factor ((O)PLS-DA cases).
#' @param parLabVc Optional character vector for the labels of observations on
#' the plot; default is NA [ie row names of 'x', if available, or indices of
#' 'x', otherwise, will be used]
#' @param parPaletteVc Optional character vector of colors to be used in the plots
#' @param parTitleL Should the titles of the plots be printed on the graphics
#' (default = TRUE); It may be convenient to set this argument to FALSE when
#' the user wishes to add specific titles a posteriori
#' @param file.pdfC Figure filename (e.g. in case of batch mode) ending with
#' '.pdf'; default is NULL (no
#' saving; displaying instead)
#' @param .sinkC Character: Name of the file for R output diversion [default =
#' NULL: no diversion]; Diversion of messages is required for the integration
#' into Galaxy
#' @param ... Currently not used.
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#' @examples
#'
#' data(sacurine)
#' attach(sacurine)
#'
#' for(typeC in c("correlation", "outlier", "overview",
#'                "permutation", "predict-train","predict-test",
#'                "summary", "x-loading", "x-score", "x-variance",
#'                "xy-score", "xy-weight")) {
#'
#'     print(typeC)
#'
#'     if(grepl("predict", typeC))
#'         subset <- "odd"
#'     else
#'         subset <- NULL
#'
#'     plsModel <- opls(dataMatrix, sampleMetadata[, "gender"],
#'                      predI = ifelse(typeC != "xy-weight", 1, 2),
#'                      orthoI = ifelse(typeC != "xy-weight", 1, 0),
#'                      permI = ifelse(typeC == "permutation", 10, 0),
#'                      subset = subset,
#'                      printL = FALSE, plotL = FALSE)
#'
#'     plot(plsModel, typeVc = typeC)
#'
#' }
#' 
#' sacPlsda <- opls(dataMatrix, sampleMetadata[, "gender"])
#' plot(sacPlsda, parPaletteVc = c("green4", "magenta"))
#'
#' detach(sacurine)
#'
#' @rdname plot
#' @export
setMethod("plot", signature(x = "opls"),
          function(x,
                   y,
                   typeVc = c("correlation",
                       "outlier",
                       "overview",
                       "permutation",
                       "predict-train",
                       "predict-test",
                       "summary",
                       "x-loading",
                       "x-score",
                       "x-variance",
                       "xy-score",
                       "xy-weight")[7],
                   parAsColFcVn = NA,
                   parCexN = 0.8,
                   parCompVi = c(1, 2),
                   parDevNewL = TRUE,
                   parEllipsesL = NA,
                   parLabVc = NA,
                   parPaletteVc = NA,
                   parTitleL = TRUE,
                   file.pdfC = NULL,

                   .sinkC = NULL,
                   ...) {


    if(!is.null(.sinkC)) ##  Diversion of messages is required for the integration into Galaxy
        sink(.sinkC, append = TRUE)


    if("summary" %in% typeVc) {
        if(!is.null(x@suppLs[["permMN"]]))
            typeVc <- c("permutation",
                        "overview",
                        "outlier",
                        "x-score")
        else
            typeVc <- c("overview",
                        "outlier",
                        "x-score",
                        "x-loading")
    }


    ## Checking arguments

    if(!all(typeVc %in% c('correlation', 'outlier', 'overview', 'permutation', 'predict-train', 'predict-test', 'x-loading', 'x-score', 'x-variance', 'xy-score', 'xy-weight')))
        stop("'typeVc' elements must be either 'correlation', 'outlier', 'overview', 'permutation', 'predict-train', 'predict-test', 'x-loading', 'x-score', 'x-variance', 'xy-score', 'xy-weight'", call. = FALSE)

    if('predict-test' %in% typeVc && length(x@subsetVi) == 0)
        stop("For the 'predict-test' graphic to be generated, 'subset' must not be kept to NULL", call. = FALSE)

    if(!any(is.na(parLabVc))) {
        if(length(x@subsetVi) > 0 && length(parLabVc) != nrow(x@suppLs[["yMCN"]])) {
            stop("When 'subset' is not NULL, 'parLabVc' vector length must be equal to the number of train + test samples (here: ", nrow(x@suppLs[["yMCN"]]), ").", call. = FALSE)
        } else if(length(parLabVc) != nrow(x@scoreMN))
            stop("'parLabVc' vector length must be equal to the number of 'x' rows")
        if(mode(parLabVc) != "character")
            stop("'parLabVc' must be of 'character' type")
    }

    if(!any(is.na(parAsColFcVn))) {
        if(length(x@subsetVi) > 0 && length(parAsColFcVn) != nrow(x@suppLs[["yMCN"]])) {
            stop("When 'subset' is not NULL, 'parAsColFcVn' vector length must be equal to the number of train + test samples (here: ", nrow(x@suppLs[["yMCN"]]), ").", call. = FALSE)
        } else if(length(parAsColFcVn) != nrow(x@scoreMN))
            stop("'parAsColFcVn' vector length must be equal to the number of 'x' rows")
        if(!(mode(parAsColFcVn) %in% c("character", "numeric")))
            stop("'parAsColFcVn' must be of 'character' or 'numeric' type")
        if(is.character(parAsColFcVn)) {
            parAsColFcVn <- factor(parAsColFcVn)
            warning("Character 'parAsColFcVn' set to a factor", call. = FALSE)
        }
    }

    if(is.null(x@suppLs[["permMN"]]) && 'permutation' %in% typeVc)
        stop("'permI' must be > 0 for 'permutation' graphic to be plotted", call. = FALSE)

    if(x@summaryDF[, "ort"] > 0)
        if(parCompVi[1] != 1) {
            parCompVi[1] <- 1
            warning("OPLS: first component to display ('parCompVi' first value) set to 1", call. = FALSE)
        }

    if("xy-weight" %in% typeVc &&
       substr(x@typeC, 1, 3) != "PLS")
       ## (is.null(yMCN) || is.na(x@summaryDF[, "ort"]) || x@summaryDF[, "ort"] > 0))
        stop("'xy-weight graphic can be displayed only for PLS(-DA) models", call. = FALSE)

    if(any(grepl('predict', typeVc)))
       if(is.null(x@suppLs[["yMCN"]]) ||
          ncol(x@suppLs[["yMCN"]]) > 1 ||
          (mode(x@suppLs[["yMCN"]]) == "character" && length(unique(drop(x@suppLs[["yMCN"]]))) > 2))
    ## if(any(grepl('predict', typeVc)) && is.matrix(x@fitted"]]) && ncol(x@fitted"]]) > 1)
        ## if(any(grepl('predict', typeVc)) && (is.null(yMCN) || ncol(yMCN) != 1))
           stop("'predict' graphics available for single response regression or binary classification only", call. = FALSE)

    if(is.na(parEllipsesL)) {
        if((all(is.na(parAsColFcVn)) && grepl("-DA$", x@typeC)) ||
           (!all(is.na(parAsColFcVn)) && is.factor(parAsColFcVn))) {
            parEllipsesL <- TRUE
        } else
            parEllipsesL <- FALSE
        ## if((x@typeC == "PCA" && !all(is.na(parAsColFcVn)) && is.factor(parAsColFcVn)) || ## PCA case
        ##    grepl("-DA$", x@typeC)) { ## (O)PLS-DA cases
        ##     parEllipsesL <- TRUE
        ## } else
        ##     parEllipsesL <- FALSE
    } else if(parEllipsesL && !grepl("-DA$", x@typeC) && (all(is.na(parAsColFcVn)) || !is.factor(parAsColFcVn)))
        stop("Ellipses can be plotted for PCA (or PLS regression) only if the 'parAsColFcVn' is a factor",
             call. = FALSE)

    if(x@summaryDF[, "pre"] + x@summaryDF[, "ort"] < 2) {

        if(!all(typeVc %in% c("permutation", "overview"))) {
            warning("Single component model: only 'overview' and 'permutation' (in case of single response (O)PLS(-DA)) plots available", call. = FALSE)
            typeVc <- "overview"
            if(!is.null(x@suppLs[["permMN"]]))
                typeVc <- c("permutation", typeVc)
        }

        tCompMN <- x@scoreMN
        pCompMN <- x@loadingMN

    } else {

        if(x@summaryDF[, "ort"] > 0) {
            if(parCompVi[2] > x@summaryDF[, "ort"] + 1)
                stop("Selected orthogonal component for plotting (ordinate) exceeds the total number of orthogonal components of the model", call. = FALSE)
            tCompMN <- cbind(x@scoreMN[, 1], x@orthoScoreMN[, parCompVi[2] - 1])
            pCompMN <- cbind(x@loadingMN[, 1], x@orthoLoadingMN[, parCompVi[2] - 1])
            colnames(pCompMN) <- colnames(tCompMN) <- c("h1", paste("o", parCompVi[2] - 1, sep = ""))
        } else {
            if(max(parCompVi) > x@summaryDF[, "pre"])
                stop("Selected component for plotting as ordinate exceeds the total number of predictive components of the model", call. = FALSE)
            tCompMN <- x@scoreMN[, parCompVi, drop = FALSE]
            pCompMN <- x@loadingMN[, parCompVi, drop = FALSE]
        }

    }

    ## if(ncol(tCompMN) > 1) {

    ##     mahInvCovMN <- solve(cov(tCompMN))

    ##     pcaResMN <- cbind(sdsVn = apply(tCompMN,
    ##                           1,
    ##                           function(x) sqrt(t(as.matrix(x)) %*% mahInvCovMN %*% as.matrix(x))),
    ##                       odsVn = apply(x@suppLs[["xModelMN"]] - tcrossprod(tCompMN, pCompMN),
    ##                           1,
    ##                           function(x) sqrt(drop(crossprod(x[complete.cases(x)])))))

    ## } else
    ##     pcaResMN <- NULL

    cxtCompMN <- cor(x@suppLs[["xModelMN"]], tCompMN,
                     use = "pairwise.complete.obs")

    if(!is.null(x@suppLs[["yModelMN"]]))
        cytCompMN <- cor(x@suppLs[["yModelMN"]], tCompMN, use = "pairwise.complete.obs")


    if(x@suppLs[["topLoadI"]] * 4 < ncol(x@suppLs[["xModelMN"]])) {

        pexVi <- integer(x@suppLs[["topLoadI"]] * ncol(pCompMN) * 2) ## 'ex'treme values

        for(k in 1:ncol(pCompMN)) {

            pkVn <-  pCompMN[, k]

            pexVi[1:(2 * x@suppLs[["topLoadI"]]) + 2 * x@suppLs[["topLoadI"]] * (k - 1)] <- c(order(pkVn)[1:x@suppLs[["topLoadI"]]],
                                                                         rev(order(pkVn, decreasing = TRUE)[1:x@suppLs[["topLoadI"]]]))

        }

    } else
        pexVi <- 1:ncol(x@suppLs[["xModelMN"]])

    pxtCompMN <- cbind(pCompMN,
                       cxtCompMN)

    if(ncol(pCompMN) == 1) {
       colnames(pxtCompMN)[2] <- paste0("cor_", colnames(pxtCompMN)[2])
    } else
        colnames(pxtCompMN)[3:4] <- paste0("cor_", colnames(pxtCompMN)[3:4])

    topLoadMN <- pxtCompMN

    topLoadMN <- topLoadMN[pexVi, , drop = FALSE]

    if(x@suppLs[["topLoadI"]] * 4 < ncol(x@suppLs[["xModelMN"]]) &&
       ncol(pCompMN) > 1) {

        topLoadMN[(2 * x@suppLs[["topLoadI"]] + 1):(4 * x@suppLs[["topLoadI"]]), c(1, 3)] <- NA
        topLoadMN[1:(2 * x@suppLs[["topLoadI"]]), c(2, 4)] <- NA

    }

    if(!any(is.na(parPaletteVc))) {

      if(mode(parPaletteVc) != "character")
        stop("'parPaletteVc' should be a vector of mode 'character'",
             call. = FALSE)
      
      palChkVl <- sapply(parPaletteVc,
                         function(colC) {
                           tryCatch(is.matrix(col2rgb(colC)), 
                                    error = function(e) FALSE)
                           
                         }) ## as proposed by Josh O'Brien on stackoverflow
      
      if(any(is.na(palChkVl)) || any(!palChkVl))
        stop("The following elements from 'parPaletteVc' could not be interpreted as colors:\n",
             paste(names(palChkVl)[is.na(palChkVl) | !palChkVl], collapse = ", "),
             call. = FALSE)
      
      if(any(is.na(palChkVl)))
        stop("The following element(s) from 'parPaletteVc' could not be interpreted as colors:\n",
             paste(names(palChkVl)[is.na(palChkVl)], collapse = ", "),
             call. = FALSE)
      
    }
    
    ## Observation and variable names and colors

    ## obsLabVc

    if(!any(is.na(parLabVc))) {
        obsLabVc <- parLabVc
    } else if(!is.null(x@suppLs[["yMCN"]]) && ncol(x@suppLs[["yMCN"]]) == 1) { ## (O)PLS of single response
        obsLabVc <- rownames(x@suppLs[["yMCN"]])
    } else { ## PCA
        if(!is.null(rownames(tCompMN))) {
            obsLabVc <- rownames(tCompMN)
        } else
            obsLabVc <- as.character(1:nrow(tCompMN))
    }

    if(length(x@subsetVi) > 0) {
        ## (O)PLS(-DA) models of a single 'y' response
        tesLabVc <- obsLabVc[-x@subsetVi]
        obsLabVc <- obsLabVc[x@subsetVi]
    } else
        tesLabVc <- ""

    ## obsColVc

    if(!any(is.na(parAsColFcVn))) {
        obsColVc <- .plotColorF(as.vector(parAsColFcVn), parPaletteVc)[["colVc"]]
        obsLegVc <- as.vector(parAsColFcVn)
    } else if(!is.null(x@suppLs[["yMCN"]]) && ncol(x@suppLs[["yMCN"]]) == 1) { ## (O)PLS of single response
        obsColVc <- .plotColorF(c(x@suppLs[["yMCN"]]), parPaletteVc)[["colVc"]]
        obsLegVc <- c(x@suppLs[["yMCN"]])
    } else { ## PCA
        obsColVc <- rep("black", nrow(tCompMN))
        obsLegVc <- NULL
    }

    if(length(x@subsetVi) > 0) {
        ## (O)PLS(-DA) models of a single 'y' response
        tesColVc <- obsColVc[-x@subsetVi]
        obsColVc <- obsColVc[x@subsetVi]
        if(!is.null(obsLegVc)) {
            tesLegVc <- obsLegVc[-x@subsetVi]
            obsLegVc <- obsLegVc[x@subsetVi]
        }
    }


    ## Layout

    if(!parDevNewL && length(typeVc) != 1)
        stop("'typeVc' must be of length 1 when 'parDevNewL' is set to FALSE", call. = FALSE)

    if(parDevNewL) {
        layRowN <- ceiling(sqrt(length(typeVc)))
        if(is.null(file.pdfC))
            dev.new()
        else
            pdf(file.pdfC)
        layout(matrix(1:layRowN^2, byrow = TRUE, nrow = layRowN))
    }

    layL <- !parDevNewL || length(typeVc) > 1


    ## Par

    if(layL) {
        marVn <- c(4.6, 4.1, 2.6, 1.6)
    } else
        marVn <- c(5.1, 4.1, 4.1, 2.1)

    par(font=2, font.axis=2, font.lab=2, lwd=2,
        mar=marVn,
        pch=18)


    ## Graph

    for(ploC in typeVc)
        .plotF(ploC,
               opl = x,
               obsColVc = obsColVc,
               obsLabVc = obsLabVc,
               obsLegVc = obsLegVc,
               layL = layL,
               parCexN = parCexN,
               parEllipsesL = parEllipsesL,
               parPaletteVc = parPaletteVc,
               parTitleL = parTitleL,
               parCompVi = parCompVi,
               typeVc = typeVc,
               tCompMN = tCompMN,
               pCompMN = pCompMN,
               cxtCompMN = cxtCompMN,
               cytCompMN = cytCompMN,
               ## pcaResMN = pcaResMN,
               topLoadMN = topLoadMN,
               pexVi = pexVi,
               tesColVc = tesColVc,
               tesLabVc = tesLabVc,
               tesLegVc = tesLegVc)

    if(layL)
        par(font=1, font.axis=1, font.lab=1, lwd=1,
            mar=c(5.1, 4.1, 4.1, 2.1),
            pch=1)

    if(!is.null(file.pdfC))
        dev.off()


    ## Closing connection

    if(!is.null(.sinkC)) ## Used in the Galaxy module
        sink()


})


####    fitted    ####

#' Fitted method for 'opls' objects
#'
#' Returns predictions of the (O)PLS(-DA) model on the training dataset
#'
#' @aliases fitted.opls fitted,opls-method
#' @param object An S4 object of class \code{opls}, created by the \code{opls}
#' function.
#' @param ... Currently not used.
#' @return Predictions (either a vector, factor, or matrix depending on the y
#' response used for training the model)
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#' @examples
#'
#' data(sacurine)
#' attach(sacurine)
#' sacurine.plsda <- opls(dataMatrix, sampleMetadata[, "gender"])
#'
#' fitted(sacurine.plsda)
#'
#' detach(sacurine)
#'
#' @rdname fitted
#' @export
setMethod("fitted", "opls",
          function(object, ...) {

              if(!is.null(object@suppLs[["yPreMN"]])) {

                  if(mode(object@suppLs[["yMCN"]]) == "character") {

                      yPredMCN <- object@suppLs[[".char2numF"]](object@suppLs[["yPreMN"]],
                                                                c2nL = FALSE)

                      if(is.vector(object@suppLs[["y"]])) {
                          fit <- c(yPredMCN)
                          names(fit) <- rownames(yPredMCN)
                      } else if(is.factor(object@suppLs[["y"]])) {
                          fit <- c(yPredMCN)
                          names(fit) <- rownames(yPredMCN)
                          fit <- factor(fit, levels = levels(object@suppLs[["y"]]))
                      } else if(is.matrix(object@suppLs[["y"]])) {
                          fit <- yPredMCN
                      } else
                          stop() ## this case should not happen

                  } else {

                      yPredMCN <- object@suppLs[["yPreMN"]]

                      if(is.vector(object@suppLs[["y"]])) {
                          fit <- c(yPredMCN)
                          names(fit) <- rownames(yPredMCN)
                      } else if(is.matrix(object@suppLs[["y"]])) {
                          fit <- yPredMCN
                      } else
                          stop() ## this case should not happen

                  }

                  return(fit)

              } else
                  return(NULL)

          }) ## fitted


####    tested    ####

#' Tested method for (O)PLS models
#'
#' Returns predictions of the (O)PLS(-DA) model on the out of the box samples
#' (when a 'subset' of samples has been selected when training the model)
#'
#' @aliases tested tested,opls-method
#' @param object An S4 object of class \code{opls}, created by \code{opls}
#' function.
#' @param ... Currently not used.
#' @return Predictions (either a vector, factor, or matrix depending on the y
#' response used for training the model)
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#' @examples
#'
#' data(sacurine)
#' attach(sacurine)
#'
#' testedorMN <- dataMatrix
#' responseFc <- sampleMetadata[, "gender"]
#'
#' sacurine.plsda <- opls(testedorMN,
#'                        responseFc,
#'                        subset = "odd")
#'
#' trainVi <- getSubsetVi(sacurine.plsda)
#'
#' table(responseFc[trainVi], fitted(sacurine.plsda))
#'
#' detach(sacurine)
#'
#' @rdname tested
#' @export
setMethod("tested", "opls",
          function(object) {

              if(!is.null(object@suppLs[["yTesMN"]])) {

                  if(mode(object@suppLs[["yMCN"]]) == "character") {

                      yTestMCN <- object@suppLs[[".char2numF"]](object@suppLs[["yTesMN"]],
                                                                             c2nL = FALSE)
                      if(is.vector(object@suppLs[["y"]])) {
                          test <- c(yTestMCN)
                          names(test) <- rownames(yTestMCN)
                      } else if(is.factor(object@suppLs[["y"]])) {
                          test <- c(yTestMCN)
                          names(test) <- rownames(yTestMCN)
                          test <- factor(test, levels = levels(object@suppLs[["y"]]))
                      } else if(is.matrix(object@suppLs[["y"]])) {
                          test <- yTestMCN
                      } else
                          stop() ## this case should not happen

                  } else {

                      yTestMCN <- object@suppLs[["yTesMN"]]

                      if(is.vector(object@suppLs[["y"]])) {
                          test <- c(yTestMCN)
                          names(test) <- rownames(yTestMCN)
                      } else if(is.matrix(object@suppLs[["y"]])) {
                          test <- yTestMCN
                      } else
                          stop() ## this case should not happen

                  }

                  return(test)

              } else
                  stop("Test results only available for (O)PLS(-DA) models", call. = FALSE)


          })


####    coef    ####

#' Coefficients method for (O)PLS models
#'
#' Coefficients of the (O)PLS(-DA) regression model
#'
#' @aliases coef.opls coef,opls-method
#' @param object An S4 object of class \code{opls}, created by \code{opls}
#' function.
#' @param ... Currently not used.
#' @return Numeric matrix of coefficients (number of rows equals the number of
#' variables, and the number of columns equals the number of responses)
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#' @examples
#'
#' data(sacurine)
#' attach(sacurine)
#'
#' sacurine.plsda <- opls(dataMatrix,
#'                        sampleMetadata[, "gender"])
#'
#' head(coef(sacurine.plsda))
#'
#' detach(sacurine)
#'
#' @rdname coef
#' @export
setMethod("coef", "opls",
          function(object, ...) {
              return(object@coefficientMN)
          })


####    residuals    ####

#' Residuals method for (O)PLS models
#'
#' Returns the residuals from the (O)PLS(-DA) regression models
#'
#' @aliases residuals.opls residuals,opls-method
#' @param object An S4 object of class \code{opls}, created by \code{opls}
#' function.
#' @param ... Currently not used.
#' @return Numeric matrix or vector (same dimensions as the modeled y
#' response); if y is a character vector or a factor (in case of
#' classification), the residuals equal 0 (predicted class identical to the
#' true class) or 1 (prediction error)
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#' @examples
#'
#' data(sacurine)
#' attach(sacurine)
#'
#' sacurine.pls <- opls(dataMatrix,
#'                      sampleMetadata[, "age"])
#'
#' head(residuals(sacurine.pls))
#'
#' detach(sacurine)
#'
#' @rdname residuals
#' @export
setMethod("residuals", "opls",
          function(object, ...) {

    if(grepl("PLS", object@typeC)) {

        fit <- fitted(object)

        if(length(object@subsetVi) == 0) {
            trainVi <- 1:length(fit)
        } else
            trainVi <- object@subsetVi

        if(mode(object@suppLs[["yMCN"]]) == "numeric") {
            y <- object@suppLs[["y"]]
            if(is.matrix(y))
                y <- y[trainVi, , drop = FALSE]
            else
                y <- y[trainVi]
            return(y - fit)
        } else
            return(as.numeric(as.character(c(object@suppLs[["y"]])[trainVi]) != as.character(c(fit))))

    } else
        stop("'residuals' defined for (O)PLS(-DA) regression models only", call. = FALSE)

}) ## residuals


####    predict    ####

#' Predict method for (O)PLS models
#'
#' Returns predictions of the (O)PLS(-DA) model on a new dataset
#'
#' @aliases predict.opls predict,opls-method
#' @param object An S4 object of class \code{opls}, created by \code{opls}
#' function.
#' @param newdata Either a data frame or a matrix, containing numeric columns
#' only, with the same number of columns (variables) as the 'x' used for model
#' training with 'opls'.
#' @param ... Currently not used.
#' @return Predictions (either a vector, factor, or matrix depending on the y
#' response used for training the model)
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#' @examples
#'
#' data(sacurine)
#' attach(sacurine)
#'
#' predictorMN <- dataMatrix
#' responseFc <- sampleMetadata[, "gender"]
#'
#' sacurine.plsda <- opls(predictorMN,
#'                        responseFc,
#'                        subset = "odd")
#'
#' trainVi <- getSubsetVi(sacurine.plsda)
#'
#' table(responseFc[trainVi], fitted(sacurine.plsda))
#'
#' table(responseFc[-trainVi],
#'       predict(sacurine.plsda, predictorMN[-trainVi, ]))
#'
#' detach(sacurine)
#'
#' @rdname predict
#' @export
setMethod("predict", "opls",
          function(object, newdata, ...) {

              if(object@typeC == "PCA")
                  stop("Predictions currently available for (O)PLS(-DA) models only (not PCA)",
                       call. = FALSE)

              if(missing(newdata)) {

                  return(fitted(object))

              } else {

                  if(is.data.frame(newdata)) {
                      if(!all(sapply(newdata, data.class) == "numeric")) {
                          stop("'newdata' data frame must contain numeric columns only", call. = FALSE)
                      } else
                          newdata <- as.matrix(newdata)
                  } else if(is.matrix(newdata)) {
                      if(mode(newdata) != "numeric")
                          stop("'newdata' matrix must be of 'numeric' mode", call. = FALSE)
                  } else
                      stop("'newdata' must be either a data.frame or a matrix", call. = FALSE)

                  if(ncol(newdata) != as.numeric(object@descriptionMC["X_variables", ])) {
                      if(length(object@xZeroVarVi) == 0) {
                          stop("'newdata' number of variables is ",
                               ncol(newdata),
                               " whereas the number of variables used for model training was ",
                               as.numeric(object@descriptionMC["X_variables", ]),
                               ".",
                               call. = FALSE)
                      } else if(ncol(newdata) - as.numeric(object@descriptionMC["X_variables", ]) ==
                                as.numeric(object@descriptionMC["near_zero_excluded_X_variables", ])) {
                          warning(as.numeric(object@descriptionMC["near_zero_excluded_X_variables", ]),
                                  " near zero variance variables excluded during the model training will be removed from 'newdata'.",
                                  call. = FALSE)
                          newdata <- newdata[, -object@xZeroVarVi, drop = FALSE]
                      } else {
                          stop("'newdata' number of variables (",
                               ncol(newdata),
                               ") does not correspond to the number of initial variables (",
                               as.numeric(object@descriptionMC["X_variables", ]),
                               ") minus the number of near zero variance variables excluded during the training (",
                               as.numeric(object@descriptionMC["near_zero_excluded_X_variables", ]),
                               ").",
                               call. = FALSE)
                      }
                  }

                  xteMN <- scale(newdata, object@xMeanVn, object@xSdVn)

                  if(object@summaryDF[, "ort"] > 0) {

                      for(noN in 1:object@summaryDF[, "ort"]) {
                          if(object@suppLs[["naxL"]]) {
                              xtoMN <- matrix(0, nrow = nrow(xteMN), ncol = 1)
                              for(i in 1:nrow(xtoMN)) {
                                  comVl <- complete.cases(xteMN[i, ])
                                  xtoMN[i, ] <- crossprod(xteMN[i, comVl], object@orthoWeightMN[comVl, noN]) / drop(crossprod(object@orthoWeightMN[comVl, noN]))
                              }
                          } else
                              xtoMN <- xteMN %*% object@orthoWeightMN[, noN]

                          xteMN <- xteMN - tcrossprod(xtoMN, object@orthoLoadingMN[, noN])
                      }

                  }

                  if(object@suppLs[["naxL"]]) {
                      yTesScaMN <- matrix(0, nrow = nrow(xteMN), ncol = ncol(object@coefficientMN),
                                          dimnames = list(rownames(xteMN), colnames(object@coefficientMN)))
                      for(j in 1:ncol(yTesScaMN))
                          for(i in 1:nrow(yTesScaMN)) {
                              comVl <- complete.cases(xteMN[i, ])
                              yTesScaMN[i, j] <- crossprod(xteMN[i, comVl], object@coefficientMN[comVl, j])
                          }
                  } else
                      yTesScaMN <- xteMN %*% object@coefficientMN

                  ## if(object@suppLs[["nayL"]])
                  ##     yTesScaMN <- yTesScaMN[!is.na(yMCN[testVi, ]), , drop = FALSE]

                  yTesMN <- scale(scale(yTesScaMN,
                                        FALSE,
                                        1 / object@ySdVn),
                                  -object@yMeanVn,
                                  FALSE)
                  attr(yTesMN, "scaled:center") <- NULL
                  attr(yTesMN, "scaled:scale") <- NULL

                  if(is.factor(fitted(object))) {

                      yTestMCN <- object@suppLs[[".char2numF"]](yTesMN,
                                                                c2nL = FALSE)
                      predMCNFcVcn <- as.character(yTestMCN)
                      names(predMCNFcVcn) <- rownames(newdata)
                      predMCNFcVcn <- factor(predMCNFcVcn, levels = levels(object@suppLs[["y"]]))

                  } else if(is.vector(fitted(object))) {

                      if(is.character(fitted(object))) {

                          yTestMCN <- object@suppLs[[".char2numF"]](yTesMN,
                                                                    c2nL = FALSE)
                          predMCNFcVcn <- as.character(yTestMCN)
                          names(predMCNFcVcn) <- rownames(newdata)

                      } else {

                          predMCNFcVcn <- as.numeric(yTesMN)
                          names(predMCNFcVcn) <- rownames(newdata)

                      }

                  } else if(is.matrix(fitted(object))) {

                      if(mode(fitted(object)) == "character") {
                          predMCNFcVcn  <- object@suppLs[[".char2numF"]](yTesMN,
                                                                         c2nL = FALSE)
                      } else
                          predMCNFcVcn <- yTesMN

                      rownames(predMCNFcVcn) <- rownames(newdata)

                  }

                  return(predMCNFcVcn)

              }

          })


####    getEset    ####

#' getEset method
#'
#' Extracts the complemented ExpressionSet when opls has been applied to an ExpressionSet
#'
#' @aliases getEset getEset, opls-method
#' @param object An S4 object of class \code{opls}, created by \code{opls}
#' function.
#' @return An S4 object of class \code{ExpressionSet} which contains the dataMatrix (t(exprs(eset))),
#' and the sampleMetadata (pData(eset)) and variableMetadata (fData(eset)) with the additional columns
#' containing the scores, predictions, loadings, VIP, coefficients etc.
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#' @examples
#'
#' data(sacurine)
#' attach(sacurine)
#' 
#' exprs <- t(dataMatrix)
#' sacSet <- Biobase::ExpressionSet(assayData = exprs, 
#'                                  phenoData = new("AnnotatedDataFrame", 
#'                                                  data = sampleMetadata), 
#'                                  featureData = new("AnnotatedDataFrame", 
#'                                                    data = variableMetadata),
#'                                  experimentData = new("MIAME", 
#'                                                       title = "sacurine"))
#'                                                       
#' sacPlsda <- opls(sacSet, "gender")
#' sacSet <- getEset(sacPlsda)
#' head(Biobase::pData(sacSet))
#' head(Biobase::fData(sacSet))
#' 
#' detach(sacurine)
#'
#' @rdname getEset
#' @export
setMethod("getEset", "opls",
          function(object) {
            return(object@eset)
          })


####    getLoadingMN    ####

#' getLoadingMN method for PCA/(O)PLS(-DA) models
#'
#' (Orthogonal) loadings of the PCA/(O)PLS(-DA) model
#'
#' @aliases getLoadingMN getLoadingMN,opls-method
#' @param object An S4 object of class \code{opls}, created by \code{opls}
#' function.
#' @param orthoL Logical: Should the orthogonal loading matrix be returned
#' @param ... Currently not used.
#' (default is FALSE and the predictive loading matrix is returned)
#' @return Numeric matrix with a number of rows equal to the number of
#' variables and a number of columns equal to the number of components
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#' @examples
#'
#' data(sacurine)
#' attach(sacurine)
#'
#' sacurine.plsda <- opls(dataMatrix,
#'                        sampleMetadata[, "gender"])
#'
#' getLoadingMN(sacurine.plsda)
#'
#' detach(sacurine)
#'
#' @rdname getLoadingMN
#' @export
setMethod("getLoadingMN", "opls",
          function(object, orthoL = FALSE) {
            if(orthoL)
              return(object@orthoLoadingMN)
            else
              return(object@loadingMN)
          })


####    getPcaVarVn    ####

#' getPcaVarVn method for PCA models
#'
#' Variance of the components (score vectors)
#'
#' @aliases getPcaVarVn getPcaVarVn,opls-method
#' @param object An S4 object of class \code{opls}, created by \code{opls}
#' function.
#' @param ... Currently not used.
#' @return Numeric vector with the same length as the number of components
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#' @examples
#'
#' data(sacurine)
#' attach(sacurine)
#'
#' sacurine.pca <- opls(dataMatrix)
#'
#' getPcaVarVn(sacurine.pca)
#'
#' detach(sacurine)
#'
#' @rdname getPcaVarVn
#' @export
setMethod("getPcaVarVn", "opls",
          function(object) {
                  return(object@pcaVarVn)
          })


####    getScoreMN    ####

#' getScoreMN method for PCA/(O)PLS(-DA) models
#'
#' (Orthogonal) scores of the (O)PLS(-DA) model
#'
#' @aliases getScoreMN getScoreMN,opls-method
#' @param object An S4 object of class \code{opls}, created by \code{opls}
#' function.
#' @param orthoL Logical: Should the orthogonal score matrix be returned
#' (default is FALSE and the predictive score matrix is returned)
#' @param ... Currently not used.
#' @return Numeric matrix with a number of rows equal to the number of samples
#' and a number of columns equal to the number of components
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#' @examples
#'
#' data(sacurine)
#' attach(sacurine)
#'
#' sacurine.plsda <- opls(dataMatrix,
#'                        sampleMetadata[, "gender"])
#'
#' getScoreMN(sacurine.plsda)
#'
#' detach(sacurine)
#'
#' @rdname getScoreMN
#' @export
setMethod("getScoreMN", "opls",
          function(object, orthoL = FALSE) {
              if(orthoL)
                  return(object@orthoScoreMN)
              else
                  return(object@scoreMN)
          })




####    getSubsetVi    ####

#' getSubsetVi method for (O)PLS(-DA) models
#'
#' Extracts the indices of the samples used for building the model (when a
#' subset argument has been specified)
#'
#' @aliases getSubsetVi getSubsetVi,opls-method
#' @param object An S4 object of class \code{opls}, created by \code{opls}
#' function.
#' @param ... Currently not used.
#' @return Integer vector with the indices of the samples used for training
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#' @examples
#'
#' data(sacurine)
#' attach(sacurine)
#'
#' predictorMN <- dataMatrix
#' responseFc <- sampleMetadata[, "gender"]
#'
#' sacurine.plsda <- opls(predictorMN,
#'                        responseFc,
#'                        subset = "odd")
#'
#' trainVi <- getSubsetVi(sacurine.plsda)
#'
#' table(responseFc[trainVi], fitted(sacurine.plsda))
#'
#' detach(sacurine)
#'
#' @rdname getSubsetVi
#' @export
setMethod("getSubsetVi", "opls",
          function(object) {
                  return(object@subsetVi)
          })


####    getSummaryDF    ####

#' getSummaryDF method for PCA/(O)PLS models
#'
#' Summary of model metrics
#'
#' @aliases getSummaryDF getSummaryDF,opls-method
#' @param object An S4 object of class \code{opls}, created by \code{opls}
#' function.
#' @param ... Currently not used.
#' @return Data frame
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#' @examples
#'
#' data(sacurine)
#' attach(sacurine)
#'
#' sacurine.plsda <- opls(dataMatrix,
#'                        sampleMetadata[, "gender"])
#'
#' getSummaryDF(sacurine.plsda)
#'
#' detach(sacurine)
#'
#' @rdname getSummaryDF
#' @export
setMethod("getSummaryDF", "opls",
          function(object) {
            return(object@summaryDF)
          })


####    getVipVn    ####

#' getVipVn method for (O)PLS(-DA) models
#'
#' (Orthogonal) VIP of the (O)PLS(-DA) model
#'
#'
#' @aliases getVipVn getVipVn,opls-method
#' @param object An S4 object of class \code{opls}, created by \code{opls}
#' function.
#' @param orthoL Logical: Should the orthogonal VIP be returned (default is
#' FALSE and the predictive VIP is returned)
#' @param ... Currently not used.
#' @return Numeric vector with a length equal to the number of variables and a
#' number of columns equal to the number of components
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#' @references Galindo-Prieto B., Eriksson L. and Trygg J. (2014). Variable
#' influence on projection (VIP) for orthogonal projections to latent
#' structures (OPLS). Journal of Chemometrics 28, 623-632.
#' @examples
#'
#' data(sacurine)
#' attach(sacurine)
#'
#' sacurine.plsda <- opls(dataMatrix,
#'                        sampleMetadata[, "gender"])
#'
#' getVipVn(sacurine.plsda)
#'
#' detach(sacurine)
#'
#' @rdname getVipVn
#' @export
setMethod("getVipVn", "opls",
          function(object, orthoL = FALSE) {
            if(orthoL)
              return(object@orthoVipVn)
            else
              return(object@vipVn)
          })


####    getWeightMN    ####

#' getWeightMN method for (O)PLS(-DA) models
#'
#' (Orthogonal) weights of the (O)PLS(-DA) model
#'
#'
#' @aliases getWeightMN getWeightMN,opls-method
#' @param object An S4 object of class \code{opls}, created by \code{opls}
#' function.
#' @param orthoL Logical: Should the orthogonal weight matrix be returned
#' @param ... Currently not used.
#' (default is FALSE and the predictive weight matrix is returned)
#' @return Numeric matrix with a number of rows equal to the number of
#' variables and a number of columns equal to the number of components
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#' @examples
#'
#' data(sacurine)
#' attach(sacurine)
#'
#' sacurine.plsda <- opls(dataMatrix,
#'                        sampleMetadata[, "gender"])
#'
#' getWeightMN(sacurine.plsda)
#'
#' detach(sacurine)
#'
#' @rdname getWeightMN
#' @export
setMethod("getWeightMN", "opls",
          function(object, orthoL = FALSE) {
            if(orthoL)
              return(object@orthoWeightMN)
            else
              return(object@weightMN)
          })


#' Checking the consistency of an ExpressionSet instance with W4M format
#'
#' @param eset An S4 object of class ExpressionSet.
#' @param ... Currently not used.
#' @return Invisible TRUE logical in case of success (otherwise generates an
#' error)
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#' @examples
#' sacSet <- fromW4M(file.path(path.package("ropls"), "extdata"))
#' print(checkW4M(sacSet))
#' @rdname checkW4M
setMethod("checkW4M", "ExpressionSet",
         function(eset, ...) {

             datMN <- t(exprs(eset))
             samDF <- pData(eset)
             varDF <- fData(eset)

             chkL <- .checkW4mFormatF(datMN, samDF, varDF)

             if(!chkL) {
                 stop("Problem with the sample or variable names in the tables to be imported from (exported to) W4M", call. = FALSE)
             } else
                 return(TRUE)
         })

#' Exporting ExpressionSet instance into 3 tabulated files.
#'
#' The 3 .tsv files are written with the indicated \code{file} prefix, and
#' '_dataMatrix.tsv', '_sampleMetadata.tsv', and '_variableMetadata.tsv'
#' suffices, respectively. Note that the \code{dataMatrix} is transposed before
#' export (e.g., the samples are written column wise in the 'dataMatrix.tsv'
#' exported file).
#'
#' @param eset An S4 object of class \code{ExpressionSet}
#' function.
#' @param filePrefixC Character: common prefix (including repository full path)
#' of the three file names: for example, the 'c:/mydata/setname' value will
#' result in writting the 'c:/mydata/setname_dataMatrix.tsv',
#' 'c:/mydata/setname_sampleMetadata.tsv', and
#' 'c:/mydata/setname_variableMetadata.tsv' files.
#' @param verboseL Logical: should comments be printed?
#' @param ... Currently not used.
#' @return No object returned.
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#' @examples
#'  sacSet <- fromW4M(file.path(path.package("ropls"), "extdata"))
#'  toW4M(sacSet)
#' @rdname toW4M
setMethod("toW4M", "ExpressionSet",
          function(eset, filePrefixC = paste0(getwd(), "/out_"), verboseL = TRUE, ...){

              if(checkW4M(eset)) {

                  datMN <- exprs(eset)
                  datDF <- cbind.data.frame(dataMatrix = rownames(datMN),
                                            as.data.frame(datMN))

                  filDatC <- paste0(filePrefixC, "dataMatrix.tsv")
                  filSamC <- paste0(filePrefixC, "sampleMetadata.tsv")
                  filVarC <- paste0(filePrefixC, "variableMetadata.tsv")

                  write.table(datDF,
                              file = filDatC,
                              quote = FALSE,
                              row.names = FALSE,
                              sep = "\t")

                  samDF <- pData(eset)
                  samDF <- cbind.data.frame(sampleMetadata = rownames(samDF),
                                            samDF)
                  write.table(samDF,
                              file = filSamC,
                              quote = FALSE,
                              row.names = FALSE,
                              sep = "\t")

                  varDF <- fData(eset)
                  varDF <- cbind.data.frame(variableMetadata = rownames(varDF),
                                            varDF)
                  write.table(varDF,
                              file = filVarC,
                              quote = FALSE,
                              row.names = FALSE,
                              sep = "\t")

                  if(verboseL) {
                      cat("The following 3 files:\n")
                      print(basename(filDatC))
                      print(basename(filSamC))
                      print(basename(filVarC))
                      cat("have been written in the following directory:\n")
                      print(dirname(filDatC))
                  }

              }

          })
