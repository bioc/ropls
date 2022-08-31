#### checkW4M ####

#' Checking the consistency of a SummarizedExperiment or ExpressionSet instance with W4M format
#'
#' @param x An S4 object of class \code{SummarizedExperiment} or \code{ExpressionSet}
#' @return Invisible TRUE logical in case of success (otherwise generates an
#' error)
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#' @examples
#' sacSet <- fromW4M(file.path(path.package("ropls"), "extdata"))
#' print(checkW4M(sacSet))
#' @rdname checkW4M
#' @export
setGeneric("checkW4M", function(x) {standardGeneric("checkW4M")})


####   getEset    ####

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
#' sacSet <- sacurine[["eset"]]
#' sacPlsda <- opls(sacSet, "gender")
#' sacSet <- getEset(sacPlsda)
#' head(Biobase::pData(sacSet))
#' head(Biobase::fData(sacSet))
#'
#' @rdname getEset
#' @export
setGeneric("getEset",
           function(object) {standardGeneric("getEset")})


####   getLoadingMN    ####

#' getLoadingMN method for PCA/(O)PLS(-DA) models
#'
#' (Orthogonal) loadings of the PCA/(O)PLS(-DA) model
#'
#' @aliases getLoadingMN getLoadingMN,opls-method
#' @param object An S4 object of class \code{opls}, created by \code{opls}
#' function.
#' @param orthoL Logical: Should the orthogonal loading matrix be returned
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
setGeneric("getLoadingMN",
           function(object, orthoL = FALSE) {standardGeneric("getLoadingMN")})


####   getMset    ####

#' getMset method
#'
#' Extracts the complemented MultiDataSet when opls has been applied to a MultiDataSet
#'
#' @aliases getMset getMset, oplsMultiDataSet-method
#' @param object An S4 object of class \code{oplsMultiDataSet}, created by \code{opls}
#' function applied to a MultiDataSet
#' @return An S4 object of class \code{MultiDataSet}.
#' @examples
#' data(NCI60)
#' nci.mds <- NCI60[["mds"]]
#' # Restricting to the 'agilent' and 'hgu95' datasets
#' nci.mds <- nci.mds[, c("agilent", "hgu95")]
#' # Restricting to the 'ME' and 'LE' cancer types
#' sampleNamesVc <- Biobase::sampleNames(nci.mds[["agilent"]])
#' cancerTypeVc <- Biobase::pData(nci.mds[["agilent"]])[, "cancer"]
#' nci.mds <- nci.mds[sampleNamesVc[cancerTypeVc %in% c("ME", "LE")], ]
#' # Principal Component Analysis of each data set
#' nci.pca <- opls(nci.mds)
#' # Getting the MultiDataSet with additional info. in pData and fData
#' nci.mds <- getMset(nci.pca)
#' @rdname getMset
#' @export
setGeneric("getMset",
           function(object) {standardGeneric("getMset")})


####   getOpls    ####

#' Getting the models from a SummarizedExperiment or a MultiAssayExperiment object
#'
#' The models are extracted as a list
#'
#' @param object An S4 object of class \code{SummarizedExperiment} or \code{MultiAssayExperiment},
#' once processed by the \code{opls} method
#' @return List of opls models contained in the SummarizedExperiment object(s)
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#' @examples
#'
#' 
#' # Getting the sacurine data set as a SummarizedExperiment
#' data(sacurine)
#' sac.se <- sacurine[["se"]]
#'
#' # Building the PLS-DA model
#' sac.se <- opls(sac.se, "gender")
#'
#' # Getting the models
#' sac_opls.ls <- getOpls(sac.se)
#' names(sac_opls.ls)
#' 
#' # Plotting the score plot from the PLS-DA model of the gender response
#' plot(sac_opls.ls[["gender_PLSDA"]], typeVc = "x-score")
#'
#' @rdname getOpls
#' @export
setGeneric("getOpls",
           function(object) {standardGeneric("getOpls")})


####   getPcaVarVn    ####

#' getPcaVarVn method for PCA models
#'
#' Variance of the components (score vectors)
#'
#' @aliases getPcaVarVn getPcaVarVn,opls-method
#' @param object An S4 object of class \code{opls}, created by \code{opls}
#' function.
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
setGeneric("getPcaVarVn",
           function(object) {standardGeneric("getPcaVarVn")})


####   getScoreMN    ####

#' getScoreMN method for PCA/(O)PLS(-DA) models
#'
#' (Orthogonal) scores of the (O)PLS(-DA) model
#'
#' @aliases getScoreMN getScoreMN,opls-method
#' @param object An S4 object of class \code{opls}, created by \code{opls}
#' function.
#' @param orthoL Logical: Should the orthogonal score matrix be returned
#' (default is FALSE and the predictive score matrix is returned)
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
setGeneric("getScoreMN",
           function(object, orthoL = FALSE) {standardGeneric("getScoreMN")})


####   getSubsetVi    ####

#' getSubsetVi method for (O)PLS(-DA) models
#'
#' Extracts the indices of the samples used for building the model (when a
#' subset argument has been specified)
#'
#' @aliases getSubsetVi getSubsetVi,opls-method
#' @param object An S4 object of class \code{opls}, created by \code{opls}
#' function.
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
setGeneric("getSubsetVi",
           function(object) {standardGeneric("getSubsetVi")})


####    getSummaryDF    ####

#' getSummaryDF method for PCA/(O)PLS models
#'
#' Summary of model metrics
#'
#' @aliases getSummaryDF getSummaryDF,opls-method
#' @param object An S4 object of class \code{opls}, created by \code{opls}
#' function.
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
setGeneric("getSummaryDF",
           function(object) {standardGeneric("getSummaryDF")})


####   getVipVn    ####

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
setGeneric("getVipVn",
           function(object, orthoL = FALSE) {standardGeneric("getVipVn")})


####   getWeightMN    ####

#' getWeightMN method for (O)PLS(-DA) models
#'
#' (Orthogonal) weights of the (O)PLS(-DA) model
#'
#'
#' @aliases getWeightMN getWeightMN,opls-method
#' @param object An S4 object of class \code{opls}, created by \code{opls}
#' function.
#' @param orthoL Logical: Should the orthogonal weight matrix be returned? (default is FALSE)
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
setGeneric("getWeightMN",
           function(object, orthoL = FALSE) {standardGeneric("getWeightMN")})

#### gg_scoreplot ####

#' PCA and (O)PLS(-DA) score plots
#'
#' Score plot visualization for PCA and (O)PLS(-DA) models in either ggplot or ggplotly formats
#' 
#' @param x An S4 object of class \code{SummarizedExperiment} (resp. \code{opls})
#' generated by the 'ropls::opls' modeling applied to a \code{SummarizedExperiment}
#' (resp. an \code{ExpressionSet})
#' @param model.c character(1): name of the model to be plotted; 
#' use 'names(ropls::getOpls(se))' to see the available models in the se object
#' @param components.vi integer(2): number of the components to display as x and y axis
#' @param label.c character(1): name of the colData (resp. pData) column 
#' to be used for the labels
#' @param color.c character(1): name of the colData (resp. pData) column 
#' to be used for the colors
#' @param title.c character(1): plot title
#' @param legend.c character(1): position of the legend (either 'bottom', 'left',
#' 'top' or 'right' [default])
#' @param palette.c character(1): name of the RColorBrewer palette (for qualitative factor)
#' @param ellipse.l logical(1): should ellipses be drawn (for qualitative factor)
#' @param plotly.l logical(1): should the ggplot be converted to an interactive
#' plotly (default: FALSE)
#' @param info.vc character(): names of the colData (resp. pData) columns
#' to be used for the plotly info; the default 'sample_names' will return 
#' the sample names as the plotly info
#' @param size.ls list: sizes for axis labels (default: 16), axis text (default: 14),
#' points (default: 3), labels (default = 5), title (default = 20), legend title (default: 15),
#' legend text (default: 15)
#' @return invisible ggplot2 (or ggplotly) object
#' @examples
#' # loading the 'sacurine' dataset from the 'ropls' package
#' data(sacurine, package = "ropls")
#' # SummarizedExperiment
#' sac.se <- sacurine[["se"]]
#' ## computing the PCA
#' sac.se <- ropls::opls(sac.se)
#' ## score plot
#' gg_scoreplot(sac.se, "PCA")
#' gg_scoreplot(sac.se, "PCA", color.c = "age")
#' gg_scoreplot(sac.se, "PCA", color.c = "gender", plotly.l = TRUE, info.vc = "all")
#' ## PLS-DA modeling
#' sac.se <- ropls::opls(sac.se, "gender")
#' gg_scoreplot(sac.se, "gender_PLSDA", color.c = "gender")
#' gg_scoreplot(sac.se, "gender_PLSDA", color.c = "gender", plotly.l = TRUE)
#' # ExpressionSet
#' sacurine.eset <- sacurine[["eset"]]
#' ## PCA
#' sacurine.pca <- ropls::opls(sacurine.eset)
#' ## score plot
#' gg_scoreplot(sacurine.pca)
#' gg_scoreplot(sacurine.pca, color.c = "age")
#' @rdname gg_scoreplot
#' @export
setGeneric("gg_scoreplot",
           function(x,
                    model.c = "",
                    components.vi = c(1, 2),
                    label.c = c("", "sample_names")[2],
                    color.c = "",
                    title.c = "",
                    palette.c = "Set1",
                    legend.c = "right",
                    ellipse.l = TRUE,
                    plotly.l = FALSE,
                    info.vc = "sample_names",
                    size.ls = list(axis_lab.i = 16,
                                   axis_text.i = 14,
                                   point.i = 3,
                                   label.i = 5,
                                   title.i = 20,
                                   legend_title.i = 15,
                                   legend_text.i = 15)) standardGeneric("gg_scoreplot"))


#### opls ####

#' PCA, PLS(-DA), and OPLS(-DA)
#'
#' PCA, PLS, and OPLS regression, classification, and cross-validation with the
#' NIPALS algorithm
#'
#' @name opls
#' @aliases opls opls,matrix-method opls,data.frame-method
#' opls,SummarizedExperiment-method opls,ExpressionSet-method
#' opls,MultiAssayExperiment-method opls,MultiDataSet-method
#' @docType methods
#' @param x Numerical matrix, (observations x variables; NAs are
#' allowed), data.frame, SummarizedExperiment or ExpressionSet object
#' @param y Response to be modelled: Either 1) 'NULL' for PCA (default) or 2) a
#' numerical vector (same length as 'x' row number) for single response (O)PLS,
#' or 3) a numerical matrix (same row number as 'x') for multiple response PLS,
#' 4) a factor (same length as 'x' row number) for (O)PLS-DA, or 5) a character
#' indicating the name of the column of the phenoData@data to be used, when x
#' is a SummarizedExperiment or an ExpressionSet object. Note that, for convenience, character vectors
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
#' @param plotSubC Character: Graphic subtitle
#' @param fig.pdfC Character: File name with '.pdf' extension for the figure;
#' if 'interactive' (default), figures will be displayed interactively; if 'none',
#' no figure will be generated
#' @param info.txtC Character: File name with '.txt' extension for the printed
#' results (call to sink()'); if 'interactive' (default), messages will be
#' printed on the screen; if 'none', no verbose will be generated
#' @return An S4 object of class 'opls' containing the slots described below; in case x is a SummarizedExperiment, a SummarizedExperiment is returned, with the 'opls' object included in the metadata:
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
#' ## PCA
#'
#' data(foods) ## see Eriksson et al. (2001); presence of 3 missing values (NA)
#' head(foods)
#' foodMN <- as.matrix(foods[, colnames(foods) != "Country"])
#' rownames(foodMN) <- foods[, "Country"]
#' head(foodMN)
#' foo.pca <- opls(foodMN)
#'
#' ## PLS with a single response
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
#' ## PLS-DA
#'
#' data(sacurine)
#' attach(sacurine)
#' sacurine.plsda <- opls(dataMatrix, sampleMetadata[, "gender"])
#'
#' ## OPLS-DA
#'
#' sacurine.oplsda <- opls(dataMatrix, sampleMetadata[, "gender"], predI = 1, orthoI = NA)
#' 
#' detach(sacurine)
#' 
#' ## Application to a SummarizedExperiment
#' 
#' sac.se <- sacurine[["se"]]
#' sac.se <- opls(sac.se, "gender")
#' SummarizedExperiment::colData(sac.se)
#' SummarizedExperiment::rowData(sac.se)
#' sac_gender.plsda <- sac.se@metadata[["opls"]][["gender_PLSDA"]]
#' plot(sac_gender.plsda, typeVc = "x-score")
#' 
#' ## Application to a MultiAssayExperiment
#' 
#' data(NCI60)
#' nci.mae <- NCI60[["mae"]]
#' # Restricting to the 'ME' and 'LE' cancer types and to the 'agilent' and 'hgu95' datasets
#' library(MultiAssayExperiment)
#' nci.mae <- nci.mae[, nci.mae$cancer %in% c("ME", "LE"), c("agilent", "hgu95")]
#' # Principal Component Analysis of each data set
#' nci.mae <- opls(nci.mae)
#' # Coloring the score plots according to cancer types
#' for (set.c in names(nci.mae))
#' plot(getOpls(nci.mae)[[set.c]][["PCA"]],
#' parAsColFcVn = MultiAssayExperiment::colData(nci.mae)[, "cancer"],
#' typeVc = "x-score",
#' plotSubC = set.c)
#' # Building PLS-DA models for the cancer type, and getting back the updated MultiDataSet
#' nci.mae <- opls(nci.mae, "cancer", predI = 2)
#' # Viewing the new variable metadata (including VIP and coefficients)
#' lapply(names(nci.mae), function(set.c) head(SummarizedExperiment::rowData(nci.mae[[set.c]])))
#' 
#' ## Application to an ExpressionSet
#' 
#' sacSet <- sacurine[["eset"]]
#'                                                       
#' sacPlsda <- opls(sacSet, "gender")
#' sacSet <- getEset(sacPlsda)
#' head(Biobase::pData(sacSet))
#' head(Biobase::fData(sacSet))
#' 
#' ## Application to a MultiDataSet
#' 
#' data(NCI60)
#' nci.mds <- NCI60[["mds"]]
#' # Restricting to the 'agilent' and 'hgu95' datasets
#' nci.mds <- nci.mds[, c("agilent", "hgu95")]
#' # Restricting to the 'ME' and 'LE' cancer types
#' sampleNamesVc <- Biobase::sampleNames(nci.mds[["agilent"]])
#' cancerTypeVc <- Biobase::pData(nci.mds[["agilent"]])[, "cancer"]
#' nci.mds <- nci.mds[sampleNamesVc[cancerTypeVc %in% c("ME", "LE")], ]  
#' # Principal Component Analysis of each data set
#' nci.pca <- opls(nci.mds)
#' # Coloring the Score plot according to cancer types
#' plot(nci.pca, parAsColFcVn = Biobase::pData(nci.mds[["agilent"]])[, "cancer"], typeVc = "x-score")
#' # Getting the updated MultiDataSet (now including scores and loadings)
#' nci.mds <- getMset(nci.pca)
#' # Building PLS-DA models for the cancer type, and getting back the updated MultiDataSet
#' nci.plsda <- opls(nci.mds, "cancer", predI = 2)
#' nci.mds <- getMset(nci.plsda)
#' # Viewing the new variable metadata (including VIP and coefficients)
#' lapply(Biobase::fData(nci.mds), head)
#' @rdname opls
#' @export
setGeneric("opls",
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
                    
                    plotSubC = NA,                   
                    fig.pdfC = c("none", "interactive", "myfile.pdf")[2],                   
                    info.txtC = c("none", "interactive", "myfile.txt")[2]) standardGeneric("opls"))


####    tested    ####

#' Tested method for (O)PLS models
#'
#' Returns predictions of the (O)PLS(-DA) model on the out of the box samples
#' (when a 'subset' of samples has been selected when training the model)
#'
#' @aliases tested tested,opls-method
#' @param object An S4 object of class \code{opls}, created by \code{opls}
#' function.
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
setGeneric("tested",
           function(object) standardGeneric("tested"))


#### view ####

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
#' @return this method has no output
#' @examples
#' library(ropls)
#' 
#' # Get the sacurine dataset
#' 
#' data(sacurine)
#' 
#' # Display the data matrix
#' 
#' view(sacurine[["dataMatrix"]])
#' view(sacurine[["dataMatrix"]][, 1:40], mainC = "'Sacurine' dataset", rowAllL = TRUE,
#' colAllL = TRUE, colTruncI = 13, colMarN = 7)
#' view(sacurine[["dataMatrix"]][, 1:40], mainC = "'Sacurine' dataset", paletteC = "ramp")
#' 
#' # Display the sample metadata (dataframe)
#' 
#' view(sacurine[["sampleMetadata"]])
#' 
#' # Display the SummarizedExperiment
#' 
#' view(sacurine[["se"]])
#' 
#' # Display the ExpressionSet
#' 
#' view(sacurine[["eset"]])
#' 
#' @rdname view
#' @export
setGeneric("view",
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
                    colMarN = 3.1,
                    colLabC = "",
                    colTruncI = 0,
                    drawScaleL = TRUE,
                    delimitReplicatesL = FALSE,
                    standardizeL = FALSE,
                    fig.pdfC = "interactive") {standardGeneric("view")})


#### toW4M ####

#' Exporting a SummarizedExperiment or ExpressionSet instance into 3 tabulated files.
#'
#' The 3 .tsv files are written with the indicated \code{file} prefix, and
#' '_dataMatrix.tsv', '_sampleMetadata.tsv', and '_variableMetadata.tsv'
#' suffices, respectively. Note that the \code{dataMatrix} is transposed before
#' export (e.g., the samples are written column wise in the 'dataMatrix.tsv'
#' exported file).
#'
#' @param x An S4 object of class \code{SummarizedExperiment} or \code{ExpressionSet}
#' function.
#' @param filePrefixC Character: common prefix (including repository full path)
#' of the three file names: for example, the 'c:/mydata/setname' value will
#' result in writting the 'c:/mydata/setname_dataMatrix.tsv',
#' 'c:/mydata/setname_sampleMetadata.tsv', and
#' 'c:/mydata/setname_variableMetadata.tsv' files.
#' @param verboseL Logical: should comments be printed?
#' @return No object returned.
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#' @examples
#'  sacSet <- fromW4M(file.path(path.package("ropls"), "extdata"))
#'  toW4M(sacSet)
#' @rdname toW4M
#' @export
setGeneric("toW4M", function(x,
                             filePrefixC = paste0(getwd(), "/out_"),
                             verboseL = TRUE) standardGeneric("toW4M"))
