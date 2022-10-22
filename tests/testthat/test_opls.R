library(ropls)

context("Testing 'ropls'")

#### PCA ####

test_that("PCA", {
  
  data(foods) ## contains 3 NA
  
  fooMN <- as.matrix(foods[, colnames(foods) != "Country"])
  rownames(fooMN) <- foods[, "Country"]
  
  foo.pca <- opls(fooMN, permI = 0, fig.pdfC = "none", info.txtC = "none")
  testthat::expect_equivalent(getPcaVarVn(foo.pca)[1],
                              6.19,
                              tolerance = 1e-3)
  testthat::expect_equivalent(getScoreMN(foo.pca)[1, 1],
                              1.309494,
                              tolerance = 1e-5)
  
})

#### PCA_sacurine ####

test_that("PCA_sacurine",  {
  
  data(sacurine)
  
  sac.pca <- opls(sacurine[["dataMatrix"]],
                  fig.pdfC = "none", info.txtC = "none")
  testthat::expect_equivalent(getSummaryDF(sac.pca)["Total", "R2X(cum)"],
                              0.501,
                              tolerance = 1e-3)
  testthat::expect_equivalent(sac.pca@modelDF["p3", "R2X"],
                              0.067,
                              tolerance = 1e-3)
  testthat::expect_equivalent(sac.pca@modelDF["p2", "R2X(cum)"],
                              0.252,
                              tolerance = 1e-3)
  
  sac.pcaCenter <- opls(sacurine[["dataMatrix"]],
                        scaleC = "center",
                        fig.pdfC = "none", info.txtC = "none")
  testthat::expect_equivalent(sac.pcaCenter@modelDF["p2", "R2X(cum)"],
                              0.223,
                              tolerance = 1e-3)
  
  sac.pcaPareto <- opls(sacurine[["dataMatrix"]],
                        scaleC = "pareto",
                        fig.pdfC = "none", info.txtC = "none")
  testthat::expect_equivalent(sac.pcaPareto@modelDF["p3", "R2X"],
                              0.063,
                              tolerance = 1e-3)
  
  sac.pcaSvd <- opls(sacurine[["dataMatrix"]],
                     algoC = "svd",
                     fig.pdfC = "none", info.txtC = "none")
  testthat::expect_equivalent(sac.pcaSvd@modelDF["p3", "R2X"],
                              0.067,
                              tolerance = 1e-3)
  testthat::expect_equivalent(getScoreMN(sac.pcaSvd)[1, 1],
                              -8.744009,
                              tolerance = 1e-5)
  
})

#### PLS_single ####

test_that("PLS_single", {
  
  data(cornell) ## see Tenenhaus, 1998
  
  corPlsLs <- opls(as.matrix(cornell[, grep("x", colnames(cornell))]),
                   matrix(cornell[, "y"], ncol = 1),
                   permI = 0, fig.pdfC = "none", info.txtC = "none")
  testthat::expect_equivalent(getScoreMN(corPlsLs)[1, 1],
                              2.0513,
                              tolerance = 1e-3)
  cornell.pls <- opls(as.matrix(cornell[, grep("x", colnames(cornell))]),
                      cornell[, "y"],
                      permI = 0, fig.pdfC = "none", info.txtC = "none")
  testthat::expect_equivalent(getVipVn(cornell.pls)[1],
                              1.125,
                              tolerance = 1e-3)
  
})

#### PLS_multiple ####

test_that("PLS_multiple", {
  
  ## lowarp
  
  data(lowarp) ## see Eriksson et al. (2001); presence of NAs
  
  lowarp.pls <- opls(as.matrix(lowarp[, c("glas", "crtp", "mica", "amtp")]),
                     as.matrix(lowarp[, grepl("^wrp", colnames(lowarp)) |
                                        grepl("^st", colnames(lowarp))]),
                     permI = 0, fig.pdfC = "none", info.txtC = "none")
  testthat::expect_equivalent(lowarp.pls@weightStarMN[2, 1],
                              -0.0861,
                              tolerance = 1e-3)
  testthat::expect_equivalent(lowarp.pls@uMN["s5", "p2"],
                              0.596598855,
                              tolerance = 1e-8)
  
  ## sacurine
  
  data(sacurine)
  
  agebmiPls <- opls(sacurine[["dataMatrix"]],
                    as.matrix(sacurine[["sampleMetadata"]][, c("age","bmi")]),
                    fig.pdfC = "none", info.txtC = "none")
  agebmiPreMN <- predict(agebmiPls)
  testthat::expect_identical(colnames(agebmiPreMN),
                             c("age", "bmi"))
  testthat::expect_equivalent(agebmiPreMN[1, 2], 20.5831139196, tol = 1e-10)
  
})


#### PLS_predict_sacurine ####

test_that("PLS_predict_sacurine", {
  
  data(sacurine)
  
  sac.opls <- opls(sacurine[["dataMatrix"]],
                   sacurine[["sampleMetadata"]][, "age"],
                   predI = 1,
                   orthoI = 1,
                   permI = 0, fig.pdfC = "none", info.txtC = "none")
  testthat::expect_equivalent(predict(sac.opls)[107],
                              42.97439,
                              tolerance = 1e-5)
  
  sac.opls <- opls(sacurine[["dataMatrix"]],
                   sacurine[["sampleMetadata"]][, "gender"],
                   predI = 1,
                   orthoI = 1,
                   permI = 0, fig.pdfC = "none", info.txtC = "none")
  pred107Fc <- "M"
  names(pred107Fc) <- "HU_125"
  pred107Fc <- factor(pred107Fc, levels = c("M", "F"))
  testthat::expect_identical(predict(sac.opls)[107],
                             pred107Fc)
  
  trainVi <- 1:floor(nrow(sacurine[["dataMatrix"]]) / 2)
  
  sac.pls.tr <- opls(sacurine[["dataMatrix"]],
                     sacurine[["sampleMetadata"]][, "age"],
                     subset = trainVi,
                     permI = 0, fig.pdfC = "none", info.txtC = "none")
  testthat::expect_equivalent(as.numeric(sac.pls.tr@descriptionMC["samples", 1]),
                              length(trainVi),
                              tolerance = 0)
  testthat::expect_equivalent(predict(sac.pls.tr)[90],
                              42.94277,
                              tolerance = 1e-5)
  testVi <- setdiff(1:nrow(sacurine[["dataMatrix"]]), trainVi)
  predVn <- predict(sac.pls.tr, sacurine[["dataMatrix"]][testVi, ])
  testthat::expect_equivalent(predVn["HU_207"],
                              36.38991,
                              tolerance = 1e-5)
  
})

#### PLSDA_sacurine ####

test_that("PLSDA_sacurine", {
  
  data(sacurine)
  
  sac.plsda <- opls(sacurine[["dataMatrix"]],
                    matrix(sacurine[["sampleMetadata"]][, "gender"], ncol = 1),
                    permI = 0, fig.pdfC = "none", info.txtC = "none")
  testthat::expect_equivalent(getSummaryDF(sac.plsda)["Total", "Q2(cum)"],
                              0.584,
                              tolerance = 1e-3)
  sac.plsda <- opls(sacurine[["dataMatrix"]],
                    sacurine[["sampleMetadata"]][, "gender"],
                    permI = 0, fig.pdfC = "none", info.txtC = "none")
  testthat::expect_equivalent(getVipVn(sac.plsda)[5],
                              0.7034828,
                              tolerance = 1e-7)
  
  ## splitting the dataset between a reference and a test subsets
  
  sac.plsdaCrv <- opls(sacurine[["dataMatrix"]],
                       sacurine[["sampleMetadata"]][, "gender"],
                       subset = setdiff(1:nrow(sacurine[["dataMatrix"]]), 1:10),
                       permI = 0, fig.pdfC = "none", info.txtC = "none")
  testthat::expect_equivalent(getSummaryDF(sac.plsdaCrv)["Total", "RMSEP"],
                              0.275,
                              tolerance = 1e-3)
  ## for subset = "odd"
  ##     R2X(cum)  R2Y(cum)   Q2(cum)     RMSEE     RMSEP pre ort
  ## h2 0.1802609 0.7666581 0.5618918 0.2446342 0.3422395   2   0
  
})



#### PLSDA_sacurine_pareto ####

test_that("PLSDA_sacurine_pareto", {
  
  data(sacurine)
  
  sac.plsdaPar <- opls(sacurine[["dataMatrix"]],
                       matrix(sacurine[["sampleMetadata"]][, "gender"], ncol = 1),
                       scaleC = "pareto",
                       permI = 0, fig.pdfC = "none", info.txtC = "none")
  testthat::expect_equivalent(getSummaryDF(sac.plsdaPar)["Total", "Q2(cum)"],
                              0.521,
                              tolerance = 1e-3)
  sac.plsdaPar <- opls(sacurine[["dataMatrix"]],
                       sacurine[["sampleMetadata"]][, "gender"],
                       scaleC = "pareto",
                       permI = 0, fig.pdfC = "none", info.txtC = "none")
  testthat::expect_equivalent(getVipVn(sac.plsdaPar)[5],
                              0.7710519,
                              tolerance = 1e-7)
  
})

#### PLSDA_multiclass ####

test_that("PLSDA_multiclass", {
  
  data(sacurine)
  ageVn <- sacurine[["sampleMetadata"]][, "age"]
  ageVc <- ifelse(ageVn < 25, "young",
                  ifelse(ageVn < 40, NA,
                         ifelse(ageVn < 45, "normal",
                                ifelse(ageVn < 57, NA,
                                       "old"))))
  
  sacMN <- sacurine[["dataMatrix"]][!is.na(ageVc), ]
  ageVc <- ageVc[!is.na(ageVc)]
  
  oplsda.mul <- opls(sacMN, ageVc, predI = 2,
                     fig.pdfC = "none", info.txtC = "none")
  
  testthat::expect_equivalent(getSummaryDF(oplsda.mul)[, "Q2(cum)"],
                              0.0394,
                              tolerance = 1e-3)
  
  testthat::expect_equivalent(length(getVipVn(oplsda.mul, orthoL = TRUE)),
                              0,
                              tolerance = 1e-10)
  
  testthat::expect_error(opls(sacMN, ageVc, predI = 1, orthoI = NA,
                              fig.pdfC = "none", info.txtC = "none"),
                         silent = TRUE)
  
  oplsda.mul.par <- opls(sacMN, ageVc,
                         predI = 2,
                         scaleC = "pareto",
                         fig.pdfC = "none", info.txtC = "none")
  
  testthat::expect_equivalent(oplsda.mul.par@modelDF["p2", "R2Y(cum)"],
                              0.442,
                              tolerance = 1e-3)
  
  testthat::expect_equivalent(getSummaryDF(oplsda.mul.par)[, "Q2(cum)"],
                              0.0531,
                              tolerance = 1e-3)
  
})

####    OPLSDA_sacurine    ####

test_that("OPLSDA_sacurine", {
  
  data(sacurine)
  
  sac.oplsda <- opls(sacurine[["dataMatrix"]],
                     matrix(sacurine[["sampleMetadata"]][, "gender"],
                            ncol = 1),
                     predI = 1,
                     orthoI = 1,
                     permI = 0, fig.pdfC = "none", info.txtC = "none")
  testthat::expect_equivalent(getScoreMN(sac.oplsda, orthoL = TRUE)[1, 1],
                              4.980604,
                              tolerance = 1e-7)
  ##     R2X(cum) R2Y(cum) Q2(cum) RMSEE pre ort
  ## sum    0.185    0.668   0.554 0.289   1   1
  
  testthat::expect_equivalent(getVipVn(sac.oplsda)[2],
                              1.551154,
                              tolerance = 1e-6)
  
  testthat::expect_equivalent(getVipVn(sac.oplsda, orthoL = TRUE)[3],
                              1.642173,
                              tolerance = 1e-6)
  
  ## permutation testing
  
  sac.oplsdaPer <- opls(sacurine[["dataMatrix"]],
                        sacurine[["sampleMetadata"]][, "gender"],
                        predI = 1,
                        orthoI = 1,
                        permI = 10,
                        fig.pdfC = "none", info.txtC = "none")
  testthat::expect_equivalent(sac.oplsdaPer@suppLs[["permMN"]][1, 1],
                              0.185,
                              tolerance = 1e-3)
  
  ## number of components
  
  sac.oplsda <- opls(sacurine[["dataMatrix"]],
                     matrix(sacurine[["sampleMetadata"]][, "gender"],
                            ncol = 1),
                     predI = 1,
                     orthoI = NA,
                     permI = 0, fig.pdfC = "none", info.txtC = "none")
  testthat::expect_equivalent(getSummaryDF(sac.oplsda)[, "ort"],
                              2,
                              tolerance = 1e-3)
  
  sac.oplsda <- opls(sacurine[["dataMatrix"]],
                     matrix(sacurine[["sampleMetadata"]][, "age"],
                            ncol = 1),
                     predI = 1,
                     orthoI = NA,
                     permI = 0, fig.pdfC = "none", info.txtC = "none")
  testthat::expect_equivalent(getSummaryDF(sac.oplsda)[, "ort"],
                              1,
                              tolerance = 1e-3)
  
})

#### OPLSDA_subset ####

test_that("OPLSDA_subset", {
  ## bug fixed following T. Souza remark
  
  tiagoMN <- c(10, 11, 12, 11, 10, 11, 10, 11, 10, 11, 10, 11, 10, 11, 10, 11, 10, 11, 10, 8, 10, 11, 10, 12, 10, 11, 10, 11, 10, 11, 11, 10, 11, 10, 11, 10, 11, 10, 12, 10, 5, 6, 5, 6, 7, 6, 5, 6, 5, 6, 15, 16, 15, 16, 15, 16, 15, 14, 15, 15, 16, 16, 16, 14, 16, 15, 16, 14, 16, 14, 2, 3, 2, 3, 22, 3, 4, 3, 2, 3, 3, 7, 3, 7, 3, 7, 4, 7, 3, 7, 3, 7, 3, 7, 3, 4, 3, 7, 3, 7, 2, 8, 2, 8, 20, 8, 2, 4, 2, 8, 7, 3, 7, 3, 4, 3, 7, 3, 7, 7, 3, 7, 3, 7, 3, 7, 3, 7, 6, 8, 25, 29, 25, 24, 25, 29, 25, 29, 25, 25, 23, 27, 23, 27, 23, 27, 23, 27, 23, 25, 20, 21, 22, 21, 20, 21, 20, 21, 20, 21)
  dim(tiagoMN) <- c(20, 8)
  dimnames(tiagoMN) <- list(1:20, letters[1:8])
  tiagoFc <- rep(c("P", "F"), each = 10)
  
  opls(tiagoMN, tiagoFc, predI = 1, orthoI = 1, subset = "odd",
       fig.pdfC = "none", info.txtC = "none")
  
  oplExc <- opls(tiagoMN, tiagoFc, subset = seq(1, 20, by = 2),
                 fig.pdfC = "none", info.txtC = "none")
  testthat::expect_error(plot(oplExc, parLabVc = paste0("s", seq(1, 20, by = 2))),
                         silent = TRUE)
  
})

#### OPLSDA_ns ####

test_that("OPLSDA_ns", {
  
  datMN <- c(76043, 412165, 44943, 27242, 436566, 173175, 242549, 57066, 559869, 3732, 339188, 471368, 262271, 127285, 451270, 212500, 79673, NA, 891129, 43907, 30689, 6877586, 52217, 3158, 10789748, 229568, 4763576, 3878773, 976436, 831937, 608298, 1605075, 72021, 442510, 1107705, 1464339, 31250, 2724553, 72900, 32742, 47259, 544877, 60885, 34582, 529874, 168264, 176500, 76457, 610110, 16262, 279156, 524468, 451573, 591487, 433529, 161069, 214392, 13781, 1580343, 39315, 357351, 1030464, 301983, 67604, 306862, 1028110, 1530493, 270027, 1378535, 289677, 808334, 1132813, 871209, 895435, 715190, 1563158, 784738, 146195, 994336, 239030, 483755, 579287, 1132413, 157113, 1577570, 1469735, 1085454, 477909, 814755, 245417, 610681, 763706, 2406336, 827531, 992508, 569605, 355321, 150259, 1334200, 271010, 2644620, 727587, 1661412, 619181, 136278, 2755434, 593863, 837865, 3526136, 2003278, 1608814, 3446611, 1941527, 113937, 3132404, 2893445, 2092753, 1034666, 1517319, 841661, 250551, 1046138, 456162, 159386, 1013302, 808657, 614370, 250403, 768004, 242085, 504108, 1014041, 1362408, 1057660, 1110050, 566050, 411886, 142233, 1992420, 284775, 560002, 771533, 575790, 392284, 888498, 785428, 645785, 591569, 960658, 910201, 639437, 1092885, 1409045, 2292023, 1246459, 1945577, 710519, 773384, 1061418, 622898, 34236, 58249, 85944, NA, 342102, 129886, 175800, 13154, 230242, NA, 440223, 315368, 10657, 419508, 48673, 28361, 514579, 23108, 867108, 73831, 1252089, 2547452, 905408, 371059, 4983588, 5140022, 2658555, 814523, 2558923, 859466, 4184204, 3865723, 3236644, 2615560, 3820724, 3577833, 2295288, 625924, 7517724, 1341900, 2569205, 26023086, 1604999, 430453, 8103558, 26222916, 257139, 675754, 59906109, 263055, 31151730, 18648127, 14989438, 1554658, 20249262, 5588731, 871010, 15920, 9120781, 44276, 747080, 13420742, 595872, 1172376, 7172632, 3143654, 4059767, 1433702, 5593888, 5402629, 2477288, 3346077, 4230072, 7621236, 8960828, 10335722, 7037373, 1574738, 3359238, 2540044, 374028, 1144386, 539206, 178517, 1046190, 959381, 605191, 310260, 1253319, 477259, 477995, 825691, 1157093, 1089284, 1411802, 1020206, 782673, 346761, 1824553, 387811, 53304, 319783, 280560, 85009, 1333877, 556003, 590779, 209285, 342532, 198512, 569970, 525240, 246282, 1140422, 542345, 1171008, 827723, 222953, 438839, 85554, 368600, 616555, 94936, 622468, 180988, 293988, 352855, 767894, 268331, 167246, 310918, 1248919, 577184, 10985, 335711, 403815, 80614, 63393, 454489, 616061)
  
  dim(datMN) <- c(20, 15)
  
  ageVn <- c(41, 41, 52, 24, 55, 46, 61, 53, 23, 50, 33, 48, 42, 35, 26, 35, 60, 42, 31, 27)
  
  data_age.oplsda <- opls(datMN,
                          ageVn,
                          predI = 1,
                          orthoI = NA,
                          permI = 0,
                          fig.pdfC = "none", info.txtC = "none")
  
  testthat::expect_equivalent(cumprod(dim(getSummaryDF(data_age.oplsda)))[2],
                              0,
                              tolerance = 1e-15)
  
  data_age_pareto.oplsda <- opls(datMN,
                                 ageVn,
                                 predI = 1,
                                 orthoI = NA,
                                 permI = 0,
                                 scaleC = "pareto",
                                 fig.pdfC = "none", info.txtC = "none")
  
  testthat::expect_equivalent(cumprod(dim(getSummaryDF(data_age_pareto.oplsda)))[2],
                              0,
                              tolerance = 1e-15)
  
})


####    SummarizedExperiment    ####

test_that("SummarizedExperiment", {
  
  data(sacurine)
  
  sac.se <- sacurine[["se"]]
  
  ## PCA
  
  sac.se <- opls(sac.se, fig.pdfC = "none", info.txtC = "none")
  sac.pca <- getOpls(sac.se)[["PCA"]]
  
  testthat::expect_equivalent(getSummaryDF(sac.pca)["Total", "R2X(cum)"],
                              0.501,
                              tolerance = 1e-3)
  
  ## PLS-DA
  
  sac.se <- opls(sac.se, "gender", fig.pdfC = "none", info.txtC = "none")
  sac_gender.plsda <- getOpls(sac.se)[["gender_PLSDA"]]
  
  testthat::expect_equivalent(getSummaryDF(sac_gender.plsda)["Total", "Q2(cum)"],
                              0.584,
                              tolerance = 1e-3)
  
  testthat::expect_equivalent(SummarizedExperiment::colData(sac.se)["HU_011", "gender_PLSDA_xscor-p2"],
                              3.396,
                              tolerance = 1e-3)
  testthat::expect_equivalent(SummarizedExperiment::rowData(sac.se)["1-Methyluric acid", "gender_PLSDA_VIP"],
                              0.994,
                              tolerance = 1e-3)
  
  ## OPLS-DA
  
  sac.se <- opls(sac.se, "gender", predI = 1, orthoI = NA,
                 fig.pdfC = "none", info.txtC = "none")
  sac_gender.oplsda <- getOpls(sac.se)[["gender_OPLSDA"]]
  
  testthat::expect_equivalent(getSummaryDF(sac_gender.oplsda)["Total", "Q2(cum)"],
                              0.602,
                              tolerance = 1e-3)
  
})


####    MultiAssayExperiment    ####

test_that("MultiAssayExperiment", {
  
  ## Application to a MultiAssayExperiment
  data(NCI60)
  nci.mae <- NCI60[["mae"]]
  # Restricting to the 'ME' and 'LE' cancer types and to the 'agilent' and 'hgu95' datasets
  nci.mae <- nci.mae[, nci.mae$cancer %in% c("ME", "LE"), c("agilent", "hgu95")]
  
  ## PCA
  
  nci.mae <- opls(nci.mae, fig.pdfC = "none", info.txtC = "none")
  
  testthat::expect_equivalent(getSummaryDF(getOpls(nci.mae[["agilent"]])[["PCA"]])["Total", "R2X(cum)"],
                              0.521,
                              tolerance = 1e-3)
  
  testthat::expect_equivalent(SummarizedExperiment::colData(nci.mae[["agilent"]])["LE.CCRF_CEM", "PCA_xscor-p1"],
                              -15.04493,
                              tolerance = 1e-5)
  
  ## PLS-DA
  
  nci.mae <- opls(nci.mae, "cancer", fig.pdfC = "none", info.txtC = "none")
  
  testthat::expect_equivalent(getSummaryDF(getOpls(nci.mae[["agilent"]])[["cancer_PLSDA"]])["Total", "Q2(cum)"],
                              0.906,
                              tolerance = 1e-3)
  
  testthat::expect_equivalent(SummarizedExperiment::rowData(nci.mae[["agilent"]])["ST8SIA1", "cancer_PLSDA_VIP"],
                              1.568032,
                              tolerance = 1e-5)
  
})


####    ExpressionSet    ####

test_that("ExpressionSet", {
  
  data(sacurine)
  
  sacSet <- sacurine[["eset"]]
  
  ## PCA
  
  sacPca <- opls(sacSet, fig.pdfC = "none", info.txtC = "none")
  
  testthat::expect_equivalent(getSummaryDF(sacPca)["Total", "R2X(cum)"],
                              0.501,
                              tolerance = 1e-3)
  
  ## PLS-DA
  
  sacPlsda <- opls(sacSet, "gender", fig.pdfC = "none", info.txtC = "none")
  
  testthat::expect_equivalent(getSummaryDF(sacPlsda)["Total", "Q2(cum)"],
                              0.584,
                              tolerance = 1e-3)
  
  sacSet <- getEset(sacPlsda)
  
  testthat::expect_equivalent(Biobase::pData(sacSet)["HU_011", "gender_PLSDA_xscor-p2"],
                              3.396,
                              tolerance = 1e-3)
  testthat::expect_equivalent(Biobase::fData(sacSet)["1-Methyluric acid", "gender_PLSDA_VIP"],
                              0.994,
                              tolerance = 1e-3)
  
  ## OPLS-DA
  
  sacOplsda <- opls(sacSet, "gender", predI = 1, orthoI = NA,
                    fig.pdfC = "none", info.txtC = "none")
  
  testthat::expect_equivalent(getSummaryDF(sacOplsda)["Total", "Q2(cum)"],
                              0.602,
                              tolerance = 1e-3)
  
})


####    MultiDataSet    ####

test_that("MultiDataSet", {
  
  ## Application to a MultiDataSet
  data(NCI60)
  nci.mds <- NCI60[["mds"]]
  # Restricting to the 'agilent' and 'hgu95' datasets
  nci.mds <- nci.mds[, c("agilent", "hgu95")]
  # Restricting to the 'ME' and 'LE' cancer types
  sampleNamesVc <- Biobase::sampleNames(nci.mds[["agilent"]])
  cancerTypeVc <- Biobase::pData(nci.mds[["agilent"]])[, "cancer"]
  nci.mds <- nci.mds[sampleNamesVc[cancerTypeVc %in% c("ME", "LE")], ]  

  # Principal Component Analysis of each data set
  nciPca <- opls(nci.mds, fig.pdfC = "none", info.txtC = "none")
  
  testthat::expect_equivalent(getSummaryDF(nciPca@oplsLs[["agilent"]])["Total", "R2X(cum)"],
                              0.521,
                              tol = 1e-3)
  
  # Getting the updated MultiDataSet (now including scores and loadings)
  nci.mds <- getMset(nciPca)
  testthat::expect_equivalent(Biobase::fData(nci.mds)[["hgu95"]]["USE1", "PCA_xload-p1"],
                              0.003122651,
                              tol = 1e-8)
  
  # Building PLS-DA models for the cancer type, and getting back the updated MultiDataSet
  nciPlsda <- opls(nci.mds, "cancer", predI = 2, fig.pdfC = "none", info.txtC = "none")
  testthat::expect_equivalent(getSummaryDF(nciPlsda@oplsLs[["hgu95"]])["Total", "Q2(cum)"],
                              0.908,
                              tol = 1e-3)  
  
})


#### plot ####

test_that("plot", {
  
  data(sacurine)
  
  for (typeC in c("correlation", "outlier", "overview",
                  "permutation", "predict-train","predict-test",
                  "summary", "x-loading", "x-score", "x-variance",
                  "xy-score", "xy-weight")) {
    
    if (grepl("predict", typeC)) {
      subset <- "odd"
    } else
      subset <- NULL
    
    opLs <- opls(sacurine[["dataMatrix"]],
                 sacurine[["sampleMetadata"]][, "gender"],
                 predI = ifelse(typeC != "xy-weight", 1, 2),
                 orthoI = ifelse(typeC != "xy-weight", 1, 0),
                 permI = ifelse(typeC == "permutation", 10, 0),
                 subset = subset,
                 fig.pdfC = "none", info.txtC = "none")
    
    plot(opLs, typeVc = typeC)
    
  }
  
  ## (O)PLS-DA: Turning display of ellipses off in case parAsColFcVn is numeric
  plot(opLs,
       parAsColFcVn = sacurine[["sampleMetadata"]][, "age"])
  
  ## Converting 'parAsColFcVn' character into a factor (with warning)
  plot(opLs,
       parAsColFcVn = as.character(sacurine[["sampleMetadata"]][, "gender"]))
  
})

#### print ####

test_that("print", {
  
  data(sacurine)
  pcaLs <- opls(sacurine[["dataMatrix"]], predI = 2,
                fig.pdfC = "none", info.txtC = "none")
  print(pcaLs)
  plsLs <- opls(sacurine[["dataMatrix"]],
                sacurine[["sampleMetadata"]][, "gender"],
                predI = 2, fig.pdfC = "none", info.txtC = "none")
  print(plsLs)
  
  
})

#### fromW4M ####

test_that("fromW4M",{
  
  sacSet <- fromW4M(system.file("extdata/", package = "ropls"))
  
  testthat::expect_equivalent(Biobase::exprs(sacSet)["Tryptophan", "HU_020"], 4.594285, tol = 1e-6)
  
  testthat::expect_equivalent(Biobase::pData(sacSet)["HU_017", "bmi"], 23.03, tol = 1e-2)
  
  testthat::expect_identical(Biobase::fData(sacSet)["Taurine", "hmdb"], "HMDB00251")
  
  testthat::expect_identical(Biobase::sampleNames(sacSet)[3], "HU_015")
  
  testthat::expect_identical(Biobase::featureNames(sacSet)[5], "X1.3.Dimethyluric.acid")
  
})




#### imageF ####

test_that("imageF", {
  data(sacurine)
  imageF(sacurine[['dataMatrix']])
})

#### strF ####

test_that("strF", {
  
  data(foods)
  
  strF(foods)
  strF(foods, border = 3)
  
  fatML <- matrix(TRUE, nrow = 1, ncol = 1000)
  strF(fatML, bigMarkC = "'")
  
  testMC <- matrix("a", nrow = 10, ncol = 10)
  strF(testMC)
  
  data(sacurine)
  
  strF(sacurine[["dataMatrix"]])
  
  strF(sacurine[["sampleMetadata"]])
  
})


