# Tests the ASSIGN preprocess function on the lung dataset
library(ASSIGN)
library(testthat)
data(testData1)
data(trainingData1)
data(geneList1)

context("Tests the ASSIGN preprocess function on the lung data")

test_that("test assign.preprocess on the lung data with multiple parameters", {
  set.seed(0)
  training_label1 <- list(control = list(bcat = 1:10, e2f3 = 1:10, myc = 1:10,
                                         ras = 1:10, src = 1:10), bcat = 11:19,
                          e2f3 = 20:28, myc = 29:38, ras = 39:48, src = 49:55)

  processed_data <- assign.preprocess(trainingData = trainingData1,
                                      testData = testData1,
                                      trainingLabel = training_label1,
                                      geneList = geneList1, n_sigGene = NULL,
                                      progress_bar = FALSE)

  line750 <- c(0.3609345, -0.0274970, 0.0000617, 1.0123866, 0.9085872)
  expect_equal(round(as.numeric(processed_data$"S_matrix"[750, ]), 7), line750)
})

test_that("exclude genes don't affect the number of genes in the signature", {
  set.seed(0)
  trainData_sub <- trainingData1[, c(grep("Control", colnames(trainingData1)),
                                     grep("B[Cc]at", colnames(trainingData1)))]
  trainLabel_sub <- list(control = list(bcat = 1:10), bcat = 11:19)
  numgenes <- 100
  processed_data <- assign.preprocess(trainingData = trainData_sub,
                                      testData = testData1,
                                      trainingLabel = trainLabel_sub,
                                      geneList = NULL,
                                      n_sigGene = numgenes,
                                      geneselect_iter = 50,
                                      geneselect_burn_in = 10,
                                      excludeGenes = list(bcat = c("201891_s_at",
                                                                   "1555340_x_at",
                                                                   "1555339_at",
                                                                   "1558048_x_at",
                                                                   "217369_at",
                                                                   "212560_at",
                                                                   "201291_s_at",
                                                                   "222227_at",
                                                                   "217356_s_at",
                                                                   "208835_s_at",
                                                                   "202832_at")),
                                      progress_bar = FALSE)
  expect_equal(nrow(processed_data$S_matrix), numgenes)
})

test_that("excluding all available genes causes an error", {
  set.seed(0)
  trainData_sub <- trainingData1[, c(grep("Control", colnames(trainingData1)),
                                     grep("B[Cc]at", colnames(trainingData1)))]
  trainLabel_sub <- list(control = list(bcat = 1:10), bcat = 11:19)
  numgenes <- 1000
  expect_error(assign.preprocess(trainingData = trainData_sub,
                                 testData = testData1,
                                 trainingLabel = trainLabel_sub,
                                 geneList = NULL,
                                 n_sigGene = numgenes,
                                 geneselect_iter = 50,
                                 geneselect_burn_in = 10,
                                 excludeGenes = list(bcat = c("201891_s_at",
                                                              "1555340_x_at",
                                                              "1555339_at",
                                                              "1558048_x_at",
                                                              "217369_at",
                                                              "212560_at",
                                                              "201291_s_at",
                                                              "222227_at",
                                                              "217356_s_at",
                                                              "208835_s_at",
                                                              "202832_at")),
                                 progress_bar = FALSE),
               "There aren't enough genes available in the data. Check excludeGenes and n_sigGene.")
})
