library(ASSIGN)
library(testthat)
data(testData1)
data(trainingData1)
data(geneList1)

context("Test error messages that were previously unclear.")

test_that("test no trainingLabels", {
  expect_error(assign.wrapper(trainingData = trainingData1,
                              trainingLabel = NULL,
                              testData = testData1,
                              geneList = geneList1, n_sigGene = NULL,
                              outputDir = tempfile()), "You must supply trainingLabels when specifying trainingData!")
})
