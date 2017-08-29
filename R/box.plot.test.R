box.plot.test <- function(coef_test, trainingLabel, testLabel, geneList=NULL){
  if (is.null(geneList)){
    nPath <- length(trainingLabel) - 1
    pathName <- names(trainingLabel)[-1]
  } else {
    nPath <- length(geneList)
    pathName <- names(geneList)
  }

  grDevices::pdf("pathway_activity_boxplot_testset.pdf")
  for(i in 1:nPath){
    graphics::boxplot(coef_test[, i] ~ as.factor(testLabel),
                      main = paste("box-plot of", pathName[i],
                                   "pathway activity in test samples",
                                   sep = " "))
    graphics::points(jitter(as.numeric(as.factor(testLabel)), amount = 0.2),
                     coef_test[, i], pch = 17, col = 4)
  }
  invisible(grDevices::dev.off())
}
