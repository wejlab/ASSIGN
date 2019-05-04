scatter.plot.train <- function(coef_train, trainingData, trainingLabel){
  nPath <- length(trainingLabel) - 1
  trainL <- character(length = length(unlist(trainingLabel)))
  for (i in 1:(nPath + 1)){
    if (i == 1){
      x <- unique(trainingLabel[[i]])
      names(x) <- paste("control", seq_len(length(x)), sep = "")
      for (j in seq_len(length(x))){
        trainL[x[[j]]] <- rep(names(x)[j], length(x[[j]]))
      }
    } else {
      trainL[trainingLabel[[i]]] <- rep(names(trainingLabel)[i],
                                        length(trainingLabel[[i]]))
    }
  }
  trainL <- trainL[trainL != ""]

  grDevices::pdf("pathway_activity_scatterplot_trainingset.pdf")
  for (i in 1:nPath){
    HMEC_samples <- seq_len(ncol(trainingData))
    Pathway_strength_HMEC <- coef_train[, i]
    graphics::plot(HMEC_samples, Pathway_strength_HMEC, col = as.factor(trainL),
                   xlab = "HMEC sample",
                   ylab = paste(names(trainingLabel)[i + 1], "pathway activity",
                                sep = " "),
                   main = paste("Cross-validation in HMEC",
                                names(trainingLabel)[i + 1], "pathway",
                                sep = " "),
                   pch = 19, cex = 0.7)
    graphics::legend("topleft", legend = unique(trainL), pch = 19, cex = 0.7,
                     col = as.numeric(as.factor(unique(trainL))))
  }
  invisible(grDevices::dev.off())
}
