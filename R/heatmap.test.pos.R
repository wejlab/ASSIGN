heatmap.test.pos <- function(testData, Delta_pos, trainingLabel, testLabel=NULL,
                             Delta_cutoff = 0.95, coef_test, geneList = NULL){
  
  if (is.null(geneList)){
    nPath <- length(trainingLabel) - 1
    pathName <- names(trainingLabel)[-1]
  } else {
    nPath <- length(geneList)
    pathName <- names(geneList)
  }
  
  diffGeneList <- vector("list")
  for (i in 1:nPath){
    ##the cutoff can be modified
    diffGeneList[[i]] <- row.names(testData)[Delta_pos[,i] >= Delta_cutoff]
  }
  if (!is.null(testLabel)){
    cc <- as.numeric(as.factor(testLabel))
  } 
  grDevices::pdf("signature_heatmap_testset_posterior.pdf")
  for (i in 1:nPath){
    tmp <- match(diffGeneList[[i]], row.names(testData))
    path <- testData[tmp, ]
    #do not print heatmaps if there is 1 or fewer genes in the heatmap
    if(nrow(path) <= 1){
      next
    }
    if (!is.null(testLabel)){
      stats::heatmap(as.matrix(path[,order(coef_test[,i])]),Colv=NA,scale="row",
                     ColSideColors=as.character(cc[order(coef_test[,i])]),
                     col=gplots::bluered(128), margins = c(10,10),
                     main=paste(pathName[i],"signature",sep=" "))
    } else {
      stats::heatmap(as.matrix(path[,order(coef_test[,i])]),Colv=NA,scale="row",
                     col=gplots::bluered(128), margins = c(10,10),
                     main=paste(pathName[i],"signature",sep=" "))
    }
  }
  invisible(grDevices::dev.off())
}
