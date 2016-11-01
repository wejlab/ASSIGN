#' Perform the second step of ComBat
#'
#' The first ComBat step (on the signatures only) has already been
#' performed. This step performs batch correction on the test data,
#' using reference batch ComBat, to prepare the test data for ASSIGN
#' analysis.
#' 
#' This function downloads the training data from the internet, so an internet
#' connection is necessary
#'
#' @param testData The input test data to batch correct
#' @param pcaPlots a logical value indicating whether or not the function
#' should create PCA plots. The default is FALSE.
#' @param combat_train the ComBat training data data frame. If you do not have
#' this, the function will attempt to download it from the internet. Please
#' contact the developers if you have any issues with access to the file.
#'
#' @return A list of data.frames is returned, including control (GFP) and
#' signature data, as well as the batch corrected test data. This data can go
#' directly into the runassign.single and runassign.multi functions, or
#' subsetted to go directly into ASSIGN.
#'
#' @export ComBat.step2
ComBat.step2 <- function(testData, pcaPlots=FALSE, combat_train=NULL) {
  if(!("ref.batch" %in% names(as.list(args(sva::ComBat))))){
    stop("Installed version of sva: ", packageVersion("sva"), " does not have ref.batch option.\n",
         "Use devtools to install the github version of sva:\ndevtools::install_github('jtleek/sva-devel')")
  }
  if(is.null(combat_train)){
    combat_train_file <- tempfile(pattern="combat_train",fileext = ".rda")
    utils::download.file("https://dl.dropboxusercontent.com/u/62447/ASSIGN/combat_train.rda",
                         combat_train_file)
    load(combat_train_file)
    unlink(combat_train_file)
    rm(combat_train_file)
  }
  
  dat <- merge_drop(combat_train,testData)
  sub <- c(6,6,12,6,6,5,6,6,9,9,9,9,ncol(testData))
  bat <- c(rep(1,ncol(combat_train)),rep(2,ncol(testData)))
  if(pcaPlots){
    grDevices::pdf("pca_refcombat_twostep.pdf")
    pcaplot(dat,sub)
  }
  combat_expr1 <- sva::ComBat(dat=dat, batch=bat, mod=NULL, ref.batch=1)
  if(pcaPlots){
    pcaplot(combat_expr1,sub)
    invisible(grDevices::dev.off())
  }
  c_gfp      <- combat_expr1[,13:24]
  c_akt      <- combat_expr1[,25:30]
  c_bad      <- combat_expr1[,31:36]
  c_her2     <- combat_expr1[,37:41]
  c_igf1r    <- combat_expr1[,42:47]
  c_raf      <- combat_expr1[,48:53]
  c_egfr_gfp <- combat_expr1[,1:6]
  c_egfr     <- combat_expr1[,7:12]
  c_kras_gfp <- combat_expr1[,54:62]
  c_kraswt   <- combat_expr1[,63:71]
  c_krasqh   <- combat_expr1[,72:80]
  c_krasgv   <- combat_expr1[,81:89]
  c_test     <- combat_expr1[,(ncol(combat_train)+1):ncol(combat_expr1)]
  results <- list(gfp=c_gfp, akt=c_akt, bad=c_bad, her2=c_her2, igf1r=c_igf1r,
                  raf=c_raf, egfr_gfp=c_egfr_gfp, egfr=c_egfr,
                  kras_gfp=c_kras_gfp, kraswt=c_kraswt, krasqh=c_krasqh,
                  krasgv=c_krasgv, test=c_test)
  return(results)
}
