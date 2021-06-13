#' ASSIGN All-in-one function
#'
#' The assign.wrapper function integrates the assign.preprocess, assign.mcmc,
#' assign.summary, assign.output, assign.cv.output functions into one wrapper
#' function.
#'
#' The assign.wrapper function is an all-in-one function which outputs the
#' necessary results for basic users. For users who need more
#' intermediate results for model diagnosis, it is better to run the
#' assign.preprocess, assign.mcmc, assign.convergence, assign.summary functions
#' separately and extract the output values from the returned list objects of
#' those functions.
#'
#' @param trainingData The genomic measure matrix of training samples (i.g.,
#' gene expression matrix). The dimension of this matrix is probe number x
#' sample number. The default is NULL.
#' @param testData The genomic measure matrix of test samples (i.g., gene
#' expression matrix). The dimension of this matrix is probe number x sample
#' number.
#' @param trainingLabel The list linking the index of each training sample to a
#' specific group it belongs to. See examples for more information.
#' @param testLabel The vector of the phenotypes/labels of the test samples.
#' The default is NULL.
#' @param geneList The list that collects the signature genes of one/multiple
#' pathways. Every component of this list contains the signature genes
#' associated with one pathway. The default is NULL.
#' @param anchorGenes A list of genes that will be included in the signature
#' even if they are not chosen during gene selection.
#' @param excludeGenes A list of genes that will be excluded from the signature
#' even if they are chosen during gene selection.
#' @param n_sigGene The vector of the signature genes to be identified for one
#' pathway. n_sigGene needs to be specified when geneList is set NULL. The
#' default is NA. See examples for more information.
#' @param adaptive_B Logicals. If TRUE, the model adapts the
#' baseline/background (B) of genomic measures for the test samples. The
#' default is TRUE.
#' @param adaptive_S Logicals. If TRUE, the model adapts the signatures (S) of
#' genomic measures for the test samples. The default is FALSE.
#' @param mixture_beta Logicals. If TRUE, elements of the pathway activation
#' matrix are modeled by a spike-and-slab mixture distribution. The default is
#' TRUE.
#' @param outputDir The path to the directory to save the output files. The
#' path needs to be quoted in double quotation marks.
#' @param p_beta p_beta is the prior probability of a pathway being activated
#' in individual test samples. The default is 0.01.
#' @param theta0 The prior probability for a gene to be significant, given that
#' the gene is NOT defined as "significant" in the signature gene lists
#' provided by the user. The default is 0.05.
#' @param theta1 The prior probability for a gene to be significant, given that
#' the gene is defined as "significant" in the signature gene lists provided by
#' the user. The default is 0.9.
#' @param iter The number of iterations in the MCMC. The default is 2000.
#' @param burn_in The number of burn-in iterations. These iterations are
#' discarded when computing the posterior means of the model parameters. The
#' default is 1000.
#' @param sigma_sZero Each element of the signature matrix (S) is modeled by a
#' spike-and-slab mixture distribution. Sigma_sZero is the variance of the
#' spike normal distribution. The default is 0.01.
#' @param sigma_sNonZero Each element of the signature matrix (S) is modeled by
#' a spike-and-slab mixture distribution. Sigma_sNonZero is the variance of
#' the slab normal distribution. The default is 1.
#' @param S_zeroPrior Logicals. If TRUE, the prior distribution of signature
#' follows a normal distribution with mean zero. The default is TRUE.
#' @param pctUp By default, ASSIGN bayesian gene selection chooses the
#' signature genes with an equal fraction of genes that increase with pathway
#' activity and genes that decrease with pathway activity. Use the pctUp
#' parameter to modify this fraction. Set pctUP to NULL to select the most
#' significant genes, regardless of direction. The default is 0.5
#' @param geneselect_iter The number of iterations for bayesian gene selection.
#' The default is 500.
#' @param geneselect_burn_in The number of burn-in iterations for bayesian gene
#' selection. The default is 100
#' @param outputSignature_convergence Create a pdf of the MCMC chain. The
#' default is FALSE.
#' @param ECM Logicals. If TRUE, ECM algorithm, rather than Gibbs sampling, is
#' applied to approximate the model parameters. The default is FALSE.
#' @param progress_bar Display a progress bar for MCMC and gene selection.
#' Default is TRUE.
#' @param override_S_matrix Replace the S_matrix created by assign.preprocess
#' with the matrix provided in override_S_matrix. This can be used to indicate
#' the expected directions of genes in a signature if training data is not
#' provided.
#' @return The assign.wrapper returns one/multiple pathway activity for each
#' individual training sample and test sample, scatter plots of pathway
#' activity for each individual pathway in the training and test data,
#' heatmap plots for gene expression signatures for each individual pathway,
#' heatmap plots for the gene expression of the prior and posterior
#' signatures (if adaptive_S equals TRUE) of each individual pathway in the test
#' data
#' @author Ying Shen and W. Evan Johnson
#' @examples
#'
#' \dontshow{
#' tempdir <- file.path(tempdir(), "assign_wrapper")
#' }
#' data(trainingData1)
#' data(testData1)
#' data(geneList1)
#'
#' trainingLabel1 <- list(control = list(bcat=1:10, e2f3=1:10, myc=1:10,
#'                                       ras=1:10, src=1:10),
#'                        bcat = 11:19, e2f3 = 20:28, myc= 29:38, ras = 39:48,
#'                        src = 49:55)
#' testLabel1 <- rep(c("subtypeA","subtypeB"), c(53,58))
#'
#' assign.wrapper(trainingData=trainingData1, testData=testData1,
#'                trainingLabel=trainingLabel1, testLabel=testLabel1,
#'                geneList=geneList1, adaptive_B=TRUE, adaptive_S=FALSE,
#'                mixture_beta=TRUE, outputDir=tempdir, p_beta=0.01,
#'                theta0=0.05, theta1=0.9, iter=20, burn_in=10)
#'
#' @export assign.wrapper
assign.wrapper <- function(trainingData = NULL, testData, trainingLabel,
                           testLabel = NULL, geneList = NULL,
                           anchorGenes = NULL, excludeGenes = NULL,
                           n_sigGene = NA, adaptive_B = TRUE,
                           adaptive_S = FALSE, mixture_beta = TRUE, outputDir,
                           p_beta = 0.01, theta0 = 0.05, theta1 = 0.9,
                           iter = 2000, burn_in = 1000, sigma_sZero = 0.01,
                           sigma_sNonZero = 1, S_zeroPrior=FALSE, pctUp=0.5,
                           geneselect_iter=500, geneselect_burn_in=100,
                           outputSignature_convergence = FALSE, ECM = FALSE,
                           progress_bar = TRUE, override_S_matrix = NULL) {

  if (!dir.exists(outputDir)) {
    dir.create(outputDir)
  }

  if (any(file.exists(
    file.path(outputDir, "parameters.yaml"),
    file.path(outputDir, "signature_gene_list_prior.csv"),
    file.path(outputDir, "pathway_activity_trainingset.csv"),
    file.path(outputDir, "pathway_activity_testset.csv"),
    file.path(outputDir, "Signature_convergence.pdf"),
    file.path(outputDir, "signature_heatmap_trainingset.pdf"),
    file.path(outputDir, "signature_heatmap_testset_prior.pdf"),
    file.path(outputDir, "signature_heatmap_testset_posterior.pdf"),
    file.path(outputDir, "posterior_delta.csv"),
    file.path(outputDir, "pathway_activity_scatterplot_trainingset.pdf"),
    file.path(outputDir, "pathway_activity_scatterplot_testset.pdf"),
    file.path(outputDir, "pathway_activity_boxplot_testset.pdf"),
    file.path(outputDir, "output.rds")))) {
    stop("Output files already exist. Delete the following files to run assign.cv.output():\n",
         file.path(outputDir, "parameters.yaml"), "\n",
         file.path(outputDir, "signature_gene_list_prior.csv"), "\n",
         file.path(outputDir, "pathway_activity_trainingset.csv"), "\n",
         file.path(outputDir, "pathway_activity_testset.csv"), "\n",
         file.path(outputDir, "Signature_convergence.pdf"), "\n",
         file.path(outputDir, "signature_heatmap_trainingset.pdf"), "\n",
         file.path(outputDir, "signature_heatmap_testset_prior.pdf"), "\n",
         file.path(outputDir, "signature_heatmap_testset_posterior.pdf"), "\n",
         file.path(outputDir, "posterior_delta.csv"), "\n",
         file.path(outputDir, "pathway_activity_scatterplot_trainingset.pdf"), "\n",
         file.path(outputDir, "pathway_activity_scatterplot_testset.pdf"), "\n",
         file.path(outputDir, "pathway_activity_boxplot_testset.pdf"), "\n",
         file.path(outputDir, "output.rds"))
  }

  if (is.null(geneList)) {
    pathName <- names(trainingLabel)[-1]
  } else {
    pathName <- names(geneList)
  }
  processed.data <- assign.preprocess(trainingData, testData, anchorGenes,
                                      excludeGenes, trainingLabel, geneList,
                                      n_sigGene, theta0, theta1, pctUp = pctUp,
                                      geneselect_iter = geneselect_iter,
                                      geneselect_burn_in = geneselect_burn_in,
                                      progress_bar = progress_bar)
  if (!is.null(trainingData)) {
    message("Estimating model parameters in the training dataset...")
    mcmc.chain.trainingData <- assign.mcmc(Y = processed.data$trainingData_sub,
                                           Bg = processed.data$B_vector,
                                           X = processed.data$S_matrix,
                                           Delta_prior_p = processed.data$Pi_matrix,
                                           iter = iter,
                                           sigma_sZero = sigma_sZero,
                                           sigma_sNonZero = sigma_sNonZero,
                                           S_zeroPrior = S_zeroPrior,
                                           adaptive_B = FALSE,
                                           adaptive_S = FALSE,
                                           mixture_beta = TRUE,
                                           ECM = ECM,
                                           progress_bar = progress_bar)
    mcmc.pos.mean.trainingData <- assign.summary(test = mcmc.chain.trainingData,
                                                 burn_in = burn_in, iter = iter,
                                                 adaptive_B = FALSE,
                                                 adaptive_S = FALSE,
                                                 mixture_beta = TRUE)

  }
  message("Estimating model parameters in the test dataset...")

  if (is.null(trainingData) & !is.null(override_S_matrix)) {
    Smat <- override_S_matrix
    message("Using override_S_matrix")
  } else {
    Smat <- processed.data$S_matrix
  }

  mcmc.chain.testData <- assign.mcmc(Y = processed.data$testData_sub,
                                     Bg = processed.data$B_vector,
                                     X = Smat,
                                     Delta_prior_p = processed.data$Pi_matrix,
                                     iter = iter, sigma_sZero = sigma_sZero,
                                     sigma_sNonZero = sigma_sNonZero,
                                     S_zeroPrior = S_zeroPrior,
                                     adaptive_B = adaptive_B,
                                     adaptive_S = adaptive_S,
                                     mixture_beta = mixture_beta,
                                     p_beta = p_beta,
                                     ECM = ECM)
  mcmc.pos.mean.testData <- assign.summary(test = mcmc.chain.testData,
                                           burn_in = burn_in, iter = iter,
                                           adaptive_B = adaptive_B,
                                           adaptive_S = adaptive_S,
                                           mixture_beta = mixture_beta)

  message("Outputing results...")
  if (mixture_beta) {
    if (!is.null(trainingData)) {
      coef_train <- mcmc.pos.mean.trainingData$kappa_pos
    }
    coef_test <- mcmc.pos.mean.testData$kappa_pos
  } else {
    if (!is.null(trainingData)) {
      coef_train <- mcmc.pos.mean.trainingData$beta_pos
    }
    coef_test <- mcmc.pos.mean.testData$beta_pos
  }

  #Save ASSIGN parameters
  param <- data.frame(ASSIGN = "ASSIGN parameters",
                      version = as.character(utils::packageVersion("ASSIGN")),
                      n_sigGene = if (is.null(n_sigGene)) {NA} else n_sigGene,
                      adaptive_B = adaptive_B,
                      adaptive_S = adaptive_S,
                      mixture_beta = mixture_beta,
                      p_beta = p_beta,
                      theta0 = theta0,
                      theta1 = theta1,
                      iter = iter,
                      burn_in = burn_in,
                      sigma_sZero = sigma_sZero,
                      sigma_sNonZero = sigma_sNonZero,
                      S_zeroPrior = S_zeroPrior,
                      pctUp = pctUp,
                      geneselect_iter = geneselect_iter,
                      geneselect_burn_in = geneselect_burn_in,
                      ECM = ECM,
                      outputDir = outputDir)
  y1 <- yaml::as.yaml(param, column.major = FALSE)
  writeLines(y1, con = file.path(outputDir, "parameters.yaml"))

  #Include the gene list and prior coefficient
  if (!is.null(trainingData)) {
    rownames(coef_train) <- colnames(processed.data$trainingData_sub)
    colnames(coef_train) <- pathName
    utils::write.csv(processed.data$S_matrix, file = file.path(outputDir, "signature_gene_list_prior.csv"))
    utils::write.csv(coef_train, file = file.path(outputDir, "pathway_activity_trainingset.csv"))
  }

  rownames(coef_test) <- colnames(processed.data$testData_sub)
  colnames(coef_test) <- pathName
  utils::write.csv(coef_test, file = file.path(outputDir, "pathway_activity_testset.csv"))
  if (!is.null(trainingData)) {
    heatmap.train(diffGeneList = processed.data$diffGeneList,
                  trainingData, trainingLabel,
                  outPath = file.path(outputDir, "signature_heatmap_trainingset.pdf"))
  }

  heatmap.test.prior(diffGeneList = processed.data$diffGeneList,
                     testData, trainingLabel, testLabel, coef_test, geneList,
                     outPath = file.path(outputDir, "signature_heatmap_testset_prior.pdf"))

  if (adaptive_S) {
    heatmap.test.pos(testData = processed.data$testData_sub,
                     Delta_pos = mcmc.pos.mean.testData$Delta_pos,
                     trainingLabel, testLabel, Delta_cutoff = 0.95,
                     coef_test, geneList,
                     outPath = file.path(outputDir, "signature_heatmap_testset_posterior.pdf"))
    # Plot signature convergence (can be a very large PDF when performing many iterations)
    if (outputSignature_convergence) {
      grDevices::pdf(file.path(outputDir, "Signature_convergence.pdf"))
      graphics::plot(mcmc.chain.testData$S_mcmc)
      graphics::abline(h = 0, col = "red")
      invisible(grDevices::dev.off())
    }

    dimnames(mcmc.pos.mean.testData$Delta_pos) <- dimnames(processed.data$S_matrix)
    deltas <- cbind(Smat, processed.data$Delta_matrix,
                    mcmc.pos.mean.testData$S_pos,
                    mcmc.pos.mean.testData$Delta_pos)

    colnames(deltas) <- c(paste("Prior change in expression",
                                pathName, sep = ":"),
                          paste("Prior probability of inclusion",
                                pathName, sep = ":"),
                          paste("Posterior change in expression",
                                pathName, sep = ":"),
                          paste("Posterior probability of inclusion",
                                pathName, sep = ":"))
    delta_in <- NULL

    for (i in seq_len(ncol(deltas))) {
      delta_in[i] <- (strsplit(colnames(deltas), ":")[[i]][2])
    }

    utils::write.csv(round(deltas[, order(delta_in)], digits = 4),
                     file.path(outputDir, "posterior_delta.csv"), quote = FALSE)
  }

  if (!is.null(trainingData)) {
    scatter.plot.train(coef_train, trainingData, trainingLabel,
                       outPath = file.path(outputDir, "pathway_activity_scatterplot_trainingset.pdf"))
  }
  scatter.plot.test(coef_test, trainingLabel, testLabel, geneList,
                    outPath = file.path(outputDir, "pathway_activity_scatterplot_testset.pdf"))
  if (!is.null(testLabel)) {
    box.plot.test(coef_test, trainingLabel, testLabel, geneList,
                  outPath = file.path(outputDir, "pathway_activity_boxplot_testset.pdf"))
  }
  if (!is.null(trainingData)) {
    output.data <- list(processed.data = processed.data,
                        mcmc.pos.mean.trainingData = mcmc.pos.mean.trainingData,
                        mcmc.pos.mean.testData = mcmc.pos.mean.testData)
  }
  else {
    output.data <- list(processed.data = processed.data,
                        mcmc.pos.mean.testData = mcmc.pos.mean.testData)
  }
  saveRDS(output.data, file = file.path(outputDir, "output.rds"))
}
