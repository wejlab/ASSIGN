#' Optimize GFRN gene lists lengths
#'
#' This function runs ASSIGN pathway prediction on gene list lengths from 5 to
#' 500 to find the optimum gene list length for the GFRN pathways by correlating
#' the ASSIGN predictions to a matrix of correlation data that you provide. This
#' function takes a long time to run because you are running ASSIGN many times
#' on many pathways, so I recommend parallelizing by pathway or running the
#' ASSIGN predictions first (long and parallelizable) and then running the
#' correlation step (quick) separately.
#'
#' @param indata The list of data frames from ComBat.step2
#' @param correlation A matrix of data to correlate ASSIGN predictions to.
#' The number of rows should be the same and in the same order as indata
#' @param correlationList A list that shows which columns of correlation should
#' be used for each pathway. See below for more details
#' @param run specifies the pathways to predict. The default list will
#' cause all eight pathways to be run in serial. Specify a pathway ("akt",
#' "bad", "egfr", etc.) or list of pathways to run those pathways only.
#' @param run_ASSIGN_only a logical value indicating if you want to run the
#' ASSIGN predictions only. Use this to parallelize ASSIGN runs across a compute
#' cluster or across compute threads
#' @param correlation_only a logical value indicating if you want to run the
#' correlation step only. The function will find the ASSIGN runs in the cwd and
#' optimize them based on the correlation data matrix.
#' @param keep_optimized_only a logical value indicating if you want to keep
#' all of the ASSIGN run results, or only the runs that provided the optimum
#' ASSIGN correlations. This will delete all directories in the current working
#' directory that match the pattern "_gene_list". The default is FALSE
#' @param pathway_lengths The gene list lengths that should be run. The default
#' is the 20 pathway lengths that were used in the paper, but this list can
#' be customized to which pathway lengths you are willing to accept
#' @param iter The number of iterations in the MCMC.
#' @param burn_in The number of burn-in iterations. These iterations are
#' discarded when computing the posterior means of the model parameters.
#'
#' @return ASSIGN runs are output to the current workingdirectory. This function
#' returns the correlation data and the optimized gene lists that you can use
#' with runassignGFRN to try these lists on other data.
#'
#' @examples
#' \dontrun{
#' testData <- read.table("https://drive.google.com/uc?authuser=0&id=1mJICN4z_aCeh4JuPzNfm8GR_lkJOhWFr&export=download",
#'                        sep='\t', row.names=1, header=1)
#' corData <- read.table("https://drive.google.com/uc?authuser=0&id=1MDWVP2jBsAAcMNcNFKE74vYl-orpo7WH&export=download",
#'                       sep='\t', row.names=1, header=1)
#' corData$negAkt <- -1 * corData$Akt
#' corData$negPDK1 <- -1 * corData$PDK1
#' corData$negPDK1p241 <- -1 * corData$PDK1p241
#'
#' corList <- list(akt=c("Akt","PDK1","PDK1p241"),
#'                 bad=c("negAkt","negPDK1","negPDK1p241"),
#'                 egfr=c("EGFR","EGFRp1068"),
#'                 her2=c("HER2","HER2p1248"),
#'                 igf1r=c("IGFR1","PDK1","PDK1p241"),
#'                 krasgv=c("EGFR","EGFRp1068"),
#'                 krasqh=c("EGFR","EGFRp1068"),
#'                 raf=c("MEK1","PKCalphap657","PKCalpha"))
#'
#' combat.data <- ComBat.step2(testData, pcaPlots = TRUE)
#'
#' optimization_results <- optimizeGFRN(combat.data, corData, corList)
#' }
#'
#' @export optimizeGFRN
#'
optimizeGFRN <- function(indata, correlation, correlationList,
                         run=c("akt", "bad", "egfr", "her2", "igf1r", "krasgv",
                               "krasqh", "raf"), run_ASSIGN_only=FALSE,
                         correlation_only=FALSE, keep_optimized_only=FALSE,
                         pathway_lengths=c(seq(5, 20, 5), seq(25, 275, 25),
                                           seq(300, 500, 50)), iter=100000,
                         burn_in=50000) {
  #check that correlationList and run list are identical
  if (!identical(names(correlationList), run)){
    stop("Make sure the run list and correlationList list names are identical and in the same order")
  }

  utils::data("gfrn_geneList", package = "ASSIGN", envir = environment())
  gfrn_geneList <- get("gfrn_geneList", envir = environment())

  # run the pathway predictions
  if (!(correlation_only)){
    for (curr_path in run){
      for (curr_len in pathway_lengths){
        geneList <- list()
        geneList[[curr_path]] <- c(gfrn_geneList[[paste(curr_path, "up", sep = "_")]][1:floor(curr_len / 2)],
                                   gfrn_geneList[[paste(curr_path, "down", sep = "_")]][1:ceiling(curr_len / 2)])
        message("Running: ", curr_path, " ", curr_len, " genes")
        runassignGFRN(indata, run = curr_path, optimized_geneList = geneList,
                      iter = iter, burn_in = burn_in)
      }
    }
  }

  #gather the results into a matrix
  results <- ASSIGN::gather_assign_results()

  #check the rownames, warnings if there are extras
  if (sum(!(rownames(results) %in% rownames(correlation))) != 0){
    warning(sum(!(rownames(results) %in% rownames(correlation))),
            " out of ", length(rownames(results)),
            " samples in input data not in correlation data:\n",
            paste(rownames(results)[!(rownames(results) %in% rownames(correlation))], collapse = ", "),
            "\nThese samples will be removed from the correlation results.")
  }
  if (sum(!(rownames(correlation) %in% rownames(results))) != 0){
    warning(sum(!(rownames(correlation) %in% rownames(results))),
            " out of ", length(rownames(correlation)),
            " samples in correlation data not in input data:\n",
            paste(rownames(correlation)[!(rownames(correlation) %in% rownames(results))], collapse = ", "),
            "\nThese samples will be removed from the correlation results.")
  }

  combined <- ASSIGN::merge_drop(results, correlation)
  results <- combined[, seq_len(length(colnames(results)))]
  correlation <- combined[, (length(colnames(results)) + 1):(length(colnames(results)) + length(colnames(correlation)))]

  #correlate the specific columns with the pathways according to the
  #correlationList, optimize the correlation and get the gene list
  correlation_results <- list()
  optimized_path <- list()
  optimizedGeneList <- list()
  for (curr_path in names(correlationList)){
    cor_column_idx <- which(colnames(correlation) %in% correlationList[[curr_path]])
    if (length(cor_column_idx) != length(correlationList[[curr_path]])){
      stop("correlationList columns do not match columns found in correlation data. Check the correlation column names")
    }
    run_column_idx <- grep(paste("^", curr_path, "_", sep = ""), colnames(results))
    #for each column in results
    cor_mat <- matrix(0,
                      nrow = length(run_column_idx),
                      ncol = length(cor_column_idx),
                      dimnames = list(colnames(results)[run_column_idx],
                                      colnames(correlation)[cor_column_idx]))
    for (cor_column in cor_column_idx){
      for (run_column in run_column_idx){
        temp <- stats::cor.test(correlation[, cor_column], results[, run_column],
                                use = "pairwise", method = "spearman")
        cor_mat[colnames(results)[run_column], colnames(correlation)[cor_column]] <- temp$estimate
      }
    }
    correlation_results[[curr_path]] <- cor_mat

    #get the result that has the optimized path
    optimized_path[[curr_path]] <- names(which.max(rowMeans(correlation_results[[curr_path]])))

    #get the gene list for the optimized path
    optimum_length <- as.numeric(strsplit(optimized_path[[curr_path]], split = "_")[[1]][2])
    optimizedGeneList[[curr_path]] <- c(gfrn_geneList[[paste(curr_path, "up", sep = "_")]][1:floor(optimum_length / 2)],
                                        gfrn_geneList[[paste(curr_path, "down", sep = "_")]][1:ceiling(optimum_length / 2)])
  }

  #delete the non-optimal outputs if keep_optimized_only is set
  if (keep_optimized_only){
    message("Deleting the results that are not optimium")
    unlink(dir(pattern = "gene_list")[!(dir(pattern = "_gene_list") %in% as.vector(unlist(optimized_path)))], recursive = TRUE)
  }

  #return the optimized gene list and the correlation matrices
  outlist <- list()
  outlist[["optimizedGeneList"]] <- optimizedGeneList
  outlist[["correlationResults"]] <- correlation_results
  return(outlist)
}
