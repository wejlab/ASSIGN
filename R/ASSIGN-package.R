#' Pathway signature gene sets
#'
#' Signature genes for 5 oncogenic pathways.
#'
#' @name geneList1
#' @docType data
#' @format List with 5 components representing each pathway. 200 signature
#' genes are selected for each pathway.
#' @source Bild et al. (2006) Oncogenic pathway signatures in human cancers as
#' a guide to targeted therapies. Nature, 439, 353-357.
#' @keywords datasets
NULL

#' Gene expression profiling from cancer patients (test dataset)
#'
#' Gene expression datasets for 111 lung cancer patient samples, including 53
#' cases of lung adenocarcinoma and 58 cases of lung squamous carcinoma.
#'
#' @name testData1
#' @docType data
#' @format Data frame with 1000 genes/probes (rows) and 111 samples (columns)
#' @source Bild et al. (2006) Oncogenic pathway signatures in human cancers as
#' a guide to targeted therapies. Nature, 439, 353-357.
#' @keywords datasets
NULL

#' Gene expression profiling from cell line perturbation experiments (training
#' dataset)
#'
#' Gene expression datasets for 5 oncogenic pathway perturbation experiments,
#' including B-Catenin, E2F3, MYC, RAS, and SRC pathways.
#'
#' @name trainingData1
#' @docType data
#' @format Data frame with 1000 genes/probes (rows) and 55 samples (columns)
#' @source Bild et al. (2006) Oncogenic pathway signatures in human cancers as
#' a guide to targeted therapies. Nature, 439, 353-357.
#' @keywords datasets
NULL

#' Pathway Signature Gene Lists
#'
#' Pathway signature gene lists have been optimized based on correlations of
#' pathway activity data and protein data. The gene lists can be used to avoid
#' the bayesian gene selection step of ASSIGN, which will decrease the amount of
#' time it takes to run ASSIGN.
#'
#' @name gfrn_geneList
#' @docType data
#' @format List of gene lists for akt, bad, egfr, her2, igf1r, krasgv, krasqh,
#' and raf
#' @source Bild et al.
#' @keywords datasets
NULL

#' Exclude Gene List
#'
#' Overexpression signatures may contain genes that are consistently
#' differentially expressed. This list was compiled based on the GFRN gene list.
#' These genes appear in at least 60% of all signatures at various lengths.
#'
#' @name excludegenes
#' @docType data
#' @format character vector of commonly differentially expressed genes
#' @source Bild et al.
#' @keywords datasets
NULL
