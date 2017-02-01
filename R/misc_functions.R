#' Combine two data frames
#'
#' @param x The first data frame to be coerced to one.
#' @param y The second data frame to be coerced to one.
#' @param by specifications of the columns used for merging. The default is
#' by row names
#' @param ... arguments to be passed to or from methods.
#'
#' @return The returned data frame is the combination of x and y, with the
#' rownames properly assigned.
#' @export merge_drop
#'
#' @examples
#'
#' \dontrun{
#' merged.df <- merge_drop(df1,df2)
#' }
merge_drop<-function(x,y,by=0,...){
  new_m<-merge(x,y,by=by,...)
  rownames(new_m)<-new_m$Row.names
  return(new_m[,2:length(colnames(new_m))])
}

#' Display a PCA Plot of the Data
#'
#' @param mat The data frame on which to perform pca.
#' @param sub The number of samples in this batch, from left to right in the
#' data frame
#' @param center a logical value indicating whether the variables should be
#' shifted to be zero centered. The default is TRUE
#' @param scale a logical value indicating whether the variables should be
#' scaled to have unit variance before the analysis takes place. The default
#' is TRUE
#'
#' @return A PCA plot is displayed
#' @export pcaplot
#'
pcaplot<-function(mat,sub,center=T,scale=T, plottitle="PCA"){
  if(length(sub)!=length(mat)){
    stop("verify the subscripts...exiting now")
  }
  else{
    pca_mat <- stats::prcomp(t(mat), center=center,scale=scale)
    pca_mat_plot <- data.frame(pca_mat$x[,1:2])
    pca_mat_plot$Group <- factor(sub)
    return(ggplot2::ggplot(pca_mat_plot, ggplot2::aes(PC1, PC2)) + 
             ggplot2::geom_point(ggplot2::aes(colour = Group), size=2) +
             ggplot2::ggtitle(plottitle))
  }
}

#' Gather the ASSIGN results in a specific directory
#'
#' @return A data frame of ASSIGN predictions from each run in the directory
#' 
#' @export gather_assign_results
#'
gather_assign_results <- function(){
  curr_files <- list.files(pattern="pathway_activity_testset.csv",
                           recursive = T)
  results_df <- data.frame()
  for (i in curr_files){
    curr <- utils::read.csv(i, header=T, row.names=1)
    colnames(curr) <- strsplit(i, split="/")[[1]][1]
    if(ncol(results_df) ==0){
      results_df <- curr
    }
    else{
      results_df <- cbind(results_df, curr)
    }
  }
  rownames(results_df) <- substr(rownames(results_df),1,14)
  return(results_df)
}
