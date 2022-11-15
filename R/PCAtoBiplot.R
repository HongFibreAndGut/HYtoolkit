#' @title PCAtoBiplot
#' @description This function returns a list containing elements required for PCA biplot. list[[1]] is the PCA axis data. list[[2]] is a vector including the PC variation explained percentages.
#' list [[3]] is the loadings
#'
#' @param df Numeric data.frame. Sample names as rownames and variables as colnames.
#' @param scale Logical. Default - T.
#' @param center Logical. Default - T.
#'
#' @return Returns a list. See description
#' @export
#'
#' @examples None.
PCAtoBiplot <- function(df,scale = T,center = T){

  fit <- prcomp(df[complete.cases(df),],center=TRUE, scale.=TRUE)
  data <- data.frame(fit$x)
  datapc <- data.frame(varnames=rownames(fit$rotation), fit$rotation)
  eigs <- fit$sdev^2
  PCvector <- c()
  for (i in 1: length(eigs)){
    PCvector <- append(PCvector, round(eigs[i]/sum(eigs) * 100, 3))
  }
  ReturnList <- list()
  ReturnList[[1]] <- data
  ReturnList[[2]] <- PCvector
  ReturnList[[3]] <- datapc
  names(ReturnList) <- c("PCA_data","PC_percent","Loadings")
  return(ReturnList)
}



