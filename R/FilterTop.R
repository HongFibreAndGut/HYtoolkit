#' @title FilterTop
#' @description This functions filters the top abundant data with a specified threshold
#' @param df A data.frame. Variables (e.g. taxa) as rows; Samples as colnames.
#' @param TaxonCol A character. The colname of the Taxon column. Should be specified if Others = T.
#' @param AbunCol A numeric/character vector, or a number/character, of the collumn number(s)/name(s), to transfer to relative abundance
#' @param Others Logical. If T, an "Others" row will be generated.
#' @param Threshold A number. The threshold that the min values across all collumns less than it will not be included in the output data.frame.
#'
#' @return Returns a data.frame with the top abundant taxa
#' @export
#'
#' @examples library(FermentKinetics)
#' df <- data.frame(Taxon = paste0("Taxon_",1:10),
#'                 x = abs(rnorm(10,2)),
#'                 y = abs(rnorm(10,1.5)),
#'                 z = abs(rnorm(10,3)))
#' df <- RelAbun(df, AbunCol = c(2:4))
#' Top5df <- FilterTop(df, TaxonCol = "Taxon",AbunCol = c(2:4), Threshold = 5)
#'

FilterTop <- function(df,
                      TaxonCol = NULL,
                      AbunCol,
                      Others = F,
                      Threshold){
  df1 <- df[apply(df[,AbunCol],1,min)>=Threshold,]
  if (Others == T){
    OtherAbun <- 100 - apply(df1[,AbunCol],2,sum)
    df1[nrow(df1)+1,AbunCol] <- OtherAbun
    df1[nrow(df1),TaxonCol] <- "Others"
  }
  return(df1)
}



