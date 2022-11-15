#' @title RelAbun
#' @description Convert an input data.frame to relative abundance %
#' @param df a data.frame, with sample name(s) as colname(s)
#' @param AbunCol A numeric/character vector, or a number/character, of the collumn number(s)/name(s), to transfer to relative abundance
#' @return Returns a relative abundace data.frame
#' @export
#' @examples df <- data.frame(Taxon = paste0("Taxon_",1:30),
#'                 x = abs(rnorm(30)),
#'                 y = abs(rnorm(30)),
#'                 z = abs(rnorm(30)))
#'                 Reldf <- RelAbun(df, c(2:4))
#'
#'
RelAbun <- function(df,
                    AbunCol){
  df[,AbunCol] <- sweep(df[,AbunCol],2,colSums(df[,AbunCol]),`/`) *100
  return(df)
}
