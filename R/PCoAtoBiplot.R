#' @title PCoAtoBiplot
#' @description This function returns a list containing elements required for PCoA biplot. list[[1]] is the PCoA axis data. list[[2]] is a vector including the PCoA axis percentage.
#' list [[3]] is the standardised covariance matrix required for the biplot arrows
#'
#' @param NormDf a normalised data.frame. Sample names as rownames and variables as colnames.
#' @param method PCoA method. This is based on vegdist function in vegan package. Please see help(vegdist). Default - "bray".
#'
#' @importFrom dplyr %>%
#' @importFrom ape pcoa
#' @importFrom vegan vegdist
#'
#' @return Returns a list. See description
#' @export
#'
#' @examples None
PCoAtoBiplot <- function(NormDf, # rowname - sample name, colnames - variable name
                       method = "bray"
){
  # generate pcoa matrix
  pcoa <- pcoa(vegdist(NormDf, method = method)) # you can choose other methods
  pcoa_data <- as.data.frame(pcoa$vectors)

  compute_arrows <-  function(given_pcoa, trait_df) {

    # Keep only quantitative or ordinal variables
    # /!\ Change this line for different dataset
    #     or select only quantitative/ordinal var. /!\

    #trait_df = trait_df[, c(4:6, 20, 21)]

    n <- nrow(trait_df)
    points.stand <- scale(given_pcoa$vectors)

    # Compute covariance of variables with all axes
    S <- cov(trait_df, points.stand)

    # Select only positive eigenvalues
    pos_eigen = given_pcoa$values$Eigenvalues[seq(ncol(S))]

    # Standardize value of covariance (see Legendre & Legendre 1998)
    U <- S %*% diag((pos_eigen/(n - 1))^(-0.5))
    colnames(U) <- colnames(given_pcoa$vectors)

    # Add values of covariances inside object
    given_pcoa$U <- U

    return(given_pcoa)
  }

  pcoa_arrows <- compute_arrows(pcoa, as.data.frame(NormDf))

  # if (ShowVari. != "all"){
  #   pcoa_arrows$U <- pcoa_arrows$U[ShowVari.,]
  # }
  #
  # mult <- ifelse(is.null(SegLengthMultiply),  min(
  #   (max(pcoa_data[,AxisVector[2]]) - min(pcoa_data[,AxisVector[2]])/(max(pcoa_arrows$U[,AxisVector[2]])-min(pcoa_arrows$U[,AxisVector[2]]))),
  #   (max(pcoa_data[,AxisVector[1]]) - min(pcoa_data[,AxisVector[1]])/(max(pcoa_arrows$U[,AxisVector[1]])-min(pcoa_arrows$U[,AxisVector[1]])))),
  #   SegLengthMultiply
  # )
  #
  # arrows_df <- transform(pcoa_arrows$U,
  #                        v1 = mult * (get(colnames(pcoa_arrows$U)[AxisVector[1]])),
  #                        v2 = mult * (get(colnames(pcoa_arrows$U)[AxisVector[2]]))
  # )

  #arrows_df <- arrows_df[,AxisVector]


  ReturnList <- list()
  ReturnList[[1]] <- pcoa_data
  ReturnList[[2]] <- pcoa$values$Relative_eig * 100
  ReturnList[[3]] <- pcoa_arrows$U

  names(ReturnList) <- c("PCoA_data","PCoA_percent","Standardised_cov")
  return(ReturnList)
}
