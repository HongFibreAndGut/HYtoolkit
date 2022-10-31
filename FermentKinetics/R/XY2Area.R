#' @title XYdf2Area
#'
#' @description This function converts coordinates (X and Y) of a spectrum to the Abs/Rel area for each molecule.
#' A dataframe (PkRange augment) showing the peak range is required as an input. Please see the examples section to generate your PkRange df.
#'
#'
#' @param df Data.frame with the coordinate numbers as two columns.
#' @param X Numeric. The horizontal axis coordicates.
#' @param Y Numeric. The vertical axis coordicates.
#' @param InRef Character. Internal reference. Default - "none". If InRef is specify, the InRef should be included in the PkRange. Also, the relative area for each molecule will be reported.
#' @param ShiftPeakTo A vecter containing the name of a molecule and a number. The spectrum will be shifted horizontally so that the horizontal coordinate of the specified molecule is the number. e.g. c("TSP", 0). Default - "none".
#' @param PkRange A data.frame which shows the peak range of each molecule.
#' @param GivePlots Logical. If T, a plot will be generated and store in the outputlist[[2]]. Default - T.
#' @param ncolplot Numeric. Number of figures in the output plot. Default - 5.
#'
#' @return This function returns a list containing a dataframe of Abs/Rel Area values, and a plot summarising the integral details
#' @export
#'
#' @examples
#' The PkRange data.frame should have 2 columns for each molecule, with form of, e.g. Formate_from, Formate_to. The first row should be the
#' horizontal axis coordinate of the _from and _to.
#'
#' Here is a sample:
#'
#' PkRange <- setNames(data.frame(matrix(ncol = 26, nrow = 1)), paste0(rep(c("TSP","Formate","Acetate",
#'                                                                       "Propionate","Butyrate","i-Butyrate",
#'                                                                       "Valerate","i-Valerate","Methanol",
#'                                                                       "Ethanol","Lactate","Succinate",
#'                                                                       "4-Methylvalerate"), each = 2), c("_from","_to")))
#' PkRange[1,] <- c(-0.01,	0.01,	8.235, 8.25,	2.07,	2.1,
#'                  1.056,	1.1,	1.6, 1.642,	2.561,	2.64,
#'                  1.29,	1.36,	2.238,	2.27,	3.333,	3.36,
#'                  1.155,	1.19,	1.402,	1.43,	2.664,	2.675,
#'                  0.86,	0.895)
#'


XYdf2Area <- function(df,
                      X,
                      Y,
                      InRef = "none",
                      ShiftPeakTo = "none", #c("TSP", 0),
                      PkRange,
                      GivePlots = T,
                      ncolplot = 5
                      ){

ReturnList <- list()

rownames(PkRange) <- PkRange[,"SampleName"]
  PkRange <- PkRange[,-1]

    GetAbsArea <- function(Metab){
    MetNm <- Metab[1]
    Fm0 <- tmpPkRange[1,which(colnames(tmpPkRange) == paste0(MetNm,"_from"))]
    To0 <- tmpPkRange[1,which(colnames(tmpPkRange) == paste0(MetNm,"_to"))]
    Fm <- min(Fm0, To0)
    To <- max(Fm0, To0)
    # range for peak
    RangePk <- tmpdf[which(tmpdf$ChemShift > Fm & tmpdf$ChemShift < To),]
    # matrix for baseline
    blmatri <- matrix(ncol = nrow(RangePk), nrow = 1, dimnames = list("1",RangePk$ChemShift))
    blmatri[1,] <- RangePk$Intensity
    # basline
    blp <- new("hyperSpec", spc = blmatri) %>% spc.rubberband(., spline = F)
    bl <- blp$spc %>% as.data.frame()

    RangeBl <- data.frame(ChemShift = RangePk$ChemShift,
                          Intensity = bl[1,] %>% as.numeric)
    # delete NA rows
    if(T %in% is.na(RangeBl$Intensity)){
      RangePk <- RangePk[-which(is.na(RangeBl$Intensity)),]
      RangeBl <- RangeBl[-which(is.na(RangeBl$Intensity)),]
    }

    integrate(approxfun(RangePk$ChemShift,
                        RangePk$Intensity - RangeBl$Intensity),
              lower = RangePk$ChemShift[1],
              upper = RangePk$ChemShift[nrow(RangePk)],
              subdivisions=2000,
              stop.on.error = F) %>% .$value %>% return()
  }

  # Setup final SampleRelArea Table
  SmpRelArea <- setNames(data.frame(matrix(ncol = ncol(PkRange)/2, nrow = 1)), c(colnames(PkRange) %>% gsub("_.*","",.) %>% unique()))

  tmpdf <- df %>% arrange(., X %>% parse(text = .) %>% eval())
  tmpPkRange <- PkRange

  # If you want to shift the whole spectrum to somewhere a peak should be 0

  if (is.vector(ShiftPeakTo)){
    thePk <- ShiftPeakTo[1]
    tmpTSP <- tmpdf[tmpdf[,X] %>% between(.,
                                          str_detect(colnames(tmpPkRange), thePk) %>% tmpPkRange[,.] %>% min(),
                                          str_detect(colnames(tmpPkRange), thePk) %>% tmpPkRange[,.] %>% max()),]
    tmpChemShCal <- tmpTSP[,X][which(tmpTSP[,Y] == max(tmpTSP[,Y]))]
    tmpdf[,X] <- tmpdf[,X] - tmpChemShCal + ShiftPeakTo[2] %>% as.numeric()
  }

  # Setup final table
  if (InRef != "none"){
    Area <- setNames(data.frame(matrix(ncol = 3, nrow = ncol(tmpPkRange)/2)),c("Metabolite","Abs_Area","Rel_Area"))
    Area$Metabolite <- colnames(tmpPkRange) %>% gsub("_.*","",.) %>% unique()
    Area$Abs_Area <- apply(Area, 1, GetAbsArea)
    Area$Rel_Area <- Area$Abs_Area / Area[Area$Metabolite == InRef,"Abs_Area"]
  } else {
    Area <- setNames(data.frame(matrix(ncol = 2, nrow = ncol(tmpPkRange)/2)),c("Metabolite","Abs_Area"))
    Area$Metabolite <- colnames(tmpPkRange) %>% gsub("_.*","",.) %>% unique()
    Area$Abs_Area <- apply(Area, 1, GetAbsArea)
  }
  ReturnList[[1]] <- Area
  names(ReturnList)[1] <- "Area"

  # To plot
  if (GivePlots == T){
    Plotlist <- sapply(Area$Metabolite,function(x) NULL)
    for (e in 1: length(Plotlist)){
      MetNm <- names(Plotlist[e])
      Fm0 <- tmpPkRange[1,which(colnames(tmpPkRange) == paste0(MetNm,"_from"))]
      To0 <- tmpPkRange[1,which(colnames(tmpPkRange) == paste0(MetNm,"_to"))]
      Fm <- min(Fm0, To0)
      To <- max(Fm0, To0)
      # range for peak
      RangePk <- tmpdf[which(tmpdf[,X] > Fm & tmpdf[,X] < To),]
      # matrix for baseline
      blmatri <- matrix(ncol = nrow(RangePk), nrow = 1, dimnames = list("1",RangePk[,X]))
      blmatri[1,] <- RangePk[,Y]
      # basline
      blp0 <- new("hyperSpec", spc = blmatri) %>% spc.rubberband(., spline = F) %>% .[[1]] %>% t() %>% as.data.frame()
      blp <- data.frame(x = as.numeric(rownames(blp0)), y = blp0$`1`)
      names(blp) <- c(X,Y)
      # plot
      Plotlist[[e]] <- ggplot() +
        geom_line(data = RangePk[complete.cases(blp),], aes_string(x = X, y = Y), color = "blue") +
        geom_line(data = blp[complete.cases(blp),], aes_string(x = X, y = Y), color = "red") +
        theme_bw() + ggtitle(MetNm)
    }
    ReturnList[[2]] <- paste0("Plotlist[[",1: length(Plotlist),"]]") %>% paste0(., collapse = "+") %>%
      paste0(., "+plot_layout(ncol = ", ncolplot,")") %>%
      parse(text = .) %>% eval()
    names(ReturnList)[2] <- "IntegralPlots"
  }
  return(ReturnList)
}


