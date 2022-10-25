#' @title Bruker2conc
#'
#' @description Convert Bruker output files to final concentration. All spectra should be at least displayed once in Topspin.
#' An EndproductInfo.xlsx file is requried to include the chemical shift range and the calculation information. Please download a demo from
#' my github repository hyfood/DeomoData/EndproductInfo.xlsx. The DemoData Directory also includes a Bruker output directory as a demo for this function.
#'
#' @param BrukerDir A Character. Path to the folder where all the Bruker out directories are stored
#' @param InRef A character. Internal standard. Default - "TSP"
#' @param InRef2 A character OR False. A second internal standard. Default = "4-Methylvalerate"
#' @param CorInRef2 Logical. If the concentration should be corrected by the second internal standard. Default - T.
#' @param ConcRef2 Numeric. Concentration of the second internal standard in nM. Only works if CorInRef2 = T. Default - 4.
#' @param CorChemshift2Zero Logical. If the chemical shift should be corrected. If T, all chemical shift will be shifting so that for InRef it is 0. Default - T.
#' @param Samp2run "All" or a character vector containing all samples names to process. Default - "All"
#' @param ncolplot Numeric. number of column in the intergral detail plots. Default - 5.
#' @param EndProductInfo Character. Path to the EndProductInfo.xlsx file. This file should include two sheets, namely ChemRange and Calculation. Default - "EndProductInfo.xlsx"
#' @param OutDir Character. Path to export resulting files. Default - getwd
#'
#' @return Returns a Directory containing the intergral detail plots, and two excel files showing the relative area and the concentration in mM.
#' @export
#'
#' @examples Omnitted
Bruker2conc <- function(BrukerDir,
                        InRef = "TSP",
                        InRef2 = "4-Methylvalerate",
                        CorInRef2 = T,
                        ConcRef2 = 4,
                        CorChemshift2Zero = T,
                        Samp2run = "All",
                        ncolplot = 5,
                        EndProductInfo = "EndProductInfo.xlsx",
                        OutDir = getwd()
                        ){
  dir.create(file.path(paste0(OutDir,"/","Bruker2conc_Output")), showWarnings = FALSE)
  dir.create(file.path(paste0(OutDir,"/Bruker2conc_Output/","IntegralDetail")), showWarnings = FALSE)
  ChemRange = read_xlsx(EndProductInfo, sheet = "ChemRange") %>% as.data.frame()
  Calculation = read_xlsx(EndProductInfo, sheet = "Calculation") %>% as.data.frame()
# Setup function(s)
GetAbsArea <- function(Metab){
  MetNm <- Metab[1]
  Fm0 <- tmpChemRange[1,which(colnames(tmpChemRange) == paste0(MetNm,"_from"))]
  To0 <- tmpChemRange[1,which(colnames(tmpChemRange) == paste0(MetNm,"_to"))]
  Fm <- min(Fm0, To0)
  To <- max(Fm0, To0)
  # range for peak
  RangePk <- Rawdata0[which(Rawdata0$ChemShift > Fm & Rawdata0$ChemShift < To),]
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

folderlist <- BrukerDir %>% list.dirs() %>% .[which(grepl("pdata/1",.))] %>% as.list()
SampleNm <- BrukerDir %>% list.dirs() %>% .[which(grepl("pdata/1",.))] %>% gsub("/pdata/1","",.) %>% gsub(".*/","",.)

if (!"All" %in% Samp2run) {
  folderlist <- BrukerDir %>% list.dirs() %>% .[which(grepl("pdata/1",.))] %>% .[which(str_detect(., Samp2run))] %>% as.list() %>% suppressWarnings()
  SampleNm <- Samp2run
}


ReadNMRxy <- function(datapath){
  xy <- readBruker(datapath, dimension = "1D")[[1]]
  return(data.frame(ChemShift = names(xy) %>% as.numeric(),
                    Intensity = xy,
                    row.names = NULL))
}

datalist <- lapply(folderlist, ReadNMRxy)
names(datalist) <- SampleNm

# Reformat ChemRange
rownames(ChemRange) <- ChemRange$SampleName
ChemRange <- ChemRange[,-1]

# Setup final SampleRelArea Table
SmpRelArea <- setNames(data.frame(matrix(ncol = ncol(ChemRange)/2+1, nrow = length(SampleNm))), c("Sample",colnames(ChemRange) %>% gsub("_.*","",.) %>% unique()))
SmpRelArea$Sample <- SampleNm
Allplot <- list()
for (f in 1: length(datalist)){
  Rawdata0 <- datalist[[f]] %>% arrange(., ChemShift)
  tmpSmp <- names(datalist[f])


  # Reformat Chemshift
  tmpChemRange <- ChemRange
  if (tmpSmp %in% rownames(tmpChemRange)){
    ChangeShift <- which(!is.na(tmpChemRange[which(rownames(tmpChemRange) == tmpSmp),]))
    tmpChemRange[1,ChangeShift] <- tmpChemRange[which(rownames(tmpChemRange) == tmpSmp),ChangeShift]
  }
  tmpChemRange <- tmpChemRange[1,]

  # Calibrate all chemical shifts
  if (CorChemshift2Zero == T){
      tmpTSP <- Rawdata0[Rawdata0$ChemShift %>% between(.,
                                                    str_detect(colnames(tmpChemRange), InRef) %>% tmpChemRange[,.] %>% min() %>% -0.02,
                                                    str_detect(colnames(tmpChemRange), InRef) %>% tmpChemRange[,.] %>% max() %>% +0.02),]
  tmpChemShCal <- tmpTSP$ChemShift[which(tmpTSP$Intensity == max(tmpTSP$Intensity))]
  Rawdata0$ChemShift <- Rawdata0$ChemShift - tmpChemShCal
  }

  # Setup final table
  RelArea <- setNames(data.frame(matrix(ncol = 3, nrow = ncol(tmpChemRange)/2)),c("Metabolite","Abs_Area","Rel_Area"))
  RelArea$Metabolite <- colnames(tmpChemRange) %>% gsub("_.*","",.) %>% unique()


  RelArea$Abs_Area <- apply(RelArea, 1, GetAbsArea)
  RelArea$Rel_Area <- RelArea$Abs_Area / RelArea[RelArea$Metabolite == InRef,"Abs_Area"]
  SmpRelArea[SmpRelArea$Sample == tmpSmp,2:ncol(SmpRelArea)] <- RelArea$Rel_Area

  # To plot
  Plotlist <- sapply(RelArea$Metabolite,function(x) NULL)
  for (e in 1: length(Plotlist)){
    MetNm <- names(Plotlist[e])
    Fm0 <- tmpChemRange[1,which(colnames(tmpChemRange) == paste0(MetNm,"_from"))]
    To0 <- tmpChemRange[1,which(colnames(tmpChemRange) == paste0(MetNm,"_to"))]
    Fm <- min(Fm0, To0)
    To <- max(Fm0, To0)
    # range for peak
    RangePk <- Rawdata0[which(Rawdata0$ChemShift > Fm & Rawdata0$ChemShift < To),]
    # matrix for baseline
    blmatri <- matrix(ncol = nrow(RangePk), nrow = 1, dimnames = list("1",RangePk$ChemShift))
    blmatri[1,] <- RangePk$Intensity
    # basline
    blp0 <- new("hyperSpec", spc = blmatri) %>% spc.rubberband(., spline = F) %>% .[[1]] %>% t() %>% as.data.frame()
    blp <- data.frame(ChemShift = as.numeric(rownames(blp0)), Intensity = blp0$`1`)
    # plot
    Plotlist[[e]] <- ggplot() +
      geom_line(data = RangePk[complete.cases(blp),], aes(x = ChemShift, y = Intensity), color = "blue") +
      geom_line(data = blp[complete.cases(blp),], aes(x = ChemShift, y = Intensity), color = "red") +
      theme_bw() + ggtitle(MetNm)
  }

  Allplot[[f]] <- paste0("Plotlist[[",1: length(Plotlist),"]]") %>% paste0(., collapse = "+") %>%
      paste0(., "+plot_layout(ncol = ", ncolplot,")") %>%
      parse(text = .) %>% eval()
}

# Setup final Concentratioin  Table
ConcTable <- SmpRelArea[,-which(colnames(SmpRelArea) == InRef)]

for (i in 2: ncol(ConcTable)){
 MetaNm <- colnames(ConcTable)[i]
 ConcTable[,MetaNm] <- SmpRelArea[,MetaNm] *
   Calculation[Calculation$Metabolite == MetaNm,"ProtonInRef"] *
   Calculation[Calculation$Metabolite == MetaNm,"MwInRef"] /
   Calculation[Calculation$Metabolite == MetaNm,"ProtonSmp"] *
   Calculation[Calculation$Metabolite == MetaNm,"DilutionFactor"] *
   Calculation[Calculation$Metabolite == MetaNm,"PkFactor"]
 if (T %in% (ConcTable$Sample %in% colnames(Calculation))){
   for (j in ConcTable$Sample[ConcTable$Sample %in% colnames(Calculation)]){
     ConcTable[ConcTable$Sample == j,-1] <- ConcTable[ConcTable$Sample == j,-1] *  Calculation[,j]
   }
 }
}

ConcTable1 <- ConcTable
if (CorInRef2 == T){
  for(j in 2: ncol(ConcTable1)){
    ConcTable1[,j] <- ConcTable[,j] * ConcRef2 /ConcTable[,InRef2]
  }
}

for (i in 1: length(Allplot)){
  png(paste0(OutDir,"/Bruker2conc_Output/","IntegralDetail","/",SampleNm[i],".png"),
    height = 240 * ceiling(length(Plotlist) / ncolplot),
    width = 240 * ncolplot,
    res = 80)
print(Allplot[[i]])
graphics.off()
}

write_xlsx(SmpRelArea, ifelse("All" %in% Samp2run, paste0(OutDir,"/Bruker2conc_Output/","AllSample_RelArea.xlsx"), paste0(OutDir,"/","SelectedSample_RelArea.xlsx")))
write_xlsx(ConcTable1,ifelse("All" %in% Samp2run, paste0(OutDir,"/Bruker2conc_Output/","AllConcentration Table.xlsx"), paste0(OutDir,"/","SelectedConcentration Table.xlsx")))
}
