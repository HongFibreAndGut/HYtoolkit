#' @title PCAbiplot
#' @description This function returns a list for information of PCA biplot, and a ggplot object (list[[4]] in the list) which can be export directly. Otherwise, you can generate your own
#' ggplot object just from the first three elements. list[[1]] is a data.frame showing the point for PCA plot and the grouping information.
#' list[[2]] is a vector of the percentages of each PC. list[[3]] is a data.frame for loadings (rotation). list[[4]] is a data.frame giving the information from which you can generate the arrows by: geom_segment(data=list[[3]], aes(x=0, y=0, xend=v1, yend=v2))
#'
#' @param df A data.frame. Can include the grouping info column(s).
#' @param Vari. Character vector. Which variables (columns) are used.
#' @param scale Logical. Default - T.
#' @param center Logical. Default - T.
#' @param Grps Character vector. Which groups (columns) are used. Default - NULL.
#' @param GrpCol Character. Which group (column) is used to define the color. Default - NULL.
#' @param GrpColLevels Character vector, setting the levels of the GrpCol. Default - NULL.
#' @param GrpShape Character. Which group (column) is used to define the shape. Default - NULL.
#' @param GrpShapeLevels Character vector, setting the levels of the GrpShape. Default - NULL.
#' @param GrpSize Character. Which group (column) is used to define the size. Default - NULL.
#' @param GrpSizeLevels Character vector, setting the levels of the GrpSize. Default - NULL.
#' @param NoGrpCol Character. Color to use if not define the GrpCol. Default - "black".
#' @param NoGrpShape Character. Shape to use if not define the GrpCol. Default - "circle".
#' @param NoGrpSize Numeric. Size to use if not define the GrpCol. Default - 3.
#' @param PointAlpha Numeric. Transparency of points. Default - 0.5.
#' @param PointText  Character. Which column is used to show the text on plot. Default - NULL.
#' @param Colmannual Character vector. Define the color. Default - NULL.
#' @param Shapemannual Character vector. Define the shape. Default - NULL.
#' @param Sizemannual Character vector. Define the size. Default - NULL.
#' @param DrawElips Logical. Draw elipse or not.Default - T.
#' @param GrpElips Character. Which group (column) is used to Draw elipse. Default - NULL.
#' @param ElipsCol Character. Which group (column) is used to define elipse color. Default - NULL.
#' @param ShowArrow Logical. Draw arrows or not. Default - T
#' @param ArrowSize Numeric. Arrow size. Default - 1.
#' @param ArrowLength Numeric. Arrow length. Default - 0.3.
#' @param SegLengthMultiply Numeric. Length factor to multiply to the segment. Default - NULL.
#' @param ArrowAlpha Numeric. Arrow transparency. Default = 1.
#' @param ArrowCol Character. Arrow color. Default - "black".
#' @param ArrowTextSize Numeric. Text size of arrow. Default - 3.
#' @param axis.title Numeric. Text size of axis title. Default - 15.
#' @param axis.text.x Numeric. Text size of x axis. Default - 15.
#' @param axis.text.y Numeric. Text size of y axis. Default - 15.
#' @param legend.key.size Numeric. Legend key size. Default - 1.
#' @param legend.title Numeric. Legend title size. Default - 16.
#' @param legend.text Numeric. Legend text size Default - 16.
#'
#' @return Returns a list. See description
#' @export
#' @import ggplot2
#' @import ggrepel
#' @importFrom dplyr %>%
#' @import tidyr
#'
#' @examples Download a demo data from my github: hyfood/DemoData/PCAbiplotDemo.txt, and read into R.
#' PCAbiplotDemo$Dnr <- as.character(PCAbiplotDemo$Dnr)
#' biplotList <- PCAbiplot(df = PCAbiplotDemo,
#' Vari. = c("A", "B", "C", "cCV", "TRMax", "RMax", "LagT"),
#' Grps = c("Substrate", "Dnr", "PWgt"),
#' GrpCol = "Dnr",
#' GrpShape = "Substrate",
#' GrpSize = "PWgt",
#' DrawElips = T,
#' GrpElips = "Dnr",
#' ElipsCol = "Dnr")

PCAbiplot <- function(# PCA info
  df,
  Vari.,
  scale = T,
  center = T,
  Grps = NULL,

  # Point - color, shape, size info
  GrpCol = NULL,
  GrpColLevels = NULL,
  GrpShape = NULL,
  GrpShapeLevels = NULL,
  GrpSize = NULL,
  GrpSizeLevels = NULL,
  NoGrpCol = "black",
  NoGrpShape = "circle",
  NoGrpSize = 3,
  PointAlpha = 0.5,
  PointText = NULL,
  Colmannual = NULL,
  Shapemannual = NULL,
  Sizemannual = NULL,

  # Elipse info
  DrawElips = F,
  GrpElips = NULL,
  ElipsCol = NULL,

  # biplot - info
  ShowArrow = T,
  ArrowSize = 1,
  ArrowLength = 0.3,
  SegLengthMultiply = NULL,
  ArrowAlpha = 1,
  ArrowCol = "black",
  ArrowTextSize = 3,

  # theme setting
  axis.title = 15,
  axis.text.x = 15,
  axis.text.y = 15,
  legend.key.size = 1,
  legend.title = 16,
  legend.text = 16
){
  fit <- prcomp(df[,Vari.][complete.cases(df[,Vari.]),],center=TRUE, scale.=TRUE)
  df2 <- df[complete.cases(df[,Vari.]),]

  # do biplot
  data <- data.frame(fit$x)
  datapc <- data.frame(varnames=rownames(fit$rotation), fit$rotation)
  datapc0 <- datapc
  mult <- ifelse(is.null(SegLengthMultiply),  min(
    (max(data[,"PC2"]) - min(data[,"PC2"])/(max(datapc[,"PC2"])-min(datapc[,"PC2"]))),
    (max(data[,"PC1"]) - min(data[,"PC1"])/(max(datapc[,"PC1"])-min(datapc[,"PC1"])))),
    SegLengthMultiply
  )
  datapc <- transform(datapc,
                      v1 = mult * (get("PC1")),
                      v2 = mult * (get("PC2"))
  )

  #data$obsnames <- to_PCA$DnrSubWgt
  data <- cbind(data, df2[,Grps])
  #data$Substrate <- factor(data$Substrate, levels = c("Pect","ACW","APart"))
  eigs <- fit$sdev^2
  PCvector <- c()
  for (i in 1: length(eigs)){
    PCvector <- append(PCvector, round(eigs[i]/sum(eigs) * 100, 2))
  }

  if (!is.null(GrpColLevels)){ data[, GrpCol] <- factor(data[, GrpCol], levels = GrpColLevels)}
  if (!is.null(GrpShapeLevels)){ data[, GrpShape] <- factor(data[, GrpShape], levels = GrpColLevels)}
  if (!is.null(GrpSizeLevels)){ data[, GrpSize] <- factor(data[, GrpSize], levels = GrpColLevels)}

  # geom_point info
  if (is.null(GrpShape) & is.null(GrpCol) & is.null(GrpSize)){ pointList <- list(geom_point(shape = NoGrpShape,color = NoGrpCol, size = NoGrpSize, alpha = PointAlpha))}
  if (is.null(GrpShape) & is.null(GrpCol) & !is.null(GrpSize)){ pointList <- list(geom_point(aes_string(size = GrpSize), shape = NoGrpShape,color = NoGrpCol, alpha = PointAlpha))}
  if (is.null(GrpShape) & !is.null(GrpCol) & is.null(GrpSize)){ pointList <- list(geom_point(aes_string(color = GrpCol), shape = NoGrpShape, size = NoGrpSize, alpha = PointAlpha))}
  if (is.null(GrpShape) & !is.null(GrpCol) & !is.null(GrpSize)){ pointList <- list(geom_point(aes_string(color = GrpCol, size = GrpSize), shape = NoGrpShape, alpha = PointAlpha))}
  if (!is.null(GrpShape) & is.null(GrpCol) & is.null(GrpSize)){ pointList <- list(geom_point(aes_string(shape = GrpShape)), color = NoGrpCol, size = NoGrpSize, alpha = PointAlpha)}
  if (!is.null(GrpShape) & is.null(GrpCol) & !is.null(GrpSize)){ pointList <- list(geom_point(aes_string(shape = GrpShape, size = GrpSize), color = NoGrpCol, alpha = PointAlpha))}
  if (!is.null(GrpShape) & !is.null(GrpCol) & !is.null(GrpSize)){ pointList <- list(geom_point(aes_string(shape = GrpShape, color = GrpCol, size = NoGrpSize, alpha = PointAlpha)))}
  if (!is.null(GrpShape) & !is.null(GrpCol) & is.null(GrpSize)){ pointList <- list(geom_point(aes_string(shape = GrpShape, color = GrpCol), size = NoGrpSize, alpha = PointAlpha))}

  if (!is.null(Colmannual)){ ColmannualList <- list(scale_color_manual(values = Colmannual))}
  if (is.null(Colmannual)){ ColmannualList <- list()}
  if (!is.null(Shapemannual)){ ShapemannualList <- list(scale_shape_manual(values = Shapemannual))}
  if (is.null(Shapemannual)){ ShapemannualList <- list()}
  if (!is.null(Sizemannual)){ SizemannualList <- list(scale_size_manual(values = Sizemannual))}
  if (is.null(Sizemannual)){ SizemannualList <- list()}

  # point text info
  if (!is.null(PointText)){ pointTextList <- list(geom_text(aes_string(x = "PC1", y = "PC2", label = PointText)))}
  if (is.null(PointText)){ pointTextList <- list()}
  PointText <- "Dnr"
  # geom_segment info
  if (ShowArrow == T){ArrowList <- list(geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), size = ArrowSize, arrow=arrow(length=unit(ArrowLength,"cm")), alpha=ArrowAlpha, color=ArrowCol),
                                        geom_text_repel(data=datapc, aes(x=v1, y=v2, label=varnames),size = ArrowTextSize, color=ArrowCol))}
  if (ShowArrow == F){ArrowList <- list()}

  # Elips info
  if (DrawElips == T & !is.null(ElipsCol)){ElipsList <- stat_ellipse(aes_string(group = GrpElips, col = ElipsCol))}
  if (DrawElips == T & is.null(ElipsCol)){ElipsList <- stat_ellipse(aes_string(group = GrpElips), col = "black")}
  if (DrawElips == F){ElipsList <- list()}

  # labs info
  labList <- list(labs(# Add Explained Variance per axis
    x = paste0("PCA 1 (", PCvector[1], "%)"),
    y = paste0("PCA 2 (", PCvector[2], "%)")))

  # theme
  themeList <- list(theme(axis.title = element_text(size = axis.title),
                          axis.text.x = element_text(size = axis.text.x),
                          axis.text.y = element_text(size = axis.text.y ),
                          legend.key.size = unit(legend.key.size, 'cm'),
                          legend.title = element_text(size = legend.title),
                          legend.text = element_text(size = legend.text)) +
                      theme(panel.border = element_rect(color="black", fill='transparent', size = 2)) +
                      theme(panel.background = element_rect(color = 'black', fill = 'transparent')))
  p <- ggplot(data = data, aes(x = PC1, y = PC2)) +
    pointList +
    ArrowList +
    ElipsList +
    pointTextList +
    ColmannualList + ShapemannualList + SizemannualList +
    labList +
    themeList

  ReturnList <- list()
  ReturnList[[1]] <- data
  ReturnList[[2]] <- PCvector
  ReturnList[[3]] <- datapc0
  ReturnList[[4]] <- datapc
  ReturnList[[4]] <- p

  names(ReturnList) <- c("PCandGrpInfo","PCpercent","Loadings","Segment4biplot","ggplotOut")

  return(ReturnList)
}



