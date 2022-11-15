#' @title doBiplot
#' @description This function returns a ggplot object for biplot.
#'
#' @param toBiplotList The output list of PCAtoBiplot or PCoAtoBiplot.
#' @param GrpInfo Grouping Info data.frame if the color, size, shape and elipse need to be defined. Same rownames(sample names) as the toBiplotList[[1]]. Default - NULL.
#' @param AxisVector Numeric vector with length of 2, wo define which two axis to do the biplot. Default - c(1,2).
#' @param Vari2Arrow Character vector to determine which variable(s) to show on the biplot as arrow(s).
#' @param AxisTitle A Character. Name of Axis title. e.g. "PCA" or "PCoA" or "Axis".
#' @param Decimal Numeric. Decimal points of axis title percentage. Default - 2.
#' @param GrpCol Character. Which group (column) is used to define the color. Default - NULL.
#' @param GrpShape Character. Which group (column) is used to define the shape. Default - NULL.
#' @param GrpSize Character. Which group (column) is used to define the size. Default - NULL.
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
#' @return This function returns a ggplot object for biplot.
#' @export
#' @import ggplot2
#' @import ggrepel
#' @importFrom dplyr %>%
#' @import tidyr

#' @examples Coming...
doBiplot <- function(# PCA info
  toBiplotList,
  GrpInfo = NULL, # Grouping Info data.frame. Same rownames(sample names) as the toBiplotList[[1]]. Default - Null.
  AxisVector = c(1,2),
  Vari2Arrow = NULL, #Character vector to determine which variable(s) to show on the biplot as arrow(s).
  AxisTitle = NA, # A Character. Name of Axis title. e.g. "PCA" or "PCoA" or "Axis".
  Decimal = 2, # Numeric. Decimal points. Default - 2.

  # Point - color, shape, size info
  GrpCol = NULL, # The colname of GrpInfo data.frame to determine the color of points
  #GrpColLevels = NULL,
  GrpShape = NULL,
  #GrpShapeLevels = NULL,
  GrpSize = NULL,
  #GrpSizeLevels = NULL,
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
  data <- toBiplotList[[1]]
  AxisPercent <- toBiplotList[[2]]
  toSeg <- toBiplotList[[3]]

  if (!is.null(Vari2Arrow)){
    toSeg <- toSeg[Vari2Arrow,]
  }

  mult <- ifelse(is.null(SegLengthMultiply),  min(
    (max(data[,colnames(toSeg)[3]]) - min(data[,colnames(toSeg)[3]])/(max(toSeg[,colnames(toSeg)[3]])-min(toSeg[,colnames(toSeg)[3]]))),
    (max(data[,colnames(toSeg)[2]]) - min(data[,colnames(toSeg)[2]])/(max(toSeg[,colnames(toSeg)[2]])-min(toSeg[,colnames(toSeg)[2]])))),
    SegLengthMultiply
  )

  toSeg <- transform(toSeg,
                     v1 = mult * (get(colnames(toSeg)[2])),
                     v2 = mult * (get(colnames(toSeg)[3]))
  )

  if (!is.null(GrpInfo)){
    data <- cbind(data, GrpInfo)
  }

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
  if (!is.null(PointText)){ pointTextList <- list(geom_text(aes_string(x = colnames(data)[AxisVector[1]], y = colnames(data)[AxisVector[2]], label = PointText)))}
  if (is.null(PointText)){ pointTextList <- list()}

  # geom_segment info
  if (ShowArrow == T){ArrowList <- list(geom_segment(data=toSeg, aes(x=0, y=0, xend=v1, yend=v2), size = ArrowSize, arrow=arrow(length=unit(ArrowLength,"cm")), alpha=ArrowAlpha, color=ArrowCol),
                                        geom_text_repel(data=toSeg, aes(x=v1, y=v2, label=varnames),size = ArrowTextSize, color=ArrowCol))}
  if (ShowArrow == F){ArrowList <- list()}

  # Elips info
  if (DrawElips == T & !is.null(ElipsCol)){ElipsList <- stat_ellipse(aes_string(group = GrpElips, col = ElipsCol))}
  if (DrawElips == T & is.null(ElipsCol)){ElipsList <- stat_ellipse(aes_string(group = GrpElips), col = "black")}
  if (DrawElips == F){ElipsList <- list()}

  # labs info
  labList <- list(labs(# Add Explained Variance per axis
    x = paste0(AxisTitle," ",AxisVector[1], "(", sprintf(paste0("%0.", Decimal,"f"),AxisPercent[1]), "%)"),
    y = paste0(AxisTitle," ",AxisVector[2], "(", sprintf(paste0("%0.", Decimal,"f"),AxisPercent[2]), "%)")))

  # theme
  themeList <- list(theme(axis.title = element_text(size = axis.title),
                          axis.text.x = element_text(size = axis.text.x),
                          axis.text.y = element_text(size = axis.text.y ),
                          legend.key.size = unit(legend.key.size, 'cm'),
                          legend.title = element_text(size = legend.title),
                          legend.text = element_text(size = legend.text)) +
                      theme(panel.border = element_rect(color="black", fill='transparent', size = 2)) +
                      theme(panel.background = element_rect(color = 'black', fill = 'transparent')))
  p <- ggplot(data = data, aes_string(x = colnames(data)[AxisVector[1]], y = colnames(data)[AxisVector[2]])) +
    pointList +
    ArrowList +
    ElipsList +
    pointTextList +
    ColmannualList + ShapemannualList + SizemannualList +
    labList +
    themeList

  return(p)
}

