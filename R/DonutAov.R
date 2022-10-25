#' @title DonutAov
#'
#' @description Making Donut charts summarising the Eta2 of two-way or three-way ANOVA based on Aov2TukeySummary OR Aov2TukeySummary output data.frame
#'
#' @param aovdf The data.frame output from Aov2TukeySummary OR Aov3TukeySummary function
#' @param Vari. A character vector contaning variables to make donut charts
#' @param single_out Logical. Whether to export every single Donut .png file. Default - F
#' @param all_out Logical. Whether to export all Donut figures in a .png file. Default - F
#' @param all_out_col Numeric. Column numbers in the all_out .png file. Default - 1
#' @param outfolder Character. Folder path to export the .png file. default - getwd. Better to setwd prior to this function
#' @param res Numeric. Resolution of the output .png file. Default - 150
#' @param height Numeric. Height of each single .png file. Default - 1000
#' @param width Numeric. Width of each single .png file. Default - 1000
#' @param VariSize Numeric. Font size of the Vari. on Donut chart. Default - 8
#' @param labelSize Numeric. Font size of the labels, i.e. Variable&Interaction&Residuals. Default - 4
#' @param DonutFill Character vector must have a length of 7 OR 10. Colors to choose for donut chart. Anti-clockwise, outer donut and then inner donut. Default - Default.
#' @return One or more .png files
#' @examples
#' Omitted
#' @export
#' @import ggplot2
#' @import reshape2
#' @import ggrepel
#' @import tidyverse
#' @import patchwork

DonutAov <- function(aovdf, Vari.,
                     single_out = F,
                     all_out = F,
                     all_out_col = 1,
                     outfolder = getwd(),
                     res = 150,
                     height = 1000,
                     width = 1000,
                     VariSize = 8,
                     labelSize = 4,
                     DonutFill = "default"){
colnames(aovdf)[1] <- "Eta2"
df0 <- aovdf[(which(aovdf$Eta2 == "Eta2")+1):nrow(aovdf),c("Eta2",Vari.)]
df1 <- melt(df0, id.vars = "Eta2")
df1$Eta2 <- gsub("Eta2_","",df1$Eta2)
df1$value <- as.numeric(df1$value) * 100
df1$Type <- ifelse(grepl(":", df1$Eta2), "Interaction", ifelse(df1$Eta2 == "Residuals", "Residuals", "Variable"))

plotlist <- list()
for (i in Vari.){
  tmpdf <- df1[df1$variable == i,]
  tmpGrp <- aggregate(value~Type, tmpdf[,c("variable","value","Type")], sum)
  tmpGrp <- arrange(tmpGrp, match(Type, c("Variable","Interaction","Residuals")))
  tmpGrp$Type <- factor(tmpGrp$Type, levels = c("Variable","Interaction","Residuals"))
  tmpdf$Eta2 <- factor(tmpdf$Eta2, levels = tmpdf$Eta2)
  # Get the positions
  dfpos <- tmpdf %>%
    mutate(csum = rev(cumsum(rev(value))),
           pos = value/2 + lead(csum, 1),
           pos = if_else(is.na(pos), value/2, pos))

  dfposGrp <- tmpGrp %>%
    mutate(csum = rev(cumsum(rev(value))),
           pos = value/2 + lead(csum, 1),
           pos = if_else(is.na(pos), value/2, pos))

  VariCir <- data.frame(x = 1,
                        y  = 0)

  if (DonutFill != "default"){
    col0 <- c(tmpdf$Eta2, tmpGrp$Type) %>% as.character() %>% unique()
    col1 <- col0 %>% sort()
    Fillcolor0 <- match(col1, col0) %>% DonutFill[.]


    plotlist[[i]] <- ggplot() +
      geom_col(aes(x, y),data = VariCir, width = 0.01) +
      geom_col(aes(x = 3, y = value, fill = Eta2), data = tmpdf, color = "black", width = 1)  +
      geom_col(aes(x = 1.5, y = value, fill = Type),alpha= 1/3, data = tmpGrp, width = 1, color = "black",) +

      geom_label_repel(data = dfpos,
                       aes(x = 3, y = pos, label = paste0(Eta2,"\n",value, "%")),
                       size = labelSize, nudge_x = 0.5, show.legend = FALSE, max.overlaps = 10000) +
      geom_text_repel(data = dfposGrp,
                      aes(x = 1.7, y = pos, label = paste0(Type,"\n",value, "%")),
                      size = labelSize, show.legend = FALSE,fontface =2, max.overlaps = 10000) +
      scale_fill_manual(values = Fillcolor0) +
      annotate("text", x = 0, y = 0, label = i, size = VariSize) +
      coord_polar("y") +
      theme_void() +
      theme(legend.position = "none",
            plot.margin = margin(0, 0.5, 0.5, 0.5, "cm"))
  } else {
    plotlist[[i]] <- ggplot() +
      geom_col(aes(x, y),data = VariCir, width = 0.01) +
      geom_col(aes(x = 3, y = value, fill = Eta2), data = tmpdf, color = "black", width = 1)  +
      geom_col(aes(x = 1.5, y = value, fill = Type),alpha= 1/3, data = tmpGrp, width = 1, color = "black",) +

      geom_label_repel(data = dfpos,
                       aes(x = 3, y = pos, label = paste0(Eta2,"\n",value, "%")),
                       size = labelSize, nudge_x = 0.5, show.legend = FALSE, max.overlaps = 10000) +
      geom_text_repel(data = dfposGrp,
                      aes(x = 1.7, y = pos, label = paste0(Type,"\n",value, "%")),
                      size = labelSize, show.legend = FALSE,fontface =2, max.overlaps = 10000) +
      annotate("text", x = 0, y = 0, label = i, size = VariSize) +
      coord_polar("y") +
      theme_void() +
      theme(legend.position = "none",
            plot.margin = margin(0, 0.5, 0.5, 0.5, "cm"))
  }
  if (single_out == T){
  png(paste0(outfolder,"/",i,".png"),
      width = width,
      height = height,
      res = res)
  print(plotlist[[i]])
  print()
  graphics.off()
}
}

if (all_out == T){
  png(paste0(outfolder,"/all_out.png"),
      width = all_out_col * width,
      height = ceiling(length(plotlist)/all_out_col) * height,
      res = res)
  (paste0("plotlist[[",1: length(plotlist),"]]") %>%
    paste0(.,collapse = "+") %>% paste0(.,"+plot_layout(ncol = ", all_out_col,")")) %>%
    parse(text = .) %>% eval() %>% print()
  graphics.off()
}
}
