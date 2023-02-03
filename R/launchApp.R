#' Launch ShinyMultiSDA.
#' @title Launch ShinySDA.
#' @description Launch ShinySDA.
#' @keywords shiny SDA
#' @export
#' @return Shiny application.
#' @import shiny
#' @import shinydashboard
#' @import shinyWidgets
#' @import shinyFiles
#' @import data.table
#' @rawNamespace import(dplyr, except = c(last, first, select, between, combine)) 
#' @import ggplot2
#' @import ggforce
#' @import ggthemes
#' @import ggthemes
#' @import ggrepel
#' @import rclipboard
#' @import grid
#' @import gridExtra
#' @import stringr
#' @import viridis
#' @import RColorBrewer
#' @import BiocParallel
#' @import clusterProfiler
#' @import AnnotationHub
#' @import org.Hs.eg.db
#' @import org.Mmu.eg.db
#' @import biomaRt
#' @import Seurat
#' 
launchShinyMultiSDA <- function(...) {
  ## runApp() does not work w shiny-server
  shinyAppDir(appDir = system.file("app", package = "ShinySDA"))
  
}
