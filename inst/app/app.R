library(shiny)
library(shinyWidgets)
library(shinydashboard)
library(shinyFiles)

library(ggplot2)
library(ggrepel)
library(viridis)
library(RColorBrewer)

library(grid)
library(gridExtra) 

library(data.table)
library(dplyr)

library(Seurat)

library(SDAtools)

library(Rlabkey)
library(Rdiscvr)

library(roxygen2)

library(rclipboard)


library("BiocParallel")
register(MulticoreParam(4))

library(ShinyMultiSDA)


if (Sys.getenv("SCRATCH_DIR") != "") {
  init.path = paste0(Sys.getenv("SCRATCH_DIR"), "/data")
}  else {
  init.path = getwd()
}

# source(system.file('app/fxs.R', package = 'ShinySDA', mustWork = TRUE), local = TRUE)

print(Sys.getenv("SCRATCH_DIR"))

ui <- dashboardPage(skin="red",
                    dashboardHeader(title = "ShinyMultiSDA"),
                    #https://rstudio.github.io/shinydashboard/appearance.html#icons
                    dashboardSidebar(
                      sidebarMenu(
                        menuItem("Folder Input", tabName = "FolderInput", icon = icon("dashboard")),
                        menuItem("Gene-loading Cor HM", tabName = "GL_cor_HM", icon = icon("wrench")),
                        menuItem("Save Out", tabName = "SaveOut", icon = icon("save")),
                        menuItem("@eisamahyari", icon = icon("heart"), 
                                 href = "https://eisascience.github.io")
                      )
                    ),
                    
                    dashboardBody(
                      # useShinyjs(),
                      tags$head(
                        tags$style(HTML("
                                        .content-wrapper {
                                        background-color: black !important;
                                        }
                                        .main-sidebar {
                                        background-color: black !important;
                                        }
                                        .multicol .shiny-options-group{
                                        -webkit-column-count: 5; /* Chrome, Safari, Opera */
                                        -moz-column-count: 5;    /* Firefox */
                                        column-count: 5;
                                        -moz-column-fill: balanced;
                                        -column-fill: balanced;
                                        }
                                        .checkbox{
                                        margin-top: 0px !important;
                                        -webkit-margin-after: 0px !important; 
                                        }
                                        "))),
                      tabItems(
                        
                        # Folder input------
                        tabItem(tabName = "FolderInput",
                                h2("Select Needed SDA models"),
                                fluidRow(
                                  valueBoxOutput("InfoBox_Folder", width = 6),
                                  
                                  box(textInput("SDAroot", "Path to SDA folders. Expects dimnames in one dir up.", 
                                                value =init.path),
                                      actionButton("loadSDAfiles", "Load all files"),
                                      actionButton("CombineSDAs", "Combine SDA Loadings"),
                                      width = 10, background = "teal"
                                      
                                  ),
                                  box(uiOutput("files_ui"),
                                      width = 10, background = "teal"
                                      
                                  ),
                                  
                                  )
                                
                                ),
                        
                       
                        
                        # Save out----------
                        tabItem(tabName = "SaveOut",
                                h2("Save the results for downstream analysis"),
                                fluidRow(
                                  box(
                                    title = "Save as Seurat Object", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    actionButton("SaveAsSerObj", "Save as Seurat Obj"),
                                    width = 5, background = "black"
                                  )))
                        
                      ) #end of tabItems
                    ) #end of body
) #end UI

# Server --------

server <- function(input, output, session) {
  
  
  ## loading local files
  shinyFileChoose(input,'file', session=session,roots=c(wd='.'))
  
  
  ## environment defaults
  envv=reactiveValues(y=NULL)
  envv$InfoBox_sub = "Load in SDA, use either the Prime-seq or Folder input tabs"
  
  envv$selected_files = NULL
  envv$input_obj_ls = list()
  
  ## controls how to handle pipe depending on data
  envv$Origin = "unk"
  
  output$files_ui <- renderUI({
    req(input$SDAroot)
    path <- input$SDAroot
    files <- list.files(path=path, pattern=".rds", full.names=TRUE, recursive=F)

    files = files[!grepl("rawData", files)]
    files = files[!grepl("SDAtools", files)]
    # envv$input_SDA_files = files
    checkboxGroupInput("files", "Available .rds files", choices=basename(files))
  })
  
  source("app_OE.R",local = TRUE)
  
  
  
  
}

shinyApp(ui, server)

