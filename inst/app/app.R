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
                        menuItem("Folder Input", tabName = "FolderInput_sda", icon = icon("dashboard")),
                        menuItem("Gene Loading HM", tabName = "GL_HM", icon = icon("wrench")),
                        menuItem("Seurat Obj Input", tabName = "FolderInput_ser", icon = icon("dashboard")),
                        menuItem("ProjX SDA on Seurat", tabName = "SDA2Ser_projx", icon = icon("wrench")),
                        
                        # menuItem("Save Out", tabName = "SaveOut", icon = icon("save")),
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
                        tabItem(tabName = "FolderInput_sda",
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
                        
                        
                        # Folder input------
                        tabItem(tabName = "FolderInput_ser",
                                h2("Load Seurat Obj rds file"),
                                fluidRow(
                                  valueBoxOutput("InfoBox_Folder_ser", width = 6),
                                  
                                  box(textInput("SerObjroot", "Path to Seurat Obj.", 
                                                value =init.path),
                                      actionButton("loadSerObj", "Load SerObj"),
                                      # actionButton("CombineSDAs", "Combine SDA Loadings"),
                                      width = 10, background = "teal"
                                      
                                  ),
                                  box(#uiOutput("files_ui"),
                                    uiOutput("select.folder_ser"),
                                    
                                      width = 10, background = "teal"
                                      
                                  ),
                                  box(#uiOutput("files_ui"),
                                    plotOutput("DimPlot_base_ser"),
                                    
                                    width = 10, background = "teal"
                                    
                                  ),
                                  
                                  
                                  
                                )
                                
                        ),
                        
                        
                        # SDA2Ser_projx------
                        tabItem(tabName = "SDA2Ser_projx",
                                h2("Project SDA Components on the Seurat Object"),
                                fluidRow(
                                  box(
                                    actionButton("PrevComp", "<< Prev Comp"),
                                    actionButton("NextComp", "Next Comp >>"),
                                   
                                    
                                    width = 5, background = "black"
                                  ),
                                  
                                  box(
                                    plotOutput("FeaturePlot_projx"),
                                    
                                      width = 10, background = "teal"
                                      
                                  ),
                                  box(title = "Pos. Top Genes", status = "info", solidHeader = TRUE, width = 4,
                                      tableOutput("packageTablePos")
                                  ),
                                  box(title = "Neg. Top Genes", status = "info", solidHeader = TRUE, width = 4,
                                      tableOutput("packageTableNeg")
                                  )
                                  
                                )
                                
                        ),
                        
                        
                       
                        # Gene loadings Heatmap ----------
                        tabItem(tabName = "GL_HM",
                                h2("Combined Gene Loadings Heatmap"),
                                fluidRow(
                                  box(
                                    title = "Threshold gene loading", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("Comb_GL_Hist_Zscore"),
                                    sliderInput(inputId = "GL_Zscore_slider", 
                                                label = "Gene Loading Z-score Value:", 
                                                min = 0, 
                                                max = 10, 
                                                value = 7, 
                                                step = 1),
                                    sliderInput(inputId = "GL_Kmeans_slider", 
                                                label = "K Cut + Kmeans", 
                                                min = 2, 
                                                max = 15, 
                                                value = 5, 
                                                step = 1),
                                    width = 10, background = "black"
                                  ),
                                  box(
                                    title = "Combine Gene Loadings (Z-score) Heatmap", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("Comb_GL_HM_Zscore"),
                                    width = 10, background = "teal"
                                  ),
                                  box(
                                    title = "Combine Gene Loadings (Z-score) Heatmap", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("Comb_GL_HM_Kmeans"),
                                    width = 10, background = "teal"
                                  ),
                                  
                                  
                                  
                                  ))#,
                        
                        
                        
                        # Save out----------
                        # tabItem(tabName = "SaveOut",
                        #         h2("Save the results for downstream analysis"),
                        #         fluidRow(
                        #           box(
                        #             title = "Save as Seurat Object", status = "primary", solidHeader = TRUE,
                        #             collapsible = TRUE,
                        #             actionButton("SaveAsSerObj", "Save as Seurat Obj"),
                        #             width = 5, background = "black"
                        #           )))
                        
                      ) #end of tabItems
                    ) #end of body
) #end UI

# Server --------

server <- function(input, output, session) {
  
  
  ## loading local files
  shinyFileChoose(input,'file_sda', session=session,roots=c(wd='.'))
  
  
  ## environment defaults
  envv=reactiveValues(y=NULL)
  envv$InfoBox_sub = "Load in SDA, use either the Prime-seq or Folder input tabs"
  
  envv$selected_files = NULL
  envv$input_obj_ls = list()
  
  ## controls how to handle pipe depending on data
  envv$Origin = "unk"
  
 
  
  #SDA files UI handler ---- 
  output$files_ui <- renderUI({
    req(input$SDAroot)
    path <- input$SDAroot
    files <- list.files(path=path, pattern=".rds", full.names=TRUE, recursive=F)

    files = files[!grepl("rawData", files)]
    files = files[!grepl("SDAtools", files)]
    files = files[!grepl("tSNE", files)]
    
    names(files) = gsub("sda.", "", gsub(".rds", "", basename(files)))
    
    if(file.exists(paste0(path, "/rawData_idDF.rds"))){
      DataMetaDF = readRDS( paste0(path, "/rawData_idDF.rds"))

      fileChoices = paste0(names(files), " - ", DataMetaDF[names(files), ]$ObjName)
    } else {
      fileChoices = basename(files)
    }

    # envv$input_SDA_files = files
    checkboxGroupInput("files_sda", "Available .rds files", choiceValues=basename(files), 
                       choiceNames = fileChoices)
  })
  
  
  #Ser Obj file UI handler ---- 
  
  output$select.folder_ser <- renderUI({
    req(input$SerObjroot)
    path <- input$SerObjroot
    files <- list.files(path=path, pattern=".rds", full.names=TRUE, recursive=F)
    
    files = files[!grepl("rawData", files)]
    files = files[!grepl("SDAtools", files)]
    files = files[!grepl("tSNE", files)]
    files = files[!grepl("DEres", files)]
    
    
    names(files) = gsub("sda.", "", gsub(".rds", "", basename(files)))

    # envv$input_SDA_files = files
    selectInput("file_ser", "Available .rds files", 
                choices = basename(files), 
                multiple = F)
  })
  
  
  # ### SDA local folder
  # output$select.folder <-
  #   renderUI(expr = selectInput(inputId = 'folder.name',
  #                               label = 'Folder Name',
  #                               choices = list.dirs(path = input$SDAroot,
  #                                                   full.names = FALSE,
  #                                                   recursive = FALSE)))
  
  
  
  source("app_Figs.R",local = TRUE)
  source("app_OE.R",local = TRUE)
  
  output$packageTablePos <- renderTable({
    
    req(envv$proc_obj_mat)
    req(envv$SDAcomp2Projx)
    
    topNfeats = 20

    sort(envv$proc_obj_mat[envv$SDAcomp2Projx, ], decreasing = T)[1:topNfeats] %>% names()
      
  
  }, digits = 1)
  
  output$packageTableNeg <- renderTable({
    
    req(envv$proc_obj_mat)
    req(envv$SDAcomp2Projx)
    
    topNfeats = 20
    
    sort(envv$proc_obj_mat[envv$SDAcomp2Projx, ], decreasing = F)[1:topNfeats] %>% names()
    
  }, digits = 1)
  
  
  
  
  output$GL_Zscore_slider <- renderUI({
    sliderInput("GL_Zscore_slider", "Gene Loading Z-score Value:", min = 0, max = 10, value = 7)
  })
  
  output$GL_Zscore_slider <- renderUI({
    sliderInput("GL_Kmeans_slider", "K Cut rows + Kmeans:", min = 2, max = 15, value = 5)
  })
  
  
  
  
}

shinyApp(ui, server)

