observeEvent(input$NextComp, {
  
  if(!is.null(envv$proc_obj_mat)){
    if(is.null(envv$SDAcomp2Projx)){
      envv$SDAcomp2Projx = 1
    } else {
      if(envv$SDAcomp2Projx > ncol(envv$proc_obj_mat)){
        envv$SDAcomp2Projx = envv$SDAcomp2Projx
      } else {
        envv$SDAcomp2Projx = envv$SDAcomp2Projx + 1
      }
    }
  }
  
  
})

observeEvent(input$PrevComp, {
  
  if(!is.null(envv$proc_obj_mat)){
    if(is.null(envv$SDAcomp2Projx)){
      envv$SDAcomp2Projx = 1
    } else {
      if(envv$SDAcomp2Projx == 1){
        envv$SDAcomp2Projx = envv$SDAcomp2Projx
      } else {
        envv$SDAcomp2Projx = envv$SDAcomp2Projx - 1
      }
    }
  }
  
  
})


observeEvent(input$loadSerObj, {
  req(input$file_ser)
  require(Seurat)
  require(SeuratDisk)
  
  full_file = paste0(input$SerObjroot, '/', input$file_ser)
  print(full_file)
  envv$selected_file_ser = full_file
  
  envv$input_obj_ser = readRDS(full_file)
  

  print("Ser Obj loaded")
  
  
})




observeEvent(input$loadSDAfiles, {

  req(input$files_sda)
  print(input$SDAroot)
  print(input$files_sda)
  
  full_files = paste0(input$SDAroot, '/', input$files_sda)
  print(full_files)
  envv$selected_files = full_files
  envv$input_obj_ls = lapply(full_files, readRDS)
  names(envv$input_obj_ls)  = gsub(".rds", "", basename(full_files))
  
  print("files loaded")


})

observeEvent(input$CombineSDAs, {
  
  req(envv$input_obj_ls)
  print(length(envv$input_obj_ls))
  
  envv$common_col_names <- Reduce(intersect, lapply(envv$input_obj_ls, function(x){
    colnames(x$loadings[[1]])
  }))
  
  envv$proc_obj_ls = lapply(envv$input_obj_ls, function(x){
    x$loadings[[1]][,envv$common_col_names]
  })
  
  print(length(envv$proc_obj_ls))
  names(envv$proc_obj_ls) = names(envv$input_obj_ls)
  print( names(envv$proc_obj_ls))
  
  envv$proc_obj_mat =  do.call(rbind, envv$proc_obj_ls)

  rownames(envv$proc_obj_mat) = unlist(lapply(1:length(envv$proc_obj_ls), function(xN){
    paste0( names(envv$proc_obj_ls)[xN], "_V", 1:nrow(envv$proc_obj_ls[[xN]]))
  }))
  
  envv$col_sums = colSums(abs((envv$proc_obj_mat)))
  
  envv$TotNComp = ncol(envv$proc_obj_mat)
  
})

# observeEvent(input$goButton, {
#   selected_value <- input$slider
#   
#   # Perform the desired calculation or plotting using the selected value
#   # ...
# })

