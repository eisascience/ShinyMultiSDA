



observeEvent(input$loadSDAfiles, {
  req(input$files)
  full_files = paste0(input$SDAroot, '/', input$files)
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
  
  

  
})

