



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
  # print(length(envv$input_obj_ls))
  
  envv$common_col_names <- Reduce(intersect, lapply(envv$input_obj_ls, function(x){
    colnames(x$loadings[[1]])
  }))
  
  
  envv$proc_obj_ls = lapply(envv$input_obj_ls, function(x){
    x$loadings[[1]][,envv$common_col_names]
  })
  names(envv$proc_obj_ls) = names(envv$input_obj_ls)
  
  envv$proc_obj_mat =  do.call(rbind, envv$proc_obj_ls)
  
  rownames(envv$proc_obj_mat) = unlist(lapply(1:length(envv$proc_obj_ls), function(xN){
    paste0( names(envv$proc_obj_ls)[xN], "_V", 1:nrow(envv$proc_obj_ls[[xN]]))
  }))
  
  envv$col_sums = colSums(abs((envv$proc_obj_mat)))
  
})

# observeEvent(input$goButton, {
#   selected_value <- input$slider
#   
#   # Perform the desired calculation or plotting using the selected value
#   # ...
# })

