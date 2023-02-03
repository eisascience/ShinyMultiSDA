
# Envv update fx ---------

#' This function to update the environment list when loading from local
#'
#' @param envv environment list associated with ShinySDA
#' @param input environment list associated with innput paramters in ShinySDA
#' @return the updated environment list
List_files_local_evv <- function(envv, input){
  
  # print(head(paste0(input$SDAroot, "/", input$folder.name)))
  
  envv$path2SDA_dyn <- paste0(input$SDAroot, "/", input$folder.name)
  
  if(file.exists(envv$path2SDA_dyn)) {
    
    print(envv$path2SDA_dyn)
    
    
    
  } else { 
    # updateTextInput(session, "loadSDAmsg", value = "File not found")
    envv$InfoBox_sub = "File not found loading error"
    
    
    print(envv$path2SDA_dyn)
    
  }
  return(envv)
}