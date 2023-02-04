

output$Comb_GL_HM_Zscore <- renderPlot({
  if(is.null(envv$proc_obj_mat)){
    plot(x=0, y=0, main="Load an SDA")
  } else {
    
    col_sums = envv$col_sums
    
    envv$Loading_thr_genes = colnames(envv$proc_obj_mat)[col_sums > input$GL_Zscore_slider]
    
    pheatmap::pheatmap(asinh(scale(envv$proc_obj_mat[,envv$Loading_thr_genes], center = T, scale = T)))
    
    
  }
})


output$Comb_GL_Hist_Zscore <- renderPlot({
  if(is.null(envv$proc_obj_mat)){
    plot(x=0, y=0, main="Load an SDA")
  } else {
    
    
    col_sums = envv$col_sums    
    
    ggplot(data.frame(col_sums), aes(col_sums)) + 
      geom_histogram(bins = 200, fill = "gray") + 
      geom_vline(xintercept = c(0.5, 1, 2, 4, 7), color = "red", linetype = 2) + 
      annotate("text", x = 0.8, y = 1500, label = paste0(">0.5:   ", length(colnames(envv$proc_obj_mat)[colSums(abs(envv$proc_obj_mat)) > 0.5])), hjust = 1, angle = 90) +
      annotate("text", x = 1.5, y = 1500, label = paste0(">1:   ", length(colnames(envv$proc_obj_mat)[colSums(abs(envv$proc_obj_mat)) > 1])), hjust = 1, angle = 90) +
      annotate("text", x = 2.5, y = 1500, label = paste0(">2:   ", length(colnames(envv$proc_obj_mat)[colSums(abs(envv$proc_obj_mat)) > 2])), hjust = 1, angle = 90) + 
      annotate("text", x = 4.5, y = 1500, label = paste0(">4:   ", length(colnames(envv$proc_obj_mat)[colSums(abs(envv$proc_obj_mat)) > 4])), hjust = 1, angle = 90) + 
      annotate("text", x = 7.5, y = 1500, label = paste0(">7:   ", length(colnames(envv$proc_obj_mat)[colSums(abs(envv$proc_obj_mat)) > 7])), hjust = 1, angle = 90) + 
      theme_bw() + ggtitle("Thresholding of Gene loading \n No. genes to the right")
    
    
  }
})
