


output$FeaturePlot_projx <- renderPlot({
  if(is.null(envv$input_obj_ser) | is.null(envv$proc_obj_mat)){
    plot(x=0, y=0, main="Load an SDA")
  } else {
    
    require(Seurat)
    
    if(is.null(envv$SDAcomp2Projx)){
      envv$SDAcomp2Projx = 1
    } 
    
    envv$SDAgenes2Projx = names(envv$proc_obj_mat[envv$SDAcomp2Projx, ])
    envv$SDAgenes2Projx = envv$SDAgenes2Projx[envv$SDAgenes2Projx %in% rownames(envv$input_obj_ser)]
    
    SDAscore = envv$proc_obj_mat[envv$SDAcomp2Projx, envv$SDAgenes2Projx] %*% envv$input_obj_ser@assays$RNA@data[envv$SDAgenes2Projx,]
    
    
    envv$input_obj_ser = AddMetaData(envv$input_obj_ser , as.data.frame(t(asinh(SDAscore))), 
                                     col.name = rownames(envv$proc_obj_mat)[envv$SDAcomp2Projx])
    
    Seurat::FeaturePlot(envv$input_obj_ser, 
                        features = c(rownames(envv$proc_obj_mat)[envv$SDAcomp2Projx]), 
                        # max.cutoff = 'q99', min.cutoff = 'q01',
                        order = T) + coord_flip()  + scale_y_reverse() +
      theme_classic(base_size = 14) + #NoLegend()+
      theme(axis.line = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank() #,plot.title = element_blank()
      )  + ggtitle(rownames(envv$proc_obj_mat)[envv$SDAcomp2Projx]) &
      ggplot2::scale_colour_gradientn(colours = c("navy", "dodgerblue", "gold", "red"))
    
    
    
    
  
  
  
    
    
    
  }
})



output$DimPlot_base_ser <- renderPlot({
  if(is.null(envv$input_obj_ser)){
    plot(x=0, y=0, main="Load an SDA")
  } else {
   require(Seurat)
    DimPlot(envv$input_obj_ser)
    
  }
})


output$Comb_GL_HM_Zscore <- renderPlot({
  if(is.null(envv$proc_obj_mat)){
    plot(x=0, y=0, main="Load an SDA")
  } else {
    
    col_sums = envv$col_sums
    
    envv$Loading_thr_genes = colnames(envv$proc_obj_mat)[col_sums > input$GL_Zscore_slider]
    
    pheatmap::pheatmap(asinh(scale(envv$proc_obj_mat[,envv$Loading_thr_genes], center = T, scale = T)), 
                       clustering_distance_rows = "euclidean",
                       clustering_distance_cols = "euclidean", 
                       clustering_method = "ward.D2",
                       cutree_rows = input$GL_Kmeans_slider, cutree_cols = NA
    )    
    
  }
})

output$Comb_GL_HM_Kmeans <- renderPlot({
  if(is.null(envv$proc_obj_mat)){
    plot(x=0, y=0, main="Load an SDA")
  } else {
    
    col_sums = envv$col_sums
    
    envv$Loading_thr_genes = colnames(envv$proc_obj_mat)[col_sums > input$GL_Zscore_slider]
    
    pheatmap::pheatmap(asinh(scale(envv$proc_obj_mat[,envv$Loading_thr_genes], center = T, scale = T)), 
                       clustering_distance_rows = "euclidean",
                       clustering_distance_cols = "euclidean", 
                       clustering_method = "ward.D2",
                       cutree_rows = input$GL_Kmeans_slider, cutree_cols = NA, kmeans_k = input$GL_Kmeans_slider
    )    
    
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
