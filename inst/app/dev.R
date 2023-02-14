library(ggplot2)
library(dplyr)

DevMode = T

path = "/Volumes/Maggie/Work/OHSU/Bimber/Expts/RIRA_manuscript/data/sda/primeseq"


envv = list()
input = list()





DataMetaDF = data.frame(objID= c(510801L, 511966L, 508076L, 509929L, 508067L, 506842L, 
                                 510762L, 510759L, 510184L, 510134L, 505660L, 505421L,
                                 505374L, 513620L, 513984L, 514034L, 514273L, 514232L, 
                                 514215L, 514213L, 514166L, 514095L, 514084L, 514079L, 
                                 514060L,514274L, 511244L, 507934L, 510764L, 506832L, 514282L),
ObjName=c("CD8Ef", "CD4Tot", "CD4TotDS", "CD8Tot", 
          "CD8TotDS", "TotEf", "CD4Ef", "CD8Ef", 
          "TotEf", "TotEf", "TotEf", "NK50", 
          "TotEff50C", "MoMacDC", "Cd8Mem",
          "NK_NKT", "TotNaives", "CD4Mem",
          "RIRA", "CD8EffwNK", "Bcells",
          "other", "CD8Naives", "TNK",
          "CD8Ef", "CD4Naives", "NK", "NK", "CD4eff", "GdT", "FACS"
          ) )

rownames(DataMetaDF)= DataMetaDF$objID

DataMetaDF = DataMetaDF[naturalsort::naturalorder(DataMetaDF$ObjName),]

# saveRDS(DataMetaDF, paste0(path, "/rawData_idDF.rds"))



files <- list.files(path=path, pattern=".rds", full.names=TRUE, recursive=F)

files = files[!grepl("rawData", files)]
files = files[!grepl("SDAtools", files)]
files = files[!grepl("tSNE", files)]

names(files) = gsub("sda.", "", gsub(".rds", "", basename(files)))

DataMetaDF = DataMetaDF[names(files),]



DataMetaDF$rowId = 1:nrow(DataMetaDF)

DataMetaDF

#select which files to load
files = files[names(files) %in% DataMetaDF$objID[c(1:16, 18, 19, 20, 22)]]



input$files = basename(files)
envv$selected_files = as.character(files)

envv$input_obj_ls = lapply(envv$selected_files, readRDS)
names(envv$input_obj_ls)  = gsub(".rds", "", basename(envv$selected_files))
names(envv$input_obj_ls) = paste0(names(envv$input_obj_ls), "_", DataMetaDF[gsub("sda.", "", names(envv$input_obj_ls)),]$ObjName)



envv$input_obj_ls = lapply(envv$input_obj_ls , function(sdaobj){
  sdaobj = ShinySDA:::AddCompStats(sdaobj)
  ShinySDA:::CalculateFilters(sdaobj, threshold = 100, plot = F, quantThr = .95)
  
})



## Gene loadings
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







#if c = T and s = T => Z-score, but if c=F and s=T, => just / sdtv

envv$proc_obj_mat_zscore = as.matrix(scale(envv$proc_obj_mat, center = T, scale = T))


KeepComps = as.logical(unlist(lapply(envv$input_obj_ls , function(sdaobj){
  !sdaobj$FailingFilters
})))

# envv$top_loaded_genes = lapply(envv$input_obj_ls, function(x){
# c(sort(x$loadings[[1]][1,], decreasing = T)[1:topNfeats] %>% names(),
#   sort(x$loadings[[1]][1,], decreasing = F)[1:topNfeats] %>% names())
# }) ) %>% unlist() %>% unique()

# cowplot::plot_grid(plotlist = lapply(1:2, function(xN){
#   plot(sort(abs(envv$proc_obj_mat_zscore[xN,]), decreasing = T), pch=20)
# }), ncol=2)


# OOSAP:::FindElbow(sort(abs(envv$proc_obj_mat_zscore[1,])))

topNfeats = 5

envv$top_loaded_genes = lapply(1:nrow(envv$proc_obj_mat_zscore), function(xN){
  x=envv$proc_obj_mat_zscore[xN,]
  c(sort(x, decreasing = T)[1:topNfeats] %>% names(),
   sort(x, decreasing = F)[1:topNfeats] %>% names())
}) %>% unlist() %>% unique()


envv$top_loaded_genes2 = lapply(1:nrow(envv$proc_obj_mat), function(xN){
  x=envv$proc_obj_mat[xN,]
  c(sort(x, decreasing = T)[1:topNfeats] %>% names(),
    sort(x, decreasing = F)[1:topNfeats] %>% names())
}) %>% unlist() %>% unique()



envv$col_sums = colSums(abs(envv$proc_obj_mat_zscore[,]))
hist(envv$col_sums, breaks=200, main="Zscore distribution")


envv$row_sums = asinh(rowSums(abs(envv$proc_obj_mat_zscore)))
hist(envv$row_sums, breaks=200, main="Zscore distribution")


envv$Loading_thr_genes = colnames(envv$proc_obj_mat[,])[envv$col_sums > 280]
envv$Loading_thr_comps = rownames(envv$proc_obj_mat)[envv$row_sums > 9.5]
round(length(envv$Loading_thr_comps)/nrow(envv$proc_obj_mat)*100, 1)


ggvenn::ggvenn(list(topLoadedGenesSc = envv$top_loaded_genes ,
                    # topLoadedGenes = envv$top_loaded_genes2 ,
                    ThreshGenes = envv$Loading_thr_genes ))



## check cell scores

# DGE_mat = envv$input_obj_ls$sda.507934_NK30$scores %*% (envv$input_obj_ls$sda.507934_NK30$loadings[[1]])
# 
# Scores_mat = DGE_mat %*% t(envv$input_obj_ls$sda.507934_NK30$loadings[[1]])
# Scores_mat2 = DGE_mat[, unique(c(envv$Loading_thr_genes, envv$top_loaded_genes))] %*% t(envv$input_obj_ls$sda.507934_NK30$loadings[[1]][, unique(c(envv$Loading_thr_genes, envv$top_loaded_genes))])
# 
# head(Scores_mat)
# head(envv$input_obj_ls$sda.507934_NK30$scores)
# 
# plot(x=envv$input_obj_ls$sda.507934_NK30$scores[,1],
#      y=(Scores_mat[,1]), pch=20, xlab="Original Score C1 from SDA model", ylab="Dot product recomputed C1 Score")
# 
# plot(x=envv$input_obj_ls$sda.507934_NK30$scores[,1],
#      y=(Scores_mat2[,1]), pch=20, xlab="Original Score C1 from SDA model", ylab="Dot product recomputed C1 Score")
# 
# plot(x=Scores_mat[,1],
#      y=(Scores_mat2[,1]), pch=20, xlab="Dot product full recomputed C1 Score", ylab="Dot product reduced recomputed C1 Score")
# 
# 
# plot(x=envv$input_obj_ls$sda.507934_NK30$scores[,1],
#      y=envv$input_obj_ls$sda.507934_NK30$scores[,2], pch=20)


if(!DevMode) pheatmap::pheatmap(t(asinh(envv$proc_obj_mat_zscore[envv$Loading_thr_comps, 
                                                                 unique(c(envv$Loading_thr_genes, envv$top_loaded_genes))])),
                                clustering_distance_rows = "euclidean",
                                clustering_distance_cols = "euclidean", 
                                clustering_method = "ward.D2",
                                color = c("navy", "dodgerblue", "grey", "grey", "gold", "red"),
                                filename = paste0(path, "/GL_heatmap_Combo.png"),
                                width = 15, height = 30,
                                cutree_rows = NA, cutree_cols = NA#, kmeans_k = 10
)




# if(!DevMode) pheatmap::pheatmap(t(asinh(envv$proc_obj_mat[, 
#                                                                  unique(c(envv$Loading_thr_genes, envv$top_loaded_genes))])),
#                                 clustering_distance_rows = "euclidean",
#                                 clustering_distance_cols = "euclidean", 
#                                 clustering_method = "ward.D2",
#                                 color = c("navy", "dodgerblue", "grey", "grey", "gold", "red"),
#                                 # filename = paste0(path, "/GL_heatmap_Combo.png"),
#                                 width = 20, height = 15,
#                                 cutree_rows = NA, cutree_cols = NA#, kmeans_k = 10
# )

if(!DevMode) pheatmap::pheatmap(t(asinh(envv$proc_obj_mat_zscore[KeepComps, 
                                                                 unique(c(envv$Loading_thr_genes, envv$top_loaded_genes))])),
                                clustering_distance_rows = "euclidean",
                                clustering_distance_cols = "euclidean", 
                                clustering_method = "ward.D2",
                                color = c("navy", "dodgerblue", "grey", "grey", "gold", "red"),
                                # filename = paste0(path, "/GL_heatmap_Combo.png"),
                                width = 20, height = 15,
                                cutree_rows = NA, cutree_cols = NA#, kmeans_k = 10
)



if(!DevMode) pheatmap::pheatmap(t(asinh(envv$proc_obj_mat_zscore[envv$Loading_thr_comps, envv$Loading_thr_genes])),
                                clustering_distance_rows = "euclidean",
                                clustering_distance_cols = "euclidean", 
                                clustering_method = "ward.D2",
                                color = c("navy", "dodgerblue", "grey", "grey", "gold", "red"),
                                # filename = paste0(path, "/GL_heatmap_thr.png"),
                                width = 30, height = 15,
                                cutree_rows = NA, cutree_cols = NA#, kmeans_k = 10
                   )


if(!DevMode) pheatmap::pheatmap(t(asinh(envv$proc_obj_mat_zscore[envv$Loading_thr_comps, envv$top_loaded_genes])), 
                                clustering_distance_rows = "euclidean",
                                clustering_distance_cols = "euclidean", 
                                clustering_method = "ward.D2",
                                color = c("navy", "dodgerblue", "grey", "grey", "gold", "red"),
                                # filename = paste0(path, "/GL_heatmap_toploaded.png"),
                                width = 30, height = 100,
                                cutree_rows = NA, cutree_cols = NA#, kmeans_k = 10
)



# if(!DevMode) pheatmap::pheatmap(asinh(envv$proc_obj_mat_zscore[,envv$Loading_thr_genes]), 
#                    clustering_distance_rows = "euclidean",
#                    clustering_distance_cols = "euclidean", 
#                    clustering_method = "ward.D2",
#                    color = c("navy", "dodgerblue", "grey", "grey", "gold", "red"),
#                    cutree_rows = NA, cutree_cols = NA, kmeans_k = 4
# )


head.path = paste0( "ComboGL_genes_" , #paste0(gsub("sda.", "", sort(names(envv$input_obj_ls))), collapse = "."),
                    "tSNE")
tSNE_n.iter = 1000
tsnepp = 20 

envv$tSNE_path = paste0(path, "/", head.path, "_it", tSNE_n.iter,"_pp",tsnepp, ".rds")


envv$tsne_GL_genes = scCustFx:::.run_tSNE_mat(myMat = 
                                   t(as.matrix(envv$proc_obj_mat_zscore[envv$Loading_thr_comps, 
                                                unique(c(envv$Loading_thr_genes, envv$top_loaded_genes))])),
                                   save_path = envv$tSNE_path, 
                                   tsnepp=tsnepp, tSNE_n.iter=tSNE_n.iter, 
                                   num_threads = 8, check_duplicates = F, 
                                   selectComps = NULL, overwrite=T )
rownames(envv$tsne_GL_genes$Y)  <- unique(c(envv$Loading_thr_genes, envv$top_loaded_genes)) #envv$top_loaded_genes



envv$umap_GL_genes = scCustFx:::RunUMAP.Matrix(DGEmat=  
                      t(as.matrix(envv$proc_obj_mat_zscore[envv$Loading_thr_comps, 
                                      unique(c(envv$Loading_thr_genes, envv$top_loaded_genes))])),
                       assay = NULL,
                       n.neighbors = 30L,
                       n.components = 2L,
                       metric = "correlation",
                       n.epochs = NULL,
                       learning.rate = 1.0,
                       min.dist = 0.3,
                       spread = 1.0,
                       set.op.mix.ratio = 1.0,
                       local.connectivity = 1L,
                       repulsion.strength = 1,
                       negative.sample.rate = 5,
                       a = NULL,
                       b = NULL,
                       seed.use = 42,
                       metric.kwds = NULL,
                       angular.rp.forest = FALSE,
                       reduction.key = 'UMAP_',
                       verbose = TRUE)



scCustFx:::plot2Dredux(plotDF = as.data.frame(envv$umap_GL_genes), title = "UMAP of SDA Gene Loadings", NfeatPerClus = 15)
scCustFx:::plot2Dredux(plotDF = as.data.frame(envv$tsne_GL_genes$Y), title = "tSNE of SDA Gene Loadings", NfeatPerClus = 15)











input$species = "human"

library(AnnotationHub) # source("https://bioconductor.org/biocLite.R"); biocLite("AnnotationHub")
library(clusterProfiler) # source("https://bioconductor.org/biocLite.R"); biocLite("clusterProfiler")

hub <- AnnotationHub()

# if(input$species == "mouse") qrhub="org.MM.eg"
if(input$species == "human") qrhub="org.Hs.eg.db"
# if(input$species == "rhesus") qrhub="org.Mmu.eg.db"

RefGenome.names <- query(hub, qrhub)#org.Hs.eg.db  org.MM.eg


# if(input$species == "mouse") qrhub.id="AH84123"
# if(input$species == "human") qrhub.id="org.Hs.eg.db"
# if(input$species == "rhesus") qrhub.id="org.Mmu.eg.db"

print(RefGenome.names$ah_id)
RefGenome <- hub[[RefGenome.names$ah_id]] 


envv$GOAnn <- list(RefGenome = RefGenome, RefGenome.names = RefGenome.names)


GO_data <- list()

for (i in 1:max(as.numeric(levels(tsneDF$km)))){
  print(i)
  ego <- enrichGO(gene = subset(tsneDF, km==i)$labs1,
                  universe = envv$common_col_names,
                  OrgDb = RefGenome,
                  keyType = 'SYMBOL',
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 1,
                  qvalueCutoff = 1)
  
  ego@result$Enrichment <- ShinySDA::frac_to_numeric(ego@result$GeneRatio)/ShinySDA::frac_to_numeric(ego@result$BgRatio)
  
  ego@result$GeneOdds <- unlist(lapply(strsplit(ego@result$GeneRatio, "/", fixed = TRUE), function(x){ x<-as.numeric(x) ;x[1] / (x[2]-x[1])}))
  ego@result$BgOdds <- unlist(lapply(strsplit(ego@result$BgRatio, "/", fixed = TRUE), function(x){ x<-as.numeric(x) ;x[1] / (x[2]-x[1])}))
  
  GO_data[[paste0("Clus",i)]] = ego@result
}



cowplot::plot_grid(plotlist = lapply(names(GO_data), function(x){
  ShinySDA:::go_volcano_plot(x=GO_data, component=x, extraTitle="") + theme_bw(base_size = 12)}),
  ncol = 3)












head.path = paste0( "ComboGL_comps_" , paste0(gsub("sda.", "", sort(names(envv$input_obj_ls))), collapse = "."), "_tSNE")
tSNE_n.iter = 3000
tsnepp = 35 

envv$tSNE_path = paste0(path, "/", head.path,tSNE_n.iter,"_pp",tsnepp, ".rds")


envv$tsne_GL_comps = ShinyMultiSDA:::.run_tSN_mat(myMat = 
                                                    (as.matrix(envv$proc_obj_mat_zscore[,envv$Loading_thr_genes])),
                                                  save_path = envv$tSNE_path, 
                                                  tsnepp=tsnepp, tSNE_n.iter=tSNE_n.iter, 
                                                  num_threads = 8, check_duplicates = F, 
                                                  selectComps = NULL, overwrite=T )

tsneDF <- as.data.frame(envv$tsne_GL_comps$Y)
nrow(tsneDF)


rownames(tsneDF)  <- rownames(as.matrix(envv$proc_obj_mat_zscore))
colnames(tsneDF) <- c("tSNE1", "tSNE2")

tsneDF$labs1 = rownames(tsneDF)

tsneDF$sums <- envv$row_sums



# tsneDF$km = factor(kmeans( t(as.matrix(envv$proc_obj_mat_zscore[,envv$Loading_thr_genes])), centers = 4)$cluster)
tsneDF$km = factor(kmeans(tsneDF[,1:2] , centers = 4)$cluster)

# group by km
keepLabs <- tsneDF %>% group_by(km) %>% 
  sample_n(10, replace = TRUE) %>% 
  pull(labs1)

tsneDF$labs2 = ""

tsneDF[keepLabs, ]$labs2 = keepLabs




ggplot(tsneDF, aes(tSNE1, tSNE2, col=sums)) +
  geom_point(size = 1) + theme_bw() + 
  scale_color_distiller(palette = "Spectral") +
  theme(legend.position = "bottom") +
  ggtitle("tSNE SDA qc Components\n Sum absolute-cell-scores normalized by its mean \n ")

ggplot(tsneDF, aes(tSNE1, tSNE2, col=km)) +
  geom_point(size = 1) + theme_bw() +
  # scale_color_distiller(palette = "Spectral") +
  theme(legend.position = "bottom") +
  ggrepel::geom_text_repel(aes(label=labs2, col="10"), size=3, max.overlaps =100, segment.color="grey50", nudge_y = 0.25) +
  ggtitle("tSNE SDA qc Components\n Sum absolute-cell-scores normalized by its mean \n ") + 
  scale_color_manual(values = scCustFx:::ColorTheme()$col_vector)

