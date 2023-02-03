
envv = list()
input = list()

path = "/Volumes/Maggie/Work/OHSU/Bimber/Expts/RIRA_manuscript/data/sda/primeseq"

files <- list.files(path=path, pattern=".rds", full.names=TRUE, recursive=F)

files = files[!grepl("rawData", files)]
files = files[!grepl("SDAtools", files)]

input$files = basename(files[1:2])
envv$selected_files = files[1:2]

envv$input_obj_ls = lapply(envv$selected_files, readRDS)
names(envv$input_obj_ls)  = gsub(".rds", "", basename(envv$selected_files))




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



hist(colSums(abs(envv$proc_obj_mat)), breaks = 200)
abline(v=0.5, col="red", lty=2)
abline(v=1, col="red", lty=2)
abline(v=2, col="red", lty=2)

library(ggplot2)

col_sums = colSums(abs((envv$proc_obj_mat)))


ggplot(data.frame(col_sums), aes(col_sums)) + 
  geom_histogram(bins = 200, fill = "gray") + 
  geom_vline(xintercept = c(0.5, 1, 2, 4, 7), color = "red", linetype = 2) + 
  annotate("text", x = 0.8, y = 1500, label = paste0(">0.5:   ", length(colnames(envv$proc_obj_mat)[colSums(abs(envv$proc_obj_mat)) > 0.5])), hjust = 1, angle = 90) +
  annotate("text", x = 1.5, y = 1500, label = paste0(">1:   ", length(colnames(envv$proc_obj_mat)[colSums(abs(envv$proc_obj_mat)) > 1])), hjust = 1, angle = 90) +
  annotate("text", x = 2.5, y = 1500, label = paste0(">2:   ", length(colnames(envv$proc_obj_mat)[colSums(abs(envv$proc_obj_mat)) > 2])), hjust = 1, angle = 90) + 
  annotate("text", x = 4.5, y = 1500, label = paste0(">4:   ", length(colnames(envv$proc_obj_mat)[colSums(abs(envv$proc_obj_mat)) > 4])), hjust = 1, angle = 90) + 
  annotate("text", x = 7.5, y = 1500, label = paste0(">4:   ", length(colnames(envv$proc_obj_mat)[colSums(abs(envv$proc_obj_mat)) > 7])), hjust = 1, angle = 90) + 
  theme_bw() + ggtitle("Thresholding of Gene loading \n No. genes to the right")




envv$Loading_thr_genes = colnames(envv$proc_obj_mat)[col_sums > 4]


#if c = T and s = T => Z-score, but if c=F and s=T, => just / sdtv
pheatmap::pheatmap(asinh(scale(envv$proc_obj_mat[,envv$Loading_thr_genes], center = T, scale = T)))


