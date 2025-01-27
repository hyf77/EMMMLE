##1. Package loading
##--------------------------------------------------------------------------------------------
library(Seurat)
library(SeuratData)
File_path0<-"D:/scGeneNet_test/Real data analysis/p500/Benchmarking on scRNA-seq data/Zheng's file/"
File_path1<-"D:/scGeneNet_test/Real data analysis/p300/Benchmarking on scRNA-seq data/Zheng's file/"
#2.1 Data loading and pre-process
load(paste(File_path0,"Data/expression_profile_use.Rdata",sep = ""))
load(paste(File_path0,"Data/gene_GRN_intop500HVG.Rdata",sep = ""))
load("D:/scGeneNet_test/Real data analysis/p500/Benchmarking on scRNA-seq data/Zheng's file/Evaluation/Obj_use/TF_list_all.Rdata")
#2.2 Choose 300 high variable genes as gene set of interest
P3se_ifnb_stim = CreateSeuratObject(counts = expression_profile_use,min.cells = 3)
P3se_ifnb_stim <- NormalizeData(P3se_ifnb_stim, normalization.method = "LogNormalize", scale.factor = 10000)
P3se_ifnb_stim <- FindVariableFeatures(P3se_ifnb_stim,nfeatures = 300)
variable_gene_stim<-Seurat::VariableFeatures(P3se_ifnb_stim)
gene_names<-unique(c(variable_gene_stim[which(variable_gene_stim %in% TF_list_all[[1]])],variable_gene_stim))
gene_GRN<-c(gene_names[which(gene_names %in% TF_list_all[[1]])],setdiff(gene_names,gene_names[which(gene_names %in% TF_list_all[[1]])]))
save(gene_GRN,file = paste(File_path1,"Data/gene_GRN_intop300HVG.Rdata",sep = ""))
##3.1 resample for calculating stability
for (i in 1:10) {
  sampled_index<-sample(ncol(expression_profile_use), size = 0.9 * ncol(expression_profile_use))
  sampled_data <- expression_profile_use[,sampled_index]
  sampled_cluster_true<-cluster_true[sampled_index]
  save(sampled_index,file = paste(File_path1,"Sampled/Data/sampled_index_rep",i,".Rdata",sep = ""))
  save(sampled_data,file = paste(File_path1,"Sampled/Data/sampled_data_rep",i,".Rdata",sep = ""))
  save(sampled_cluster_true,file = paste(File_path1,"Sampled/Data/sampled_cluster_true_rep",i,".Rdata",sep = ""))
}