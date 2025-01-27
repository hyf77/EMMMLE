##########################################
rm(list = ls())
#1. Loading package 
##--------------------------------------------------------------------------------------------
library(Matrix)
library(MASS)
library(glasso)
library(mclust)
library(mltools)
library(data.table)
library(doParallel)
library(PLNmodels)
library(scGeneNet)
library(cluster)
library(Seurat)
library(pryr)
library(philentropy)
library(umap)
library(igraph)
library(stats)
library(fpc)
library(XMRF)
library(CVXR)
library(orthopolynom)
library(Rcpp)
library(RcppArmadillo)
library(PLNetNew)
###
library(EMMMLE)
library(foreach)
method_name<-'Glasso'
num_core<-10
cl <- makeCluster(num_core, outfile ="debug_glasso.txt")
registerDoParallel(cl)
##--------------------------------------------------------------------------------------------
File_path<-"/home/ruibinxi_pkuhpc/lustre1/huangyf/R_project/scGeneNet/run_true/Benchmarking_on_scRNA_data/Zheng_file"
##2.2 Path for load used data
initial_load_path<-paste(File_path,"/Method/EMPLN/Output_Obj/EMPLN_init_4group.Rdata",sep = "")
gene_GRN_path<-paste(File_path,"/Data/gene_GRN_intop300HVG_4group.Rdata",sep = "")

##2.3 Path for save result
result_save_path<-paste(File_path,"/Method/",method_name,"/Output_Obj",sep = "")
##--------------------------------------------------------------------------------------------

##3. Load the data
##--------------------------------------------------------------------------------------------
load(initial_load_path)
load(gene_GRN_path)
##--------------------------------------------------------------------------------------------
##4. Set the parameter
##--------------------------------------------------------------------------------------------
p_TF<-8
celltypes_num<-4
p_GRN<-length(gene_GRN)
p_nonTF<-p_GRN - p_TF
##
run_fun <- function(g){
  ##2.3 Path for save result
  Glasso_save_path<-paste(File_path,"/Method/Glasso/Output_Obj",sep = "")
  print(g)
  #############################
  ########pre function
  #tr
  tr<-function(mat_1){
    return(sum(diag(as.matrix(mat_1))))
  }
  ##
  log_det<-function(matrix_use){
    eigen_res<-eigen(matrix_use)
    eig_value<-eigen_res$values
    eig_value<-ifelse(eig_value<=0,1e-8,eig_value)
    return(sum(log(eig_value)))
  }

  BIC_glasso<-function(data_1,pre_mat_array){
    #data 1: n*p
    #pre_mat_array: p*p*L
    sample_size<-nrow(data_1)
    cov_mat<-cov(data_1)
    BIC_vec<-c()
    for(l in 1:dim(pre_mat_array)[3]){
      Omega_use<-(pre_mat_array[,,l] + t(pre_mat_array[,,l]))/2
      BIC_vec<-c(BIC_vec,sample_size*(-1*log_det(Omega_use) + tr(cov_mat %*% Omega_use))+
                   log(sample_size)*sum(ifelse(Omega_use != 0,1,0))/2)
    }
    return(BIC_vec)
  }

  ##exact_zeroindex
  exact_zeroindex<-function(support = NULL){
    zeroindex_res<-NULL
    if(!is.null(support)){
      zeroindex_res<-matrix(NA,nrow = 0,ncol = 2)
      ##
      zero_mat<-ifelse(support==0,1,0)
      for(i in 1:(nrow(zero_mat) - 1)){
        for(j in (i+1): nrow(zero_mat) ){
          if(zero_mat[i,j] == 1){
            zeroindex_res<-rbind(zeroindex_res,c(i,j))
          }
        }
      }
    }
    ##
    return(zeroindex_res)
  }
  #########
  ##
  load(file = initial_load_path)
  ##
  locator_vec<-EMPLN_init$celltypes_label
  ##
  data_obs<-as.matrix(EMPLN_init$obs_mat)
  exp_count_GRN<-data_obs[,1:p_GRN]
  S_depth<-EMPLN_init$S_depth
  rm(EMPLN_init);gc()
  zero_mat<-diag(p_GRN)
  zero_mat[1:p_TF,]<-1
  zero_mat[,1:p_TF]<-1
  zero_mat_use<-exact_zeroindex(support = zero_mat)
  ##
  exp_nor<-log((as.matrix(exp_count_GRN)+1)/as.vector(S_depth))
  lambda_vec<-exp(seq(from=(log(1e-3)),to=(log(10)),length.out = 80))
  ####
  exp_nor_group<-exp_nor[which(locator_vec == g),]
  cov_X = cov(exp_nor_group)
  cov_X = cov_X + 1e-3 * diag(dim(cov_X)[1])
  ##
  graph_lasso_array<-array(NA,dim = c(p_GRN,p_GRN,length(lambda_vec)))
  for(l in 1:length(lambda_vec)){
    print(paste("lambda: ",l,sep = ""))
    glasso_res<-glasso(s = cov_X, rho = lambda_vec[l], penalize.diagonal = TRUE,thr = 1e-3,maxit = 1000,
                       zero = zero_mat_use)
    graph_lasso_array[,,l]<-glasso_res$wi
  }
  BIC_glasso_vec<-BIC_glasso(exp_nor_group,graph_lasso_array)

  ##
  save(BIC_glasso_vec,file = paste(Glasso_save_path,"/BIC_Glasso_vec_pen_top300GRN_g",g,".Rdata",sep = ""))
  save(graph_lasso_array,file = paste(Glasso_save_path,"/graph_lasso_array_pen_top300GRN_g",g,".Rdata",sep = ""))
}

repres_glasso<-foreach (
  g = 1:celltypes_num,
  .inorder = TRUE,
  .export = ls(.GlobalEnv),
  .packages = c('MASS','glasso',
                'mclust','mltools','data.table','PLNmodels','dplyr',
                "cluster","pryr","philentropy","umap","Seurat","igraph","Rcpp","scGeneNet","Matrix","stats","GMPR","peakRAM")
)%dopar%{
  ##--------------------------------------------------------------------------------------------
  run_fun(g)
}
##--------------------------------------------------------------------------------------------
