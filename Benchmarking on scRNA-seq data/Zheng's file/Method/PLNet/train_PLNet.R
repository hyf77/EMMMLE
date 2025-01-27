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
method_name<-'PLNet'
num_core<-10
cl <- makeCluster(num_core, outfile =paste("debug_",method_name,".txt",sep = ""))
registerDoParallel(cl)

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

mle_newton_new<-function(data_use,S_depth,core_num = 1,k_max = 10){
  mlemu<-list()
  mlesigmahat<-list()
  mu_grad<-list()
  sigma_grad<-list()
  early_stop_iternum<-1
  ##
  ##1. Perparsion
  ##-------------------------------------------------------------------------------
  dim_use<-ncol(data_use)
  sample_size<-nrow(data_use)
  ##-------------------------------------------------------------------------------
  
  ##2. Initialized estiamation of mu and sigma
  ##-------------------------------------------------------------------------------
  data_use_nor<-data_use/as.vector(S_depth)
  temp1<-colMeans(data_use/as.vector(S_depth))
  # temp1 <- replace(temp1, temp1 == 0, eps)
  log_vec1<-as.vector(log(temp1))
  temp11 <- (t(data_use_nor) %*% data_use_nor) / nrow(data_use_nor)
  # temp11 <- replace(temp11, temp11 == 0, eps^2)
  sigmahat<-t(log(temp11) - log_vec1) - log_vec1
  ##
  temp12 <-colMeans((data_use * (data_use - 1)) / as.vector(S_depth^2))
  # temp12 <- replace(temp12, temp12 == 0, eps^2)
  diag(sigmahat)<-(log(temp12) - 2 * log(temp1))
  ##judge if the diagonal of sigmahat is all positive
  neg_index<-which(diag(sigmahat)<=0)
  if(length(neg_index)>0){
    diag(sigmahat)[neg_index]<-min(diag(sigmahat)[setdiff(1:dim(sigmahat)[1],neg_index)])
  }
  diag(sigmahat)[is.na(diag(sigmahat))]<-min(diag(sigmahat)[!is.na(diag(sigmahat))])
  sigma_me1<-sigmahat
  isinfinite_mat<-ifelse(is.infinite(sigma_me1)| is.na(sigma_me1),1,0)
  isinfinite_mat[lower.tri(isinfinite_mat)]<-0
  min_vec<-rep(NA,ncol(sigma_me1))
  for(i in 1:ncol(sigma_me1)){
    min_vec[i] <- min(sigma_me1[i,-i][is.finite(c(sigma_me1[i,-i]))])
  }
  for(i in 1:(ncol(sigma_me1)-1)){
    index_select<-which(isinfinite_mat[i,]==1)
    index_select<-index_select[which(index_select>i)]
    if(length(index_select)>0){
      for(j in 1:length(index_select)){
        isinfinite_mat[i,index_select[j]]<-min(min_vec[i],min_vec[index_select[j]])
      }
    }
  }
  isinfinite_mat<-isinfinite_mat + t(isinfinite_mat)
  diag(isinfinite_mat)<-0
  
  sigma_me1<-ifelse(is.infinite(sigma_me1)| is.na(sigma_me1),0,sigma_me1)
  sigma_me1<-sigma_me1 + isinfinite_mat
  sigmahat<-sigma_me1
  temp_mu0<-colMeans(data_use/as.vector(S_depth))
  mu<-log(temp_mu0) - diag(sigmahat)/2
  temp_num_mu<-min(mu[is.finite(mu)])
  mu[is.infinite(mu)| is.na(mu)]<-temp_num_mu
  ##
  mlemu[[1]]<-mu
  mlesigmahat[[1]]<-sigmahat
  ##
  ##-------------------------------------------------------------------------------
  
  ##3. Use orthopolynom package for solving the Hermite polynomials
  ##-------------------------------------------------------------------------------
  t.root<-polynomial.roots(monic.polynomial.recurrences(hermite.h.recurrences(10, normalized=FALSE)))
  omega.root<-2^9*factorial(10)*sqrt(pi)/100/polynomial.values(hermite.h.polynomials(10,normalized=FALSE), t.root[[11]])[[10]]^2
  ##-------------------------------------------------------------------------------
  
  ##4. Calculate the gradiant, hessian and delta of mu and sigma
  ##-------------------------------------------------------------------------------
  gradiant_iter_mat<-matrix(NA,nrow = k_max,ncol = ncol(sigma_me1))
  for (k in 2:(k_max + 1)){
    ##re-initial
    sigmahat<-mlesigmahat[[k-1]]
    mu<-mlemu[[k-1]]
    ##
    sigma_diag_max<-max(diag(sigmahat))
    z.root_mat<-t(t(matrix(t.root[[11]]*sqrt(2),ncol = 1) %*% matrix(diag(sqrt(sigmahat)),nrow = 1)) + mu)
    m.root_mat<-matrix( omega.root*exp(t.root[[11]]^2)*sqrt(2), ncol = 1) %*% matrix(diag(sqrt(sigmahat)),nrow = 1)
    ##
    share12_mat<-t(t(z.root_mat) - mu)
    share11_mat<-exp(t((-1) * (((t(share12_mat)))^2)/(2 * as.vector(diag(sigmahat)))))
    ##
    root1_array<-array(0,dim = c(10,dim(data_use)))
    for(root_index in 1:10){
      root1_array[root_index,,]<-t(t(data_use) * as.vector(z.root_mat[root_index,])) - matrix(S_depth, ncol = 1) %*% matrix(exp(z.root_mat[root_index,]),nrow = 1)
    }
    max_mat<-apply(root1_array,MARGIN = c(2,3),max)
    share1_array<-array(0,dim = c(10,dim(data_use)))
    for(root_index in 1:10){
      share1_array[root_index,,]<-exp(root1_array[root_index,,] - max_mat)
    }
    ##
    sigma_diag<-diag(sigmahat)
    gz.up.mu_array<-sapply(X = 1:sample_size,FUN = function(sample_index){
      return(t(t(share1_array[,sample_index,] * share12_mat) / as.vector(sigma_diag)) * share11_mat)},simplify = "array")
    
    array_temp<-sapply(X = 1:sample_size,FUN = function(sample_index){
      return(t(t(share1_array[,sample_index,] * (share12_mat)^2) / as.vector(sigma_diag^2) ))},simplify = "array")
    
    gz.up.mu.2_array<-sapply(X = 1:sample_size,FUN = function(sample_index){
      return(t(t(share1_array[,sample_index,] * (share12_mat)^2) / as.vector(sigma_diag^2)) * share11_mat)},simplify = "array")
    
    gz.up.sigma_array<-sapply(X = 1:sample_size,FUN = function(sample_index){
      return(gz.up.mu.2_array[,,sample_index]/2)},simplify = "array")
    
    mid_mat1<-t(t(share12_mat^4) /(4 * (as.vector(sigma_diag))^4)) - t(t(share12_mat^2) /((as.vector(sigma_diag))^3))
    gz.up.sigma.2_array<-sapply(X = 1:sample_size,FUN = function(sample_index){
      return(share1_array[,sample_index,] *
               (mid_mat1)*
               share11_mat)},simplify = "array")
    
    gz.down_array<-sapply(X = 1:sample_size,FUN = function(sample_index){
      return(share1_array[,sample_index,] * share11_mat)},simplify = "array")
    
    mid_mat2<-t(t(share12_mat^3) /(2 * (as.vector(sigma_diag))^3)) - t(t(share12_mat) /((as.vector(sigma_diag))^2))
    gz.up.inter_array<-sapply(X = 1:sample_size,FUN = function(sample_index){
      return(share1_array[,sample_index,] *
               (mid_mat2) *
               share11_mat)},simplify = "array")
    ##
    
    mat_temp1<-t(sapply(X = 1:sample_size,FUN = function(sample_index){colSums(gz.down_array[,,sample_index] * m.root_mat)
    }))
    
    share2_mat<-t(sapply(X = 1:sample_size,FUN = function(sample_index){
      colSums(gz.up.mu_array[,,sample_index] * m.root_mat)/mat_temp1[sample_index,]
    }))
    share3_mat<-t(sapply(X = 1:sample_size,FUN = function(sample_index){
      colSums(gz.up.sigma_array[,,sample_index] * m.root_mat)/mat_temp1[sample_index,]
    }))
    share4_mat<-t(sapply(X = 1:sample_size,FUN = function(sample_index){
      colSums(gz.up.sigma.2_array[,,sample_index] * m.root_mat)/mat_temp1[sample_index,]
    }))
    share5_mat<-t(sapply(X = 1:sample_size,FUN = function(sample_index){
      colSums(gz.up.inter_array[,,sample_index] * m.root_mat)/mat_temp1[sample_index,]
    }))
    share6_mat<-t(sapply(X = 1:sample_size,FUN = function(sample_index){
      colSums(gz.up.mu.2_array[,,sample_index] * m.root_mat)/mat_temp1[sample_index,]
    }))
    
    share2_mat<-ifelse(is.na(share2_mat),0,share2_mat)
    share3_mat<-ifelse(is.na(share3_mat),0,share3_mat)
    share4_mat<-ifelse(is.na(share4_mat),0,share4_mat)
    share5_mat<-ifelse(is.na(share5_mat),0,share5_mat)
    share6_mat<-ifelse(is.na(share6_mat),0,share6_mat)
    ##
    gradiant.mu<-colMeans(share2_mat)
    gradiant.sigma<-colMeans(share3_mat)-sigma_diag^(-1)/2
    Hessian.mu<-colMeans(share6_mat)- sigma_diag^(-1) - colMeans((share2_mat)^2)
    Hessian.sigma<- colMeans(share4_mat) - colMeans((share3_mat)^2) + (0.5) * (sigma_diag^(-2))
    Hessian.int<-colMeans(share5_mat) - colMeans(share2_mat * share3_mat)
    ##
    mu_grad[[k-1]]<-gradiant.mu
    sigma_grad[[k-1]]<-gradiant.sigma
    ##
    delta.mu<-c()
    delta.sigma<-c()
    for(dim_index in 1:dim_use){
      Hessian.matrix<-matrix(c(Hessian.mu[dim_index],Hessian.int[dim_index],Hessian.int[dim_index],Hessian.sigma[dim_index]),2,2)
      gradiant<-matrix(c(gradiant.mu[dim_index],gradiant.sigma[dim_index]),2,1)
      gradiant_iter_mat[k-1,dim_index]<-gradiant[1,1]
      # ##
      Hessian.matrix_solve<-solve(Hessian.matrix)
      mat_temp<-Hessian.matrix_solve%*%gradiant
      delta.mu[dim_index]<-mat_temp[1]
      delta.sigma[dim_index]<-mat_temp[2]
      ##
      if (sigmahat[dim_index,dim_index]-delta.sigma[dim_index]<=0){
        delta.mu[dim_index]<-gradiant.mu[dim_index]/Hessian.mu[dim_index]
        delta.sigma[dim_index]<-gradiant.sigma[dim_index]/Hessian.sigma[dim_index]
      }
      if (abs(delta.sigma[dim_index]) > sigma_diag_max){
        delta.mu[dim_index]<-0
        delta.sigma[dim_index]<-0
      }
    }
    
    mlemu1<-mu-delta.mu
    ##
    sigmahat<-mlesigmahat[[k-1]]
    ##judge if the diagonal of sigmahat is all positive
    neg_index<-which(diag(sigmahat)<=0)
    if(length(neg_index)>0){
      diag(sigmahat)[neg_index]<-min(diag(sigmahat)[setdiff(1:dim(sigmahat)[1],neg_index)])
    }
    ##
    base_diag<-delta.sigma
    reduce_diag<-diag(sigmahat) - base_diag
    adjust_index<-which(reduce_diag<=0)
    if(length(adjust_index)>0){
      #adjust the step size
      reduce_mat<-diag(delta.sigma)
      diag(reduce_mat)[adjust_index]<-0
      mlesigmahat10<-sigmahat - reduce_mat
      diag(mlesigmahat10)[adjust_index]<-min(diag(mlesigmahat10)[-adjust_index])
      mlesigmahat1<-mlesigmahat10
    }else{
      mlesigmahat1<-sigmahat-diag(delta.sigma)
    }
    ##
    mlesigmahat1_diag<-diag(mlesigmahat1)
    ##
    ##
    mlemu[[k]]<-mlemu1
    mlesigmahat[[k]]<-mlesigmahat1
  }
  if_convergence_feature<-ifelse(abs(sigma_grad[[k_max]])<1e-3,TRUE,FALSE)
  ##
  update_iter_mat<-matrix(NA,nrow = nrow(gradiant_iter_mat),ncol = ncol(gradiant_iter_mat))
  for(dim_index in 1:ncol(gradiant_iter_mat)){
    vec_use<-abs(gradiant_iter_mat[,dim_index])
    aaa_vec<-which(diff(vec_use)>0)
    if(length(aaa_vec)>0){
      update_iter_mat[min(min(aaa_vec) + 1,k_max) :k_max,dim_index]<-FALSE
      update_iter_mat[1 :min(aaa_vec),dim_index]<-TRUE
      
    }else{
      update_iter_mat[,dim_index]<-TRUE
    }
  }
  ##--------------------------------------------------------------------------
  gradiant_iter_mat<-matrix(NA,nrow = k_max,ncol = ncol(sigma_me1))
  if_update_feature<-rep(TRUE,ncol(sigma_me1))
  for (k in 2:(k_max + 1)){
    ##re-initial
    sigmahat<-mlesigmahat[[k-1]]
    mu<-mlemu[[k-1]]
    ##
    sigma_diag_max<-max(diag(sigmahat))
    z.root_mat<-t(t(matrix(t.root[[11]]*sqrt(2),ncol = 1) %*% matrix(diag(sqrt(sigmahat)),nrow = 1)) + mu)
    m.root_mat<-matrix( omega.root*exp(t.root[[11]]^2)*sqrt(2), ncol = 1) %*% matrix(diag(sqrt(sigmahat)),nrow = 1)
    ##
    share12_mat<-t(t(z.root_mat) - mu)
    share11_mat<-exp(t((-1) * (((t(share12_mat)))^2)/(2 * as.vector(diag(sigmahat)))))
    ##
    root1_array<-array(0,dim = c(10,dim(data_use)))
    for(root_index in 1:10){
      root1_array[root_index,,]<-t(t(data_use) * as.vector(z.root_mat[root_index,])) - matrix(S_depth, ncol = 1) %*% matrix(exp(z.root_mat[root_index,]),nrow = 1)
    }
    max_mat<-apply(root1_array,MARGIN = c(2,3),max)
    share1_array<-array(0,dim = c(10,dim(data_use)))
    for(root_index in 1:10){
      share1_array[root_index,,]<-exp(root1_array[root_index,,] - max_mat)
    }
    ##
    sigma_diag<-diag(sigmahat)
    gz.up.mu_array<-sapply(X = 1:sample_size,FUN = function(sample_index){
      return(t(t(share1_array[,sample_index,] * share12_mat) / as.vector(sigma_diag)) * share11_mat)},simplify = "array")
    
    array_temp<-sapply(X = 1:sample_size,FUN = function(sample_index){
      return(t(t(share1_array[,sample_index,] * (share12_mat)^2) / as.vector(sigma_diag^2) ))},simplify = "array")
    
    gz.up.mu.2_array<-sapply(X = 1:sample_size,FUN = function(sample_index){
      return(t(t(share1_array[,sample_index,] * (share12_mat)^2) / as.vector(sigma_diag^2)) * share11_mat)},simplify = "array")
    
    gz.up.sigma_array<-sapply(X = 1:sample_size,FUN = function(sample_index){
      return(gz.up.mu.2_array[,,sample_index]/2)},simplify = "array")
    
    mid_mat1<-t(t(share12_mat^4) /(4 * (as.vector(sigma_diag))^4)) - t(t(share12_mat^2) /((as.vector(sigma_diag))^3))
    gz.up.sigma.2_array<-sapply(X = 1:sample_size,FUN = function(sample_index){
      return(share1_array[,sample_index,] *
               (mid_mat1)*
               share11_mat)},simplify = "array")
    
    gz.down_array<-sapply(X = 1:sample_size,FUN = function(sample_index){
      return(share1_array[,sample_index,] * share11_mat)},simplify = "array")
    
    mid_mat2<-t(t(share12_mat^3) /(2 * (as.vector(sigma_diag))^3)) - t(t(share12_mat) /((as.vector(sigma_diag))^2))
    gz.up.inter_array<-sapply(X = 1:sample_size,FUN = function(sample_index){
      return(share1_array[,sample_index,] *
               (mid_mat2) *
               share11_mat)},simplify = "array")
    ##
    
    mat_temp1<-t(sapply(X = 1:sample_size,FUN = function(sample_index){colSums(gz.down_array[,,sample_index] * m.root_mat)
    }))
    
    share2_mat<-t(sapply(X = 1:sample_size,FUN = function(sample_index){
      colSums(gz.up.mu_array[,,sample_index] * m.root_mat)/mat_temp1[sample_index,]
    }))
    share3_mat<-t(sapply(X = 1:sample_size,FUN = function(sample_index){
      colSums(gz.up.sigma_array[,,sample_index] * m.root_mat)/mat_temp1[sample_index,]
    }))
    share4_mat<-t(sapply(X = 1:sample_size,FUN = function(sample_index){
      colSums(gz.up.sigma.2_array[,,sample_index] * m.root_mat)/mat_temp1[sample_index,]
    }))
    share5_mat<-t(sapply(X = 1:sample_size,FUN = function(sample_index){
      colSums(gz.up.inter_array[,,sample_index] * m.root_mat)/mat_temp1[sample_index,]
    }))
    share6_mat<-t(sapply(X = 1:sample_size,FUN = function(sample_index){
      colSums(gz.up.mu.2_array[,,sample_index] * m.root_mat)/mat_temp1[sample_index,]
    }))
    
    share2_mat<-ifelse(is.na(share2_mat),0,share2_mat)
    share3_mat<-ifelse(is.na(share3_mat),0,share3_mat)
    share4_mat<-ifelse(is.na(share4_mat),0,share4_mat)
    share5_mat<-ifelse(is.na(share5_mat),0,share5_mat)
    share6_mat<-ifelse(is.na(share6_mat),0,share6_mat)
    ##
    gradiant.mu<-colMeans(share2_mat)
    gradiant.sigma<-colMeans(share3_mat)-sigma_diag^(-1)/2
    Hessian.mu<-colMeans(share6_mat)- sigma_diag^(-1) - colMeans((share2_mat)^2)
    Hessian.sigma<- colMeans(share4_mat) - colMeans((share3_mat)^2) + (0.5) * (sigma_diag^(-2))
    Hessian.int<-colMeans(share5_mat) - colMeans(share2_mat * share3_mat)
    ##
    mu_grad[[k-1]]<-gradiant.mu
    sigma_grad[[k-1]]<-gradiant.sigma
    ##
    delta.mu<-c()
    delta.sigma<-c()
    for(dim_index in 1:dim_use){
      Hessian.matrix<-matrix(c(Hessian.mu[dim_index],Hessian.int[dim_index],Hessian.int[dim_index],Hessian.sigma[dim_index]),2,2)
      gradiant<-matrix(c(gradiant.mu[dim_index],gradiant.sigma[dim_index]),2,1)
      gradiant_iter_mat[k-1,dim_index]<-gradiant[1,1]
      
      # ##
      Hessian.matrix_solve<-solve(Hessian.matrix)
      mat_temp<-Hessian.matrix_solve%*%gradiant
      delta.mu[dim_index]<-mat_temp[1]
      delta.sigma[dim_index]<-mat_temp[2]
      ##
      if((if_convergence_feature[dim_index] == FALSE)){
        
        if((update_iter_mat[(k-1),dim_index] == FALSE) | (if_update_feature[dim_index] == FALSE)){
          if_update_feature[dim_index]<-FALSE
          delta.mu[dim_index]<-0
          delta.sigma[dim_index]<-0
        }else{
          if (sigmahat[dim_index,dim_index]-delta.sigma[dim_index]<=0){
            delta.mu[dim_index]<-gradiant.mu[dim_index]/Hessian.mu[dim_index]
            delta.sigma[dim_index]<-gradiant.sigma[dim_index]/Hessian.sigma[dim_index]
          }
          if (abs(delta.sigma[dim_index]) > sigma_diag_max){
            delta.mu[dim_index]<-0
            delta.sigma[dim_index]<-0
          }
        }
      }else{
        if(if_update_feature[dim_index] == TRUE){
          if (sigmahat[dim_index,dim_index]-delta.sigma[dim_index]<=0){
            delta.mu[dim_index]<-gradiant.mu[dim_index]/Hessian.mu[dim_index]
            delta.sigma[dim_index]<-gradiant.sigma[dim_index]/Hessian.sigma[dim_index]
            
          }
          if (abs(delta.sigma[dim_index]) > sigma_diag_max){
            delta.mu[dim_index]<-0
            delta.sigma[dim_index]<-0
          }
        }
      }
      ##
    }
    #####
    
    mlemu1<-mu-delta.mu
    ##
    sigmahat<-mlesigmahat[[k-1]]
    ##judge if the diagonal of sigmahat is all positive
    neg_index<-which(diag(sigmahat)<=0)
    if(length(neg_index)>0){
      diag(sigmahat)[neg_index]<-min(diag(sigmahat)[setdiff(1:dim(sigmahat)[1],neg_index)])
    }
    ##
    base_diag<-delta.sigma
    reduce_diag<-diag(sigmahat) - base_diag
    adjust_index<-which(reduce_diag<=0)
    if(length(adjust_index)>0){
      #adjust the step size
      reduce_mat<-diag(delta.sigma)
      diag(reduce_mat)[adjust_index]<-0
      mlesigmahat10<-sigmahat - reduce_mat
      diag(mlesigmahat10)[adjust_index]<-min(diag(mlesigmahat10)[-adjust_index])
      mlesigmahat1<-mlesigmahat10
    }else{
      mlesigmahat1<-sigmahat-diag(delta.sigma)
    }
    ##
    mlesigmahat1_diag<-diag(mlesigmahat1)
    
    if(k == (k_max+1)){
      integrate_list<-PLNet::integrated_fun(data_use = data_use, S_depth = S_depth, mlemu1 = mlemu1, mlesigmahat1 = mlesigmahat1,
                                            t_root_vec = t.root[[11]], omega_root_vec = omega.root,
                                            core_num = core_num)
      gradiant.int<-integrate_list$gradiant_int
      Hessian.int<-integrate_list$Hessian_int
      ##
      gradiant.int<-ifelse(is.na(gradiant.int),0,gradiant.int)
      Hessian.int<-ifelse(is.na(Hessian.int),1,Hessian.int)
      ##
      mlesigmahat1<-basic_fun( data_use = data_use,  S_depth = S_depth,  mlemu1 = mlemu1,  mlesigmahat1 = mlesigmahat1,
                               gradiant_int = gradiant.int,  Hessian_int = Hessian.int,
                               core_num = core_num)
      ##
      mlesigmahat1[lower.tri(mlesigmahat1)]<-0
      mlesigmahat1<-mlesigmahat1 + t(mlesigmahat1)
      diag(mlesigmahat1)<-mlesigmahat1_diag
    }
    ##
    mlemu[[k]]<-mlemu1
    mlesigmahat[[k]]<-mlesigmahat1
    
  }
  
  return(list(mlemu = mlemu,
              mlesigmahat = mlesigmahat))
  
}

##--------------------------------------------------------------------------------------------
repres_glasso<-foreach (
  g = 1:celltypes_num,
  .inorder = TRUE,
  .export = ls(.GlobalEnv),
  .packages = c('MASS', 'glasso','PLNmodels','doParallel','CVXR','igraph',"PLNetNew","orthopolynom","Rcpp","RcppArmadillo")
)%dopar%{
  ##--------------------------------------------------------------------------------------------
  g2<-g
  zero_mat<-diag(p_GRN)
  zero_mat[1:p_TF,]<-1
  zero_mat[,1:p_TF]<-1
  zero_mat_use<-exact_zeroindex(support = zero_mat)
  n_lambda<-100
  BIC_mat_shrink<-c()
  ##
  print(paste("group(",g2,") task of PLNet."))
  time_PLNet_start<-Sys.time()
  ##load the data
  locator_base<-EMPLN_init$celltypes_label
  data_obs<-as.matrix(EMPLN_init$obs_mat)
  exp_count_GRN<-data_obs[,1:p_GRN]
  Y<-exp_count_GRN[locator_base==g2,]
  ##Estimate the library size by total sum scaling.
  S_depth<-EMPLN_init$S_depth[locator_base==g2]
  k_max<-10 ## The maximal iter for MLE calculation by newton-raphson algorithm
  ##MLE for convariance matrix
  cov_input<-mle_newton_new(data_use = Y,
                            S_depth = S_depth,
                            core_num = 1,
                            k_max = k_max)
  save(cov_input,file = paste(result_save_path,"/cov_input_g",g2,"_4group.Rdata",sep = ""))
  ##
  print(paste("group(",g2,") task of P."))
  k_index<-k_max
  ##Use dtrace-loss for precision matrix's estimation
  PLNet_res<-PLNet_main(obs_mat = Y,
                        Sd_est = S_depth,
                        n_lambda = n_lambda,
                        penalize.diagonal = FALSE,
                        cov_input = cov_input$mlesigmahat[[1 + k_index]],
                        weight_mat = NULL,zero_mat = NULL,
                        core_num = 1)
  lambda_vec<-PLNet_res$lambda_vec
  Theta_mat_array_hat_all<-array(NA,dim = c(p_GRN,p_GRN,n_lambda))
  for (l in 1:length(lambda_vec)){
    Theta_mat_array_hat_all[,,l]<-as.matrix(PLNet_res$Omega_est[[l]])
    BIC_mat_shrink <- c(BIC_mat_shrink,PLNet_res$BIC_vec[l])
  }
  ##
  save(Theta_mat_array_hat_all,file = paste(result_save_path,"/Theta_mat_array_hat_all_g",g2,"_4group.Rdata",sep = ""))
  save(BIC_mat_shrink,file = paste(result_save_path,"/BIC_mat_shrink_g",g2,"_4group.Rdata",sep = ""))
  save(lambda_vec,file = paste(result_save_path,"/lambda_vec_g",g2,"_4group.Rdata",sep = ""))
}
##


