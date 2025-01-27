##1. Package required 
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
library(EMMMLE)

time_whole_start<-Sys.time()
##--------------------------------------------------------------------------------------------

##2. Define some function
##--------------------------------------------------------------------------------------------
#tr
tr<-function(mat_1){
  return(sum(diag(as.matrix(mat_1))))
}

#log of det of matrix
log_det<-function(matrix_use){
  eigen_res<-eigen((matrix_use + t(matrix_use))/2)
  eig_value<-eigen_res$values
  eig_value<-ifelse(eig_value<=0,1e-8,eig_value)
  return(sum(log(eig_value)))
}

#BIC_glasso
BIC_glasso<-function(data_1,pre_mat_array){
  #data 1: n*p
  #pre_mat_array: L*p*p
  sample_size<-nrow(data_1)
  cov_mat<-cov(data_1)
  BIC_vec<-c()
  for(l in 1:dim(pre_mat_array)[1]){
    BIC_vec<-c(BIC_vec,sample_size*(-1*log_det(pre_mat_array[l,,])+tr(cov_mat%*%pre_mat_array[l,,]))+
                 log(sample_size)*sum(ifelse(pre_mat_array[l,,]!=0,1,0))/2)  
  }
  return(BIC_vec)
}

##prior tran
prior_tran<-function(p_GRN,group_num,gene_GRN_name,prior_nw_list = NULL, prior_TF = NULL,weight_prior_nw = 0.3,penalize_diagonal = TRUE){
  ################
  ##weighted network
  prior_nw_fill = list()
  if(!is.null(prior_nw_list)){
    # prior_nw_fill = prior_nw_list
    prior_nw_fill = lapply(X = 1:length(prior_nw_list),FUN = function(g){return(prior_nw_list[[g]][1:p_GRN,1:p_GRN])})
    for(g in 1:group_num){
      if(is.null(prior_nw_fill[[g]])){
        prior_nw_fill[[g]] <- matrix(0,nrow = p_GRN, ncol = p_GRN)
      }
    }
  }else{
    for(g in 1:group_num){
      prior_nw_fill[[g]] <- matrix(0,nrow = p_GRN, ncol = p_GRN)
    }
  }
  
  weight_mat_list<-list()
  for(g in 1:group_num){
    weight_mat0<-prior_nw_fill[[g]] * weight_prior_nw + (1 - prior_nw_fill[[g]])
    if(penalize_diagonal == TRUE){
      diag(weight_mat0)<-1
    }else{
      diag(weight_mat0)<-0
    }
    weight_mat_list[[g]]<-weight_mat0
    ##
    rm(weight_mat0);gc()
  }
  rm(prior_nw_fill)
  
  ##zero mat
  prior_adjoint_mat = matrix(1,ncol = p_GRN,nrow = p_GRN)
  if(!is.null(prior_TF)){
    index <- which(gene_GRN_name %in% prior_TF)
    index = setdiff(1:p_GRN, index) 
    prior_adjoint_mat[index,index] = 0
  }
  
  zero_mat<-NULL
  zero_num<-length(which(prior_adjoint_mat[upper.tri(prior_adjoint_mat)] == 0))
  if(zero_num > 0){
    zero_mat<-matrix(NA,nrow = 0,ncol = 2)
    for(i in 1:(p_GRN - 1)){
      for(j in (i+1):p_GRN){
        if(prior_adjoint_mat[i,j] == 0){
          zero_mat<-rbind(zero_mat,c(i,j))
        }
      }
    }
  }
  
  return(list(weight_mat_list = weight_mat_list,
              zero_mat = zero_mat))
}

##evaluator_2
evaluator_2<-function(input_array_array,
                      Theta_array_true,
                      lambda_vec1,
                      cluster_predict_true_vec
){
  #####
  ##input_array_array: dim_obs_data,dim_obs_data,group_num,n_lambda
  ##Theta_array_true: group_num * p * p
  ## performance criterion: AUC AUPR
  #####
  group_num<-length(cluster_predict_true_vec)
  res_mat<-matrix(NA,nrow = group_num,ncol =2)
  
  for(g in 1:group_num){
    
    Pre_mat<-Theta_array_true[cluster_predict_true_vec[g],,]
    
    ##Calculate the AUC and AUPR
    ROCXY<-matrix(0,nrow = 2,ncol = length(lambda_vec1))
    PRXY<-matrix(0,nrow = 2,ncol = length(lambda_vec1))
    
    for (i in 1:length(lambda_vec1)) {
      table1<-matrix(0,nrow = 2,ncol = 2)
      bb<-as.matrix(input_array_array[,,g,i])
      table1[1,1]<-sum((bb[upper.tri(bb)]==0) & (Pre_mat[upper.tri(Pre_mat)]==0))
      table1[1,2]<-sum((bb[upper.tri(bb)]==0) & (Pre_mat[upper.tri(Pre_mat)]!=0))
      table1[2,1]<-sum((bb[upper.tri(bb)]!=0) & (Pre_mat[upper.tri(Pre_mat)]==0))
      table1[2,2]<-sum((bb[upper.tri(bb)]!=0) & (Pre_mat[upper.tri(Pre_mat)]!=0))
      ROCXY[,i]<-table1[2,]/colSums(table1)
      PRXY[1,i]<-table1[2,2]/sum(table1[,2])
      PRXY[2,i]<-ifelse(is.nan((table1[,2]/rowSums(table1))[2]),1,(table1[,2]/rowSums(table1))[2])
      rm(bb)
    }
    
    ROCXY_2<-cbind(matrix(c(1,1),nrow = 2),ROCXY)
    ROCXY_1<-as.data.frame(t(rbind(sort(ROCXY_2[1,]),ROCXY_2[2,order(ROCXY_2[1,])])))
    names(ROCXY_1)<-c("FPR","TPR")
    ROCXY_1<-rbind(c(0,0),ROCXY_1,c(1,1))
    
    PRXY_2<-cbind(matrix(1,colSums(table1)[2]/colSums(table1)[1],nrow = 2),PRXY)
    PRXY_1<-as.data.frame(t(rbind(sort(PRXY_2[1,]),PRXY_2[2,order(PRXY_2[1,])])))
    names(PRXY_1)<-c("Recall","Precison")
    PRXY_1<-rbind(c(0,1),PRXY_1,c(1,0))
    PRXY_1<-PRXY_1[which(rowSums(PRXY_1)!=0),]
    
    #AUC
    AUC_Dtrace<-sum(diff(ROCXY_1[,1])*(ROCXY_1[-nrow(ROCXY_1),2]+ROCXY_1[-1,2])/2)
    #AUPR
    AUPR_Dtrace<-sum(diff(PRXY_1[,1])*(PRXY_1[-nrow(PRXY_1),2]+PRXY_1[-1,2])/2)
    
    res_mat[g,1]<-AUC_Dtrace
    res_mat[g,2]<-AUPR_Dtrace
    
  }
  ##
  return(res_mat)
}

MMLE_init<-function(data,res_init){
  iter_max_t<-50
  ELBO_threshold<-1e-4
  #####################################################################################################
  data_obs <-data$obs_mat
  Theta_mat_array_TRUE <-data[["group_specific_precision"]]
  cluster_true <-data[["locator_vec"]]
  dim_obs_data<-ncol(data_obs)
  num_obs_data<-nrow(data_obs)
  group_num<-nrow(data$mu_mat)
  t.root<-polynomial.roots(monic.polynomial.recurrences(hermite.h.recurrences(20, normalized=FALSE)))
  t_root_vec<-t.root[[21]]
  omega_root_vec<-2^19*factorial(20)*sqrt(pi)/400/polynomial.values(hermite.h.polynomials(20,normalized=FALSE), t.root[[21]])[[20]]^2
  ## data initialization-----------------------------------------------------
  locator_base <- res_init$celltypes_label
  ARI_orgin <- adjustedRandIndex(cluster_true,locator_base)
  S_depth<-res_init$S_depth
  ## record the change of ARI
  ARI_iter<-c()
  ARI_iter<-c(ARI_iter,ARI_orgin)
  ##initialization for EM step--------------------------------------------------------------------------
  diff_ELBO <- 1
  t<-1
  # mlemu_init<-matrix(NA,nrow = group_num,ncol = dim_obs_data)
  # mlesigmahat_init<-list()
  if_NaN_1<-0
  ELBO_new<-(-1)*1e10
  for (g in 1:group_num) {
    eigen_cov_Y_temp<-eigen(res_init$sigma_mat_list[[g]])
    if(min(eigen_cov_Y_temp$values)<0){
      if_NaN_1<-1
      cov_Y_temp<-sigma_positive(cov_input = res_init$sigma_mat_list[[g]])
      if(det(cov_Y_temp)<0){
        if_NaN_1<-2
        while(det(cov_Y_temp)<0){
          diag(cov_Y_temp)<-diag(cov_Y_temp)+1e-4
        }
      }
      res_init$sigma_mat_list[[g]]<-cov_Y_temp
    }
  }
  time_EMPLN_start1<-Sys.time()
  init_Zg<-cond_pro_Xdata_Zg(Xdata = data_obs,pi=res_init$pi_vec,mu=res_init$mu_mat,sigma=res_init$sigma_mat_list,scale = matrix(rep(S_depth,dim_obs_data),num_obs_data,dim_obs_data),obs_data_sp=res_init$obs_mat)
  # s<-list()
  # for (g in 1:3) {
  #   s[[g]]<-data$group_specific_precision[g,,]
  # }
  # data_obs<-data$obs_mat
  #
  # init_Zg<-cond_pro_Xdata_Zg(Xdata = data_obs,pi=c(1/3,1/3,1/3),mu=data$mu_mat,sigma=s,scale = matrix(rep(data$ls_vec,dim_obs_data),num_obs_data,dim_obs_data),obs_data_sp=res_init$obs_mat)
  #
  #
  pi_iter<-init_Zg$pi.new
  post_dist_Z_iter<-init_Zg$post_dist_Z
  ELBO_new<-sum(init_Zg$elbo_res)
  if(is.infinite(ELBO_new))
  {
    ELBO_new<-(-1)*1e10
  }
  locator_base1<-apply(post_dist_Z_iter,MARGIN = 1,which.max)
  cluster_predict_true_cur1<-as.vector(apply(table(locator_base1,cluster_true),MARGIN = 1,which.max))
  table_reshape<-matrix(NA,nrow = group_num,ncol = group_num)
  for(g1 in 1:group_num){
    for(g2 in 1:group_num){
      table_reshape[g1,g2]<-length(which(locator_base1 == g1 & cluster_true == cluster_predict_true_cur1[g2]))
    }
  }
  ARI_orgin11 <- adjustedRandIndex(cluster_true,locator_base1)
  ARI_iter<-c(ARI_iter,ARI_orgin11)
  mlemu_init<-matrix(NA,nrow = group_num,ncol = dim_obs_data)
  mlesigmahat_init<-list()
  for(g in 1:group_num){
    #select_Z<-ifelse(parameter_update$post_dist_Z[,g] > 0.9, 1, 0)
    data_obs_group<-data_obs[locator_base1==g,]
    Y<-data_obs_group
    data_use<-Y
    S_depth_group<-res_init$S_depth[locator_base1==g]
    ##
    ##1. Perparsion
    ##-------------------------------------------------------------------------------
    dim_use<-ncol(data_use)
    sample_size<-nrow(data_use)
    ##-------------------------------------------------------------------------------
    
    ##2. Initialized estiamation of mu and sigma
    ##-------------------------------------------------------------------------------
    data_use_nor<-data_use/as.vector(S_depth_group)
    log_vec1<-as.vector(log(colMeans(data_use_nor)))
    sigmahat<-t(log((t(data_use_nor) %*% data_use_nor)/nrow(data_use_nor)) - log_vec1) - log_vec1
    ##
    diag(sigmahat)<-(log(colMeans((data_use * (data_use - 1)) / as.vector(S_depth_group^2))) - 2 * log(colMeans(data_use/as.vector(S_depth_group))))
    mu<-log(colMeans(data_use/as.vector(S_depth_group))) - diag(sigmahat)/2
    ##judge if the diagonal of sigmahat is all positive
    neg_index<-which(diag(sigmahat)<=0)
    if(length(neg_index)>0){
      diag(sigmahat)[neg_index]<-min(diag(sigmahat)[setdiff(1:dim(sigmahat)[1],neg_index)])
    }
    ##
    ##adjust
    allzero_set<-which(colSums(data_use) == 0)
    zero_plus_one_set<-setdiff(which(colSums(ifelse(data_use>=2,1,0)) == 0),allzero_set)
    choose_index<-setdiff(1:ncol(data_use),allzero_set)
    sigma_me1<-sigmahat[choose_index,choose_index]
    diag(sigma_me1)[which(choose_index %in% zero_plus_one_set)]<-0
    ##
    isinfinite_mat<-ifelse(is.infinite(sigma_me1),1,0)
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
    sigma_me1<-ifelse(is.infinite(sigma_me1),0,sigma_me1)
    sigma_me1<-sigma_me1 + isinfinite_mat
    sigmahat<-sigma_me1
    # cov_Y1<-sigmahat
    # eigen_cov_Y1<-eigen(cov_Y1)
    # if(min(eigen_cov_Y1$values)<0){
    #   cov_Y1<-sigma_positive(cov_input = cov_Y1)
    # }
    # if(det(cov_Y1)<0){
    #   while(det(cov_Y1)<0){
    #     diag(cov_Y1)<-diag(cov_Y1)+1e-4
    #   }
    #   if_NaN_1<-TRUE
    # }
    # sigmahat<-cov_Y1
    ##
    mlemu_init[g,]<-mu
    mlesigmahat_init[[g]]<-sigmahat
    rm(mu);gc()
    rm(sigmahat);gc()
  }
  updated_mu<-matrix(0,nrow = dim(data$mu_mat)[1],ncol = dim_obs_data)
  updated_sigma<-list()
  updated_nonpos_sigma<-list()
  positions <- expand.grid(row = 1:dim_obs_data, col = 1:dim_obs_data)
  position_input <- positions[positions$row > positions$col, ]
  position_input<-as.matrix(position_input)
  if_NaN_2<-FALSE
  for(g in 1:group_num){
    k_max<-10
    cond_Z_input<-post_dist_Z_iter[,g]
    #data_obs_group<-data_obs
    res_init_mu<-mlemu_init[g,]
    res_init_sigma<-mlesigmahat_init[[g]]
    mmle_newton_res<-mmle_newton_fun(data_obs,cond_Z_input,S_depth,res_init_mu,res_init_sigma,t_root_vec,omega_root_vec,k_max,20,position_input,50)
    updated_mu[g,]<-mmle_newton_res$mlemu
    updated_nonpos_sigma[[g]]<-mmle_newton_res$mlesigmahat
    updated_sigma[[g]]<-mmle_newton_res$mlesigmahat
    cov_Y1<-updated_sigma[[g]]
    eigen_cov_Y1<-eigen(cov_Y1)
    if(min(eigen_cov_Y1$values)<0){
      cov_Y1<-sigma_positive(cov_input = cov_Y1)
    }
    if(det(cov_Y1)<0){
      while(det(cov_Y1)<0){
        diag(cov_Y1)<-diag(cov_Y1)+1e-4
      }
      if_NaN_2<-TRUE
    }
    updated_sigma[[g]]<-cov_Y1
  }
  # parameter_init$mu_mat<-mlemu_init
  # parameter_init$sigma_mat_list<-mlesigmahat_init
  parameter_init<-list(pi_vec=pi_iter,mu_mat=updated_mu,sigma_mat_list=updated_sigma,nonpos_sigma=updated_nonpos_sigma)
  parameter_update<-list()
  time_cost<-c()
  if_NaN_2<-FALSE
  ARI_old<-ARI_orgin
  cond_pro_Z<-cond_pro_Xdata_Zg(Xdata = data_obs,pi=parameter_init$pi_vec,mu=parameter_init$mu_mat,sigma=parameter_init$sigma_mat_list,scale = matrix(rep(S_depth,dim_obs_data),num_obs_data,dim_obs_data),obs_data_sp=res_init$obs_mat)
  
  pi_iter<-cond_pro_Z$pi.new
  post_dist_Z_iter<-cond_pro_Z$post_dist_Z
  
  while(diff_ELBO>ELBO_threshold && t <= iter_max_t){
    ELBO_old<-ELBO_new
    ## update pi and the posterior distribution of Z
    #cond_pro_Z<-cond_pro_Xdata_Zg(Xdata = data_obs,pi=parameter_init$pi_vec,mu=parameter_init$mu_mat,sigma=parameter_init$sigma_mat_list,scale = matrix(rep(S_depth,dim_obs_data),num_obs_data,dim_obs_data))
    
    locator_base_temp<-apply(post_dist_Z_iter,MARGIN = 1,which.max)
    ARI_temp <- adjustedRandIndex(cluster_true,locator_base_temp)
    ARI_iter<-c(ARI_iter,ARI_temp)
    if(t>1 && ARI_temp<ARI_old)
    {
      print("ARI decrease")
      break
    }
    ARI_old<-ARI_temp
    parameter_update$pi_vec<-pi_iter
    parameter_update$post_dist_Z<-post_dist_Z_iter
    ## update mu and sigma
    time3_1<-Sys.time()
    ##-------------------------------------------------------------------------------
    # mlemu_init<-matrix(NA,nrow = group_num,ncol = dim_obs_data)
    # mlesigmahat_init<-list()
    # for(g in 1:group_num){
    #   #select_Z<-ifelse(parameter_update$post_dist_Z[,g] > 0.9, 1, 0)
    #   data_obs_group<-data_obs[locator_base1==g,]
    #   Y<-data_obs_group
    #   data_use<-Y
    #   S_depth_group<-res_init$S_depth[locator_base1==g]
    #   ##
    #   ##1. Perparsion
    #   dim_use<-ncol(data_use)
    #   sample_size<-nrow(data_use)
    #   ##2. Initialized estiamation of mu and sigma
    #   data_use_nor<-data_use/as.vector(S_depth_group)
    #   log_vec1<-as.vector(log(colMeans(data_use_nor)))
    #   sigmahat<-t(log((t(data_use_nor) %*% data_use_nor)/nrow(data_use_nor)) - log_vec1) - log_vec1
    #   ##
    #   diag(sigmahat)<-(log(colMeans((data_use * (data_use - 1)) / as.vector(S_depth_group^2))) - 2 * log(colMeans(data_use/as.vector(S_depth_group))))
    #   mu<-log(colMeans(data_use/as.vector(S_depth_group))) - diag(sigmahat)/2
    #   ##judge if the diagonal of sigmahat is all positive
    #   neg_index<-which(diag(sigmahat)<=0)
    #   if(length(neg_index)>0){
    #     diag(sigmahat)[neg_index]<-min(diag(sigmahat)[setdiff(1:dim(sigmahat)[1],neg_index)])
    #   }
    #   ##
    #   ##adjust
    #   allzero_set<-which(colSums(data_use) == 0)
    #   zero_plus_one_set<-setdiff(which(colSums(ifelse(data_use>=2,1,0)) == 0),allzero_set)
    #   choose_index<-setdiff(1:ncol(data_use),allzero_set)
    #   sigma_me1<-sigmahat[choose_index,choose_index]
    #   diag(sigma_me1)[which(choose_index %in% zero_plus_one_set)]<-0
    #   ##
    #   isinfinite_mat<-ifelse(is.infinite(sigma_me1),1,0)
    #   isinfinite_mat[lower.tri(isinfinite_mat)]<-0
    #   min_vec<-rep(NA,ncol(sigma_me1))
    #   for(i in 1:ncol(sigma_me1)){
    #     min_vec[i] <- min(sigma_me1[i,-i][is.finite(c(sigma_me1[i,-i]))])
    #   }
    #   for(i in 1:(ncol(sigma_me1)-1)){
    #     index_select<-which(isinfinite_mat[i,]==1)
    #     index_select<-index_select[which(index_select>i)]
    #     if(length(index_select)>0){
    #       for(j in 1:length(index_select)){
    #         isinfinite_mat[i,index_select[j]]<-min(min_vec[i],min_vec[index_select[j]])
    #       }
    #     }
    #   }
    #   isinfinite_mat<-isinfinite_mat + t(isinfinite_mat)
    #   diag(isinfinite_mat)<-0
    #   sigma_me1<-ifelse(is.infinite(sigma_me1),0,sigma_me1)
    #   sigma_me1<-sigma_me1 + isinfinite_mat
    #   sigmahat<-sigma_me1
    #   # cov_Y1<-sigmahat
    #   # eigen_cov_Y1<-eigen(cov_Y1)
    #   # if(min(eigen_cov_Y1$values)<0){
    #   #   cov_Y1<-sigma_positive(cov_input = cov_Y1)
    #   # }
    #   # if(det(cov_Y1)<0){
    #   #   while(det(cov_Y1)<0){
    #   #     diag(cov_Y1)<-diag(cov_Y1)+1e-4
    #   #   }
    #   #   if_NaN_1<-TRUE
    #   # }
    #   # sigmahat<-cov_Y1
    #   ##
    #   mlemu_init[g,]<-mu
    #   mlesigmahat_init[[g]]<-sigmahat
    #   rm(mu);gc()
    #   rm(sigmahat);gc()
    # }
    ##-------------------------------------------------------------------------------
    mlemu_init<-parameter_init$mu_mat
    mlesigmahat_init<-parameter_init$sigma_mat_list
    updated_mu<-matrix(0,nrow = dim(data$mu_mat)[1],ncol = dim_obs_data)
    updated_sigma<-list()
    updated_nonpos_sigma<-list()
    for(g in 1:group_num){
      k_max<-10
      cond_Z_input<-parameter_update$post_dist_Z[,g]
      res_init_mu<-mlemu_init[g,]
      res_init_sigma<-mlesigmahat_init[[g]]
      mmle_newton_res<-mmle_newton_fun(data_obs,cond_Z_input,S_depth,res_init_mu,res_init_sigma,t_root_vec,omega_root_vec,k_max,20,position_input,50)
      updated_mu[g,]<-mmle_newton_res$mlemu
      updated_nonpos_sigma[[g]]<-mmle_newton_res$mlesigmahat
      updated_sigma[[g]]<-mmle_newton_res$mlesigmahat
      cov_Y1<-updated_sigma[[g]]
      eigen_cov_Y1<-eigen(cov_Y1)
      if(min(eigen_cov_Y1$values)<0){
        cov_Y1<-sigma_positive(cov_input = cov_Y1)
      }
      if(det(cov_Y1)<0){
        while(det(cov_Y1)<0){
          diag(cov_Y1)<-diag(cov_Y1)+1e-4
        }
        if_NaN_2<-TRUE
      }
      updated_sigma[[g]]<-cov_Y1
    }
    time3_2<-Sys.time()
    time3<-as.double(difftime(time3_2,time3_1,units = "hours"))
    parameter_update$mu_mat<-updated_mu
    parameter_update$sigma_mat_list<-updated_sigma
    parameter_update$nonpos_sigma<-updated_nonpos_sigma
    ## calculate the ELBO
    sigma_input_list<-parameter_update$sigma_mat_list
    ELBO_update<-cond_pro_Xdata_Zg(Xdata = data_obs,pi=pi_iter,mu=parameter_update$mu_mat,
                                   sigma=sigma_input_list,scale = matrix(rep(S_depth,dim_obs_data),num_obs_data,dim_obs_data),obs_data_sp=res_init$obs_mat)
    ELBO_iter<-ELBO_update$elbo_res
    print(ELBO_iter)
    if(length(which(is.infinite(ELBO_iter)==TRUE))>0)
    {
      print("ELBO_iter has NaN")
      break
    }
    pi_iter<-ELBO_update$pi.new
    post_dist_Z_iter<-ELBO_update$post_dist_Z
    if(t>1 && sum(ELBO_iter)<ELBO_old)
    {
      print(abs(sum(ELBO_iter)-ELBO_old)/(abs(ELBO_old)))
      parameter_update<-parameter_init
      #parameter_update$sigma_mat_list<-sigma_input_old
      break
    }
    ELBO_new<-sum(ELBO_iter)
    diff_ELBO<-abs(ELBO_new-ELBO_old)/(abs(ELBO_old))
    print(diff_ELBO)
    #sigma_input_old<-sigma_input_list
    parameter_init<-parameter_update
    t<-t+1
  }
  time_EMPLN_end1<-Sys.time()
  time_cost<-c(time_cost,as.double(difftime(time_EMPLN_end1,time_EMPLN_start1,units = "hours")))
  locator_cur<-apply(parameter_update$post_dist_Z,MARGIN = 1,which.max)
  cluster_predict_true_cur<-as.vector(apply(table(locator_cur,cluster_true),MARGIN = 1,which.max))
  sigma_error<-c()
  for(g in 1:group_num){
    temp_matrix0<-parameter_update$sigma_mat_list[[g]]-solve(data$group_specific_precision[cluster_predict_true_cur[g],,])
    sigma_error<-c(sigma_error,norm(temp_matrix0, type ="F"))
  }
  mu_sigma_update<-list(mu_mat=parameter_update$mu_mat,
                        sigma_list=parameter_update$sigma_mat_list,
                        post_dist_Z=parameter_update$post_dist_Z,
                        pi_vec=parameter_update$pi_vec,
                        nonpos_sigma_list=parameter_update$nonpos_sigma,
                        cluster_predict_true_cur=cluster_predict_true_cur,
                        sigma_error=sigma_error,
                        time_cost=time_cost,
                        ARI_iter=ARI_iter,
                        if_NaN_1=if_NaN_1,
                        if_NaN_2=if_NaN_2)
  return(mu_sigma_update)
}

##prior_generator
prior_generator<-function(adjoint_mat,eta=0.3){
  edge_index<-which(adjoint_mat[upper.tri(adjoint_mat)]==1)
  edge_index_sel<-sample(edge_index,floor(length(edge_index) * eta),replace = FALSE)
  prior_adjoint<-diag(nrow(adjoint_mat))
  prior_adjoint[upper.tri(prior_adjoint)][edge_index_sel]<-1
  prior_adjoint<-prior_adjoint + t(prior_adjoint) - diag(nrow(adjoint_mat))
  ##
  return(prior_adjoint)
}

sub_diag<-function(p,diag_num){
  AAA<-matrix(0,nrow = p,ncol = p)
  for(j in 1:(p-1)){
    AAA[j,]<-(c(rep(0,j),rep(1,diag_num),rep(0,ifelse(p-j-diag_num>0,p-j-diag_num,0))))[1:p]
  }
  AAA_1<-AAA+t(AAA)+diag(p)
  return(AAA_1)
}

##data_generator_new
data_generator<-function(n,dim_use,group_num,
                         network_type=c("ER-graph","AFF-graph","PA-graph","HUB-graph","banded-graph","banded-random-graph"),
                         densy_degree = 0.1,v = 0.3,p_v = 0.5,
                         num_de = 15,ub = 1,lb = -3.5, ub_non = 1, lb_non = -1,
                         l_mu = log(10),l_sd = 0.05,sigma_scale = 1,
                         de_ratio_bet = 0.5,alpha_com = 0.01,
                         if_fix = FALSE,
                         fix_element = NULL){
  if(if_fix == FALSE){
    #Generate Group-specific adjoint matrix & precision
    group_specific_adjoint<-array(0,dim=c(group_num,dim_use,dim_use))
    group_specific_precision<-array(0,dim=c(group_num,dim_use,dim_use))
    hub_index<-NULL
    if(network_type == "HUB-graph"){
      hub_index<-sample(1:dim_use, 2 * dim_use * densy_degree, replace = FALSE)
    }
    ##adjoint matrix
    ##adjoint
    if(network_type == "ER-graph"){
      ##common_adjoint
      common_adjoint<-matrix(0,nrow = dim_use,ncol = dim_use)
      common_adjoint[upper.tri(common_adjoint)][sample(1:(dim_use * (dim_use - 1)/2),floor(dim_use * (dim_use - 1)/2 * densy_degree * alpha_com),replace = FALSE)]<-1
      common_adjoint <- common_adjoint + t(common_adjoint) + diag(dim_use)
      avil_edge<-which(common_adjoint[upper.tri(common_adjoint)]==0)
      for(r in 1:group_num){
        adjoin_mat_choose<-matrix(0,nrow = dim_use,ncol = dim_use)
        adjoin_mat_choose[upper.tri(adjoin_mat_choose)][sample(avil_edge,floor(dim_use * (dim_use - 1)/2 * densy_degree * (1-alpha_com)),replace = FALSE)]<-1
        adjoin_mat_choose<-adjoin_mat_choose+t(adjoin_mat_choose)
        
        group_specific_adjoint[r,,]<- adjoin_mat_choose + common_adjoint
      }
    }
    
    if(network_type == "AFF-graph"){
      community_label<-as.vector(matrix(rep(1:4,each = dim_use/4),ncol = 4,byrow = FALSE))
      adjoin_mat_com<-diag(dim_use)
      for(k in 1:4){
        adjoin_mat_com[which(community_label == k),which(community_label == k)]<-1
      }
      edge_index<-which(adjoin_mat_com[upper.tri(adjoin_mat_com)]==1)
      ##common_adjoint
      common_adjoint<-matrix(0,nrow = dim_use,ncol = dim_use)
      common_adjoint[upper.tri(common_adjoint)][sample(edge_index,floor(dim_use * (dim_use - 1)/2 * densy_degree * alpha_com),replace = FALSE)]<-1
      common_adjoint <- common_adjoint + t(common_adjoint) + diag(dim_use)
      avil_edge<-setdiff(edge_index,which(common_adjoint[upper.tri(common_adjoint)] == 1))
      for(r in 1:group_num){
        adjoin_mat_choose<-matrix(0,nrow = dim_use,ncol = dim_use)
        adjoin_mat_choose[upper.tri(adjoin_mat_choose)][sample(avil_edge,floor(dim_use * (dim_use - 1)/2 * densy_degree * (1 - alpha_com)),replace = FALSE)]<-1
        adjoin_mat_choose<-adjoin_mat_choose+t(adjoin_mat_choose)
        group_specific_adjoint[r,,]<-adjoin_mat_choose + common_adjoint 
      }
    }
    if(network_type == "PA-graph"){
      #calculate the densy of single network
      g1 <- sample_pa_age(dim_use, pa.exp=0.5, aging.exp=-1, aging.bin=1000)
      adjoin_mat<-as_adjacency_matrix(g1)
      adjoin_mat<-as.matrix(adjoin_mat) + t(as.matrix(adjoin_mat)) + diag(dim_use)
      ##
      num_common<-floor(densy_degree * alpha_com/as.numeric(prop.table(table(adjoin_mat[upper.tri(adjoin_mat)]))[2]))
      num_noncommon<-floor(densy_degree/as.numeric(prop.table(table(adjoin_mat[upper.tri(adjoin_mat)]))[2])) - num_common
      
      ##common_adjoint
      common_adjoint<-matrix(0,nrow = dim_use,ncol = dim_use)
      for(j in 1:num_common){
        g1 <- sample_pa_age(dim_use, pa.exp=0.5, aging.exp=-1, aging.bin=1000)
        adjoin_mat<-as_adjacency_matrix(g1)
        adjoin_mat<-as.matrix(adjoin_mat) + t(as.matrix(adjoin_mat))
        common_adjoint <- common_adjoint + adjoin_mat
      }
      common_adjoint<-ifelse(common_adjoint!=0,1,0) + diag(dim_use)
      
      for(r in 1:group_num){
        adjoin_mat_group<-matrix(0,nrow = dim_use,ncol = dim_use)
        for(j in 1:num_noncommon){
          g1 <- sample_pa_age(dim_use, pa.exp=0.5, aging.exp=-1, aging.bin=1000)
          adjoin_mat<-as_adjacency_matrix(g1)
          adjoin_mat<-as.matrix(adjoin_mat) + t(as.matrix(adjoin_mat))
          adjoin_mat_group <- adjoin_mat_group + adjoin_mat
        }
        group_specific_adjoint[r,,]<-ifelse(adjoin_mat_group + common_adjoint!=0,1,0)
      }
      
    }
    if(network_type == "HUB-graph"){
      adjoin_mat_choose<-matrix(1,nrow = dim_use,ncol = dim_use)
      adjoin_mat_choose[-hub_index,-hub_index]<-0
      diag(adjoin_mat_choose)<-1
      edge_choose<-which(adjoin_mat_choose[upper.tri(adjoin_mat_choose)]==1)
      ##common_adjoint
      common_adjoint<-matrix(0,nrow = dim_use,ncol = dim_use)
      common_adjoint[upper.tri(common_adjoint)][sample(edge_choose,floor(dim_use * (dim_use - 1)/2 * densy_degree * (alpha_com)),replace = FALSE)]<-1
      common_adjoint <- common_adjoint + t(common_adjoint) + diag(dim_use)
      avil_edge<-setdiff(edge_choose,which(common_adjoint[upper.tri(common_adjoint)] == 1))
      
      for(r in 1:group_num){
        edge_choose1<-sample(avil_edge,floor(dim_use * (dim_use - 1)/2 * densy_degree * (1-alpha_com)),replace = FALSE)
        adjoin_mat<-matrix(0,nrow = dim_use,ncol = dim_use)
        adjoin_mat[upper.tri(adjoin_mat)][edge_choose1]<-1
        adjoin_mat<-adjoin_mat + t(adjoin_mat)
        group_specific_adjoint[r,,]<-adjoin_mat + common_adjoint
      }
    }
    if(network_type == "banded-graph"){
      ##common_adjoint
      pre_matrix<-sub_diag(dim_use,2)
      common_adjoint<-matrix(0,nrow = dim_use,ncol = dim_use)
      pre_edg<-length(which(pre_matrix[upper.tri(pre_matrix)]==1))
      common_adjoint[upper.tri(common_adjoint)][sample(1:pre_edg,floor(pre_edg *alpha_com),replace = FALSE)]<-1
      common_adjoint <- common_adjoint + t(common_adjoint) + diag(dim_use)
      for(r in 1:group_num){
        group_specific_adjoint[r,,]<- sub_diag(dim_use,2) 
      }
    }
    if(network_type == "banded-random-graph"){
      ##common_adjoint
      pre_matrix<-sub_diag(dim_use,floor(dim_use*0.06))
      edge_choose<-which(pre_matrix[upper.tri(pre_matrix)]==1)
      common_adjoint<-matrix(0,nrow = dim_use,ncol = dim_use)
      common_adjoint[upper.tri(common_adjoint)][sample(edge_choose,floor(dim_use * (dim_use - 1)/2 * densy_degree * (alpha_com)),replace = FALSE)]<-1
      common_adjoint <- common_adjoint + t(common_adjoint) + diag(dim_use)
      avil_edge<-setdiff(edge_choose,which(common_adjoint[upper.tri(common_adjoint)] == 1))
      for(r in 1:group_num){
        adjoin_mat_choose<-matrix(0,nrow = dim_use,ncol = dim_use)
        adjoin_mat_choose[upper.tri(adjoin_mat_choose)][sample(avil_edge,floor(dim_use * (dim_use - 1)/2 * densy_degree * (1-alpha_com)),replace = FALSE)]<-1
        adjoin_mat_choose<-adjoin_mat_choose+t(adjoin_mat_choose)
        
        group_specific_adjoint[r,,]<- adjoin_mat_choose + common_adjoint
      }
    }
    #precision
    ##1. common part
    # rbinom_mat<-diag(dim_use)
    # rbinom_mat[upper.tri(rbinom_mat)]<-(rbinom(dim_use * (dim_use-1)/2,1,1-p_v)-0.5)*2*v
    # rbinom_mat<-rbinom_mat + t(rbinom_mat) - diag(dim_use)
    # pre_common <- common_adjoint * rbinom_mat
    # common_index<-which(pre_common[upper.tri(pre_common)]!=0)
    for(r in 1:group_num){
      # diag_value<-1
      # group_specific_precision0<-matrix(0,nrow = dim_use,ncol = dim_use)
      # edge_allow<-setdiff(which(group_specific_adjoint[r,,][upper.tri(group_specific_adjoint[r,,])]==1),common_index)
      # group_specific_precision0[upper.tri(group_specific_precision0)][edge_allow]<-(rbinom(length(edge_allow),1,1 - p_v)-0.5)*2*v
      # group_specific_precision0<-group_specific_precision0 + t(group_specific_precision0)
      # group_specific_precision0<-group_specific_precision0 + pre_common
      # con_pos<-TRUE
      # while (con_pos) {
      #   group_specific_precision1<-group_specific_precision0
      #   diag(group_specific_precision1)<-diag_value
      #   
      #   eigen_value<-eigen(group_specific_precision1)$values
      #   if(min(eigen_value)>0){
      #     con_pos<-FALSE
      #     group_specific_precision[r,,]<-group_specific_precision1 * sigma_scale
      #   }else{
      #     diag_value<-diag_value * 1.2
      #   }
      # }
      #############################
      diag_value<-1
      group_specific_precision0<-diag(dim_use)
      group_specific_precision0[upper.tri(group_specific_precision0)][which(group_specific_adjoint[r,,][upper.tri(group_specific_adjoint[r,,])]==1)]<-(rbinom(length(which(group_specific_adjoint[r,,][upper.tri(group_specific_adjoint[r,,])]==1))
                                                                                                                                                              ,1,1 - p_v)-0.5)*2*v
      group_specific_precision0<-group_specific_precision0 + t(group_specific_precision0) - diag(dim_use)
      eigen_value<-eigen(group_specific_precision0)$values
      if(min(eigen_value)<0){
        group_specific_precision0<-group_specific_precision0 + diag(dim_use) * (abs(min(eigen_value))+1e-1)
      }else{
        group_specific_precision0<-group_specific_precision0 + diag(dim_use) * (1e-1)
      }
      
      group_specific_precision[r,,]<-group_specific_precision0 * sigma_scale
      # ##
      # con_pos<-TRUE
      # while (con_pos) {
      #   group_specific_precision1<-group_specific_precision0
      #   diag(group_specific_precision1)<-diag_value
      # 
      #   eigen_value<-eigen(group_specific_precision1)$values
      #   if(min(eigen_value)>0){
      #     con_pos<-FALSE
      #     group_specific_precision[r,,]<-group_specific_precision1 * sigma_scale
      #   }else{
      #     diag_value<-diag_value * 1.2
      #   }
      # }
      # ##
      
    }
    
    #Generate the sc-RNA data
    #1. mu_mat
    differential_gene_index<-sample(1:dim_use,num_de,replace = FALSE)
    nondifferential_gene_index<-setdiff(1:dim_use,differential_gene_index)
    mu_mat<-matrix(0,nrow = group_num,ncol = dim_use)
    mu_nondifferent<-runif(length(nondifferential_gene_index),min = lb_non,max = ub_non)
    for(r in 1:group_num){
      mu_mat[r,nondifferential_gene_index]<-mu_nondifferent
    }
    for(j in 1:length(differential_gene_index)){
      mu_mat[,differential_gene_index[j]]<-sample(c(lb,ub,(lb+ub)/2),group_num,replace = TRUE)
    }
    # group_num_true = group_num
    # # group_num = max(5,group_num_true)
    # rep_con<-TRUE
    # while (rep_con) {
    #   AAA<-matrix(NA,nrow = group_num,ncol = num_de)
    #   for(j in 1:num_de){
    #     aaa_vec<-sample(c(1:3),group_num,replace = TRUE)
    #     while(var(aaa_vec) == 0){
    #       aaa_vec[sample(1:group_num,1)]<-sample(c(1:3),1,replace = TRUE)
    #     }
    #     
    #     AAA[,j]<-aaa_vec
    #   }
    #   con_1<-TRUE
    #   count<-0
    #   while(con_1 & (count<20)){
    #     count<-count+1
    #     con_vec<-c()
    #     for(g in 1:group_num){
    #       if(min(rowSums(sign(abs(t(t(AAA[-g,]) - as.vector(AAA[g,]))))))<(de_ratio_bet * num_de)){
    #         choose_index<-sample(which(sign(abs(t(t(AAA[-g,]) - as.vector(AAA[g,]))))[which.min((rowSums(sign(abs(t(t(AAA[-g,]) - as.vector(AAA[g,]))))))),]==0),num_de * de_ratio_bet - min(rowSums(sign(abs(t(t(AAA[-g,]) - as.vector(AAA[g,])))))))
    #         for(i in 1:length(choose_index)){
    #           AAA[-g,][which.min((rowSums(sign(abs(t(t(AAA[-g,]) - as.vector(AAA[g,]))))))),choose_index[i]]<-sample(setdiff(c(1:3),AAA[-g,][which.min((rowSums(sign(abs(t(t(AAA[-g,]) - as.vector(AAA[g,]))))))),choose_index[i]]),1)
    #         }
    #         con_vec<-c(con_vec,1)
    #       }else{
    #         con_vec<-c(con_vec,0)
    #       }
    #     }
    #     if(length(which(con_vec==1))==0){
    #       con_1<-FALSE
    #     }
    #   }
    #   if(con_1 == TRUE){
    #     rep_con<-TRUE
    #   }else{
    #     rep_con<-FALSE
    #   }
    # }
    # AAA = AAA[1:group_num_true,]
    # group_num = group_num_true
    # value_choose<-seq(from = lb, to = ub,length.out = 3)
    # mu_matde<-ifelse(AAA==1,value_choose[1],ifelse(AAA==2,value_choose[2],
    #                                                value_choose[3]))
    # 
    # mu_mat[,differential_gene_index]<-mu_matde
  }else{
    group_specific_precision = fix_element$group_specific_precision
    common_adjoint = fix_element$common_adjoint
    group_specific_adjoint = fix_element$group_specific_adjoint
    mu_mat = fix_element$mu_mat
    hub_gene = fix_element$hub_gene
    ################
    n<-nrow(fix_element$obs_mat)
    group_num<-nrow(mu_mat)
    dim_use<-ncol(fix_element$obs_mat)
    hub_index<-NULL
    if(network_type == "HUB-graph"){
      hub_index<-sample(1:dim_use, 2 * dim_use * densy_degree, replace = FALSE)
    }
  }
  
  
  #5.2 library size
  ls_vec<-exp(rnorm(n,mean = log(l_mu),sd=l_sd))
  
  #5.3 pi
  pi_vec<-rep(1/group_num,group_num)
  # print("parameter design complete.")
  #
  #5.4 cluster_lable & obs_mat generation
  locator_vec<-c()
  obs_mat1<-matrix(0,nrow = n,ncol = dim_use)
  inv_group_pre<-array(NA,dim = dim(group_specific_precision))
  for(g in 1:group_num){
    inv_group_pre[g,,]<-solve(group_specific_precision[g,,])
  }
  ###########################
  locator_vec<-sample(1:group_num,n,prob = pi_vec,replace = TRUE)
  log_a_mat<-matrix(NA,nrow = n,ncol = dim_use)
  for(g in 1:group_num){
    log_a_mat[which(locator_vec == g),]<-mvrnorm(n = length(which(locator_vec == g)), (mu_mat[g,]), inv_group_pre[g,,])
  }
  log_a_mat<-log_a_mat + log(ls_vec)
  
  log_a_mat<-ifelse(log_a_mat>10,10,log_a_mat)
  a_mat<-exp(log_a_mat)
  ###########################
  for(i in 1:n){
    for(j in 1:dim_use){
      obs_mat1[i,j]<-rpois(1,a_mat[i,j])
    }
  }
  colnames(obs_mat1)<-paste("Gene",1:ncol(obs_mat1),sep = "")
  rownames(obs_mat1)<-paste("Cell",1:nrow(obs_mat1),sep = "")
  hub_gene<-colnames(obs_mat1)[hub_index]
  return(list(obs_mat = obs_mat1,locator_vec = locator_vec,
              group_specific_precision = group_specific_precision,
              common_adjoint = common_adjoint,
              group_specific_adjoint = group_specific_adjoint,
              ls_vec = ls_vec,mu_mat = mu_mat,hub_gene = hub_gene))
  
}

##data_simulation
data_simulation<-function(run_index,reptime,
                          parameter_setting_fix=list(p_v=0.5,v=0.3,
                                                     clustermethod_i=1,
                                                     densy_degree=0.1,
                                                     generator_type="MPLN"),
                          dim_use_vec,
                          group_num,
                          diff_gene_num=10,
                          sample_size_mat,
                          cluster_method_vec,
                          mu_design,
                          network_type,
                          ARI_deisgn,
                          if_fix=FALSE){
  s_i<-run_index[1]
  p_i<-run_index[2]
  m_i<-run_index[3]
  net<-run_index[4]
  ARI_i<-run_index[5]
  ##1 Set up the current parameter settings
  #####################################################################
  print(paste("s_i:",s_i,"p_i:",p_i," m_i:",m_i," net:",net," ARI_i:",ARI_i,sep = ""))
  dim_use <- dim_use_vec[p_i]
  sample_size_vec<-rep(sample_size_mat[p_i,s_i],group_num)
  n <- sum(sample_size_vec)
  mu_design_vec <- mu_design[m_i,]
  nwtype <- network_type[net]
  cluster_method<-cluster_method_vec[parameter_setting_fix$clustermethod_i]
  ARI_low<-ARI_deisgn[ARI_i,1]
  ARI_upper<-ARI_deisgn[ARI_i,2]
  #####################################################################
  
  ##2 Generate the simulation data
  #####################################################################
  
  for(rep in 1:reptime){
    ##Generate data ------------------
    print(rep)
    con_use<-TRUE
    iter0<-1
    while (con_use & iter0<500) {
      iter0<-iter0+1
      set.seed(123*rep + iter0)
      if((if_fix == TRUE) & (rep>1)){
        if(parameter_setting_fix$generator_type == "MPLN"){
          data = data_generator(n,dim_use,group_num,network_type = nwtype,
                                densy_degree =parameter_setting_fix$densy_degree,v = parameter_setting_fix$v,p_v = parameter_setting_fix$p_v,
                                num_de = diff_gene_num, ub =mu_design_vec[1],
                                lb = mu_design_vec[2],ub_non = mu_design_vec[3],
                                lb_non = mu_design_vec[4],
                                l_mu = log(10),sigma_scale = 1,alpha_com = 0.01,
                                if_fix = TRUE,fix_element = data_fix) 
        }else{
          data = data_generator_mulitnomial(n,dim_use,group_num,network_type = nwtype,
                                            densy_degree =parameter_setting_fix$densy_degree,v = parameter_setting_fix$v,p_v = parameter_setting_fix$p_v,
                                            num_de =diff_gene_num, ub =mu_design_vec[1],
                                            lb = mu_design_vec[2],ub_non = mu_design_vec[3],
                                            lb_non = mu_design_vec[4],
                                            l_mu = log(10),sigma_scale = 1,alpha_com = 0.01,
                                            if_fix = TRUE,fix_element = data_fix)
        }
      }else{
        if(parameter_setting_fix$generator_type == "MPLN"){
          data = data_generator(n,dim_use,group_num,network_type = nwtype,
                                densy_degree =  parameter_setting_fix$densy_degree,v = parameter_setting_fix$v,p_v = parameter_setting_fix$p_v,
                                num_de = diff_gene_num, ub =mu_design_vec[1],
                                lb = mu_design_vec[2],ub_non = mu_design_vec[3],
                                lb_non = mu_design_vec[4],
                                l_mu = log(10),sigma_scale = 1,alpha_com = 0.01,
                                if_fix = FALSE) 
        }else{
          data = data_generator_mulitnomial(n,dim_use,group_num,network_type = nwtype,
                                            densy_degree =  parameter_setting_fix$densy_degree,v = parameter_setting_fix$v,p_v = parameter_setting_fix$p_v,
                                            num_de =diff_gene_num, ub =mu_design_vec[1],
                                            lb = mu_design_vec[2],ub_non = mu_design_vec[3],
                                            lb_non = mu_design_vec[4],
                                            l_mu = log(10),sigma_scale = 1,alpha_com = 0.01,
                                            if_fix = FALSE)
        }
      }
      
      obs_mat = data[["obs_mat"]]
      colnames(obs_mat)<-paste("Gene",1:ncol(obs_mat),sep = "")
      rownames(obs_mat)<-paste("Cell",1:nrow(obs_mat),sep = "")
      Theta_mat_array_TRUE =data[["group_specific_precision"]]
      cluster_true = data[["locator_vec"]]
      # hub_gene<-data[["hub_gene"]]
      if(length(data[["hub_gene"]])>0){
        hub_gene<-data[["hub_gene"]]
      }else{
        hub_gene<-NULL
      }
      
      n<-dim(obs_mat)[1]
      p<-dim(obs_mat)[2]
      dr = length(which( obs_mat == 0 )) / (n*p)
      
      # 1.0 initial part
      time_init_start = Sys.time()
      gene_GRN_use<-rownames(as(t(obs_mat),"sparseMatrix"))
      suppressWarnings(res_init<-cluster_init(expression_profile = as(t(obs_mat),"sparseMatrix"),
                                               celltypes_num = group_num,
                                               celltypes_ref = NULL,
                                               ls_est = "TSS",
                                               gene_GRN = gene_GRN_use,
                                               HVG_model_num = 0,
                                               zero_GRN = NULL,
                                               preprocess_Control = list(HVG_num = length(gene_GRN_use),npc = 50,
                                                                         run_umap = FALSE,label_umap = NULL,
                                                                         cluster_method = "Kmeans",resolution = 0.8),
                                               core_num = 1))
      
      
      time_init_end<-Sys.time()
      time_init<-as.double(difftime(time_init_end,time_init_start,units = "hours"))
      
      locator = res_init$celltypes_label
      ARI_orgin = adjustedRandIndex(cluster_true,locator)
      # the right order of res_init
      cluster_predict_true_base<-as.vector(apply(table(locator,cluster_true),MARGIN = 1,which.max))
      
      #
      if(ARI_orgin<ARI_low){
        diff_gene_num<-diff_gene_num + 1
        diff_gene_num<-max(diff_gene_num,1)
      }else{
        if(ARI_orgin>ARI_upper){
          diff_gene_num<-diff_gene_num - 1
          diff_gene_num<-max(diff_gene_num,1)
        }else{
          #check whether have paired component predicted
          if(length(unique(cluster_predict_true_base)) == group_num){
            con_use<-FALSE
          }
        }
      }
      print(paste("diff_gene_num: ",diff_gene_num," con_use: ",con_use," ARI_orgin: ",ARI_orgin,sep = ""))
      
    }
    ##
    if((if_fix == TRUE) & (rep == 1)){
      data_fix<-data
    }
    ##add noise to diff_gene_num
    diff_gene_num<-diff_gene_num+sample(c(-1:1),1)
    diff_gene_num<-max(diff_gene_num,1)
    print("success")
    #save
    save(data,file = paste("/home/ruibinxi_pkuhpc/lustre1/huangyf/R_project/emplnet/simulation_right/data_1115/data626_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                           "_nettype",network_type[net],"_ARI_new",ARI_i,"_density",parameter_setting_fix$densy_degree,"_EMMMLE_omega_",EMMMLE_omega,"_rep",rep,".Rdata",sep = ""))
    save(res_init,file = paste("/home/ruibinxi_pkuhpc/lustre1/huangyf/R_project/emplnet/simulation_right/data_1115/res_init626_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                               "_nettype",network_type[net],"_ARI_new",ARI_i,"_density",parameter_setting_fix$densy_degree,"_EMMMLE_omega_",EMMMLE_omega,"_rep",rep,".Rdata",sep = ""))
  }
  print("Data generation complete.")
}

##3. Define the evaluator function for each method
##--------------------------------------------------------------------------------------------
##Evaluator for EMMMLE
EMMMLE_weight_evaluator<-function(penalize_diagonal = FALSE,
                                           if_stability_eval = FALSE,
                                           diagonal_if=FALSE,omega=0.5, initial_method="kmeans",iter_emmle=3
){
  library_size_est = "TSS"
  time_EMMMLE_start1<-Sys.time()
  repres_EMMMLE<-foreach (
    rep = 1:reptime,
    .combine = cfun,
    .inorder = TRUE,
    .export = ls(.GlobalEnv),
    .packages = c('MASS','CVXR','mclust','mltools','orthopolynom','Rcpp','RcppArmadillo','doParallel',"Seurat","EMMMLE")
  )%dopar%{
    ## load data-------------
    load(paste("/home/ruibinxi_pkuhpc/lustre1/huangyf/R_project/emplnet/simulation_right/data_1115/data626_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
               "_nettype",network_type[net],"_ARI_new",ARI_i,"_density",parameter_setting_fix$densy_degree,"_EMMMLE_omega_",EMMMLE_omega,"_rep",rep,".Rdata",sep = ""))
    load(paste("/home/ruibinxi_pkuhpc/lustre1/huangyf/R_project/emplnet/simulation_right/data_1115/res_init626_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
               "_nettype",network_type[net],"_ARI_new",ARI_i,"_density",parameter_setting_fix$densy_degree,"_EMMMLE_omega_",EMMMLE_omega,"_rep",rep,".Rdata",sep = ""))
    #####################################################################################################
    Theta_mat_array_TRUE <-data[["group_specific_precision"]]
    cluster_true <-data[["locator_vec"]]
    data_use<-data$obs_mat
    S_depth<-res_init$S_depth
    dim_use<-ncol(data_use)
    sample_size<-nrow(data_use)
    group_num<-nrow(data$mu_mat)
    bic_object_list<-list()
    ## data initialization-----------------------------------------------------
    locator_base <- res_init$celltypes_label
    ARI_orgin <- adjustedRandIndex(cluster_true,locator_base)
    cluster_predict_true_base<-as.vector(apply(table(locator_base,cluster_true),MARGIN = 1,which.max))
    table_reshape<-matrix(NA,nrow = group_num,ncol = group_num)
    for(g1 in 1:group_num){
      for(g2 in 1:group_num){
        table_reshape[g1,g2]<-length(which(locator_base == g1 & cluster_true == cluster_predict_true_base[g2]))
      }
    }
    ##
    purity_orgin<-sum(diag(table_reshape))/sum(table_reshape)
    ##
    bic_object_list[["cluster_predicted_sample"]]<-rowSums(table_reshape)
    bic_object_list[["cluster_predicted_TDR"]]<-diag(table_reshape)/rowSums(table_reshape)
    bic_object_list[["ARI_orgin"]]<-ARI_orgin
    bic_object_list[["purity_orgin"]]<-purity_orgin
    ##
    rm(cluster_predict_true_base);gc()
    rm(table_reshape);gc()
    ## initial
    EMMMLE_init_res<-MMLE_init(data,res_init)
    post_dist_Z_mat_max<-ifelse(EMMMLE_init_res$post_dist_Z > 0.9, 1, 0)
    soft_cs<-mean(as.vector(abs(EMMMLE_init_res$post_dist_Z - post_dist_Z_mat_max)))
    time_cost<-c()
    time_cost<-c(time_cost,EMMMLE_init_res$time_cost)
    init_update_res<-EMMMLE_init_res
    ########################################################################
    locator_cur<-apply(init_update_res$post_dist_Z,MARGIN = 1,which.max)
    cluster_predict_true_cur<-as.vector(apply(table(locator_cur,cluster_true),MARGIN = 1,which.max))
    t.root<-polynomial.roots(monic.polynomial.recurrences(hermite.h.recurrences(20, normalized=FALSE)))
    t_root_vec<-t.root[[21]]
    omega_root_vec<-2^19*factorial(20)*sqrt(pi)/400/polynomial.values(hermite.h.polynomials(20,normalized=FALSE), t.root[[21]])[[20]]^2
    data_obs<-data$obs_mat
    positions <- expand.grid(row = 1:dim_use, col = 1:dim_use)
    position_input <- positions[positions$row > positions$col, ]
    position_input<-as.matrix(position_input)
    if(initial_method=="EMPLN"){
      mu_initial_all<-init_update_res$mu_mat
      sigma_inital_all<-array(unlist(init_update_res$nonpos_sigma_list), dim = c(dim_use, dim_use, group_num))
      pi_init<-init_update_res$pi_vec
    }
    if(initial_method=="kmeans"){
      mu_initial_all<-res_init$mu_mat
      sigma_inital_all<-array(unlist(res_init$sigma_mat_list), dim = c(dim_use, dim_use, group_num))
      pi_init<-res_init$pi_vec
    }
    p_p_post_dist_Z<-init_update_res$post_dist_Z
    k_max<-10
    ########################################################################
    EMMMLE_result_list<-EMMMLE(data_obs,pi_init,mu_initial_all,sigma_inital_all,S_depth, t_root_vec, omega_root_vec, p_p_post_dist_Z,k_max, position_input,20,omega,iter_emmle,50)
    save(EMMMLE_result_list,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/huangyf/R_project/emplnet/simulation_right/result_all_1115/EMMMLE/EMMMLE_result_list",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                               "_nettype",network_type[net],"_ARI_new_",ARI_i,"_density",parameter_setting_fix$densy_degree,"_omega_",omega,"_initial_method",initial_method,"_iter_emmle",iter_emmle,"_rep",rep,".Rdata",sep = "")))
    if(diagonal_if==FALSE){
      sigma_res<-list()
      mu_temp<-list()
      sigma_diag_res<-list()
      for(g in 1:group_num) {
        sigma_matrix12<-matrix(0,dim_use,dim_use)
        sigma_diag<-matrix(0,dim_use,dim_use)
        mu_matrix<-matrix(0,dim_use,dim_use)
        for (i in 1:nrow(position_input)) {
          row <- position_input[i, 2]
          col <- position_input[i, 1]
          sigma_matrix12[row, col] <- EMMMLE_result_list$res_sigma12_all[g,i]
          sigma_diag[row,col]<-EMMMLE_result_list$res_sigma11_all[g,i]
          sigma_diag[col,row]<-EMMMLE_result_list$res_sigma22_all[g,i]
          mu_matrix[row,col]<-EMMMLE_result_list$res_mu_all[g,1,i]
          mu_matrix[col,row]<-EMMMLE_result_list$res_mu_all[g,2,i]
        }
        sigma_diag_res[[g]]<-sigma_diag
        mu_temp[[g]]<-mu_matrix
        sigma_res[[g]]<-sigma_matrix12+t(sigma_matrix12)
      }
      mu_res<-matrix(0,nrow = group_num,ncol=dim_use)
      mu_res<-init_update_res$mu_mat
      for (g in 1:group_num) {
        diag(sigma_res[[g]])<-diag(init_update_res$sigma_list[[g]])
      }
    }
    if(diagonal_if==TRUE){
      sigma_res<-list()
      mu_temp<-list()
      sigma_diag_res<-list()
      for(g in 1:group_num) {
        sigma_matrix12<-matrix(0,dim_use,dim_use)
        sigma_diag<-matrix(0,dim_use,dim_use)
        mu_matrix<-matrix(0,dim_use,dim_use)
        for (i in 1:nrow(position_input)) {
          row <- position_input[i, 2]
          col <- position_input[i, 1]
          sigma_matrix12[row, col] <- EMMMLE_result_list$res_sigma12_all[g,i]
          sigma_diag[row,col]<-EMMMLE_result_list$res_sigma11_all[g,i]
          sigma_diag[col,row]<-EMMMLE_result_list$res_sigma22_all[g,i]
          mu_matrix[row,col]<-EMMMLE_result_list$res_mu_all[g,1,i]
          mu_matrix[col,row]<-EMMMLE_result_list$res_mu_all[g,2,i]
        }
        sigma_diag_res[[g]]<-sigma_diag
        mu_temp[[g]]<-mu_matrix
        sigma_res[[g]]<-sigma_matrix12+t(sigma_matrix12)
      }
      mu_res<-matrix(0,nrow = group_num,ncol=dim_use)
      for (g in 1:group_num) {
        mu_res[g,]<-rowSums(mu_temp[[g]])/(dim_use-1)
        diag(sigma_res[[g]])<-rowSums(sigma_diag_res[[g]])/(dim_use-1)
      }
    }
    sigma_errorF<-c()
    for(g in 1:group_num){
      temp_sigma_error2<-sigma_res[[g]]-solve(data$group_specific_precision[cluster_predict_true_cur[g],,])
      sigma_errorF<-c(sigma_errorF,norm(temp_sigma_error2, type ="F"))
    }
    ##
    n_lambda<-100
    Theta_mat_array_hat_all<-array(NA,dim = c(dim_use,dim_use,group_num,n_lambda))
    BIC_mat_shrink <- matrix(NA,nrow = group_num,ncol = n_lambda)
    locator_cur<-apply(init_update_res$post_dist_Z,MARGIN = 1,which.max)
    cluster_predict_true_cur<-as.vector(apply(table(locator_cur,cluster_true),MARGIN = 1,which.max))
    table_reshape<-matrix(NA,nrow = group_num,ncol = group_num)
    for(g1 in 1:group_num){
      for(g2 in 1:group_num){
        table_reshape[g1,g2]<-length(which(locator_cur == g1 & cluster_true == cluster_predict_true_cur[g2]))
      }
    }
    purity_cur<-sum(diag(table_reshape))/sum(table_reshape)
    bic_object_list[["purity_cur"]]<-purity_cur
    
    ##
    cluster_acc<-adjustedRandIndex(locator_cur,cluster_true)
    bic_object_list[["cluster_acc"]]<-cluster_acc
    print(paste("ARI: ",round(cluster_acc[length(cluster_acc)],3),sep = ""))
    print(paste("purity: ",round(purity_cur,3),sep = ""))
    time_cost<-c()
    cluster_predict_true_vec<-as.vector(apply(table(locator_cur,cluster_true),MARGIN = 1,which.max))
    densy_mat<-matrix(NA,nrow = group_num,ncol = n_lambda)
    lambda_mat<-matrix(NA,nrow = group_num,ncol = n_lambda)
    for(g in 1:group_num){
      k_max<-10
      data_obs_group<-data_use[which(locator_cur==g),]
      S_depth_group<-res_init$S_depth[which(locator_cur==g)]
      cov_input_equal<-sigma_res[[g]]
      k_index<-k_max
      ##
      ##Use dtrace-loss for precision matrix's estimation
      EMMMLE_res<-EMMMLE_Network(obs_mat = data_obs_group,
                               Sd_est = S_depth_group,
                               n_lambda = n_lambda,
                               penalize.diagonal = FALSE,
                               cov_input =cov_input_equal ,
                               weight_mat = NULL,zero_mat = NULL,
                               core_num = 1)
      time_cost_EQUAL<-EMMMLE_res$time_EQUAL
      time_cost<-c(time_cost,time_cost_EQUAL)
      lambda_vec1<-EMMMLE_res$lambda_vec
      lambda_mat[g,]<-lambda_vec1
      for (l in 1:length(lambda_vec1)){
        Theta_mat_array_hat_all[,,g,l]<-as.matrix(EMMMLE_res$Omega_est[[l]])
        BIC_mat_shrink[g,l]<-EMMMLE_res$BIC_vec[l]
      }
    }
    time_use<-sum(time_cost)
    cluster_predict_true_vec<-as.vector(apply(table(locator_cur,cluster_true),MARGIN = 1,which.max))
    bic_object_list[["cluster_predict_true_vec"]]<-cluster_predict_true_vec
    ##ARI & purity
    ARI_end<-cluster_acc
    purity_end<-purity_cur
    ###
    
    ##save
    save(Theta_mat_array_hat_all,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/huangyf/R_project/emplnet/simulation_right/result_all_1115/EMMMLE/Theta_mat_array_hat_all",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                                    "_nettype",network_type[net],"_ARI_new_",ARI_i,"_density",parameter_setting_fix$densy_degree,"_omega_",omega,"_diagonal_if_",diagonal_if,"_initial_method",initial_method,"_iter_emmle",iter_emmle,"_rep",rep,".Rdata",sep = "")))
    save(cluster_predict_true_vec,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/huangyf/R_project/emplnet/simulation_right/result_all_1115/EMMMLE/cluster_predict_true_vec_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                                     "_nettype",network_type[net],"_ARI_new_",ARI_i,"_density",parameter_setting_fix$densy_degree,"_omega_",omega,"_diagonal_if_",diagonal_if,"_initial_method",initial_method,"_iter_emmle",iter_emmle,"_rep",rep,".Rdata",sep = "")))
    ##
    evaluator_2_res<-evaluator_2(input_array_array=Theta_mat_array_hat_all,
                                 Theta_array_true=Theta_mat_array_TRUE,
                                 lambda_vec1=lambda_vec1,
                                 cluster_predict_true_vec= cluster_predict_true_vec)
    
    edge_true_num<-c()
    for(g in 1:group_num){
      edge_true_num<-c(edge_true_num,length(which((Theta_mat_array_TRUE[g,,])[upper.tri((Theta_mat_array_TRUE[g,,]))] != 0)))
    }
    
    res_vec<-c(
      ARI_orgin,
      ARI_end,
      ##1:2
      
      purity_orgin,
      purity_end,
      ##3:4
      
      soft_cs,time_use,
      ##5:6
      
      as.vector(sigma_errorF),
      ##7:9
      
      ##BIC part
      as.vector(evaluator_2_res),
      ##10:15
      
      edge_true_num
    )
    
    save(res_vec,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/huangyf/R_project/emplnet/simulation_right/result_all_1115/EMMMLE/res_vec_1015_emplnet",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                    "_nettype",network_type[net],"_ARI_new_",ARI_i,"_density",parameter_setting_fix$densy_degree,"_omega_",omega,"_diagonal_if_",diagonal_if,"_initial_method",initial_method,"_iter_emmle",iter_emmle,"_rep",rep,".Rdata",sep = "")))
    ##
    return(res_vec)
  }
  time_EMMMLE_end1<-Sys.time()
  time_cost_EMMMLE<-as.double(difftime(time_EMMMLE_end1,time_EMMMLE_start1,units = "hours"))
  ##stability
  # early_stability_BIC<-NULL
  # early_stability_density<-NULL
  # early_stability_density2<-NULL
  # early_stability_BIC_signed<-NULL
  # early_stability_density_signed<-NULL
  # early_stability_density2_signed<-NULL
  
  return(list(repres_EMMMLE = repres_EMMMLE,
              time_cost_EMMMLE=time_cost_EMMMLE))
}

##Evaluator for PLNet
PLNet_evaluator<-function(penalize_diagonal = FALSE,if_stability_eval = FALSE){
  library_size_est = "TSS"
  ##
  time_PLNet_start1<-Sys.time()
  repres_PLNet<-foreach (
    rep = 1:reptime,
    .combine = cfun,
    .inorder = TRUE,
    .export = ls(.GlobalEnv),
    .packages = c('MASS', 'glasso','PLNmodels','doParallel','CVXR','igraph',"PLNetNew","orthopolynom","Rcpp","RcppArmadillo")
  )%dopar%{
    ## load data-------------
    load(paste("/home/ruibinxi_pkuhpc/lustre1/huangyf/R_project/emplnet/simulation_right/data_1115/data626_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
               "_nettype",network_type[net],"_ARI_new",ARI_i,"_density",parameter_setting_fix$densy_degree,"_EMMMLE_omega_",EMMMLE_omega,"_rep",rep,".Rdata",sep = ""))
    load(paste("/home/ruibinxi_pkuhpc/lustre1/huangyf/R_project/emplnet/simulation_right/data_1115/res_init626_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
               "_nettype",network_type[net],"_ARI_new",ARI_i,"_density",parameter_setting_fix$densy_degree,"_EMMMLE_omega_",EMMMLE_omega,"_rep",rep,".Rdata",sep = ""))
    
    #
    bic_object_list<-list()
    ##
    obs_mat = data[["obs_mat"]]
    colnames(obs_mat)<-paste("Gene",1:ncol(obs_mat),sep = "")
    rownames(obs_mat)<-paste("Cell",1:nrow(obs_mat),sep = "")
    Theta_mat_array_TRUE =data[["group_specific_precision"]]
    cluster_true = data[["locator_vec"]]
    # hub_gene<-data[["hub_gene"]]
    if(length(data[["hub_gene"]])>0){
      hub_gene<-data[["hub_gene"]]
    }else{
      hub_gene<-NULL
    }
    
    n<-dim(obs_mat)[1]
    p<-dim(obs_mat)[2]
    dr = length(which( obs_mat == 0 )) / (n*p)
    
    ls_vec<-res_init$ls
    ##
    
    ##
    locator_base <- res_init$celltypes_label 
    ARI_orgin = adjustedRandIndex(cluster_true,locator_base)
    cluster_predict_true_base<-as.vector(apply(table(locator_base,cluster_true),MARGIN = 1,which.max))
    
    ###
    table_reshape<-matrix(NA,nrow = group_num,ncol = group_num)
    for(g1 in 1:group_num){
      for(g2 in 1:group_num){
        table_reshape[g1,g2]<-length(which(locator_base == g1 & cluster_true == cluster_predict_true_base[g2]))
      }
    }
    bic_object_list[["cluster_predicted_sample"]]<-rowSums(table_reshape)
    bic_object_list[["cluster_predicted_TDR"]]<-diag(table_reshape)/rowSums(table_reshape)
    ##
    diag(table_reshape)<-0
    rm(table_reshape);gc()
    ##################
    
    print("PLNet")
    
    ##################
    n_lambda<-100
    Theta_mat_array_hat_all<-array(NA,dim = c(p,p,group_num,n_lambda))
    BIC_mat_shrink <- matrix(NA,nrow = group_num,ncol = n_lambda)
    ##
    densy_mat<-matrix(NA,nrow = group_num,ncol = n_lambda)
    lambda_mat<-matrix(NA,nrow = group_num,ncol = n_lambda)
    cluster_predict_true_vec<-as.vector(apply(table(res_init$celltypes_label,cluster_true),MARGIN = 1,which.max))
    ##
    #PLNet_data = obs_mat[locator_base==g2,]
    ##
    time_cost<-c()
    sigma_error<-c()
    for (g2 in 1:group_num) {
      print(paste("group(",g2,") task of PLNet."))
      time_PLNet_start<-Sys.time()
      
      ##load the data
      Y<-obs_mat[locator_base==g2,]
      
      ##Estimate the library size by total sum scaling.
      #S_depth<-rowMeans(Y) 
      S_depth<-res_init$S_depth[locator_base==g2]
      
      k_max<-10 ## The maximal iter for MLE calculation by newton-raphson algorithm
      
      ##MLE for convariance matrix
      cov_input<-mle_newton(data_use = Y,
                            S_depth = S_depth,
                            core_num = 1,
                            k_max = k_max)
      ##
      k_index<-k_max
      temp_matrix0<-cov_input$mlesigmahat[[1 + k_index]]-solve(data$group_specific_precision[cluster_predict_true_vec[g2],,])
      sigma_error<-c(sigma_error,norm(temp_matrix0, type ="2"))
      ##
      ##Use dtrace-loss for precision matrix's estimation
      PLNet_res<-PLNet_main(obs_mat = Y,
                            Sd_est = S_depth,
                            n_lambda = n_lambda,
                            penalize.diagonal = FALSE,
                            cov_input = cov_input$mlesigmahat[[1 + k_index]],
                            weight_mat = NULL,zero_mat = NULL,
                            core_num = 1)
      lambda_vec1<-PLNet_res$lambda_vec
      lambda_mat[g2,]<-lambda_vec1
      for (l in 1:length(lambda_vec1)){
        Theta_mat_array_hat_all[,,g2,l]<-as.matrix(PLNet_res$Omega_est[[l]])
        BIC_mat_shrink[g2,l]<-PLNet_res$BIC_vec[l]
      }
      ##
      
      time_PLNet_end<-Sys.time()
      time_cost<-c(time_cost,as.double(difftime(time_PLNet_end,time_PLNet_start,units = "hours")))
    }
    
    time_use<-sum(time_cost)
    
    ######################
    ###
    
    cluster_predict_true_vec<-as.vector(apply(table(res_init$celltypes_label,cluster_true),MARGIN = 1,which.max))
    
    bic_object_list[["cluster_predict_true_vec"]]<-cluster_predict_true_vec
    
    ##
    save(Theta_mat_array_hat_all,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/huangyf/R_project/emplnet/simulation_right/result_all_1115/PLNet/Theta_array_hat_BIC_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                                    "_nettype",network_type[net],"_ARI_new_",ARI_i,"_density",parameter_setting_fix$densy_degree,"_EMMMLE_omega_",EMMMLE_omega,"_rep",rep,".Rdata",sep = "")))
    save(cluster_predict_true_vec,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/huangyf/R_project/emplnet/simulation_right/result_all_1115/PLNet/cluster_predict_true_vec_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                                     "_nettype",network_type[net],"_ARI_new_",ARI_i,"_density",parameter_setting_fix$densy_degree,"_EMMMLE_omega_",EMMMLE_omega,"_rep",rep,".Rdata",sep = "")))
    ##
    evaluator_2_res<-evaluator_2(input_array_array=Theta_mat_array_hat_all,
                                 Theta_array_true=Theta_mat_array_TRUE,
                                 lambda_vec1=lambda_vec1,
                                 cluster_predict_true_vec= cluster_predict_true_vec)
    
    #############################
    edge_true_num<-c()
    for(g in 1:group_num){
      edge_true_num<-c(edge_true_num,length(which((Theta_mat_array_TRUE[g,,])[upper.tri((Theta_mat_array_TRUE[g,,]))] != 0)))
    }
    ##
    res_vec<-c(
      time_use,
      ##1
      
      ##BIC part
      as.vector(sigma_error),
      
      ##density part
      as.vector(evaluator_2_res),
      
      edge_true_num
    )
    save(res_vec,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/huangyf/R_project/emplnet/simulation_right/result_all_1115/PLNet/res_vec_1015_emplnet",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                    "_nettype",network_type[net],"_ARI_new_",ARI_i,"_density",parameter_setting_fix$densy_degree,"_EMMMLE_omega_",EMMMLE_omega,"_rep",rep,".Rdata",sep = "")))
    
    #####
    return(res_vec)
    
  }
  time_PLNet_end1<-Sys.time()
  time_cost_PLNet<-as.double(difftime(time_PLNet_end1,time_PLNet_start1,units = "hours"))
  
  ##stability
  # early_stability_BIC<-NULL
  # early_stability_density<-NULL
  # early_stability_density2<-NULL
  # early_stability_BIC_signed<-NULL
  # early_stability_density_signed<-NULL
  # early_stability_density2_signed<-NULL
  ##
  
  return(list(repres_PLNet = repres_PLNet,
              time_cost_PLNet=time_cost_PLNet))
}

##Evaluator for VMPLN
VMPLN_evaluator<-function(penalize_diagonal = FALSE,
                              if_stability_eval = FALSE
){
  library_size_est = "TSS"
  ##
  time_scGeneNet_start1<-Sys.time()
  repres_scGeneNet<-foreach (
    rep = 1:reptime,
    .combine = cfun,
    .inorder = TRUE,
    .export = ls(.GlobalEnv),
    .packages = c('MASS','glasso',
                  'mclust','mltools','data.table','PLNmodels','dplyr',
                  "cluster","pryr","philentropy","umap","Seurat","igraph","Rcpp","scGeneNet","Matrix","stats","fpc","GMPR")
  )%dopar%{
    ## load data-------------
    load(paste("/home/ruibinxi_pkuhpc/lustre1/huangyf/R_project/emplnet/simulation_right/data_1115/data626_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
               "_nettype",network_type[net],"_ARI_new",ARI_i,"_density",parameter_setting_fix$densy_degree,"_EMMMLE_omega_",EMMMLE_omega,"_rep",rep,".Rdata",sep = ""))
    
    
    
    ##data & res_init
    obs_mat = data[["obs_mat"]]
    colnames(obs_mat)<-paste("Gene",1:ncol(obs_mat),sep = "")
    rownames(obs_mat)<-paste("Cell",1:nrow(obs_mat),sep = "")
    Theta_mat_array_TRUE =data[["group_specific_precision"]]
    cluster_true = data[["locator_vec"]]
    if(length(data[["hub_gene"]])>0){
      hub_gene<-data[["hub_gene"]]
    }else{
      hub_gene<-NULL
    }
    
    n<-dim(obs_mat)[1]
    p<-dim(obs_mat)[2]
    
    group_num<-nrow(data$mu_mat)
    bic_object_list<-list()
    gene_GRN_use<-rownames(as(t(obs_mat),"sparseMatrix"))
    suppressWarnings(res_init<-scGeneNet_init(expression_profile = as(t(obs_mat),"sparseMatrix"),
                                              celltypes_num = group_num,
                                              celltypes_ref = NULL,
                                              ls_est = "TSS",
                                              gene_GRN = gene_GRN_use,
                                              HVG_model_num = 0,
                                              zero_GRN = NULL,
                                              preprocess_Control = list(HVG_num = length(gene_GRN_use),npc = 50,
                                                                        run_umap = FALSE,label_umap = NULL,
                                                                        cluster_method = "Kmeans",resolution = 0.8),
                                              core_num = 1))
    ## data initialization------------------
    dr = length(which( obs_mat == 0 )) / (n*p)
    
    locator_base = apply(res_init$U_mat, 1, which.max)
    ARI_orgin = adjustedRandIndex(cluster_true,locator_base)
    cluster_predict_true_base<-as.vector(apply(table(locator_base,cluster_true),MARGIN = 1,which.max))
    table_reshape<-matrix(NA,nrow = group_num,ncol = group_num)
    for(g1 in 1:group_num){
      for(g2 in 1:group_num){
        table_reshape[g1,g2]<-length(which(locator_base == g1 & cluster_true == cluster_predict_true_base[g2]))
      }
    }
    ##
    purity_orgin<-sum(diag(table_reshape))/sum(table_reshape)
    ##
    bic_object_list[["cluster_predicted_sample"]]<-rowSums(table_reshape)
    bic_object_list[["cluster_predicted_TDR"]]<-diag(table_reshape)/rowSums(table_reshape)
    bic_object_list[["ARI_orgin"]]<-ARI_orgin
    bic_object_list[["purity_orgin"]]<-purity_orgin
    ##
    rm(cluster_predict_true_base);gc()
    rm(table_reshape);gc()
    
    ## method run-------------
    ######################################
    ## 1. prior information generation
    prior_nw_list<-list()
    for(g in 1:group_num){
      prior_nw_list[[g]] <- prior_generator(data$common_adjoint,eta = 0)
    }
    prior_list<-list(prior_nw_list = prior_nw_list,
                     prior_TF = NULL)
    rm(prior_nw_list);gc()
    ##
    prior_tran_res<-prior_tran(p_GRN = length(res_init$gene_GRN_index_use),
                               group_num = length(unique(res_init$celltypes_label)),
                               gene_GRN_name = res_init$gene_GRN_vec,
                               prior_nw_list = prior_list$prior_nw_list, prior_TF = prior_list$prior_TF,
                               weight_prior_nw = 0.3,penalize_diagonal = penalize_diagonal)
    weight_mat_list<-prior_tran_res$weight_mat_list
    zero_mat<-prior_tran_res$zero_mat
    rm(prior_tran_res)
    
    ## one step update for lambda = 1e-6 (U_fix = true)
    res_init<-scGeneNet_main(scGeneNet_list = res_init,lambda_use = 1e-6,
                             Theta_Control = list(penalize_diagonal = FALSE),
                             U_fix = TRUE,
                             verbose = FALSE,
                             core_num = 1)
    ##
    ######################################
    lambda_max_res <- exact_lambdamax(scGeneNet_list = res_init,
                                      Global_Control = list(ELBO_threshold = 1e-4, minit = 1, maxit = 50,
                                                            maxit_nonGRN = 10,MS_update_threshold = 1e-6),
                                      M_Control = list(ADMM_max_step = 1000, ADMM_threshold = 1e-4),
                                      S_Control = list(Newton_max_step = 1000, Newton_threshold = 1e-04),
                                      Theta_Control = list(penalize_diagonal = penalize_diagonal,
                                                           Theta_threshold = 1e-04),
                                      U_fix = FALSE,
                                      core_num = 1)
    ##
    lambda_min_res<-1e-6
    lambda_max_uni<-sort(unique(lambda_max_res))
    lambda_length0<-70
    lambda_length1<-40
    if(length(lambda_max_uni) == 1){
      
      lambda_vec<-exp(seq(from=log(lambda_min_res),to=log(lambda_max_uni[1]),length.out = lambda_length0))
      
      
    }else{
      lambda_max_uni<-c(lambda_min_res,lambda_max_uni)
      lambda_vec<-c()
      
      for(k in 1:(length(lambda_max_uni) - 1)){
        
        if(k == 1){
          lambda_vec<-c(lambda_vec,(seq(from=(lambda_max_uni[k]),to=(lambda_max_uni[k+1]),length.out = lambda_length0)))
        }else{
          lambda_vec<-c(lambda_vec,(seq(from=(lambda_max_uni[k]),to=(lambda_max_uni[k+1]),length.out = lambda_length1)))
          
        }
      }
    }
    print("lambda_vec has been determine.")
    
    ## 3. model run
    
    BIC_mat_shrink = matrix(NA,nrow = group_num,ncol = length(lambda_vec))
    BIC_mat_shrink1 = matrix(NA,nrow = group_num,ncol = length(lambda_vec))
    cluster_acc = NULL
    U_mat_array<-array(NA,dim = c(length(lambda_vec),nrow(obs_mat),group_num))
    Theta_mat_array_hat_all<-array(NA,dim = c(p,p,group_num,length(lambda_vec)))
    densy_mat<-matrix(NA,nrow = group_num,ncol = length(lambda_vec))
    soft_cs_vec<-c()
    purity_vec<-c()
    gc()
    ##
    
    
    time_cost<-c()
    res_scGeneNet_NN<-NULL
    for (l in 1:length(lambda_vec)) {
      print(paste("lambda(",l,") = ",lambda_vec[l],sep = ""))
      res_init0<-res_init
      ##
      rm(res_scGeneNet_NN);gc()
      ##
      time_MPLN_start<-Sys.time()
      res_scGeneNet_NN<-scGeneNet_main(scGeneNet_list = res_init0,lambda_use = lambda_vec[l],
                                       Global_Control = list(ELBO_threshold = 1e-4, minit = 1, maxit = 50,maxit_nonGRN = 10, MS_update_threshold = 1e-6),
                                       M_Control = list(ADMM_max_step = 1000, ADMM_threshold = 1e-4),
                                       S_Control = list(Newton_max_step = 1000, Newton_threshold = 1e-4),
                                       Theta_Control = list(penalize_diagonal = FALSE,Theta_threshold = 1e-4),
                                       U_fix = FALSE,
                                       verbose = TRUE,
                                       core_num = 1)
      
      rm(res_init0);gc()
      time_MPLN_end<-Sys.time()
      time_cost<-c(time_cost,as.double(difftime(time_MPLN_end,time_MPLN_start,units = "hours")))
      
      BIC_mat_shrink[,l]<-res_scGeneNet_NN$scGeneNet_bic_VMICL
      BIC_mat_shrink1[,l]<-res_scGeneNet_NN$scGeneNet_bic_VICL
      ##
      
      for(g in 1:group_num){
        Theta_mat_array_hat_all[,,g,l]<-as.matrix(res_scGeneNet_NN$Theta_mat_list[[g]])
      }
      ##
      for(g in 1:group_num){
        densy_mat[g,l]<- length(which((Theta_mat_array_hat_all[,,g,l])[upper.tri(Theta_mat_array_hat_all[,,g,l])]!=0))/((nrow(Theta_mat_array_hat_all[,,g,l])*(nrow(Theta_mat_array_hat_all[,,g,l])-1))/2)
      }
      ##
      U_mat_array[l,,]<-res_scGeneNet_NN$U_mat
      ## Calculate the soft clustering strength
      U_mat_max<-matrix(0,nrow = nrow(res_scGeneNet_NN$U_mat),ncol = ncol(res_scGeneNet_NN$U_mat))
      for(i in 1:nrow(U_mat_max)){
        U_mat_max[i,which.max(res_scGeneNet_NN$U_mat[i,])]<-1
      }
      soft_cs<-mean(as.vector(abs(res_scGeneNet_NN$U_mat - U_mat_max)))
      soft_cs_vec<-c(soft_cs_vec,soft_cs)
      ##purity
      locator_cur<-apply(res_scGeneNet_NN$U_mat,MARGIN = 1,which.max)
      cluster_predict_true_cur<-as.vector(apply(table(locator_cur,cluster_true),MARGIN = 1,which.max))
      table_reshape<-matrix(NA,nrow = group_num,ncol = group_num)
      for(g1 in 1:group_num){
        for(g2 in 1:group_num){
          table_reshape[g1,g2]<-length(which(locator_cur == g1 & cluster_true == cluster_predict_true_cur[g2]))
        }
      }
      ##
      purity_cur<-sum(diag(table_reshape))/sum(table_reshape)
      purity_vec<-c(purity_vec,purity_cur)
      ##
      
      ##
      cluster_acc<-c(cluster_acc,adjustedRandIndex(apply(res_scGeneNet_NN$U_mat,MARGIN = 1,which.max),cluster_true))
      print(paste("ARI: ",round(cluster_acc[length(cluster_acc)],3),sep = ""))
      print(paste("purity: ",round(purity_cur,3),sep = ""))
    }
    ###################################
    ##save
    bic_object_list[["BIC_VMICL"]]<-BIC_mat_shrink
    bic_object_list[["BIC_VICL"]]<-BIC_mat_shrink1
    bic_object_list[["soft_cs_vec"]]<-soft_cs_vec
    bic_object_list[["purity_vec"]]<-purity_vec
    ##
    soft_cs_mean<-mean(soft_cs_vec)
    ##
    bic_object_list[["cluster_acc"]]<-cluster_acc
    ##
    rm(obs_mat,prior_list,res_init0);gc()
    time_use<-sum(time_cost)
    print("model has been run.")
    #############
    ##BIC_shrink_choose
    BIC_vec_shrink<-colSums(BIC_mat_shrink)
    BIC_sum_choose_shrink<-which.min(BIC_vec_shrink)
    BIC_sep_choose_shrink<-apply(BIC_mat_shrink,MARGIN = 1,which.min)
    ##
    BIC_vec_shrink1<-colSums(BIC_mat_shrink1)
    BIC_sum_choose_shrink1<-which.min(BIC_vec_shrink1)
    BIC_sep_choose_shrink1<-apply(BIC_mat_shrink1,MARGIN = 1,which.min)
    ##
    bic_object_list[["BIC_sep_choose_VMICL"]]<-BIC_sep_choose_shrink
    bic_object_list[["BIC_sum_choose_VMICL"]]<-BIC_sum_choose_shrink
    bic_object_list[["BIC_sep_choose_VICL"]]<-BIC_sep_choose_shrink1
    bic_object_list[["BIC_sum_choose_VICL"]]<-BIC_sum_choose_shrink1
    ###
    ##ARI & purity
    ARI_mean<-mean(cluster_acc)
    ARI_max<-max(cluster_acc)
    ARI_min<-min(cluster_acc)
    ARI_BIC_sum_shrink<-cluster_acc[BIC_sum_choose_shrink]
    
    purity_mean<-mean(purity_vec)
    purity_max<-max(purity_vec)
    purity_min<-min(purity_vec)
    purity_BIC_sum_shrink<-purity_vec[BIC_sum_choose_shrink]
    ###
    ####################
    cluster_predict_true_vec<-as.vector(apply(table(apply(U_mat_array[1,,],MARGIN = 1,which.max),cluster_true),MARGIN = 1,which.max))
    evaluator_2_res<-evaluator_2(input_array_array=Theta_mat_array_hat_all,
                                 Theta_array_true=Theta_mat_array_TRUE,
                                 lambda_vec1=lambda_vec,
                                 cluster_predict_true_vec= cluster_predict_true_vec)
    
    ##save
    save(Theta_mat_array_hat_all,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/huangyf/R_project/emplnet/simulation_right/result_all_1115/scGeneNet/Theta_array_hat_VICL_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                                    "_nettype",network_type[net],"_ARI_new_",ARI_i,"_density",parameter_setting_fix$densy_degree,"_EMMMLE_omega_",EMMMLE_omega,"_rep",rep,".Rdata",sep = "")))
    save(cluster_predict_true_vec,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/huangyf/R_project/emplnet/simulation_right/result_all_1115/scGeneNet/cluster_predict_true_vec_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                                     "_nettype",network_type[net],"_ARI_new_",ARI_i,"_density",parameter_setting_fix$densy_degree,"_EMMMLE_omega_",EMMMLE_omega,"_rep",rep,".Rdata",sep = "")))
    
    #######################################################
    ##
    edge_true_num<-c()
    for(g in 1:group_num){
      edge_true_num<-c(edge_true_num,length(which((Theta_mat_array_TRUE[g,,])[upper.tri((Theta_mat_array_TRUE[g,,]))] != 0)))
    }
    ##
    res_vec<-c(
      ARI_orgin,
      ARI_mean,
      ARI_max,
      ARI_min,
      ARI_BIC_sum_shrink,
      ##1:5
      
      purity_orgin,
      purity_mean,
      purity_max,
      purity_min,
      purity_BIC_sum_shrink,
      ##6:10
      
      soft_cs_mean,time_use,dr,
      ##11:13
      
      ##non-signed
      ##VMICL part
      as.vector(evaluator_2_res),
      ##14:19
      edge_true_num
      ##20:22
    )
    save(res_vec,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/huangyf/R_project/emplnet/simulation_right/result_all_1115/scGeneNet/res_vec_1015_emplnet",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                    "_nettype",network_type[net],"_ARI_new_",ARI_i,"_density",parameter_setting_fix$densy_degree,"_EMMMLE_omega_",EMMMLE_omega,"_rep",rep,".Rdata",sep = "")))
    
    ##
    return(res_vec)
  }
  time_scGeneNet_end1<-Sys.time()
  time_cost_scGeneNet<-as.double(difftime(time_scGeneNet_end1,time_scGeneNet_start1,units = "hours"))
  
  ##stability
  # early_stability_VMICL<-NULL
  # early_stability_VICL<-NULL
  # early_stability_density<-NULL
  # early_stability_density2<-NULL
  # early_stability_VMICL_signed<-NULL
  # early_stability_VICL_signed<-NULL
  # early_stability_density_signed<-NULL
  # early_stability_density2_signed<-NULL
  ##
  
  return(list(repres_scGeneNet = repres_scGeneNet,
              time_cost_scGeneNet=time_cost_scGeneNet))
}

##Evaluator for Glasso
Glasso_evaluator<-function(penalize_diagonal = FALSE,if_stability_eval = FALSE){
  library_size_est = "TSS"
  ##
  time_Glasso_start1<-Sys.time()
  repres_Glasso<-foreach (
    rep = 1:reptime,
    .combine = cfun,
    .inorder = TRUE,
    .export = ls(.GlobalEnv),
    .packages = c('MASS','glasso',
                  'mclust','mltools','data.table','PLNmodels','dplyr',
                  "cluster","pryr","philentropy","umap","Seurat","igraph","Rcpp","scGeneNet","Matrix","stats","GMPR")
  )%dopar%{
    ## load data-------------
    load(paste("/home/ruibinxi_pkuhpc/lustre1/huangyf/R_project/emplnet/simulation_right/data_1115/data626_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
               "_nettype",network_type[net],"_ARI_new",ARI_i,"_density",parameter_setting_fix$densy_degree,"_EMMMLE_omega_",EMMMLE_omega,"_rep",rep,".Rdata",sep = ""))
    load(paste("/home/ruibinxi_pkuhpc/lustre1/huangyf/R_project/emplnet/simulation_right/data_1115/res_init626_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
               "_nettype",network_type[net],"_ARI_new",ARI_i,"_density",parameter_setting_fix$densy_degree,"_EMMMLE_omega_",EMMMLE_omega,"_rep",rep,".Rdata",sep = ""))
    ##
    bic_object_list<-list()
    ##
    obs_mat = data[["obs_mat"]]
    colnames(obs_mat)<-paste("Gene",1:ncol(obs_mat),sep = "")
    rownames(obs_mat)<-paste("Cell",1:nrow(obs_mat),sep = "")
    Theta_mat_array_TRUE =data[["group_specific_precision"]]
    cluster_true = data[["locator_vec"]]
    # hub_gene<-data[["hub_gene"]]
    if(length(data[["hub_gene"]])>0){
      hub_gene<-data[["hub_gene"]]
    }else{
      hub_gene<-NULL
    }
    
    n<-dim(obs_mat)[1]
    p<-dim(obs_mat)[2]
    dr = length(which( obs_mat == 0 )) / (n*p)
    
    ls_vec<-res_init$S_depth
    ##
    if(library_size_est == "none"){
      ls_vec<-rep(1,length(res_init$ls))
    }
    if(library_size_est == "TSS"){
      ls_vec<-rowSums(obs_mat)/10000
    }
    if(library_size_est == "GMPR"){
      ls_vec<-as.vector(GMPR(obs_mat))
    }
    ##
    locator_base <- res_init$celltypes_label 
    ARI_orgin = adjustedRandIndex(cluster_true,locator_base)
    cluster_predict_true_base<-as.vector(apply(table(locator_base,cluster_true),MARGIN = 1,which.max))
    
    ###
    table_reshape<-matrix(NA,nrow = group_num,ncol = group_num)
    for(g1 in 1:group_num){
      for(g2 in 1:group_num){
        table_reshape[g1,g2]<-length(which(locator_base == g1 & cluster_true == cluster_predict_true_base[g2]))
      }
    }
    bic_object_list[["cluster_predicted_sample"]]<-rowSums(table_reshape)
    bic_object_list[["cluster_predicted_TDR"]]<-diag(table_reshape)/rowSums(table_reshape)
    ##
    diag(table_reshape)<-0
    rm(table_reshape);gc()
    ##################
    
    print("Glasso")
    
    ##################
    Theta_mat_array_hat_all<-array(NA,dim = c(p,p,group_num,100))
    ##
    densy_mat<-matrix(NA,nrow = group_num,ncol = 100)
    BIC_mat<-matrix(NA,nrow = group_num,ncol = 100)
    BIC_sep_choose<-c()
    
    ##
    glasso_data = log((1+obs_mat)/ls_vec)
    rm(obs_mat);gc()
    lambda_vec<-exp(seq(from=log((1e-6)),to=log((50)),length.out = 100))
    
    
    ##
    time_cost<-c()
    for (g2 in 1:group_num) {
      print(paste("group(",g2,") task of Glasso."))
      time_Glasso_start<-Sys.time()
      graph_lasso_array<-array(0,dim = c(length(lambda_vec),ncol(glasso_data),ncol(glasso_data)))
      cov_X = cov(glasso_data[locator_base==g2,])
      dim_p = dim(glasso_data)[2]
      cov_X = cov_X + 1e-3 * diag(dim_p)
      for (l in 1:length(lambda_vec)) {
        print(l)
        # graph_glasso<-graph_glasso_list[,,l]
        graph_glasso<-glasso(s=cov_X,rho = lambda_vec[l],penalize.diagonal = penalize_diagonal,thr = 1e-3)$wi
        ##
        graph_lasso_array[l,,]<-graph_glasso
        Theta_mat_array_hat_all[,,g2,l]<-graph_glasso
        densy_mat[g2,l]<-length(which(graph_glasso[upper.tri(graph_glasso)]!=0))/(nrow(graph_glasso) * (nrow(graph_glasso)-1)/2)
        rm(graph_glasso);gc()
        
        gc()
      }
      rm(cov_X);gc()
      BIC_glasso_group<-BIC_glasso(glasso_data[locator_base==g2,],graph_lasso_array)
      ##
      BIC_mat[g2,]<-BIC_glasso_group
      ##
      rm(graph_lasso_array);gc()
      BIC_sep_choose<-c(BIC_sep_choose,which.min(BIC_glasso_group))
      
      time_Glasso_end<-Sys.time()
      time_cost<-c(time_cost,as.double(difftime(time_Glasso_end,time_Glasso_start,units = "hours")))
    }
    
    time_use<-sum(time_cost)
    
    ######################
    
    cluster_predict_true_vec<-as.vector(apply(table(res_init$celltypes_label,cluster_true),MARGIN = 1,which.max))
    
    bic_object_list[["cluster_predict_true_vec"]]<-cluster_predict_true_vec
    
    ##
    save(Theta_mat_array_hat_all,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/huangyf/R_project/emplnet/simulation_right/result_all_1115/Glasso/Theta_array_hat_BIC_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                                    "_nettype",network_type[net],"_ARI_new_",ARI_i,"_density",parameter_setting_fix$densy_degree,"_EMMMLE_omega_",EMMMLE_omega,"_rep",rep,".Rdata",sep = "")))
    save(cluster_predict_true_vec,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/huangyf/R_project/emplnet/simulation_right/result_all_1115/Glasso/cluster_predict_true_vec_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                                     "_nettype",network_type[net],"_ARI_new_",ARI_i,"_density",parameter_setting_fix$densy_degree,"_EMMMLE_omega_",EMMMLE_omega,"_rep",rep,".Rdata",sep = "")))
    ##
    
    evaluator_2_res<-evaluator_2(input_array_array=Theta_mat_array_hat_all,
                                 Theta_array_true=Theta_mat_array_TRUE,
                                 lambda_vec1=lambda_vec,
                                 cluster_predict_true_vec= cluster_predict_true_vec)
    
    ############################
    edge_true_num<-c()
    for(g in 1:group_num){
      edge_true_num<-c(edge_true_num,length(which((Theta_mat_array_TRUE[g,,])[upper.tri((Theta_mat_array_TRUE[g,,]))] != 0)))
    }
    ##
    res_vec<-c(
      time_use,
      ##1
      
      ##BIC part
      as.vector(evaluator_2_res),
      #2:7
      edge_true_num
    )
    save(res_vec,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/huangyf/R_project/emplnet/simulation_right/result_all_1115/Glasso/res_vec_1015_emplnet",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                    "_nettype",network_type[net],"_ARI_new_",ARI_i,"_density",parameter_setting_fix$densy_degree,"_EMMMLE_omega_",EMMMLE_omega,"_rep",rep,".Rdata",sep = "")))
    
    #####
    return(res_vec)
    
    
  }
  time_Glasso_end1<-Sys.time()
  time_cost_Glasso<-as.double(difftime(time_Glasso_end1,time_Glasso_start1,units = "hours"))
  
  # ##stability
  # early_stability_BIC<-NULL
  # early_stability_density<-NULL
  # early_stability_density2<-NULL
  # early_stability_BIC_signed<-NULL
  # early_stability_density_signed<-NULL
  # early_stability_density2_signed<-NULL
  # ##
  
  return(list(repres_Glasso = repres_Glasso,
              time_cost_Glasso=time_cost_Glasso))
}
##4. Simulation setting
##--------------------------------------------------------------------------------------------

##4.1 overall setting
#####################################################################
##generate way
generator_type_vec<-c("MPLN","compositional")
##dimension
dim_use_vec<-c(50,100,200,300,400)

##sample size for each cell type
sample_size_mat<-t(sapply(X=1:length(dim_use_vec),FUN = function(i){
  return(c(1000,600))
}))

# represent different dropout levels: 20% 30% 40% 50% 60%
mu_design <- matrix(c(1.9,-0.6,0.4,-0.6,
                      1.4,-1.1,-0.1,-1.1,
                      1.1,-1.4,-0.4,-1.4,
                      0.6,-1.9,-0.9,-1.9,
                      0.2,-2.3,-1.3,-2.3),ncol = 4,byrow = T)

#network type
network_type=c("ER-graph","AFF-graph","PA-graph","banded-graph")

##different mixing degree of component evaluated by ARI
ARI_deisgn<-matrix(c(0.9,1,
                     0.75,0.85,
                     0.6,0.7),ncol = 2,byrow = TRUE)

##Initialized clustering method
cluster_method_vec<-c("Kmeans","SNN")

## The density of network and the element of precision matrix
#densy_degree<-0.1
v<-0.3
##number of group
group_num<-3

#time of experiment
reptime<-20
#library size estimation
library_size_est_method<-"TSS"
#####################################################################

##4.2 choose the specific setting for testing
#####################################################################
EMMMLE_omega<-0.8

args <- commandArgs(trailingOnly = TRUE)
param1 <- as.numeric(args[1])
param2 <- as.numeric(args[2])
param3 <- as.numeric(args[3])
param4 <- as.numeric(args[4])
param5 <- as.numeric(args[5])
run_index<-c(param1,param2,param3,param4,param5)
s_i<-run_index[1]
p_i<-run_index[2]
m_i<-run_index[3]
net<-run_index[4]
ARI_i<-run_index[5]
clustermethod_i<-1

dim_use<-dim_use_vec[p_i]
densy_degree<-0.016
#densy_degree<-(3*dim_use -6)/(dim_use * (dim_use - 1)/2)
#####################################################################

##4.3 parallel setting (Optional)
#####################################################################
num_core<-16
cl <- makeCluster(num_core, outfile = paste("debug724_",paste(run_index,collapse = ""),".txt",sep = ""))
registerDoParallel(cl)
#combine function
cfun<-function(a,b){
  rbind(a,b)
}
#####################################################################


##5. Run
##--------------------------------------------------------------------------------------------

##5.1 simulation for data
#####################################################################
parameter_setting_fix<-list(p_v=0.5,v=0.3,
                            clustermethod_i=1,
                            densy_degree=densy_degree,
                            generator_type="MPLN")

data_simulation(run_index=run_index,reptime=reptime,
                parameter_setting_fix=parameter_setting_fix,
                dim_use_vec=dim_use_vec,
                group_num=group_num,
                sample_size_mat=sample_size_mat,
                cluster_method_vec=cluster_method_vec,
                mu_design=mu_design,
                network_type=network_type,
                ARI_deisgn=ARI_deisgn)
##5.2 evaluate the performance for each method
#####################################################################
repres_all<-list()

##5.2.1 EMMMLE
print("EMMMLE start.")
repres_EMMMLE<-EMMMLE_weight_evaluator(penalize_diagonal = FALSE,if_stability_eval = FALSE,diagonal_if = FALSE,omega=EMMMLE_omega,initial_method="EMPLN",iter_emmle=2)
repres_all[[1]]<-repres_EMMMLE

##5.2.2 PLNet
print("PLNet start.")
repres_PLNet<-PLNet_evaluator(penalize_diagonal = FALSE,if_stability_eval = FALSE)
repres_all[[2]]<-repres_PLNet

##5.2.3 VMPLN
print("VMPLN start.")
repres_VMPLN<-VMPLN_evaluator(penalize_diagonal = FALSE,if_stability_eval = FALSE)
repres_all[[3]]<-repres_VMPLN

##5.2.4 Glasso
print("Glasso start.")
repres_Glasso<-Glasso_evaluator(penalize_diagonal = FALSE,if_stability_eval = FALSE)
repres_all[[4]]<-repres_Glasso

##5.3 save the result
#####################################################################
save(repres_all,file = paste("/home/ruibinxi_pkuhpc/lustre1/huangyf/R_project/emplnet/simulation_right/result_whole_1115/repres_all_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                             "_nettype",network_type[net],"_ARI_new",ARI_i,"_density",parameter_setting_fix$densy_degree,"_EMMMLE_omega_",EMMMLE_omega,
                             ".Rdata",sep = ""))
#####################################################################
time_whole_end<-Sys.time()
time_cost_whole<-as.double(difftime(time_whole_end,time_whole_start,units = "hours"))
save(time_cost_whole,file = paste("/home/ruibinxi_pkuhpc/lustre1/huangyf/R_project/emplnet/simulation_right/result_whole_1115/time_cost_whole",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                  "_nettype",network_type[net],"_ARI_new",ARI_i,"_density",parameter_setting_fix$densy_degree,"_EMMMLE_omega_",EMMMLE_omega,
                                  ".Rdata",sep = ""))
##--------------------------------------------------------------------------------------------
