#1. Loading package 
library(Seurat)
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

cfun<-function(a,b){
  rbind(a,b)
}
args <- commandArgs(trailingOnly = TRUE)
param1 <- as.numeric(args[1])
param2 <- as.numeric(args[2])
param3 <- as.numeric(args[3])
param4 <- as.numeric(args[4])
iter_emmle<-param1
omega<-param2
if(param3==1){
  diagonal_if <-TRUE
}else{
  diagonal_if<-FALSE
}
cut_num<-param4
method_name<-'EMMMLE'
num_core<-20
cl <- makeCluster(num_core, outfile =paste("debug_4group_",method_name,"_iter",iter_emmle,"_omega",omega,"_diagonal_if",diagonal_if,"_cut",cut_num,".txt",sep = ""))
registerDoParallel(cl)

File_path<-"/home/ruibinxi_pkuhpc/lustre1/huangyf/R_project/scGeneNet/run_true/Benchmarking_on_scRNA_data/Zheng_file"
##2.2 Path for load used data
expression_profile_path<-paste(File_path,"/Data/expression_profile_use_4group.Rdata",sep = "")
cluster_true_path<-paste(File_path,"/Data/cluster_true_4group.Rdata",sep = "")
gene_GRN_path<-paste(File_path,"/Data/gene_GRN_intop300HVG_4group.Rdata",sep = "")

##2.3 Path for save result
initial_save_path<-paste(File_path,"/Method/EMPLN/Output_Obj/EMPLN_init_4group.Rdata",sep = "")
EMMMLE_init_save_path<-paste(File_path,"/Method/EMPLN/Output_Obj/EMMMLE_init_4group_cut",cut_num,".Rdata",sep = "")
EMMMLE_result_list_save_path<-paste(File_path,"/Method/EMMMLE/Output_Obj/EMMMLE_result_list_iter",iter_emmle,"_omega",omega,"_cut",cut_num,"_4group.Rdata",sep = "")
result_save_path<-paste(File_path,"/Method/",method_name,"/Output_Obj/result_all_list_iter",iter_emmle,"_omega",omega,"_diagonal_if",diagonal_if,"_cut",cut_num,"_4group.Rdata",sep = "")
print(result_save_path)
##--------------------------------------------------------------------------------------------

##3. Load the data
##--------------------------------------------------------------------------------------------
load(expression_profile_path)
load(cluster_true_path)
load(gene_GRN_path)
##--------------------------------------------------------------------------------------------

##4. Set the parameter
##--------------------------------------------------------------------------------------------
penalize_diagonal_use<-TRUE
p_TF<-8
celltypes_num<-4
p_GRN<-length(gene_GRN)
support_use<-diag(p_GRN)
support_use[1:p_TF,]<-1
support_use[,1:p_TF]<-1
zero_mat<-exact_zeroindex(support = support_use)

##5. Initialization
EMPLN_init<-cluster_init(expression_profile = expression_profile_use,
                          celltypes_num = celltypes_num,
                          celltypes_ref = NULL,
                          ls_est = "TSS",
                          gene_GRN = gene_GRN,
                          HVG_model_num = p_GRN,
                          zero_GRN = zero_mat,
                          preprocess_Control = list(HVG_num = p_GRN,npc = 50,
                                                    run_umap = TRUE,label_umap = cluster_true,
                                                    cluster_method = "SNN",resolution = 0.8),
                          core_num = 20)

save(EMPLN_init,file = initial_save_path)

EMMMLE_init_fun<-function(res_init){
  iter_max_t<-10
  ELBO_threshold<-1e-4
  data_obs <-as.matrix(res_init$obs_mat)
  dim_obs_data<-ncol(data_obs)
  num_obs_data<-nrow(data_obs)
  group_num<-celltypes_num
  t.root<-polynomial.roots(monic.polynomial.recurrences(hermite.h.recurrences(20, normalized=FALSE)))
  t_root_vec<-t.root[[21]]
  omega_root_vec<-2^19*factorial(20)*sqrt(pi)/400/polynomial.values(hermite.h.polynomials(20,normalized=FALSE), t.root[[21]])[[20]]^2
  ## data initialization
  locator_base <- res_init$celltypes_label
  S_depth<-res_init$S_depth
  t<-1
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
  init_Zg<-cond_pro_Xdata_Zg(Xdata = data_obs,pi=res_init$pi_vec,mu=res_init$mu_mat,sigma=res_init$sigma_mat_list,scale = matrix(rep(S_depth,dim_obs_data),num_obs_data,dim_obs_data),obs_data_sp=res_init$obs_mat)
  pi_iter<-init_Zg$pi.new
  post_dist_Z_iter<-init_Zg$post_dist_Z
  locator_base1<-apply(post_dist_Z_iter,MARGIN = 1,which.max)
  mlemu_init<-matrix(NA,nrow = group_num,ncol = dim_obs_data)
  mlesigmahat_init<-list()
  for(g in 1:group_num){
    data_obs_group<-data_obs[locator_base1==g,]
    Y<-data_obs_group
    data_use<-Y
    S_depth_group<-res_init$S_depth[locator_base1==g]
    dim_use<-ncol(data_use)
    sample_size<-nrow(data_use)
    data_use_nor<-data_use/as.vector(S_depth_group)
    temp1<-colMeans(data_use/as.vector(S_depth_group))
    log_vec1<-as.vector(log(temp1))
    temp11 <- (t(data_use_nor) %*% data_use_nor) / nrow(data_use_nor)
    sigmahat<-t(log(temp11) - log_vec1) - log_vec1
    temp12 <-colMeans((data_use * (data_use - 1)) / as.vector(S_depth_group^2))
    diag(sigmahat)<-(log(temp12) - 2 * log(temp1))
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
    temp_mu0<-colMeans(data_use/as.vector(S_depth_group))
    mu<-log(temp_mu0) - diag(sigmahat)/2
    temp_num_mu<-min(mu[is.finite(mu)])
    mu[is.infinite(mu)| is.na(mu)]<-temp_num_mu
    mlemu_init[g,]<-mu
    mlesigmahat_init[[g]]<-sigmahat
    rm(mu);gc()
    rm(sigmahat);gc()
  }
  updated_mu<-matrix(0,nrow = group_num,ncol = dim_obs_data)
  updated_sigma<-list()
  if_NaN_2<-FALSE
  positions <- expand.grid(row = 1:dim_obs_data, col = 1:dim_obs_data)
  position_input <- positions[positions$row > positions$col, ]
  position_input<-as.matrix(position_input)
  nonpositive_sigma<-list()
  results <- foreach(g = 1:group_num, .export = ls(.GlobalEnv),.packages = c('MASS','CVXR','mclust','mltools','orthopolynom','Rcpp','RcppArmadillo','doParallel',"Seurat","EMMMLE")) %dopar% {
    print(paste("group:",g))
    k_max<-10
    cond_Z_input<-post_dist_Z_iter[,g]
    res_init_mu<-mlemu_init[g,]
    res_init_sigma<-mlesigmahat_init[[g]]
    mmle_newton_res<-mmle_newton_fun(data_obs,cond_Z_input,S_depth,res_init_mu,res_init_sigma,t_root_vec,omega_root_vec,k_max,20,position_input,cut_num)
    updated_mu_temp<-mmle_newton_res$mlemu
    updated_sigma_temp<-mmle_newton_res$mlesigmahat
    list(mu = updated_mu_temp, sigma = updated_sigma_temp)
  }
  
  for(g in 1:group_num){
    updated_mu[g,]<-as.vector(results[[g]]$mu)
    updated_sigma[[g]]<-results[[g]]$sigma
    nonpositive_sigma[[g]]<-results[[g]]$sigma
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
  
  parameter_init_save<-list(pi_vec=pi_iter,mu_mat=updated_mu,sigma_mat_list=updated_sigma,post_dist_Z_iter=post_dist_Z_iter,nonpositive_sigma=nonpositive_sigma)
  save(parameter_init_save,file=EMMMLE_init_save_path)
  
}
EMMMLE_init_fun(EMPLN_init)

##6. Run
load(file=EMMMLE_init_save_path)
data_use<-as.matrix(EMPLN_init$obs_mat)
res_init<-EMPLN_init
cluster_true<-as.numeric(as.factor(cluster_true))
S_depth<-res_init$S_depth
dim_use<-ncol(data_use)
sample_size<-nrow(data_use)
group_num<-celltypes_num
## data initialization-----------------------------------------------------
locator_base <- res_init$celltypes_label
ARI_orgin <- adjustedRandIndex(cluster_true,locator_base)
EMPLN_moment_update_res<-list(pi_vec=parameter_init_save$pi_vec,mu_mat=parameter_init_save$mu_mat,post_dist_Z=parameter_init_save$post_dist_Z_iter,sigma_list=parameter_init_save$nonpositive_sigma)
post_dist_Z_mat_max<-ifelse(EMPLN_moment_update_res$post_dist_Z > 0.9, 1, 0)
soft_cs<-mean(as.vector(abs(EMPLN_moment_update_res$post_dist_Z - post_dist_Z_mat_max)))
time_cost<-c()
init_update_res<-EMPLN_moment_update_res
#omega<-EmPLNet2_omega
##############################################
locator_cur<-apply(init_update_res$post_dist_Z,MARGIN = 1,which.max)
cluster_predict_true_cur<-as.vector(apply(table(locator_cur,cluster_true),MARGIN = 1,which.max))
t.root<-polynomial.roots(monic.polynomial.recurrences(hermite.h.recurrences(20, normalized=FALSE)))
t_root_vec<-t.root[[21]]
omega_root_vec<-2^19*factorial(20)*sqrt(pi)/400/polynomial.values(hermite.h.polynomials(20,normalized=FALSE), t.root[[21]])[[20]]^2
data_obs<-as.matrix(EMPLN_init$obs_mat)
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
EMMMLE_result_list<-EMMMLE(data_obs,pi_init,mu_initial_all,sigma_inital_all,S_depth, t_root_vec, omega_root_vec, p_p_post_dist_Z,k_max, position_input,20,omega,iter_emmle,cut_num)
save(EMMMLE_result_list,file = EMMMLE_result_list_save_path)
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
##
n_lambda<-150
Theta_mat_array_hat_all<-array(NA,dim = c(dim_use,dim_use,group_num,n_lambda))
BIC_mat_shrink <- matrix(NA,nrow = group_num,ncol = n_lambda)
locator_cur<-apply(init_update_res$post_dist_Z,MARGIN = 1,which.max)
cluster_predict_true_cur<-as.vector(apply(table(locator_cur,cluster_true),MARGIN = 1,which.max))
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
                           penalize.diagonal = penalize_diagonal,
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

result_all_list = list(BIC_mat_shrink = BIC_mat_shrink,
                       lambda_mat = lambda_mat,
                       U_mat=EMPLN_moment_update_res$post_dist_Z,
                       cluster_predict_true_vec = cluster_predict_true_vec,
                       Theta_mat_array_hat_all = Theta_mat_array_hat_all,
                       time_use=time_use)
save(result_all_list,file=result_save_path)