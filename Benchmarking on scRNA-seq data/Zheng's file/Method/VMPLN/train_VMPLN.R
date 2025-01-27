#1. Loading package 
##--------------------------------------------------------------------------------------------
library(Rcpp)
library(Seurat)
library(glasso)
library(mclust)
library(scGeneNet)
library(doParallel)
method_name<-'scGeneNet'
num_core<-20
cl <- makeCluster(num_core, outfile =paste("debug_",method_name,".txt",sep = ""))
registerDoParallel(cl)
rm.outlier <- function(x,percent){
  tmp <- (x>=quantile(x,percent))
  Q13 <- quantile(x[tmp],c(0.25,0.75))
  if (Q13[2]==Q13[1]) Q13[2] <- Q13[1]+1
  upper <- max(Q13[2]*3-2*Q13[1],1)
  lower <- max(Q13[1]*3-2*Q13[2],0)
  outlier.index <- tmp & (x>upper | x<lower)
  mean.noout <- mean(x[tmp & (!outlier.index)])
  x[outlier.index] <- round(mean.noout)
  cat(sum(outlier.index),"\n")
  x
}
##
##--------------------------------------------------------------------------------------------

#2. Set up the path for loading and saving data
##--------------------------------------------------------------------------------------------
##2.1 Overall file path
File_path<-"/home/ruibinxi_pkuhpc/lustre1/huangyf/R_project/scGeneNet/run_true/Benchmarking_on_scRNA_data/Zheng_file"

##2.2 Path for load used data
expression_profile_path<-paste(File_path,"/Data/expression_profile_use_4group.Rdata",sep = "")
cluster_true_path<-paste(File_path,"/Data/cluster_true_4group.Rdata",sep = "")
gene_GRN_path<-paste(File_path,"/Data/gene_GRN_intop300HVG_4group.Rdata",sep = "")

##2.3 Path for save result
initial_save_path<-paste(File_path,"/Method/VMPLN/Output_Obj/scGeneNet_list_initial_4group.Rdata",sep = "")
initial_optional_save_path<-paste(File_path,"/Method/VMPLN/Output_Obj/scGeneNet_list_initial_optional_4group.Rdata",sep = "")
lambda_vec_path<-paste(File_path,"/Method/VMPLN/Output_Obj/lambda_vec_4group.Rdata",sep = "")
result_save_path<-paste(File_path,"/Method/VMPLN/Output_Obj/scGeneNet_result_4group.Rdata",sep = "")

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
##--------------------------------------------------------------------------------------------

##5. Initialization
##--------------------------------------------------------------------------------------------
##5.1 Basic initialization
scGeneNet_list<-scGeneNet_init(expression_profile = expression_profile_use,
                       celltypes_num = celltypes_num,
                       celltypes_ref = NULL,
                       ls_est = "TSS",
                       gene_GRN = gene_GRN,
                       HVG_model_num = p_GRN,
                       zero_GRN = zero_mat,
                       preprocess_Control = list(HVG_num = p_GRN,npc = 50,run_umap = TRUE,label_umap = cluster_true, cluster_method = "SNN"),
                       core_num = 20)

save(scGeneNet_list,file = initial_save_path)

##6. Set the path of hyper-parameter lambda
##--------------------------------------------------------------------------------------------
lambda_max_res <- exact_lambdamax(scGeneNet_list = scGeneNet_list,
                                  Theta_Control = list(penalize_diagonal = penalize_diagonal_use),
                                  core_num = 20)
##
lambda_max_uni<-sort(unique(lambda_max_res))
lambda_min_res<-1e-6
lambda_length<-20
if(length(lambda_max_uni) == 1){
  lambda_vec<-c(seq(from=(lambda_min_res),to=(lambda_max_uni[1]/2),length.out = lambda_length * 0.6),
                seq(from=(lambda_max_uni[1]/2),to=(lambda_max_uni[1]),length.out = lambda_length * 0.4))
}else{
  lambda_max_uni<-c(lambda_min_res,lambda_max_uni)
  lambda_vec<-c()
  for(k in 1:(length(lambda_max_uni) - 1)){
    if(k == 1){
      lambda_vec<-c(lambda_vec,seq(from=(lambda_max_uni[k]),to=(lambda_max_uni[k+1]/2),length.out = lambda_length * 0.6),
                    seq(from=(lambda_max_uni[k+1]/2),to=(lambda_max_uni[k+1]),length.out = lambda_length * 0.4))
    }else{
      lambda_vec<-c(lambda_vec,seq(from=(lambda_max_uni[k]),to=(lambda_max_uni[k+1]/2),length.out = lambda_length * 0.6),
                    seq(from=(lambda_max_uni[k+1]/2),to=(lambda_max_uni[k+1]),length.out = lambda_length * 0.4))
    }
  }
}
print(paste("The length of lambda_path: ",length(lambda_vec),sep = ""))
save(lambda_vec,file = lambda_vec_path)
##--------------------------------------------------------------------------------------------

##7. Run the scGeneNet for each specific lambda
##--------------------------------------------------------------------------------------------
p_GRN<-length(scGeneNet_list$gene_GRN)
##
BIC_mat_shrink = matrix(NA,nrow = celltypes_num,ncol = length(lambda_vec))
BIC_mat_shrink1 = matrix(NA,nrow = celltypes_num,ncol = length(lambda_vec))
cluster_acc = rep(NA,length(lambda_vec))
U_mat_array<-array(NA,dim = c(length(lambda_vec),nrow(scGeneNet_list$obs_mat),celltypes_num))
Theta_mat_array_hat_all<-array(NA,dim = c(p_GRN,p_GRN,celltypes_num,length(lambda_vec)))
##
if_separate<-FALSE
if(if_separate){
  group_lambda_index<-1
  ##
  num_lambda<-10
  lambda_vec_index<-c(((group_lambda_index - 1) * num_lambda + 1) : ((group_lambda_index) * num_lambda)) 
}else{
  lambda_vec_index<-1:length(lambda_vec)
}
##
for (l in lambda_vec_index) {
  print(paste("lambda(",l,") = ",lambda_vec[l],sep = ""))
  scGeneNet_list_cur<-scGeneNet_main(scGeneNet_list = scGeneNet_list,lambda_use = lambda_vec[l],
                             Global_Control = list(ELBO_threshold = 1e-4),
                             Theta_Control = list(penalize_diagonal = penalize_diagonal_use),
                             verbose = TRUE,
                             U_fix = FALSE,
                             core_num = 20)
  #########################################################
  ## Some quantities of interest
  ## BIC
  #########################################################
  BIC_mat_shrink[,l]<-scGeneNet_list_cur$scGeneNet_bic_VMICL
  BIC_mat_shrink1[,l]<-scGeneNet_list_cur$scGeneNet_bic_VICL
  #########################################################
  
  ## Theta_mat
  #########################################################
  for(g in 1:celltypes_num){
    if(penalize_diagonal_use == FALSE){
      Theta_mat_array_hat_all[,,g,l]<-as.matrix(scGeneNet_list_cur$Theta_mat_list[[g]])[1:p_GRN,1:p_GRN]
    }else{
      Theta_mat_array_hat_all[,,g,l]<-as.matrix(scGeneNet_list_cur$Theta_mat_list_pd[[g]])[1:p_GRN,1:p_GRN]
    }
  }
  #########################################################
  
  # U_mat
  #########################################################
  U_mat_array[l,,]<-scGeneNet_list_cur$U_mat
  #########################################################
  
  ##ARI (Optional, if "cluster_true" is given.)
  #########################################################
  if(!is.null(cluster_true)){
    cluster_acc0<-adjustedRandIndex(scGeneNet_list_cur$celltypes_label,cluster_true)
    cluster_acc[l]<-cluster_acc0
    print(paste("ARI: ",round(cluster_acc0,3),sep = ""))
  }
  #########################################################
  
}

##
result_all_list = list(BIC_mat_shrink = BIC_mat_shrink,
                       BIC_mat_shrink1 = BIC_mat_shrink1,
                       cluster_acc = cluster_acc,
                       U_mat_array = U_mat_array,
                       Theta_mat_array_hat_all = Theta_mat_array_hat_all)
##--------------------------------------------------------------------------------------------

##8. Save the result
##--------------------------------------------------------------------------------------------
if(if_separate){
  save(result_all_list,file = paste(File_path,"/Method/VMPLN/Output_Obj/result_all_list_",group_lambda_index,".Rdata",sep = ""))
}else{
  save(result_all_list,file = result_save_path)
}
##--------------------------------------------------------------------------------------------


