##Density-specific GRN selection for each method
rm(list = ls())
#1. Set up the overall file path
##--------------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------------
data_source<-"Zheng"
dim_data<-300
#2. Set up the overall file path
File_path<-paste("./Benchmarking on scRNA-seq data/",data_source,"'s file",sep = "")
load(paste(File_path,"/Data/gene_GRN_intop",dim_data,"HVG_4group.Rdata",sep = ""))
##--------------------------------------------------------------------------------------------

##2. Description of data
##--------------------------------------------------------------------------------------------
load(paste(File_path,"/Evaluation/Obj_use/description_data_4group.Rdata",sep = ""))
celltypes_num<-as.numeric(description_data["celltypes_num"])
p_TF<-as.numeric(description_data["p_TF"])
p_TF<-8
p_nonTF<-as.numeric(description_data["p_nonTF"])
p_GRN<-as.numeric(description_data["p_GRN"])
num_support<-p_TF * (p_TF-1)/2 + p_TF * p_nonTF
##--------------------------------------------------------------------------------------------

##3. set the selected density
##-------------------------------------------------------------------------------------------
density_use<-0.05
densy_choose_vec<-rep(density_use,celltypes_num)
##--------------------------------------------------------------------------------------------
index_gene0<-1:p_TF
##4. Exact the confidence matrix of each cell-type for each method

# ##4.1 EMMMLE
load(paste(File_path,"/Method/EMMMLE/Output_Obj/result_all_list_iter3_omega0.8_diagonal_ifFALSE_cut50_4group.Rdata",sep = ""))
load(paste(File_path,"/Method/EMPLN/Output_Obj/EMPLN_init_4group.Rdata",sep = ""))
result_EMPLN_p_list<-result_all_list
Theta_EMPLN_p<-array(NA,dim = c(p_GRN,p_GRN,celltypes_num,dim(result_EMPLN_p_list$Theta_mat_array_hat_all)[4]))
##4.1.1 BIC for each cell-type
BIC_mat_shrink<-result_EMPLN_p_list$BIC_mat_shrink
BIC_choose_EMPLN_p<-apply(BIC_mat_shrink,MARGIN = 1,which.min)

for (g in 1:celltypes_num) {
  dim_index<-which(colnames(EMPLN_init$obs_mat) %in% gene_GRN)
  Theta_EMPLN_p[,,g,]<-result_EMPLN_p_list$Theta_mat_array_hat_all[dim_index,dim_index,g,]
}
##4.1.2 Summary the density of each candiated GRN
densy_mat_EMPLN_p<-matrix(NA,nrow = celltypes_num,ncol = dim(Theta_EMPLN_p)[4])
for(g in 1:celltypes_num){
  for(l in 1:dim(Theta_EMPLN_p)[4]){
    densy_mat_EMPLN_p[g,l] <- length(which((c((Theta_EMPLN_p[index_gene0,index_gene0,g,l])[upper.tri(Theta_EMPLN_p[index_gene0,index_gene0,g,l])],as.vector(Theta_EMPLN_p[index_gene0,-(index_gene0),g,l]))) !=0))/num_support
  }
}

##4.1.3 choose the lambda with specific density
densy_choose_EMPLN_p<-c()
for(g in 1:celltypes_num){
  densy_choose_EMPLN_p<-c(densy_choose_EMPLN_p,(which(densy_mat_EMPLN_p[g,] - densy_choose_vec[g]>0))[which.min((abs(densy_mat_EMPLN_p[g,] - densy_choose_vec[g]))[which(densy_mat_EMPLN_p[g,] - densy_choose_vec[g]>0)])])
}
##4.1.4 The result density-specific confidence matrix
weight_mat_EMMMLE_densy<-array(NA,dim = c(celltypes_num,p_GRN,p_GRN))
num_edge_choose<-floor(num_support * density_use)
for(g in 1:celltypes_num){
  weight_densy_temp<-(-1) * (Theta_EMPLN_p[,,g,densy_choose_EMPLN_p[g]])/(matrix(sqrt(diag((Theta_EMPLN_p[,,g,densy_choose_EMPLN_p[g]]))),ncol = 1) %*% matrix(sqrt(diag((Theta_EMPLN_p[,,g,densy_choose_EMPLN_p[g]]))),nrow = 1))
  diag(weight_densy_temp)<-1
  weight_densy<-diag(ncol(weight_densy_temp))
  weight_densy[index_gene0,index_gene0][upper.tri(weight_densy[index_gene0,index_gene0])]<-weight_densy_temp[index_gene0,index_gene0][upper.tri(weight_densy_temp[index_gene0,index_gene0])]
  weight_densy[index_gene0,-(index_gene0)]<-weight_densy_temp[index_gene0,-(index_gene0)]
  choose_mat<-diag(ncol(weight_densy))
  choose_mat[upper.tri(choose_mat)][order(abs(weight_densy[upper.tri(weight_densy)]),decreasing = TRUE)[1:num_edge_choose]]<-1
  choose_mat<-choose_mat + t(choose_mat) - diag(ncol(weight_densy))
  ##
  weight_densy_choose<-weight_densy * choose_mat
  weight_mat_EMMMLE_densy[g,,]<-weight_densy_choose
}

##4.1.5 Save
save(weight_mat_EMMMLE_densy,file = paste(File_path,"/Evaluation/Obj_use/weight_mat_EMMMLE_densy_4group.Rdata",sep = ""))

##4.2 PLNet
num_edge_choose<-floor(num_support * density_use)
weight_mat_PLNet_densy<-array(NA,dim = c(celltypes_num,p_GRN,p_GRN))
for(g in 1:celltypes_num){
  # load(paste(File_path,"/Method/PLNet/Output_Obj/BIC_mat_shrink_g",g,".Rdata",sep = ""))
  load(paste(File_path,"/Method/PLNet/Output_Obj/Theta_mat_array_hat_all_g",g,"_4group.Rdata",sep = ""))
  densy_vec_PLNet<-c()
  for(l in 1:dim(Theta_mat_array_hat_all)[3]){
    densy_vec_PLNet <- c(densy_vec_PLNet,length(which((c((Theta_mat_array_hat_all[index_gene0,index_gene0,l])[upper.tri(Theta_mat_array_hat_all[index_gene0,index_gene0,l])],as.vector(Theta_mat_array_hat_all[index_gene0,-(index_gene0),l]))) !=0))/length((c((Theta_mat_array_hat_all[index_gene0,index_gene0,l])[upper.tri(Theta_mat_array_hat_all[index_gene0,index_gene0,l])],as.vector(Theta_mat_array_hat_all[index_gene0,-(index_gene0),l])))))
  }
  densy_choose_PLNet<-(which(densy_vec_PLNet - densy_choose_vec[g]>0))[which.min((abs(densy_vec_PLNet - densy_choose_vec[g]))[which(densy_vec_PLNet - densy_choose_vec[g]>0)])]
  ##
  weight_densy_temp<-(-1) * Theta_mat_array_hat_all[,,densy_choose_PLNet]/(matrix(sqrt(diag(Theta_mat_array_hat_all[,,densy_choose_PLNet])),ncol = 1) %*% matrix(sqrt(diag(Theta_mat_array_hat_all[,,densy_choose_PLNet])),nrow = 1))
  diag(weight_densy_temp)<-1
  
  weight_densy<-diag(ncol(weight_densy_temp))
  weight_densy[index_gene0,index_gene0][upper.tri(weight_densy[index_gene0,index_gene0])]<-weight_densy_temp[index_gene0,index_gene0][upper.tri(weight_densy_temp[index_gene0,index_gene0])]
  weight_densy[index_gene0,-(index_gene0)]<-weight_densy_temp[index_gene0,-(index_gene0)]
  ##
  choose_mat<-diag(ncol(weight_densy))
  choose_mat[upper.tri(choose_mat)][order(abs(weight_densy[upper.tri(weight_densy)]),decreasing = TRUE)[1:num_edge_choose]]<-1
  
  choose_mat[upper.tri(choose_mat)][order(abs(weight_densy[upper.tri(weight_densy)]),decreasing = TRUE)[1:num_edge_choose]]<-1
  choose_mat<-choose_mat + t(choose_mat) - diag(ncol(weight_densy))
  ##
  weight_densy_choose<-weight_densy * choose_mat
  weight_mat_PLNet_densy[g,,]<-weight_densy_choose
}
save(weight_mat_PLNet_densy,file = paste(File_path,"/Evaluation/Obj_use/weight_mat_PLNet_densy_4group.Rdata",sep = ""))

##--------------------------------------------------------------------------------------------
##4.1 VMPLN
load(paste(File_path,"/Method/VMPLN/Output_Obj/VMPLN_result_4group.Rdata",sep = ""))
result_VMPLN_list<-result_all_list
##4.1.1 BIC for each cell-type
BIC_mat_shrink<-result_VMPLN_list$BIC_mat_shrink1
BIC_choose_VMPLN<-apply(BIC_mat_shrink,MARGIN = 1,which.min)

##4.1.2 Summary the density of each candiated GRN
densy_mat_VMPLN<-matrix(NA,nrow = celltypes_num,ncol = dim(result_VMPLN_list$Theta_mat_array_hat_all)[4])
for(g in 1:celltypes_num){
  for(l in 1:dim(result_VMPLN_list$Theta_mat_array_hat_all)[4]){
    densy_mat_VMPLN[g,l] <- length(which((c((result_VMPLN_list$Theta_mat_array_hat_all[index_gene0,index_gene0,g,l])[upper.tri(result_VMPLN_list$Theta_mat_array_hat_all[index_gene0,index_gene0,g,l])],as.vector(result_VMPLN_list$Theta_mat_array_hat_all[index_gene0,-(index_gene0),g,l]))) !=0))/num_support
  }
}

##4.1.3 choose the lambda with specific density
densy_choose_VMPLN<-c()
for(g in 1:celltypes_num){
  densy_choose_VMPLN<-c(densy_choose_VMPLN,(which(densy_mat_VMPLN[g,] - densy_choose_vec[g]>0))[which.min((abs(densy_mat_VMPLN[g,] - densy_choose_vec[g]))[which(densy_mat_VMPLN[g,] - densy_choose_vec[g]>0)])])
}

##4.1.4 The result density-specific confidence matrix
weight_mat_VMPLN_densy<-array(NA,dim = c(celltypes_num,p_GRN,p_GRN))
num_edge_choose<-floor(num_support * density_use)
for(g in 1:celltypes_num){
  weight_densy_temp<-(-1) * (result_VMPLN_list$Theta_mat_array_hat_all[,,g,densy_choose_VMPLN[g]])/(matrix(sqrt(diag((result_VMPLN_list$Theta_mat_array_hat_all[,,g,densy_choose_VMPLN[g]]))),ncol = 1) %*% matrix(sqrt(diag((result_VMPLN_list$Theta_mat_array_hat_all[,,g,densy_choose_VMPLN[g]]))),nrow = 1))
  diag(weight_densy_temp)<-1
  weight_densy<-diag(ncol(weight_densy_temp))
  weight_densy[index_gene0,index_gene0][upper.tri(weight_densy[index_gene0,index_gene0])]<-weight_densy_temp[index_gene0,index_gene0][upper.tri(weight_densy_temp[index_gene0,index_gene0])]
  weight_densy[index_gene0,-(index_gene0)]<-weight_densy_temp[index_gene0,-(index_gene0)]
  ##
  choose_mat<-diag(ncol(weight_densy))
  choose_mat[upper.tri(choose_mat)][order(abs(weight_densy[upper.tri(weight_densy)]),decreasing = TRUE)[1:num_edge_choose]]<-1
  choose_mat<-choose_mat + t(choose_mat) - diag(ncol(weight_densy))
  ##
  weight_densy_choose<-weight_densy * choose_mat
  weight_mat_VMPLN_densy[g,,]<-weight_densy_choose
}

##4.1.5 Save
save(weight_mat_VMPLN_densy,file = paste(File_path,"/Evaluation/Obj_use/weight_mat_VMPLN_densy_4group.Rdata",sep = ""))

##--------------------------------------------------------------------------------------------
##4.3 Glasso
num_edge_choose<-floor(num_support * density_use)
weight_mat_Glasso_densy<-array(NA,dim = c(celltypes_num,p_GRN,p_GRN))
for(g in 1:celltypes_num){
  load(paste(File_path,"/Method/Glasso/Output_Obj/BIC_Glasso_vec_pen_top300GRN_g",g,".Rdata",sep = ""))
  load(paste(File_path,"/Method/Glasso/Output_Obj/graph_lasso_array_pen_top300GRN_g",g,".Rdata",sep = ""))
  densy_vec_Glasso<-c()
  for(l in 1:dim(graph_lasso_array)[3]){
    densy_vec_Glasso <- c(densy_vec_Glasso,length(which((c((graph_lasso_array[index_gene0,index_gene0,l])[upper.tri(graph_lasso_array[index_gene0,index_gene0,l])],as.vector(graph_lasso_array[index_gene0,-(index_gene0),l]))) !=0))/length((c((graph_lasso_array[index_gene0,index_gene0,l])[upper.tri(graph_lasso_array[index_gene0,index_gene0,l])],as.vector(graph_lasso_array[index_gene0,-(index_gene0),l])))))
  }
  densy_choose_Glasso<-(which(densy_vec_Glasso - densy_choose_vec[g]>0))[which.min((abs(densy_vec_Glasso - densy_choose_vec[g]))[which(densy_vec_Glasso - densy_choose_vec[g]>0)])]
  ##
  weight_densy_temp<-(-1) * graph_lasso_array[,,densy_choose_Glasso]/(matrix(sqrt(diag(graph_lasso_array[,,densy_choose_Glasso])),ncol = 1) %*% matrix(sqrt(diag(graph_lasso_array[,,densy_choose_Glasso])),nrow = 1))
  diag(weight_densy_temp)<-1
  
  weight_densy<-diag(ncol(weight_densy_temp))
  weight_densy[index_gene0,index_gene0][upper.tri(weight_densy[index_gene0,index_gene0])]<-weight_densy_temp[index_gene0,index_gene0][upper.tri(weight_densy_temp[index_gene0,index_gene0])]
  weight_densy[index_gene0,-(index_gene0)]<-weight_densy_temp[index_gene0,-(index_gene0)]
  ##
  choose_mat<-diag(ncol(weight_densy))
  choose_mat[upper.tri(choose_mat)][order(abs(weight_densy[upper.tri(weight_densy)]),decreasing = TRUE)[1:num_edge_choose]]<-1
  
  choose_mat[upper.tri(choose_mat)][order(abs(weight_densy[upper.tri(weight_densy)]),decreasing = TRUE)[1:num_edge_choose]]<-1
  choose_mat<-choose_mat + t(choose_mat) - diag(ncol(weight_densy))
  ##
  weight_densy_choose<-weight_densy * choose_mat
  weight_mat_Glasso_densy[g,,]<-weight_densy_choose
}
save(weight_mat_Glasso_densy,file = paste(File_path,"/Evaluation/Obj_use/weight_mat_Glasso_densy_4group.Rdata",sep = ""))

# ##--------------------------------------------------------------------------------------------
##4.4 GENIE3
num_edge_choose<-floor(num_support * density_use)
load(paste(File_path,"/Data/gene_GRN_intop300HVG_4group.Rdata",sep = ""))
weight_mat_GENIE3_densy<-array(NA,dim = c(celltypes_num,p_GRN,p_GRN))
for(g in 1:celltypes_num){
  load(paste(File_path,"/Method/GENIE3/Output_Obj/weightMatnor_res_use_top500GRN_g",g,".Rdata",sep = ""))
  ##reorder
  weight_GENIE3<-matrix(0,nrow = p_GRN,ncol = p_GRN)
  for(i in 1:p_TF){
    for(j in 1:p_GRN){
      weight_GENIE3[i,j]<-weightMat[which(rownames(weightMat) == gene_GRN[i]),which(colnames(weightMat) == gene_GRN[j])]
    }
  }
  weight_GENIE3<-ifelse(is.nan(weight_GENIE3),0,weight_GENIE3)
  weight_GENIE3<-ifelse(weight_GENIE3>t(weight_GENIE3),weight_GENIE3,t(weight_GENIE3))

  weight_GENIE3_trans<-weight_GENIE3
  weight_GENIE3_trans[-(1:p_TF),-(1:p_TF)]<-0
  weight_GENIE3_choose<-diag(p_GRN)
  index_choose<-order(weight_GENIE3_trans[upper.tri(weight_GENIE3_trans)],decreasing = TRUE)[1:num_edge_choose]
  weight_GENIE3_choose[upper.tri(weight_GENIE3_choose)][index_choose]<-1
  weight_GENIE3_choose<-weight_GENIE3_choose + t(weight_GENIE3_choose) - diag(p_GRN)
  #
  weight_mat_GENIE3_densy[g,,]<-weight_GENIE3_choose * weight_GENIE3
}
save(weight_mat_GENIE3_densy,file = paste(File_path,"/Evaluation/Obj_use/weight_mat_GENIE3nor_densy.Rdata",sep = ""))


##4.5 PPCOR
weight_mat_ppcor_densy<-array(NA,dim = c(celltypes_num,p_GRN,p_GRN))
for(g in 1:celltypes_num){
  load(paste(File_path,"/Method/PPCOR/Output_Obj/pearson_pcor_res_use_top500GRN_g",g,".Rdata",sep = ""))
  
  weight_ppcor_trans<-abs(pcor_res_use)
  weight_ppcor_trans[-(1:p_TF),-(1:p_TF)]<-0
  weight_ppcor_choose<-diag(p_GRN)
  index_choose<-order(weight_ppcor_trans[upper.tri(weight_ppcor_trans)],decreasing = TRUE)[1:num_edge_choose]
  weight_ppcor_choose[upper.tri(weight_ppcor_choose)][index_choose]<-1
  weight_ppcor_choose<-weight_ppcor_choose + t(weight_ppcor_choose) - diag(p_GRN)
  #
  weight_mat_ppcor_densy[g,,]<-weight_ppcor_choose * pcor_res_use
}
save(weight_mat_ppcor_densy,file = paste(File_path,"/Evaluation/Obj_use/weight_mat_ppcor_densy.Rdata",sep = ""))

