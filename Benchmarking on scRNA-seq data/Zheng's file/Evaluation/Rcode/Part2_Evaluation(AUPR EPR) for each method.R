##Evaluation for each method
rm(list = ls())
#1. Set up the overall file path
##--------------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------------
data_source<-"Zheng"
dim_data<-300
#2. Set up the overall file path
##--------------------------------------------------------------------------------------------
File_path<-paste("./Benchmarking on scRNA-seq data/",data_source,"'s file",sep = "")
##--------------------------------------------------------------------------------------------

#2. Load data
##--------------------------------------------------------------------------------------------
load(paste(File_path,"/Evaluation/Obj_use/weight_mat_EMMMLE_densy_4group.Rdata",sep = ""))
load(paste(File_path,"/Evaluation/Obj_use/weight_mat_PLNet_densy_4group.Rdata",sep = ""))
load(paste(File_path,"/Evaluation/Obj_use/weight_mat_VMPLN_densy_4group.Rdata",sep = ""))
load(paste(File_path,"/Evaluation/Obj_use/weight_mat_Glasso_densy_4group.Rdata",sep = ""))
load(paste(File_path,"/Evaluation/Obj_use/weight_mat_GENIE3nor_densy.Rdata",sep = ""))
load(paste(File_path,"/Evaluation/Obj_use/weight_mat_ppcor_densy.Rdata",sep = ""))
load(paste(File_path,"/Evaluation/Obj_use/Silver_Standard_GRN.Rdata",sep = ""))
load(paste(File_path,"/Method/EMPLN/Output_Obj/EMPLN_init_4group.Rdata",sep = ""))

load(paste(File_path,"/Data/cluster_true_4group.Rdata",sep = ""))
##--------------------------------------------------------------------------------------------

##3. source some functions
#--------------------------------------------------------------------------------------------
source(paste(File_path,"/Evaluation/prepare_function/part_of_interest.R",sep = ""))
source(paste(File_path,"/Evaluation/prepare_function/Eval_fun.R",sep = ""))
##--------------------------------------------------------------------------------------------

##4. Description of data
##--------------------------------------------------------------------------------------------
load(paste(File_path,"/Evaluation/Obj_use/description_data_4group.Rdata",sep = ""))
celltypes_num<-as.numeric(description_data["celltypes_num"])
p_TF<-8
p_nonTF<-as.numeric(description_data["p_nonTF"])
p_GRN<-as.numeric(description_data["p_GRN"])
celltype_order<-colnames(table(EMPLN_init$celltypes_label,cluster_true))[apply(table(EMPLN_init$celltypes_label,cluster_true),MARGIN = 1,which.max)]

##5. Evaluation for each method
##--------------------------------------------------------------------------------------------
AUCAUPR_EMMMLE_densy<-Eval_fun(weight_mat_array = weight_mat_EMMMLE_densy,adjoint_true_array = Silver_Standard_GRN,p_TF = p_TF)
AUCAUPR_PLNet_densy<-Eval_fun(weight_mat_array = weight_mat_PLNet_densy,adjoint_true_array = Silver_Standard_GRN,p_TF = p_TF)
AUCAUPR_VMPLN_densy<-Eval_fun(weight_mat_array = weight_mat_VMPLN_densy,adjoint_true_array = Silver_Standard_GRN,p_TF = p_TF)
AUCAUPR_Glasso_densy<-Eval_fun(weight_mat_array = weight_mat_Glasso_densy,adjoint_true_array = Silver_Standard_GRN,p_TF = p_TF)
AUCAUPR_GENIE3_densy<-Eval_fun(weight_mat_array = weight_mat_GENIE3_densy,adjoint_true_array = Silver_Standard_GRN,p_TF = p_TF)
AUCAUPR_ppcor_densy<-Eval_fun(weight_mat_array = weight_mat_ppcor_densy,adjoint_true_array = Silver_Standard_GRN,p_TF = p_TF)

##AUPRC
AUPR_ratio_res_densy<-as.matrix(rbind(AUCAUPR_EMMMLE_densy$AUPR_ratio_vec,
                                      AUCAUPR_ppcor_densy$AUPR_ratio_vec,
                                      AUCAUPR_PLNet_densy$AUPR_ratio_vec,
                                      AUCAUPR_VMPLN_densy$AUPR_ratio_vec,
                                      AUCAUPR_Glasso_densy$AUPR_ratio_vec,
                                      AUCAUPR_GENIE3_densy$AUPR_ratio_vec))
AUPR_ratio_res_densy<-ifelse(is.nan(AUPR_ratio_res_densy),0,AUPR_ratio_res_densy)
rownames(AUPR_ratio_res_densy)<-c("EM-MMLE","PPCOR","PLNet","VMPLN","Glasso","GENIE3")
celltype_order<-c("CD4+ T cells","Monocytes", "B cells","CD8+ T cells")
colnames(AUPR_ratio_res_densy)<-celltype_order
AUPR_ratio_res_densy

#Sort the methods by their median values in decreasing order
row_medians <- apply(AUPR_ratio_res_densy, 1, median)
sorted_methods <- order(row_medians, decreasing = TRUE)
AUPR_ratio_res_densy_sorted <- AUPR_ratio_res_densy[sorted_methods, ]
print(AUPR_ratio_res_densy_sorted)

##--------------------------------------------------------------------------------------------
##order
celltype_order_new<-c("B cells","Monocytes","CD8+ T cells","CD4+ T cells")
AUPR_ratio_res_densy_sorted<-AUPR_ratio_res_densy_sorted[,celltype_order_new]

save(AUPR_ratio_res_densy_sorted,file = paste(File_path,"/Evaluation/Obj_use/AUPR_ratio_res_densy_sorted_right.Rdata",sep = ""))
