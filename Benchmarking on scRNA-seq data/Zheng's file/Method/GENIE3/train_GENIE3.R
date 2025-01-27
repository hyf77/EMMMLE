
#1. Loading package 
##--------------------------------------------------------------------------------------------
library(GENIE3)
##--------------------------------------------------------------------------------------------

#2. Set up the path for loading and saving data
##--------------------------------------------------------------------------------------------
##2.1 Overall file path
File_path<-"/home/ruibinxi_pkuhpc/lustre1/huangyf/R_project/scGeneNet/run_true/Benchmarking_on_scRNA_data/Zheng_file"
##2.2 Path for load used data
initial_load_path<-paste(File_path,"/Method/EMPLN/Output_Obj/EMPLN_init_4group.Rdata",sep = "")
gene_GRN_path<-paste(File_path,"/Data/gene_GRN_intop300HVG_4group.Rdata",sep = "")

load(initial_load_path)
load(gene_GRN_path)
##2.3 Path for save result
GENIE3_save_path<-paste(File_path,"/Method/GENIE3/Output_Obj",sep = "")
##--------------------------------------------------------------------------------------------


##3. Run the model
##--------------------------------------------------------------------------------------------
p_TF<-8
celltypes_num<-4
p_GRN<-300
p_nonTF<-p_GRN - p_TF
run_fun <- function(g){
  print(g)
  ############################
  ##
  locator_vec<-EMPLN_init$celltypes_label
  ##
  data_obs<-as.matrix(EMPLN_init$obs_mat)
  exp_count_GRN<-data_obs[,1:p_GRN]
  S_depth<-EMPLN_init$S_depth
  rm(EMPLN_init);gc()
  ##
  exp_count_GRN_group<-exp_count_GRN[which(locator_vec == g),]/S_depth[which(locator_vec == g)]
  ##
  set.seed(123) # For reproducibility of results
  regulators <- c(1:p_TF)
  weightMat <- GENIE3(t(as.matrix(exp_count_GRN_group)),regulators = regulators,nCores = 1)
  ##
  ##
  save(weightMat,file = paste(GENIE3_save_path,"/weightMatnor_res_use_top500GRN_g",g,".Rdata",sep = ""))
}
for (g in 1:celltypes_num) {
  run_fun(g)
}
##--------------------------------------------------------------------------------------------