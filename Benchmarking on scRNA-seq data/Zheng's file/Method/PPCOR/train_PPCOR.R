#1. Loading package 
##--------------------------------------------------------------------------------------------
library(ppcor)
##--------------------------------------------------------------------------------------------

#2. Set up the path for loading and saving data
##--------------------------------------------------------------------------------------------
##2.1 Overall file path
File_path<-"/home/ruibinxi_pkuhpc/lustre1/huangyf/R_project/scGeneNet/run_true/Benchmarking_on_scRNA_data/Zheng_file"
initial_load_path<-paste(File_path,"/Method/EMPLN/Output_Obj/EMPLN_init_4group.Rdata",sep = "")
load(initial_load_path)
##2.2 Path for load used data

##2.3 Path for save result
PPCOR_save_path<-paste(File_path,"/Method/PPCOR/Output_Obj",sep = "")

##--------------------------------------------------------------------------------------------

##3. Run the model
##--------------------------------------------------------------------------------------------
p_TF<-8
p_GRN<-300
p_nonTF<-p_GRN - p_TF

run_fun <- function(g){
  print(g)
  ###########################
  
  ##
  obs_mat_oral<-as.matrix(EMPLN_init$obs_mat)
  locator_vec<-EMPLN_init$celltypes_label
  
  ##
  exp_count_GRN<-obs_mat_oral[,1:p_GRN]
  S_depth<-EMPLN_init$S_depth
  ##
  exp_nor<-log((as.matrix(exp_count_GRN)+1)/as.vector(S_depth))
  ####
  exp_nor_group<-exp_nor[which(locator_vec == g),]
  ##
  ##pearson
  pcor_res<-pcor(exp_nor_group,method = "pearson")
  ##
  pcor_res_use<-pcor_res$estimate
  pvalue_res_use<-pcor_res$p.value
  
  ##
  save(pcor_res_use,file = paste(PPCOR_save_path,"/pearson_pcor_res_use_top500GRN_g",g,".Rdata",sep = ""))
  save(pvalue_res_use,file = paste(PPCOR_save_path,"/pearson_pvalue_res_use_top500GRN_g",g,".Rdata",sep = ""))
  ##
  
}
for (gorup_i in 1:4) {
  run_fun(gorup_i)
}
##--------------------------------------------------------------------------------------------