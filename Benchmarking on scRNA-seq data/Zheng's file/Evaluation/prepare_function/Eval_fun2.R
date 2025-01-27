part_of_interest<-function(mat_use,p_TF){
  return(c((mat_use[2:p_TF,2:p_TF])[upper.tri((mat_use[2:p_TF,2:p_TF]))],as.vector(mat_use[2:p_TF,-(2:p_TF)])))
}

##Evaluation for each method and each celltype
Eval_fun<-function(weight_mat_array,
                      adjoint_true_array,
                      p_TF){
  
  #########################
  AUPR_vec<-rep(NA,dim(weight_mat_array)[1])
  AUPR_ratio_vec<-rep(NA,dim(weight_mat_array)[1])
  AUC_vec<-rep(NA,dim(weight_mat_array)[1])
  
  EPR_vec<-rep(NA,dim(weight_mat_array)[1])
  EPRR_vec<-rep(NA,dim(weight_mat_array)[1])
  ETPR_vec<-rep(NA,dim(weight_mat_array)[1])
  
  OPR_vec<-rep(NA,dim(weight_mat_array)[1])
  OPRR_vec<-rep(NA,dim(weight_mat_array)[1])
  OTPR_vec<-rep(NA,dim(weight_mat_array)[1])
  
  for(g in 1:length(AUPR_vec)){
    weight_pre0<-abs(weight_mat_array[g,,])
    
    weight_pre<-ifelse(weight_pre0>t(weight_pre0),weight_pre0,t(weight_pre0))
    weight_pre_upper<-part_of_interest(weight_pre,p_TF)
    ##
    adjoin_true<-adjoint_true_array[g,,]
    adjoin_true_upper<-part_of_interest(adjoin_true,p_TF)
    ##
    edge_num_pre<-length(which(weight_pre_upper!=0))
    edge_num_true<-length(which(adjoin_true_upper!=0))
    edge_density_true<-edge_num_true/length(adjoin_true_upper)
    edge_order_pre<-order(weight_pre_upper,decreasing = TRUE)[1:edge_num_pre]
    edge_true<-which(adjoin_true_upper!=0)
    ##
    ##FPR TPR PRE
    FPR_vec<-c()
    TPR_vec<-c()
    PRE_vec<-c()
    for(l in 1:length(edge_order_pre)){
      PRE_vec<-c(PRE_vec,length(which(adjoin_true_upper[edge_order_pre[1:l]] != 0))/l)
      TPR_vec<-c(TPR_vec,length(which(adjoin_true_upper[edge_order_pre[1:l]] != 0))/edge_num_true)
      FPR_vec<-c(FPR_vec,length(which(adjoin_true_upper[edge_order_pre[1:l]] == 0))/length(which(adjoin_true_upper == 0)))
    }
    
    AUROC<-sum(diff(FPR_vec) * (TPR_vec[-1] + TPR_vec[-length(TPR_vec)])/2)
    AUC_vec[g]<-AUROC
    
    AUPRC<-sum(diff(TPR_vec) * PRE_vec[-1])
    AUPR_vec[g]<-AUPRC
    
    AUPR_ratio<-AUPRC/((max(TPR_vec) - min(TPR_vec)) * edge_density_true)
    AUPR_ratio_vec[g]<-AUPR_ratio
    
    ##################
    edge_num_choose<-min(edge_num_true,edge_num_pre)
    EP<-length(which(adjoin_true_upper[edge_order_pre[1:edge_num_choose]] !=0))/edge_num_choose
    EPR_vec[g]<-EP
    EP_random<-edge_density_true
    EP_ratio<-EP/EP_random
    EPRR_vec[g]<-EP_ratio
    ETPR<-length(which(adjoin_true_upper[edge_order_pre[1:edge_num_choose]] != 0))/edge_num_true
    ETPR_vec[g]<-ETPR
    
    ################
    edge_num_choose<-edge_num_pre
    OP<-length(which(adjoin_true_upper[edge_order_pre[1:edge_num_choose]] !=0))/edge_num_choose
    OPR_vec[g]<-OP
    OP_random<-edge_density_true
    OP_ratio<-OP/OP_random
    OPRR_vec[g]<-OP_ratio
    OTPR<-length(which(adjoin_true_upper[edge_order_pre[1:edge_num_choose]] != 0))/edge_num_true
    OTPR_vec[g]<-OTPR
  }
  
  #############
  return(list(AUC_vec = AUC_vec,
              AUPR_vec = AUPR_vec,
              AUPR_ratio_vec = AUPR_ratio_vec,

              EPR_vec = EPR_vec,
              EPRR_vec = EPRR_vec,
              ETPR_vec = ETPR_vec,
              
              OPR_vec = OPR_vec,
              OPRR_vec = OPRR_vec,
              OTPR_vec = OTPR_vec))
  
}