##1. Package required 
##--------------------------------------------------------------------------------------------
library(MASS)
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
library(lava)
library(latentcor)
library(mvtnorm)
library(cubature)
library(doParallel)
library(PLNetNew)
library(Matrix)
library(glasso)

##2. Define some function
##--------------------------------------------------------------------------------------------
gaussian_density <- function(x, mean, sigma) {
  dmvnorm(x, mean = mean, sigma = sigma)
}
identity_Limit<-function(observe,delta_input){
  lowerLimit<-rep(0,length(observe))
  upperLimit<-rep(0,length(observe))
  for (i in 1:length(observe)) {
    if(observe[i]==0){
      lowerLimit[i]<--Inf
      upperLimit[i]<-delta_input[i]
    }else{
      lowerLimit[i]<-delta_input[i]
      upperLimit[i]<-Inf
    }
  }
  res<-list(lowerLimit=lowerLimit,upperLimit=upperLimit)
  return(res)
}
network_generator<-function(dim_use,edge_num,t_cons = 0.15){
  Pre_mat<-diag(dim_use)
  nondiag_len<-dim_use*(dim_use-1)/2
  z1 <-matrix(runif(nondiag_len*2, 0, 1),nrow = 2,ncol=nondiag_len)
  z2 <-matrix(runif(nondiag_len*2, 0, 1),nrow = 2,ncol=nondiag_len)
  con_use<-TRUE
  iter0<-1
  while (con_use & iter0<500) {
    print(iter0)
    iter0<-iter0+1
    set.seed(iter0)
    c1<-1
    norms_temp <- apply(z1 - z2, 2, function(col) sqrt(sum(col^2)))
    success_pro<- (edge_num/nondiag_len) * exp(norms_temp/ (2 * c1))
    ber_samples <- rbinom(length(success_pro), 1, success_pro)
    edge_temp<-length(ber_samples[ber_samples==1])
    if((edge_temp>(edge_num-30))&&(edge_temp<(edge_num+30))){
      print(edge_temp)
      con_use<-FALSE
    }else{
      if(edge_temp<(edge_num-30)){
        c1<-c1+2
      }else{
        c1<-c1-2
      }
    }
  }
  pre_omega<-matrix(rep(0,dim_use*dim_use),nrow = dim_use)
  pre_omega[upper.tri(pre_omega)]<-ber_samples*t_cons
  Pre_mat<-t(pre_omega)+pre_omega
  diag(Pre_mat)<-1
  return(Pre_mat)
}
delta_est<- function(X){
  zratios <- colMeans(as.matrix(X == 0), na.rm = TRUE)
  delta_res<-stats::qnorm(zratios)
  return(delta_res)
}
data_simulation<-function(sample_num,dim_use,reptime,seed_num=123,edge_num=100,t=0.15){
  ## Generate the simulation data
  for(rep in 1:reptime){
    ##Generate data ------------------
    print(rep)
    set.seed(seed_num+rep)
    
    Pre_mat0<-network_generator(dim_use,edge_num,t)
    Pre_mat<-Pre_mat0
    mu_vec<-rep(0, dim_use)
    Sigma_scale<-diag(1/sqrt(diag(solve(Pre_mat))))%*%solve(Pre_mat)%*%diag(1/sqrt(diag(solve(Pre_mat))))
    Z_mat<-mvrnorm(n=sample_num, mu_vec, Sigma_scale)
    C_mat<-matrix(rep(runif(dim_use, -1, 1),sample_num),nrow=sample_num,ncol=dim_use,byrow = TRUE)
    X_mat <- ifelse(Z_mat > C_mat, 1, 0)
    data<-list(X_mat=X_mat,C_mat=C_mat,Z_mat=Z_mat,Sigma=Sigma_scale,mu=mu_vec,Precision=Pre_mat0)
    print("success")
    #save
    save(data,file = paste("/home/ruibinxi_pkuhpc/lustre1/huangyf/R_project/EM-MMLE/result/data_sam",sample_num,"_dim",dim_use,"_seed",seed_num,"_rep",rep,".Rdata",sep = ""))
  }
  print("Data generation complete.")
}
sigma_positive2<-function(cov_input){
  ##1. Calculate the sample covariance matrix by moment method----------------
  ################################
  
  ##3. Optimatization by minimize the finite norm--------------------
  S<-Variable(dim(cov_input)[1], dim(cov_input)[1], PSD = TRUE)
  
  obj<-Minimize(max(abs(S-cov_input)))
  
  prob <- CVXR::Problem(obj)
  
  result <- CVXR::solve(prob,solver="SCS",verbose=FALSE)
  sigma_me2 <- result$getValue(S)
  #sigma_me2<-max(abs(sigma_me2-cov_input))* diag(ncol(cov_input)) + sigma_me2
  ########################################
  
  ##
  return(sigma_me2)
}
Theta_main1<-function(obs_mat, lambda_vec = NULL, n_lambda =100,
                      penalize.diagonal = TRUE,cov_input = NULL,
                      weight_mat = NULL,zero_mat = NULL,
                      core_num = 1){
  
  ##names
  lambda_valuemax0 = 10
  lambda_valuemin0 = 1e-5
  ##
  Cellname<-rownames(obs_mat)
  Genename<-colnames(obs_mat)
  if(is.null(weight_mat)){
    weight_mat<-matrix(1,nrow = ncol(obs_mat),ncol = ncol(obs_mat))
  }
  if(is.null(zero_mat)){
    zero_mat<-matrix(0,nrow = ncol(obs_mat),ncol = ncol(obs_mat))
  }
  #################
  
  ##1. Determine if sigma is positive------------------------
  cov_Y1<-cov_input
  eigen_cov_Y1<-eigen(cov_Y1)
  if(min(eigen_cov_Y1$values)<0){
    print("12")
    cov_Y1<-sigma_positive2(cov_input = cov_Y1)
  }
  message("1. The estimation of the mean vector and covariance matrix by moment method has been completed.")
  ###################################
  
  ##2. Determine the lambda_vec--------------------------
  eigres<-eigen(cov_Y1)
  eigvalue_vec<-ifelse(eigres$values<0,0,eigres$values)
  obs_mat_s<-t(eigres$vectors%*%diag(sqrt(eigvalue_vec))*sqrt(ncol(obs_mat)-1))
  rm(eigres);
  time_1<-Sys.time()
  if(is.null(lambda_vec)){
    #find the suitable lambda_max and lambda_min
    lambda_valuemax<-lambda_valuemax0
    con_1<-TRUE
    while (con_1) {
      # print("lambda_max ok.")
      # print(lambda_valuemax)
      time_1<-Sys.time()
      aa<-as.matrix((EQUAL(X = obs_mat_s, lambda = lambda_valuemax,
                           sdiag = penalize.diagonal,
                           weight_mat = weight_mat,
                           zero_mat = zero_mat,
                           core_num = core_num,
                           err=10^(-4),
                           maxIter =100))$Omega[[1]])
      time_2<-Sys.time()
      time_2 - time_1
      nonzero_num <- length(which(aa[upper.tri(aa)]!=0))
      if(nonzero_num==0){
        lambda_valuemax<-lambda_valuemax*0.8
      }else{
        lambda_valuemax<-lambda_valuemax * 1.2
        con_1<-FALSE
      }
      rm(aa);
    }
    
    lambda_valuemin<-min(lambda_valuemin0,lambda_valuemax * (1e-2))
    con_min<-TRUE
    con_2<-TRUE
    while (con_2 | con_min) {
      # print("lambda_min ok.")
      # print(lambda_valuemin)
      aa<-as.matrix((EQUAL(X = obs_mat_s,lambda = lambda_valuemin,
                           sdiag = penalize.diagonal,
                           weight_mat = weight_mat,
                           zero_mat = zero_mat,
                           core_num = core_num,
                           err=10^(-4),
                           maxIter =100))$Omega[[1]])
      zero_num <- length(which(aa[upper.tri(aa)]==0))
      zero_num <- zero_num - sum(zero_mat[upper.tri(zero_mat)])
      if(con_2){
        if(zero_num!=0){
          lambda_valuemin<-lambda_valuemin*0.5
        }else{
          lambda_valuemin_temp<-lambda_valuemin
          lambda_valuemin<-lambda_valuemin*2
          con_2<-FALSE
        } 
      }else{
        if(zero_num!=0){
          con_min<-FALSE
        }else{
          lambda_valuemin_temp<-lambda_valuemin
          lambda_valuemin<-lambda_valuemin*2
        }
      }
      rm(aa)
    }
    lambda_valuemin<-lambda_valuemin_temp
    # lambda_vec <- exp(seq(from = log(lambda_valuemin),to=log(lambda_valuemax),length.out = n_lambda))
    lambda_vec <- seq(from = lambda_valuemin,to=lambda_valuemax,length.out = n_lambda)
  }
  time_2<-Sys.time()
  time_2 - time_1
  message("2. Lambda_vec determination has been completed.")
  ######################################
  
  ##3. Run EQUAL for sparse precision matrix estimation-----------------------
  time_EQUAL1<-Sys.time()
  EQUAL_res<-EQUAL(X = obs_mat_s,lambda = lambda_vec,
                   sdiag = penalize.diagonal,
                   weight_mat = weight_mat,
                   zero_mat = zero_mat,
                   core_num = core_num,
                   err=10^(-5),
                   maxIter =1000)
  time_EQUAL2<-Sys.time()
  time_EQUAL<-as.double(difftime(time_EQUAL2,time_EQUAL1,units = "hours"))
  rm(obs_mat_s)
  lambda_vec<-EQUAL_res$lambda
  message("3. The estimation of precision matrix has been completed.")
  ##########################################
  
  ##5. Hyper-parameter selection by BIC criterion and AIC criterion----------------------
  ##5.1 BIC criterion
  BIC_vec<-c()
  for (l in 1:length(lambda_vec)) {
    BIC_value<-nrow(obs_mat) * sqrt(sum(((EQUAL_res$Omega[[l]] %*% cov_Y1 + cov_Y1 %*% EQUAL_res$Omega[[l]])/2 - diag(nrow(cov_Y1)))^2)) + log(nrow(obs_mat)) *
      (sum(ifelse(EQUAL_res$Omega[[l]]==0,0,1)))/2
    BIC_vec<-c(BIC_vec,BIC_value)
  }
  lambda_chooseB<-lambda_vec[which.min(BIC_vec)]
  Omega_chooseB<-EQUAL_res$Omega[[which.min(BIC_vec)]]
  
  ##5.2 AIC criterion
  AIC_vec<-c()
  for (l in 1:length(lambda_vec)) {
    AIC_value<-nrow(obs_mat) * sqrt(sum(((EQUAL_res$Omega[[l]] %*% cov_Y1 + cov_Y1 %*% EQUAL_res$Omega[[l]])/2 - diag(nrow(cov_Y1)))^2)) + 2 *
      (sum(ifelse(EQUAL_res$Omega[[l]]==0,0,1)))/2
    AIC_vec<-c(AIC_vec,AIC_value)
  }
  lambda_chooseA<-lambda_vec[which.min(AIC_vec)]
  Omega_chooseA<-EQUAL_res$Omega[[which.min(AIC_vec)]]
  
  ##5.3 HBIC criterion
  HBIC_vec<-c()
  for (l in 1:length(lambda_vec)) {
    HBIC_value<-log(log(nrow(obs_mat)))*(sum(ifelse(EQUAL_res$Omega[[l]]==0,0,1)))*log(ncol(obs_mat))/(nrow(obs_mat))+tr(EQUAL_res$Omega[[l]] %*% cov_Y1) -log(det(EQUAL_res$Omega[[l]]))
    HBIC_vec<-c(HBIC_vec,HBIC_value)
  }
  lambda_chooseH<-lambda_vec[which.min(HBIC_vec)]
  Omega_chooseH<-EQUAL_res$Omega[[which.min(HBIC_vec)]]
  
  message("5. Hyper-parameter selection by BIC criterion and AIC criterion have been completed.")
  ############################################
  
  ##names
  for(l in 1:length(EQUAL_res$lambda)){
    rownames(EQUAL_res$Omega[[l]])<-Genename
    colnames(EQUAL_res$Omega[[l]])<-Genename
  }
  rownames(Omega_chooseB)<-Genename
  colnames(Omega_chooseB)<-Genename
  rownames(Omega_chooseA)<-Genename
  colnames(Omega_chooseA)<-Genename
  rownames(Omega_chooseH)<-Genename
  colnames(Omega_chooseH)<-Genename
  # rownames(Omega_est)<-Omega_est = EQUAL_res$Omega
  List_res<-list(Omega_est = EQUAL_res$Omega,
                 lambda_vec = EQUAL_res$lambda, 
                 BIC_vec = BIC_vec, Omega_chooseB = Omega_chooseB,
                 AIC_vec = AIC_vec, Omega_chooseA = Omega_chooseA,
                 HBIC_vec = HBIC_vec, Omega_chooseH = Omega_chooseH,
                 time = time_EQUAL)
  return(List_res)
}
AMLE_fun <- function(X) {
  # Number of samples
  n <- nrow(X)
  # Sample mean
  X_bar <- colMeans(X)
  # Center the data
  X_centered <- sweep(X, 2, X_bar)
  # Compute the modified sample covariance matrix
  modified_cov <- (1/n) * (t(X_centered) %*% X_centered) + (1/3)
  ##
  lambda_max0<-10
  lambda_valuemax<-lambda_max0
  con_1<-TRUE
  while (con_1) {
    aa<-as.matrix(glasso(modified_cov, rho=lambda_valuemax)$wi)
    nonzero_num <- length(which(aa[upper.tri(aa)]!=0))
    if(nonzero_num==0){
      lambda_valuemax<-lambda_valuemax*0.8
    }else{
      lambda_valuemax<-lambda_valuemax * 1.2
      con_1<-FALSE
    }
  }
  lambda_valuemin<-lambda_valuemax
  con_2<-TRUE
  while (con_2) {
    aa<-as.matrix(glasso(modified_cov, rho=lambda_valuemin)$wi)
    zero_num <- length(which(aa[upper.tri(aa)]==0))
    if(zero_num!=0){
      lambda_valuemin<-lambda_valuemin*0.8
    }else{
      con_2<-FALSE
    }
  }
  lambda_vec <- seq(from = lambda_valuemin,to=lambda_valuemax,length.out = 100)
  sample_size<-n
  lambda_vec1<-lambda_vec
  BIC_vecG<-c()
  AIC_vec<-c()
  HBIC_vec<-c()
  Pre_mat_listG<-list()
  time_1<-Sys.time()
  for (l in 1:length(lambda_vec1)) {
    model_13<-glasso(modified_cov, rho=lambda_vec1[l])
    pre_G<-model_13$wi
    Pre_mat_listG[[l]]<-pre_G
    BIC_vecG<-c(BIC_vecG,sample_size*sqrt(sum(((pre_G %*% modified_cov + modified_cov %*% pre_G)/2 - diag(nrow(modified_cov)))^2))+ log(sample_size) *
                  (sum(ifelse(pre_G==0,0,1)))/2)
    AIC_value<-sample_size * sqrt(sum(((pre_G %*% modified_cov + modified_cov %*% pre_G)/2 - diag(nrow(modified_cov)))^2)) + 2 *
      (sum(ifelse(pre_G==0,0,1)))/2
    AIC_vec<-c(AIC_vec,AIC_value)
    HBIC_value<-log(log(sample_size))*(sum(ifelse(pre_G==0,0,1)))*log(sample_size)/(sample_size)+tr(pre_G %*% modified_cov) -log(det(pre_G))
    HBIC_vec<-c(HBIC_vec,HBIC_value)
  }
  time_2<-Sys.time()
  time_glasso<-as.double(difftime(time_2,time_1,units = "mins"))
  lambda_chooseB<-lambda_vec[which.min(BIC_vecG)]
  Omega_chooseB<-Pre_mat_listG[[which.min(BIC_vecG)]]
  lambda_chooseA<-lambda_vec[which.min(AIC_vec)]
  Omega_chooseA<-Pre_mat_listG[[which.min(AIC_vec)]]
  lambda_chooseH<-lambda_vec[which.min(HBIC_vec)]
  Omega_chooseH<-Pre_mat_listG[[which.min(HBIC_vec)]]
  
  List_res<-list(Omega_est = Pre_mat_listG,
                 lambda_vec = lambda_vec, 
                 BIC_vec = BIC_vecG,Omega_chooseB = Omega_chooseB,
                 AIC_vec = AIC_vec, Omega_chooseA = Omega_chooseA,
                 HBIC_vec = HBIC_vec, Omega_chooseH = Omega_chooseH,
                 time = time_glasso)
  return(List_res)
}

##3. Simulation setting
##--------------------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
param1 <- as.numeric(args[1])
param2 <- as.numeric(args[2])
param3 <- as.numeric(args[3])
sample_vec<-c(500,1000,3000)
dimm_vec<-c(200,500,1000)
seed_num_vec<-c(1,2)
sample_num<-sample_vec[param1]
dim_use<-dimm_vec[param2]
reptime<-50
seed_num<-seed_num_vec[param3]
cfun<-function(a,b){
  rbind(a,b)
}
run_index<-paste("sample_",sample_num,"_dim_",dim_use,"_seed",seed_num,sep = "")
num_core<-20
cl <- makeCluster(num_core, outfile = paste("debug_",paste(run_index,collapse = ""),".txt",sep = ""))
registerDoParallel(cl)
data_simulation(sample_num=sample_num,dim_use = dim_use,reptime=reptime,seed_num=seed_num)

##3. Evaluate the performance
##--------------------------------------------------------------------------------------------
for (rep in 1:reptime) {
  load(file = paste("/home/ruibinxi_pkuhpc/lustre1/huangyf/R_project/EM-MMLE/result/data_sam",sample_num,"_dim",dim_use,"_seed",seed_num,"_rep",rep,".Rdata",sep = ""))
  X_mat<-data$X_mat
  ## 3.1 Covariance estimation
  #####################################################################
  ## 3.1.1 A unified rank-based approach from Fan et al.
  time_fan1<-Sys.time()
  R_all<-latentcor(X = X_mat, types = c("bin"), method = "original",use.nearPD=FALSE)$R
  time_fan2<-Sys.time()
  time_fan_est<-as.double(difftime(time_fan2,time_fan1,units = "mins"))
  save(R_all,file = paste("/home/ruibinxi_pkuhpc/lustre1/huangyf/R_project/EM-MMLE/result/R_all_sam",sample_num,"_dim",dim_use,"_seed",seed_num,"_rep",rep,".Rdata",sep = ""))
  save(time_fan_est,file = paste("/home/ruibinxi_pkuhpc/lustre1/huangyf/R_project/EM-MMLE/result/time_fan_est_sam",sample_num,"_dim",dim_use,"_seed",seed_num,"_rep",rep,".Rdata",sep = ""))
  
  ## 3.1.2 Estimate the covariance matrix by MMLE
  delta<-delta_est(X_mat)
  res_i<-foreach (
    i= 1:(dim_use-1),
    .combine = cfun,
    .inorder = TRUE,
    .packages = c('MASS','CVXR','mclust','mltools','doParallel','latentcor',"mvtnorm","lava","cubature")
  )%dopar%{
    res_mat<-matrix(0,1,dim_use)
    res_j<-foreach (
      j= (i+1):dim_use,
      .combine = c,
      .inorder = TRUE,
      .packages = c('MASS','CVXR',
                    'mclust','mltools','orthopolynom','Rcpp','RcppArmadillo','doParallel','latentcor',"lava",
                    "EmPLNet","Seurat","mvtnorm","cubature")
    )%dopar%{
      ## load data-------------
      print(paste("iter for",i,j))
      data_jk<-X_mat[,c(i,j)]
      delta_jk<-delta[c(i,j)]
      sample_use<-nrow(data_jk)
      sample_use<-nrow(data_jk)
      sigmahat_new<-R_all[i,j]
      diff_sigma<-1
      t<-1
      while (diff_sigma>1e-5 && t < 7) {
        sigmahat<-sigmahat_new
        b_mat<-matrix(c(1,0),2,1)
        a_mat<-matrix(c(0,1),2,1)
        sigma_init<-matrix(c(1,sigmahat,sigmahat,1),2,2)
        sigma_inv0<-solve(sigma_init)
        cov_matrix<-sigma_init
        gz_up_sigma<-function(x) {
          gaussian_density(x, mean = rep(0, 2), sigma = cov_matrix) * (sigma_inv0 %*% (t(t(x))) %*% (t(x)) %*% sigma_inv0)[1,2]
        }
        gz_up_sigma_2<-function(x) {
          bbb_2 <- (sigma_inv0 %*% (t(t(x))) %*% (t(x)) %*% sigma_inv0)[1,2]
          second_part11<-(t(t(x))) %*% (t(x)) %*% sigma_inv0 %*% ((a_mat)%*%t(b_mat)+(b_mat)%*%t(a_mat))
          second_part12<- sigma_inv0 %*% ((a_mat)%*%t(b_mat)+(b_mat)%*%t(a_mat)) %*% sigma_inv0
          bbb_3 <- tr(second_part11 %*% second_part12)
          gaussian_density(x, mean = rep(0, 2), sigma = cov_matrix) * ((-1) * bbb_3 + bbb_2^2)
        }
        gz_down_int<-function(x) {
          gaussian_density(x, mean = rep(0, 2), sigma = cov_matrix) 
        }
        gradiant_up_sigma_int<-rep(NA, sample_use)
        gradiant_down_int<-rep(NA, sample_use)
        Hessian_up_sigma_int<-rep(NA, sample_use)
        case<-c()
        gradiant_up_sigma_case<-c()
        gradiant_down_case<-c()
        Hessian_up_sigma_case<-c()
        for (q in c(0,1)) {
          for (w in c(0,1)) {
            case_1<-c(q,w)
            case<-rbind(case,case_1)
            Limit_input1<-identity_Limit(case_1,delta_jk)
            gradiant_up_sigma_case<-c(gradiant_up_sigma_case,adaptIntegrate(gz_up_sigma, lowerLimit= Limit_input1$lowerLimit, upperLimit = Limit_input1$upperLimit)$integral)
            gradiant_down_case<-c(gradiant_down_case,adaptIntegrate(gz_down_int, lowerLimit= Limit_input1$lowerLimit, upperLimit =Limit_input1$upperLimit)$integral)
            Hessian_up_sigma_case<-c(Hessian_up_sigma_case,adaptIntegrate(gz_up_sigma_2, lowerLimit=Limit_input1$lowerLimit, upperLimit = Limit_input1$upperLimit)$integral)
            
          }
        }
        for (sample_index in 1:sample_use) {
          index00<-which(apply(case, 1, function(row) all(row == data_jk[sample_index,])))
          gradiant_up_sigma_int[sample_index] <- gradiant_up_sigma_case[index00]
          gradiant_down_int[sample_index]<- gradiant_down_case[index00]
          Hessian_up_sigma_int[sample_index]<- Hessian_up_sigma_case[index00]
        }
        gradiant_int <- sum(gradiant_up_sigma_int/gradiant_down_int-sigma_inv0[1,2])
        Hessian_int<-sum(Hessian_up_sigma_int /gradiant_down_int-(gradiant_up_sigma_int/gradiant_down_int)^2+sigma_inv0[1,1] * sigma_inv0[2,2] + sigma_inv0[1,2]^2)
        
        if(abs(sigmahat-gradiant_int/Hessian_int)>=1){
          break
        }else{
          sigmahat_new<-sigmahat-gradiant_int/Hessian_int
        }
        diff_sigma<-abs(sigmahat_new-sigmahat)
        print(paste("diff_sigma is",diff_sigma))
        t<-t+1
      }
      
      return(sigmahat)
    }
    res_mat[1,(i+1):dim_use]<-res_j
    return(res_mat)
  }
  res_mat0<-cfun(res_i,matrix(0,1,dim_use))
  MMLE_res<-res_mat0+t(res_mat0)
  diag(MMLE_res)<-1
  save(MMLE_res,file = paste("/home/ruibinxi_pkuhpc/lustre1/huangyf/R_project/EM-MMLE/result/MMLE_res_sam",sample_num,"_dim",dim_use,"_seed",seed_num,"_rep",rep,".Rdata",sep = ""))
  
  ## 3.2 Theta estimation
  #####################################################################
  
  ## MMLE
  MMLE_Theta_res<-Theta_main1(obs_mat = X_mat,
                          n_lambda = 100,
                          penalize.diagonal = FALSE,
                          cov_input = MMLE_res,
                          weight_mat = NULL,zero_mat = NULL,
                          core_num = 1)
  
  ## Fan et al.
  Fan_res<-Theta_main1(obs_mat = X_mat,
                        n_lambda = 100,
                        penalize.diagonal = FALSE,
                        cov_input = R_all,
                        weight_mat = NULL,zero_mat = NULL,
                        core_num = 1)
  
  ## Oracle
  covZ_res<-Theta_main1(obs_mat = X_mat,
                          n_lambda = 100,
                          penalize.diagonal = FALSE,
                          cov_input = cov(data$Z_mat),
                          weight_mat = NULL,zero_mat = NULL,
                          core_num = 1)
  
  ## AMLE
  AMLE_res<-AMLE_fun(X_mat)
  
  ## save
  save(MMLE_Theta_res,file = paste("/home/ruibinxi_pkuhpc/lustre1/huangyf/R_project/EM-MMLE/result/MMLE_Theta_res_sam",sample_num,"_dim",dim_use,"_seed",seed_num,"_rep",rep,".Rdata",sep = ""))
  save(Fan_res,file = paste("/home/ruibinxi_pkuhpc/lustre1/huangyf/R_project/EM-MMLE/result/Fan_res_sam",sample_num,"_dim",dim_use,"_seed",seed_num,"_rep",rep,".Rdata",sep = ""))
  save(covZ_res,file = paste("/home/ruibinxi_pkuhpc/lustre1/huangyf/R_project/EM-MMLE/result/covZ_res_sam",sample_num,"_dim",dim_use,"_seed",seed_num,"_rep",rep,".Rdata",sep = ""))
  save(AMLE_res,file = paste("/home/ruibinxi_pkuhpc/lustre1/huangyf/R_project/EM-MMLE/result/AMLE_res_sam",sample_num,"_dim",dim_use,"_seed",seed_num,"_rep",rep,".Rdata",sep = ""))

}



























#MMLE-Dtrace-BIC
