# EMMMLE_Network
EMMMLE_Network<-function(obs_mat, Sd_est = "TSS", lambda_vec = NULL, n_lambda =100,
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
  ##1. Determine the S_depth--------------------
  if(length(Sd_est) == 1){
    S_depth<-compute_offset(obs_mat,offset = Sd_est)
  }else{
    S_depth <- Sd_est
  }
  message("1. library size determination has been completed.")
  #################

  ##2. Estimate the mean vector and covariance matrix by moment method------------------------
  ##time cost 1
  time_sigmamoment1<-Sys.time()
  if(is.null(cov_input)){
    cov_Y1<-sigma_moment(data_use = obs_mat,S_depth = S_depth)
  }else{
    cov_Y1<-cov_input
    eigen_cov_Y1<-eigen(cov_Y1)
    if(min(eigen_cov_Y1$values)<0){
      cov_Y1<-sigma_positive(cov_input = cov_Y1)
    }
  }
  mean_est<-log(colMeans(obs_mat/as.vector(S_depth))) - diag(cov_Y1)/2
  time_sigmamoment2<-Sys.time()
  time_sigmamoment<-as.double(difftime(time_sigmamoment2,time_sigmamoment1,units = "hours"))
  message("2. The estimation of the mean vector and covariance matrix by moment method has been completed.")
  ###################################

  ##3. Determine the lambda_vec--------------------------
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
  message("3. Lambda_vec determination has been completed.")
  ######################################

  ##4. Run EQUAL for sparse precision matrix estimation-----------------------
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
  message("4. The estimation of precision matrix has been completed.")
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
  message("5. Hyper-parameter selection by BIC criterion and AIC criterion have been completed.")
  ############################################

  ##names
  for(l in 1:length(EQUAL_res$lambda)){
    rownames(EQUAL_res$Omega[[l]])<-Genename
    colnames(EQUAL_res$Omega[[l]])<-Genename
  }
  names(S_depth)<-Cellname
  rownames(Omega_chooseB)<-Genename
  colnames(Omega_chooseB)<-Genename
  rownames(Omega_chooseA)<-Genename
  colnames(Omega_chooseA)<-Genename

  # rownames(Omega_est)<-Omega_est = EQUAL_res$Omega
  PLNet_res<-list(Omega_est = EQUAL_res$Omega,
                  mean_est = mean_est,
                  lambda_vec = EQUAL_res$lambda, S_depth = S_depth,
                  BIC_vec = BIC_vec, Omega_chooseB = Omega_chooseB,
                  AIC_vec = AIC_vec, Omega_chooseA = Omega_chooseA,
                  time_sigmamoment = time_sigmamoment,time_EQUAL = time_EQUAL)
  return(PLNet_res)
}
