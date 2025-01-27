# cond_pro_Xdata_Zg
cond_pro_Xdata_Zg <- function(Xdata, pi, mu, sigma, scale=NULL, obs_data_sp){
  log_Y_prod<-log_cul_mat_fun1(obs_data_sp)
  log_Y_prod<- as.matrix(log_Y_prod)
  log_Y_prod_rows<-rowSums(log_Y_prod)
  # define the values of scale K
  if(is.null(scale)){ scale = Xdata*0+1 }
  dim_use<-ncol(Xdata)
  sample_size<-nrow(Xdata)
  num_submodel<-dim(mu)[1]
  #the conditional probability of X, Gxn matrix
  log_cond_pro_Xdata<-matrix(0,nrow = num_submodel,ncol = sample_size)
  log_cond_pro_Xdata_elbo<-matrix(0,nrow = num_submodel,ncol = sample_size)
  for(g in 1:num_submodel){
    mu_g<-mu[g,]
    sigma_g<-sigma[[g]]
    precision_mat_g<-solve(sigma_g)
    A<-precision_mat_g
    # compute the conditional probability of X
    log_cond_pro_Xdata[g,]<-Update_log_cond_prob_X(mu=mu_g,precision_mat=precision_mat_g,Xdata=Xdata, scale=scale,core_num = 20)
    ##if(g==3){print(log_cond_pro_Xdata[g,])}
    log_cond_pro_Xdata_elbo[g,]<-log_cond_pro_Xdata[g,]+rowSums(log(scale)*Xdata)-log_Y_prod_rows
  }
  myfun<-function(x){x-max(x)}
  log_cond_pro_Xdata0<-apply(log_cond_pro_Xdata,2,myfun)
  cond_pro_Xdata0<-exp(log_cond_pro_Xdata0)
  pi.old<-pi
  # update pi
  pi.new<-c()
  post_dist_Z<-matrix(0,nrow = sample_size,ncol =num_submodel )
  for(i in 1:sample_size){
    post_dist_Z[i,]<-pi.old*cond_pro_Xdata0[,i]/(sum(pi.old*cond_pro_Xdata0[,i]))
  }
  pi.new<-colMeans(post_dist_Z)
  elbo_res<-c()
  for(g in 1:num_submodel){
    elbo_res[g]<-sum(pi[g]*log_cond_pro_Xdata_elbo[g,])
  }
  reture_item<-list(pi.new = pi.new, post_dist_Z = post_dist_Z,elbo_res=elbo_res)
  return(reture_item)
}

