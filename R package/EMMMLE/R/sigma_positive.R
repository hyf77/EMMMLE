# sigma_positive
sigma_positive<-function(cov_input){
  ##1. Calculate the sample covariance matrix by moment method----------------
  S<-Variable(dim(cov_input)[1], dim(cov_input)[1], PSD = TRUE)

  obj<-Minimize(max(abs(S-cov_input)))

  prob <- CVXR::Problem(obj)

  result <- CVXR::solve(prob,solver="SCS",verbose=FALSE)
  sigma_me2 <- result$getValue(S)
  sigma_me2<-max(abs(sigma_me2-cov_input))* diag(ncol(cov_input)) + sigma_me2

  ##
  return(sigma_me2)
}
