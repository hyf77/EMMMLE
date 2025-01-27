##part_of_interest
part_of_interest<-function(mat_use,p_TF){
  return(c((mat_use[1:p_TF,1:p_TF])[upper.tri((mat_use[1:p_TF,1:p_TF]))],as.vector(mat_use[1:p_TF,-(1:p_TF)])))
}