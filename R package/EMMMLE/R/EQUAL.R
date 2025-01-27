# EQUAL
EQUAL<-function(X,type=TRUE,sdiag=FALSE,lambda=NULL,lambda.min=sqrt(log(ncol(X))/nrow(X)),nlambda=50,err=10^(-5),maxIter =1000,rho=1,
                weight_mat,
                zero_mat,
                core_num)
{p=ncol(X);
n=nrow(X);
Sn<-cov(X);
A<-abs(Sn-diag(diag(Sn)));
if (is.null(lambda)){lambda=exp(seq(log(lambda.min),0,length.out =nlambda))*max(A)}
lambda<-sort(lambda,decreasing =TRUE)

if (type){obj<-equal2(X,lambda =lambda,diag=as.numeric(sdiag),err=err,maxIter =maxIter,rho=rho,
                      weight_mat = weight_mat,
                      zero_mat = zero_mat,
                      core_num = core_num)}
else {obj<-equal1(X,lambda =lambda,diag=as.numeric(sdiag),err=err,maxIter =maxIter,rho=rho,
                  weight_mat = weight_mat,
                  zero_mat = zero_mat,
                  core_num = core_num)}
return(obj)
}
