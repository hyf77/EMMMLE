#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <omp.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace arma;

static omp_lock_t lock;

// [[Rcpp::export]]
arma::mat soft(arma::mat A, arma::mat a_mat, int diag=0){
  a_mat=abs(a_mat);
  arma::mat C=(A>=a_mat)%(A-a_mat)+(A<=(-a_mat))%(A+a_mat);
if (diag==0){
  arma::mat B=A-arma::diagmat(arma::diagvec(A));
  B=(B>=a_mat)%(B-a_mat)+(B<=(-a_mat))%(B+a_mat);
  B=B+arma::diagmat(arma::diagvec(A));
  C=B;}
  return C;}

// [[Rcpp::export]]
Rcpp::List equal1(arma::mat X,arma::vec lambda,arma::mat weight_mat, arma::mat zero_mat,
                  double err=10^(-5),int maxIter=10^3,double rho=1, int diag=0,
                  int core_num = 1){
  int n=X.n_rows;
  int p=X.n_cols;
  int m=(p>n)*n+(p<=n)*p;
  int nlambda=lambda.size();
  //
  arma::mat zero_mat_choose = 1 - zero_mat;
  /*Centering int*/
  arma::mat dX=(arma::eye(n,n)-arma::ones(n,n)/n)*X/sqrt(n-1);
  arma::mat U;
  arma::mat Uv;
  arma::vec eigd;
  svd(U, eigd, Uv, dX.t()); 
  U=U.cols(0,m-1);
  eigd=eigd%eigd;
  arma::mat D=U*arma::diagmat(eigd/(eigd+rho));
  Rcpp::List Omega_all(nlambda);
  for (int k=0;k<nlambda;++k) {
    arma::mat temp_mat(p,p,fill::zeros);
    Omega_all(k) = temp_mat;
  }
  arma::vec niter=lambda;
  /*Intialization*/
  arma::mat aZ=arma::eye(p,p);
  arma::mat aU=arma::zeros(p,p);
  arma::mat aX;
  arma::mat L;
  arma::mat Z1;
  //
  omp_init_lock(&lock);
  omp_set_lock(&lock);
  
  omp_set_num_threads(core_num);
#pragma omp parallel for
  for (int k=0;k<nlambda;++k) {
    double lam=lambda(k);
    int i=0;
    arma::mat lambda_mat = weight_mat * (lam/rho);
    double ee=1;
    //
    while (((i<maxIter)&&(ee>err))||(i==0))
    { Z1=aZ;
      L=arma::eye(p,p)/rho+aZ-aU;
      aX=L-D*(U.t()*L);
      //
      aZ=soft(aX+aU, lambda_mat, diag);
      //
      aZ = aZ % zero_mat_choose;
      //
      aU=aU+aX-aZ;
      ee=mean(mean(abs(aZ-Z1)));
      i=i+1;
    }
    Omega_all(k)=arma::sp_mat((abs(aZ)<abs(aZ.t()))%aZ+(abs(aZ)>=abs(aZ.t()))%aZ.t());
    niter(k)=i;
  }
  omp_destroy_lock(&lock);
  //
  return Rcpp::List::create(Rcpp::Named("Omega") =Omega_all,
                            Rcpp::Named("lambda") =lambda,
                            Rcpp::Named("niter") =niter); }

// [[Rcpp::export]]
Rcpp::List equal2(arma::mat X,arma::vec lambda,arma::Mat<double> weight_mat, arma::mat zero_mat,
                  double err=10^(-5),int maxIter=10^3,double rho=1,int diag=0,
                  int core_num = 1){
  int n=X.n_rows;
  int p=X.n_cols;
  int m=(p>n)*n+(p<=n)*p;
  int nlambda=lambda.size();
  //
  arma::mat zero_mat_choose = 1 - zero_mat;
  /*Centering int*/
  arma::mat dX=(arma::eye(n,n)-arma::ones(n,n)/n)*X/sqrt(n-1);
  arma::mat U;
  arma::mat Uv;
  arma::vec eigd;
  svd(U, eigd, Uv, dX.t());
  eigd=eigd%eigd;
  arma::mat U1=U.cols(0,m-1)*arma::diagmat(sqrt(eigd/(eigd+2*rho)));
 arma::vec ld=exp(eigd);
  arma::mat D=2*rho/(log(ld*ld.t())+2*rho)+1;
 Rcpp::List Omega_all(nlambda);
 arma::vec niter=lambda;
 /*Intialization*/
 arma::mat aZ=arma::eye(p,p);
 arma::mat aU=arma::zeros(p,p);
 arma::mat aX;
 arma::mat L;
 arma::mat L1;
 arma::mat L2;
 arma::mat Z1;
 double lam;
 double ee=1;
 //
 omp_init_lock(&lock);
 omp_set_lock(&lock);
 
 omp_set_num_threads(core_num);
#pragma omp parallel for
 for (int k=0;k<nlambda;++k) {
   lam=lambda(k);
   int i=0;
   arma::mat lambda_mat = weight_mat * (lam/rho);
   //
   while (((i<maxIter)&&(ee>err))||(i==0))
   { Z1=aZ;
     L=arma::eye(p,p)/rho+aZ-aU;
     L=(L+L.t())/2;
     L1=L*U1;
     L2=L1*U1.t();
     aX=L-L2-L2.t()+U1*(D%(U1.t()*L1))*U1.t();
     aZ=soft(aX+aU,lambda_mat,diag);
     //
     aZ = aZ % zero_mat_choose;
     //
     aU=aU+aX-aZ;
     ee=mean(mean(abs(aZ-Z1)));
     i=i+1;
   }
   Omega_all(k)=arma::sp_mat(aZ);
   niter(k)=i;
 }
 omp_destroy_lock(&lock);
 //
 return Rcpp::List::create(Rcpp::Named("Omega") =Omega_all,
                           Rcpp::Named("lambda") =lambda,
                           Rcpp::Named("niter") =niter); }
