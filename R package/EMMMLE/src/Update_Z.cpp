// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]
#include <omp.h>
#include <iomanip>
#include <RcppEigen.h>
#include <Rcpp.h>
#include <Eigen/Sparse>
#include <math.h>

using namespace Rcpp;
using namespace std;
using namespace Eigen;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::ArrayXd;
using Eigen::ArrayXf;
using std::cout;
using std::endl;
static omp_lock_t emlock;

// [[Rcpp::export]]
Eigen::VectorXd Update_log_cond_prob_X(Eigen::VectorXd mu, Eigen::MatrixXd precision_mat, Eigen::MatrixXd Xdata, Eigen::MatrixXd scale, double eps = 1e-6, int max_iter = 100,int core_num=1) {
  int sample_size=scale.rows();
  Eigen::VectorXd Xdata_Zg(sample_size);
  //
  omp_init_lock(&emlock);
  omp_set_lock(&emlock);

  omp_set_num_threads(core_num);
#pragma omp parallel for
  for(int i=0;i<(sample_size);i++){
    //cout<< "i"<< i <<endl;
    Eigen::VectorXd scale_n=scale.row(i);
    Eigen::VectorXd Xn=Xdata.row(i);
    Eigen::MatrixXd A = precision_mat;
    int p = Xn.size();
    //Eigen::VectorXd eta_int = (Xn.array()/scale_n.array()+0.1).log();
    Eigen::VectorXd eta_int = (Xn.array()+1).log();
    double diff = 1;
    int iter = 1;
    Eigen::VectorXd  eta_new = eta_int;
    Eigen::VectorXd gradiant_old(p);
    double diff_grad=0;
    while(diff > eps && iter <= max_iter) {
      Eigen::VectorXd eta = eta_new;
      Eigen::VectorXd gradiant(p);
      Eigen::MatrixXd Hessian0(p, p);
      Eigen::MatrixXd Hessian_inv(p, p);
      gradiant = Xn-(scale_n.array() * eta.array().exp()).matrix()-A*(eta - mu);
      Eigen::VectorXd temp_vec=scale_n.array() * eta.array().exp();
      Hessian0=temp_vec.asDiagonal();
      Hessian_inv=(Hessian0+A).inverse();
      eta_new = eta + Hessian_inv * gradiant;
      diff = (eta - eta_new).array().abs().maxCoeff();
      if(iter >1){
        diff_grad=(gradiant.array() -gradiant_old.array()).abs().maxCoeff();
        if(diff_grad<1e-4)
        {
          break;
        }
      }
      gradiant_old=gradiant;
      iter++;
    }
    Eigen::VectorXd temp_vec0=scale_n.array() * eta_new.array().exp();
    Eigen::MatrixXd A_star =temp_vec0.asDiagonal();
    A_star=A_star+A;
    double A_det_sqrt_log=log(sqrt(A.determinant()));
    int xx=2;
    while((!isfinite(A_det_sqrt_log)) ){
      A_det_sqrt_log=log(sqrt((A/xx).determinant()))+0.5*p*log(xx);
      xx=xx+2;
    }
    double A_star_det_log=log(A_star.determinant());
    int yy=2;
    while((!isfinite(A_star_det_log)) ){
      A_star_det_log=log((A_star/yy).determinant())+p*log(yy);
      yy=yy+2;
    }

    double log_p_i = A_det_sqrt_log+(-1) * 0.5 *(eta_new - mu).transpose()* A * (eta_new - mu)- (scale_n.array() * eta_new.array().exp()).sum();
    log_p_i=log_p_i+(Xn.array() * eta_new.array()).sum() - 0.5 * A_star_det_log;

    Xdata_Zg(i)=log_p_i;
  }
  omp_destroy_lock(&emlock);

  return Xdata_Zg;
}
