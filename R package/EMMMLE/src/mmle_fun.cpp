#include <RcppArmadillo.h>
#include <omp.h>
#include <RcppEigen.h>
#include <Eigen/Sparse>
#include <Rcpp.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace std;
using namespace arma;
using namespace Eigen;
using Eigen::MatrixXd;
using Eigen::VectorXd;

static omp_lock_t lock_temp;

// [[Rcpp::export]]
Eigen::SparseMatrix<double> log_cul_mat_fun1(Eigen::SparseMatrix<int> mat){
  int num_row = mat.rows();
  int num_col = mat.cols();
  Eigen::SparseMatrix<double> mat_res(num_row,num_col);
  
  for(int k=0;k<mat.outerSize();++k){
    for(SparseMatrix<int>::InnerIterator it(mat,k);it;++it){
      Eigen::VectorXd vec_cur;
      int it_value = it.value();
      vec_cur.setZero(it_value);
      for(int j = 0;j<it_value;j++){
        vec_cur[j] = j+1;
      }
      
      mat_res.insert(it.row(),it.col()) = vec_cur.array().log().sum();
    }
  }
  //
  mat_res.makeCompressed();
  
  return mat_res;
}



// [[Rcpp::export]]
Rcpp::List mmle_weighted_step1(arma::mat data_use, arma::vec cond_Z_input, arma::vec S_depth, arma::vec mu_init, arma::mat sigma_init,
                               arma::vec t_root_vec, arma::vec omega_root_vec,int k_max, int core_num,int cut_num){
  int dim_use = data_use.n_cols;
  int sample_size = data_use.n_rows;
  int length = t_root_vec.n_elem;
  arma::mat mlemu(1,dim_use);
  arma::mat mlesigmahat(1,dim_use);
  arma::mat gradiant_iter_mat(k_max,dim_use);
  arma::mat update_iter_mat(k_max,dim_use);
  arma::vec t_root_vec1 = t_root_vec * sqrt(2);
  arma::vec oemga_troot_vec = omega_root_vec % exp(pow(t_root_vec,2)) * sqrt(2) ;
  arma::vec if_convergence_feature(dim_use);
  arma::vec cond_Z=cond_Z_input;
  //
  arma::mat sigmahat_input=sigma_init;
  arma::vec mu_input=mu_init;
  int Newstep=1;
  while (Newstep<(k_max+1)) {
    int k=Newstep;
    double sigma_diag_max = sigmahat_input.diag().max();
    arma::vec delta_mu(dim_use);
    arma::vec delta_sigma(dim_use);
    omp_init_lock(&lock_temp);
    omp_set_lock(&lock_temp);
    omp_set_num_threads(core_num);
#pragma omp parallel for
    for(int i=0;i<dim_use;i++){
      arma::vec data_use_dim=data_use.col(i);
      double mu=mu_input(i);
      double sigmahat=sigmahat_input(i,i);
      //z.root_mat(10)
      arma::vec z_root = t_root_vec1 * sqrt(sigmahat) + mu;
      //m.root_mat(10)
      arma::vec m1_root = oemga_troot_vec * sqrt(sigmahat);
      arma::vec share12_vec=z_root-mu;
      arma::vec share11_vec=exp((-1)*pow(share12_vec,2) /(2* sigmahat));
      arma::mat root1_mat(length,sample_size);
      for(int u = 0;u<length;++u){
        arma::vec res_mat =data_use_dim*z_root(u)-S_depth * exp(z_root(u));
        root1_mat.row(u) = res_mat.t();
      }
      arma::vec max_vec(sample_size);
      for(int sample_index = 0;sample_index < sample_size; ++sample_index){
        arma::vec mat_temp1 = root1_mat.col(sample_index);
        max_vec(sample_index) = mat_temp1.max();
      }
      arma::mat share1_mat(length,sample_size);
      for(int u = 0;u<length;++u){
        share1_mat.row(u)=exp(root1_mat.row(u)-max_vec.t());
      }
      //std::cout << "share1_mat"<<share1_mat << std::endl;
      arma::mat gz_up_mu_mat(length,sample_size);
      arma::mat gz_up_mu2_mat(length,sample_size);
      arma::mat gz_up_sigma_mat(length,sample_size);
      arma::vec mid_vec1(length);
      arma::vec mid_vec2(length);
      arma::mat gz_up_sigma2_mat(length,sample_size);
      //arma::mat gz_down_mat(length,sample_size);
      //arma::mat gz_up_inter_mat(length,sample_size);
      mid_vec1=pow(share12_vec,4)/(4*pow(sigmahat,4))-pow(share12_vec,2)/(pow(sigmahat,3));
      mid_vec2=pow(share12_vec,3)/(2*pow(sigmahat,3))-share12_vec/(pow(sigmahat,2));
      arma::vec vec_temp1(sample_size);
      arma::vec share2_vec(sample_size);
      arma::vec share3_vec(sample_size);
      arma::vec share4_vec(sample_size);
      arma::vec share5_vec(sample_size);
      arma::vec share6_vec(sample_size);
      for(int sample_index = 0;sample_index < sample_size; ++sample_index){
        gz_up_mu_mat.col(sample_index) = (share1_mat.col(sample_index)%share12_vec/sigmahat)%share11_vec;
        gz_up_mu2_mat.col(sample_index) = (share1_mat.col(sample_index)% pow(share12_vec,2) /pow(sigmahat,2))%share11_vec;
        gz_up_sigma_mat.col(sample_index) = gz_up_mu2_mat.col(sample_index)/2;
        gz_up_sigma2_mat.col(sample_index) = share1_mat.col(sample_index)%mid_vec1%share11_vec;
        //gz_down_mat.col(sample_index) = share1_mat.col(sample_index)%share11_vec;
        //gz_up_inter_mat.col(sample_index) = share1_mat.col(sample_index)%mid_vec2%share11_vec;
        vec_temp1(sample_index)=sum(share1_mat.col(sample_index)%share11_vec%m1_root);
        share2_vec(sample_index)=sum(gz_up_mu_mat.col(sample_index)%m1_root)/vec_temp1(sample_index);
        share3_vec(sample_index)=sum(gz_up_sigma_mat.col(sample_index)%m1_root)/vec_temp1(sample_index);
        share4_vec(sample_index)=sum(gz_up_sigma2_mat.col(sample_index)%m1_root)/vec_temp1(sample_index);
        share5_vec(sample_index)=sum(share1_mat.col(sample_index)%mid_vec2%share11_vec%m1_root)/vec_temp1(sample_index);
        share6_vec(sample_index)=sum(gz_up_mu2_mat.col(sample_index)%m1_root)/vec_temp1(sample_index);
      }
      share2_vec.elem(arma::find_nonfinite(share2_vec)).zeros();
      share3_vec.elem(arma::find_nonfinite(share3_vec)).zeros();
      share4_vec.elem(arma::find_nonfinite(share4_vec)).zeros();
      share5_vec.elem(arma::find_nonfinite(share5_vec)).zeros();
      share6_vec.elem(arma::find_nonfinite(share6_vec)).zeros();
      double weight_temp=sum(cond_Z)/(double(sample_size));
      double gradiant_mu_k=sum(share2_vec%cond_Z)/(double(sample_size));
      double gradiant_sigma_k=sum(share3_vec%cond_Z)/(double(sample_size))-(0.5*weight_temp)/sigmahat;
      arma::vec temp_share2_2=pow(share2_vec,2);
      arma::vec temp_share3_2=pow(share3_vec,2);
      double hessian_mu_k=sum(share6_vec%cond_Z)/(double(sample_size))-sum(temp_share2_2%cond_Z)/(double(sample_size))-(weight_temp)/sigmahat;
      double hessian_sigma_k=sum(share4_vec%cond_Z)/(double(sample_size))-sum(temp_share3_2%cond_Z)/(double(sample_size))+(0.5*weight_temp)/pow(sigmahat,2);
      double hessian_int_k=sum(share5_vec%cond_Z)/(double(sample_size))-sum(share2_vec%share3_vec%cond_Z)/(double(sample_size));
      //
      arma::mat hessian_mat(2,2);
      arma::vec gradiant(2);
      gradiant(0)=gradiant_mu_k;
      gradiant(1)=gradiant_sigma_k;
      hessian_mat(0,0)=hessian_mu_k;
      hessian_mat(0,1)=hessian_int_k;
      hessian_mat(1,0)=hessian_int_k;
      hessian_mat(1,1)=hessian_sigma_k;
      arma::mat hessian_inv=inv(hessian_mat);
      arma::vec delta=hessian_inv*gradiant;
      gradiant_iter_mat(k-1,i)=gradiant(0);
      if(k==k_max){
        if(abs(gradiant(1))<1e-3){if_convergence_feature(i)=true;}
        else{if_convergence_feature(i)=false;}
      }
      delta_mu(i)=delta(0);
      delta_sigma(i)=delta(1);
      if (sigmahat-delta_sigma(i)<= 0) {
        delta_mu(i) = gradiant_mu_k / hessian_mu_k;
        delta_sigma(i)= gradiant_sigma_k / hessian_sigma_k;
      }
      if (abs(delta_sigma(i)) >sigma_diag_max) {
        delta_mu(i)= 0;
        delta_sigma(i) = 0;
      }
    }
    omp_destroy_lock(&lock_temp);
    arma::vec mlemu1=mu_input-delta_mu;
    arma::mat sigmahat_tempp=sigmahat_input;
    arma::vec diag_elements = sigmahat_tempp.diag();
    // Find indices where the diagonal elements are non-positive
    arma::uvec neg_index = arma::find(diag_elements <= 0);
    if (!neg_index.is_empty()) {
      // Find indices where the diagonal elements are positive
      arma::uvec pos_index = arma::find(diag_elements > 0);
      double min_pos_diag = diag_elements(pos_index).min();
      diag_elements(neg_index).fill(min_pos_diag);
      sigmahat_tempp.diag() = diag_elements;
    }
    // Adjust diagonal based on delta_sigma
    arma::vec reduce_diag = sigmahat_tempp.diag() - delta_sigma;
    arma::uvec adjust_index = arma::find(reduce_diag <= 0);
    arma::mat mlesigmahat1;
    if (!adjust_index.is_empty()) {
      std::cout << "Adjust diagonal" << std::endl;
      arma::mat reduce_mat = arma::diagmat(delta_sigma);
      arma::vec reduce_mat_diag=reduce_mat.diag();
      reduce_mat_diag(adjust_index).zeros();
      reduce_mat.diag()=reduce_mat_diag;
      arma::mat mlesigmahat10 = sigmahat_tempp - reduce_mat;
      arma::uvec non_adjust_index = arma::find(reduce_diag > 0);
      arma::vec mlesigmahat10_diag=mlesigmahat10.diag();
      double min_value =mlesigmahat10_diag(non_adjust_index).min();
      mlesigmahat10_diag(adjust_index).fill(min_value);
      mlesigmahat10.diag()=mlesigmahat10_diag;
      mlesigmahat1 = mlesigmahat10;
    } else {
      mlesigmahat1 = sigmahat_tempp - arma::diagmat(delta_sigma);
    }
    arma::vec set_diag = mlesigmahat1.diag();
    for (int i = 0; i < set_diag.n_elem; ++i) {
      if (set_diag[i] > cut_num) {
        set_diag(i) = cut_num;
      }
      if (set_diag(i)< ((-1)*cut_num)) {
        set_diag(i) = (-1)*cut_num;
      }
    }
    mlesigmahat1.diag()=set_diag;
    arma::vec mlesigmahat1_diag=mlesigmahat1.diag();
    mu_input=mlemu1;
    sigmahat_input=mlesigmahat1;
    Newstep=Newstep+1;
  }
  omp_init_lock(&lock_temp);
  omp_set_lock(&lock_temp);
  omp_set_num_threads(core_num);
#pragma omp parallel for
  for(int i=0;i<dim_use;i++){
    //put here
    arma::vec vecUse = abs(gradiant_iter_mat.col(i));
    arma::vec aaaVec = diff(vecUse);
    arma::uvec aVec=find(aaaVec>0);
    if (aVec.n_elem > 0) {
      int minAAA = arma::min(aVec) + 1;
      //std::cout << minAAA << std::endl;
      int minVal = std::min(minAAA, k_max-1);
      update_iter_mat.submat(minVal, i, k_max-1, i).fill(false);
      int minbbb=arma::min(aVec);
      update_iter_mat.submat(0, i, minbbb, i).fill(true);
      //std::cout <<update_iter_mat.col(i) << std::endl;
    } else {
      update_iter_mat.col(i).fill(true);
    }
  }
  omp_destroy_lock(&lock_temp);
  //
  return Rcpp::List::create(Rcpp::Named("if_convergence_feature") = if_convergence_feature,
                            Rcpp::Named("update_iter_mat")=update_iter_mat);

}

// [[Rcpp::export]]
Rcpp::List mmle_weighted_step2(arma::mat data_use, arma::vec cond_Z_input, arma::vec S_depth,arma::vec convergence_vec, arma::mat update_mat, arma::vec mu_init, arma::mat sigma_init,
                               arma::vec t_root_vec, arma::vec omega_root_vec,int k_max, int core_num,int cut_num){
  int dim_use = data_use.n_cols;
  int sample_size = data_use.n_rows;
  int length = t_root_vec.n_elem;
  arma::mat mlemu(1,dim_use);
  arma::mat mlesigmahat(1,dim_use);
  arma::vec mlemu1(dim_use);
  arma::mat mlesigmahat1(dim_use,dim_use);
  arma::mat update_iter_mat=update_mat;
  arma::vec t_root_vec1 = t_root_vec * sqrt(2);
  arma::vec oemga_troot_vec = omega_root_vec % exp(pow(t_root_vec,2)) * sqrt(2) ;
  arma::vec if_convergence_feature=convergence_vec;
  arma::vec cond_Z=cond_Z_input;
  //
  arma::vec if_update_feature(dim_use, arma::fill::ones);
  //
  arma::mat sigmahat_input=sigma_init;
  arma::vec mu_input=mu_init;
  int Newstep=1;
  while(Newstep<k_max+1){
    int k=Newstep;
    double sigma_diag_max = sigmahat_input.diag().max();
    arma::vec delta_mu(dim_use);
    arma::vec delta_sigma(dim_use);
    omp_init_lock(&lock_temp);
    omp_set_lock(&lock_temp);
    omp_set_num_threads(core_num);
#pragma omp parallel for
    for(int i=0;i<dim_use;i++){
      arma::vec data_use_dim=data_use.col(i);
      double mu=mu_input(i);
      double sigmahat=sigmahat_input(i,i);
      //z.root_mat(10)
      arma::vec z_root = t_root_vec1 * sqrt(sigmahat) + mu;
      //m.root_mat(10)
      arma::vec m1_root = oemga_troot_vec * sqrt(sigmahat);
      arma::vec share12_vec=z_root-mu;
      arma::vec share11_vec=exp((-1)*pow(share12_vec,2) /(2* sigmahat));
      arma::mat root1_mat(length,sample_size);
      for(int u = 0;u<length;++u){
        arma::vec res_mat =data_use_dim*z_root(u)-S_depth * exp(z_root(u));
        root1_mat.row(u) = res_mat.t();
      }
      arma::vec max_vec(sample_size);
      for(int sample_index = 0;sample_index < sample_size; ++sample_index){
        arma::vec mat_temp1 = root1_mat.col(sample_index);
        max_vec(sample_index) = mat_temp1.max();
      }
      arma::mat share1_mat(length,sample_size);
      for(int u = 0;u<length;++u){
        share1_mat.row(u)=exp(root1_mat.row(u)-max_vec.t());
      }
      arma::mat gz_up_mu_mat(length,sample_size);
      arma::mat gz_up_mu2_mat(length,sample_size);
      arma::mat gz_up_sigma_mat(length,sample_size);
      arma::vec mid_vec1(length);
      arma::vec mid_vec2(length);
      arma::mat gz_up_sigma2_mat(length,sample_size);
      //arma::mat gz_down_mat(length,sample_size);
      //arma::mat gz_up_inter_mat(length,sample_size);
      mid_vec1=pow(share12_vec,4)/(4*pow(sigmahat,4))-pow(share12_vec,2)/(pow(sigmahat,3));
      mid_vec2=pow(share12_vec,3)/(2*pow(sigmahat,3))-share12_vec/(pow(sigmahat,2));
      arma::vec vec_temp1(sample_size);
      arma::vec share2_vec(sample_size);
      arma::vec share3_vec(sample_size);
      arma::vec share4_vec(sample_size);
      arma::vec share5_vec(sample_size);
      arma::vec share6_vec(sample_size);
      for(int sample_index = 0;sample_index < sample_size; ++sample_index){
        gz_up_mu_mat.col(sample_index) = (share1_mat.col(sample_index)%share12_vec/sigmahat)%share11_vec;
        gz_up_mu2_mat.col(sample_index) = (share1_mat.col(sample_index)% pow(share12_vec,2) /pow(sigmahat,2))%share11_vec;
        gz_up_sigma_mat.col(sample_index) = gz_up_mu2_mat.col(sample_index)/2;
        gz_up_sigma2_mat.col(sample_index) = share1_mat.col(sample_index)%mid_vec1%share11_vec;
        //gz_down_mat.col(sample_index) = share1_mat.col(sample_index)%share11_vec;
        //gz_up_inter_mat.col(sample_index) = share1_mat.col(sample_index)%mid_vec2%share11_vec;
        vec_temp1(sample_index)=sum(share1_mat.col(sample_index)%share11_vec%m1_root);
        share2_vec(sample_index)=sum(gz_up_mu_mat.col(sample_index)%m1_root)/vec_temp1(sample_index);
        share3_vec(sample_index)=sum(gz_up_sigma_mat.col(sample_index)%m1_root)/vec_temp1(sample_index);
        share4_vec(sample_index)=sum(gz_up_sigma2_mat.col(sample_index)%m1_root)/vec_temp1(sample_index);
        share5_vec(sample_index)=sum(share1_mat.col(sample_index)%mid_vec2%share11_vec%m1_root)/vec_temp1(sample_index);
        share6_vec(sample_index)=sum(gz_up_mu2_mat.col(sample_index)%m1_root)/vec_temp1(sample_index);
      }
      share2_vec.elem(arma::find_nonfinite(share2_vec)).zeros();
      share3_vec.elem(arma::find_nonfinite(share3_vec)).zeros();
      share4_vec.elem(arma::find_nonfinite(share4_vec)).zeros();
      share5_vec.elem(arma::find_nonfinite(share5_vec)).zeros();
      share6_vec.elem(arma::find_nonfinite(share6_vec)).zeros();
      double weight_temp=sum(cond_Z)/(double(sample_size));
      double gradiant_mu_k=sum(share2_vec%cond_Z)/(double(sample_size));
      double gradiant_sigma_k=sum(share3_vec%cond_Z)/(double(sample_size))-(0.5*weight_temp)/sigmahat;
      arma::vec temp_share2_2=pow(share2_vec,2);
      arma::vec temp_share3_2=pow(share3_vec,2);
      double hessian_mu_k=sum(share6_vec%cond_Z)/(double(sample_size))-sum(temp_share2_2%cond_Z)/(double(sample_size))-(weight_temp)/sigmahat;
      double hessian_sigma_k=sum(share4_vec%cond_Z)/(double(sample_size))-sum(temp_share3_2%cond_Z)/(double(sample_size))+(0.5*weight_temp)/pow(sigmahat,2);
      double hessian_int_k=sum(share5_vec%cond_Z)/(double(sample_size))-sum(share2_vec%share3_vec%cond_Z)/(double(sample_size));
      //
      arma::mat hessian_mat(2,2);
      arma::vec gradiant(2);
      gradiant(0)=gradiant_mu_k;
      gradiant(1)=gradiant_sigma_k;
      hessian_mat(0,0)=hessian_mu_k;
      hessian_mat(0,1)=hessian_int_k;
      hessian_mat(1,0)=hessian_int_k;
      hessian_mat(1,1)=hessian_sigma_k;
      arma::mat hessian_inv=inv(hessian_mat);
      arma::vec delta=hessian_inv*gradiant;

      if (if_convergence_feature(i)==0) {
        if ((update_iter_mat(k-1,i)==0) || (if_update_feature(i)==0)) {
          if_update_feature(i)=0;
          delta(0) = 0;
          delta(1) = 0;
        } else {
          if (sigmahat-delta(1)<= 0) {
            delta(0) = gradiant_mu_k / hessian_mu_k;
            delta(1)= gradiant_sigma_k / hessian_sigma_k;
          }
          if (abs(delta(1)) >sigma_diag_max) {
            delta(0)= 0;
            delta(1) = 0;
          }
        }
      } else {
        if (if_update_feature(i)==1) {
          if (sigmahat-delta(1)<= 0) {
            delta(0) = gradiant_mu_k / hessian_mu_k;
            delta(1)= gradiant_sigma_k / hessian_sigma_k;
          }
          if (abs(delta(1)) >sigma_diag_max) {
            delta(0)= 0;
            delta(1) = 0;
          }
        }
      }
      delta_mu(i)=delta(0);
      delta_sigma(i)=delta(1);
      //mlemu_i(k)=mu-delta(0);
      //mlesigma_i(k)=sigmahat-delta(1);
    }
    omp_destroy_lock(&lock_temp);
    mlemu1=mu_input-delta_mu;
    arma::mat sigmahat_tempp=sigmahat_input;
    arma::vec diag_elements = sigmahat_tempp.diag();
    // Find indices where the diagonal elements are non-positive
    arma::uvec neg_index = arma::find(diag_elements <= 0);
    if (!neg_index.is_empty()) {
      // Find indices where the diagonal elements are positive
      arma::uvec pos_index = arma::find(diag_elements > 0);
      double min_pos_diag = diag_elements(pos_index).min();
      diag_elements(neg_index).fill(min_pos_diag);
      sigmahat_tempp.diag() = diag_elements;
    }
    // Adjust diagonal based on delta_sigma
    arma::vec reduce_diag = sigmahat_tempp.diag() - delta_sigma;
    arma::uvec adjust_index = arma::find(reduce_diag <= 0);
    if (!adjust_index.is_empty()) {
      arma::mat reduce_mat = arma::diagmat(delta_sigma);
      arma::vec reduce_mat_diag=reduce_mat.diag();
      reduce_mat_diag(adjust_index).zeros();
      reduce_mat.diag()=reduce_mat_diag;
      arma::mat mlesigmahat10 = sigmahat_tempp - reduce_mat;
      arma::uvec non_adjust_index = arma::find(reduce_diag > 0);
      arma::vec mlesigmahat10_diag=mlesigmahat10.diag();
      double min_value =mlesigmahat10_diag(non_adjust_index).min();
      mlesigmahat10_diag(adjust_index).fill(min_value);
      mlesigmahat10.diag()=mlesigmahat10_diag;
      mlesigmahat1 = mlesigmahat10;
    } else {
      mlesigmahat1 = sigmahat_tempp - arma::diagmat(delta_sigma);
    }
    arma::vec set_diag = mlesigmahat1.diag();
    for (int i = 0; i < set_diag.n_elem; ++i) {
      if (set_diag[i] > cut_num) {
        set_diag(i) = cut_num;
      }
      if (set_diag(i)< ((-1)*cut_num)) {
        set_diag(i) = (-1)*cut_num;
      }
    }
    mlesigmahat1.diag()=set_diag;
    arma::vec mlesigmahat1_diag=mlesigmahat1.diag();
    mu_input=mlemu1;
    sigmahat_input=mlesigmahat1;
    Newstep=Newstep+1;
  }
  arma::vec sigmadiag_end=mlesigmahat1.diag();
  //
  return Rcpp::List::create(Rcpp::Named("mlemu") = mlemu1,
                            Rcpp::Named("update_iter_mat")=update_iter_mat,
                            Rcpp::Named("mlesigmahat") = sigmadiag_end);

}

// [[Rcpp::export]]
Rcpp::List Update_sigma12_fun(arma::mat data_use, arma::vec cond_Z_input, arma::vec S_depth, arma::vec mlemu1, arma::mat mlesigmahat1,
                              arma::vec t_root_vec, arma::vec omega_root_vec,int core_num,arma::mat position_input){
  //def
  arma::mat a_mat("0 ; 1");
  arma::mat b_mat("1 ; 0");
  arma::vec cond_Z=cond_Z_input;
  int dim_use = data_use.n_cols;
  int sample_size_use = data_use.n_rows;
  int position_input_rows=position_input.n_rows;
  arma::vec t_root_vec1 = t_root_vec * sqrt(2);
  arma::vec oemga_troot_vec = omega_root_vec;
  int length1 = t_root_vec1.n_elem;
  int length2 = oemga_troot_vec.n_elem;
  //
  arma::mat gradiant_int(dim_use,dim_use,arma::fill::zeros);
  arma::mat Hessian_int(dim_use,dim_use,arma::fill::zeros);
  //
  omp_init_lock(&lock_temp);
  omp_set_lock(&lock_temp);
  //
  //
  omp_set_num_threads(core_num);
#pragma omp parallel for
  for (int element = 0; element < position_input_rows; ++element){
    int i = (position_input(element, 1))-1;
    int j = (position_input(element, 0))-1;
    if (element % 100 == 0) {
      std::cout << "Iter for " << i << "," << j << std::endl;
    }
    arma::vec data_1 = data_use.col(i);
    double mu1 = mlemu1(i);
    double sigma1 = mlesigmahat1(i,i);
    arma::vec data_2 = data_use.col(j);
    double mu2 = mlemu1(j);
    double sigma2 = mlesigmahat1(j,j);
    //
    arma::vec mu_vec(2);
    mu_vec(0) = mu1;
    mu_vec(1) = mu2;
    //
    double sigmaint0 = mlesigmahat1(i,j);
    arma::mat sigma_matrix0(2,2);
    sigma_matrix0(0,0) = sigma1;
    sigma_matrix0(1,1) = sigma2;
    sigma_matrix0(1,0) = sigmaint0;
    sigma_matrix0(0,1) = sigmaint0;

    //
    double det_sigma_matrix0 = det(sigma_matrix0);
    if(det_sigma_matrix0 < 0){
      sigmaint0 = arma::sign(sigmaint0) * (sqrt(sigma1 * sigma2) * 0.99);
      sigma_matrix0(1,0) = sigmaint0;
      sigma_matrix0(0,1) = sigmaint0;
    }
    //double sigma_new1=sqrt(sigma1-(sigmaint0*sigmaint0/sigma2));
    //double sigma_new2=sqrt(sigma2-(sigmaint0*sigmaint0/sigma1));
    //
    arma::mat sigma_inv0 = inv_sympd(sigma_matrix0);
    arma::mat sqrt_sigma(2,2);
    sqrt_sigma(0,0) = sqrt(sigma1);
    sqrt_sigma(1,1) = sqrt((sigma1*sigma2-sigmaint0*sigmaint0)/sigma1);
    sqrt_sigma(1,0) = sigmaint0/sqrt(sigma1);
    sqrt_sigma(0,1) = 0;
    //z.root_mat(10)
    //arma::vec z1_root = t_root_vec1 * sigma_new1 + mu1;
    //arma::vec z2_root = t_root_vec1 * sigma_new2 + mu2;
    arma::vec z1_root = t_root_vec1;
    arma::vec z2_root = t_root_vec1;
    double det_sigma_matrix1 = det(sigma_matrix0);
    double sqrt_det_sigma=sqrt(det_sigma_matrix1);
    arma::vec m1_root = oemga_troot_vec;
    arma::vec m2_root = oemga_troot_vec;
    //
    arma::cube gz_up_sigma_int(sample_size_use,length1,length2);
    arma::cube gz_down_int(sample_size_use,length1,length2);
    arma::cube gz_up_sigma_2_int(sample_size_use,length1,length2);
    //
    arma::cube kkk1_array(sample_size_use,length1,length2);
    arma::vec col_index(2);
    col_index(0) = i;
    col_index(1) = j;


    for(int u = 0;u<length1;++u){
      arma::mat z_root_vec(2, length2);
      for(int v = 0;v<length2;++v){
        z_root_vec(0,v) = z1_root(u);
        z_root_vec(1,v) = z2_root(v);
      }
      z_root_vec=sqrt_sigma*z_root_vec;
      arma::mat z_root_vec_center = z_root_vec;
      z_root_vec.row(0) = z_root_vec.row(0) + mu_vec(0);
      z_root_vec.row(1) = z_root_vec.row(1) + mu_vec(1);
      //z_root_vec.each_col() += mu_vec;
      //
      arma::mat data_use_sub(sample_size_use,2);
      data_use_sub.col(0) = data_use.col(i);
      data_use_sub.col(1) = data_use.col(j);
      arma::vec temp_z_root=(z_root_vec.row(1)).t();
      //
      arma::mat temp_mat1(1,length2,arma::fill::ones);
      arma::mat res_mat = data_use_sub * z_root_vec - (S_depth * exp(z_root_vec(0,0))) * temp_mat1 - S_depth * (exp(temp_z_root)).t();
      //root1_array
      kkk1_array.subcube(0,u,0, (sample_size_use - 1),u,(length2 - 1)) = res_mat;
    }
    //max_mat
    arma::vec kkk1(sample_size_use);
    for(int sample_index = 0;sample_index < sample_size_use; ++sample_index){
      arma::mat mat_temp1 = kkk1_array.subcube(sample_index,0,0, sample_index,(length1 - 1),(length2 - 1));
      kkk1(sample_index) = mat_temp1.max();
    }
    for(int u = 0;u<length1;++u){
      arma::mat z_root_vec(2, length2);
      for(int v = 0;v<length2;++v){
        z_root_vec(0,v) = z1_root(u);
        z_root_vec(1,v) = z2_root(v);
      }
      z_root_vec=sqrt_sigma*z_root_vec;
      arma::mat z_root_vec_center = z_root_vec;
      z_root_vec.row(0) = z_root_vec.row(0) + mu_vec(0);
      z_root_vec.row(1) = z_root_vec.row(1) + mu_vec(1);
      //
      arma::vec bbb_2(length2);
      arma::vec bbb_3(length2);
      for(int v = 0;v<length2;++v){
        arma::mat mat_test = sigma_inv0 * (z_root_vec_center.col(v)) * (z_root_vec_center.col(v)).t() * sigma_inv0;
        bbb_2(v) = mat_test(0,1);
        arma::mat second_part11=(z_root_vec_center.col(v)) * (z_root_vec_center.col(v)).t() * sigma_inv0 * ((a_mat) *(b_mat).t() + (b_mat)* (a_mat).t());
        arma::mat second_part12=sigma_inv0 * ((a_mat)* (b_mat).t()+(b_mat) * (a_mat).t())*sigma_inv0;
        double second_part1 = trace(second_part11 *second_part12);
        bbb_3(v)=second_part1;
      }
      //
      arma::mat kkk2_0 = kkk1_array.subcube(0,u,0, (sample_size_use - 1),u,(length2 - 1));
      //kkk2_0.n_cols is 10
      for(int sample_index = 0; sample_index<sample_size_use;++sample_index){
        kkk2_0.row(sample_index) = kkk2_0.row(sample_index) - kkk1(sample_index);
      }
      arma::mat kkk2 = (exp(kkk2_0)).t();
      //kkk2.n_cols is n
      //
      arma::mat gz_up_sigma_int_res = kkk2;
      arma::mat gz_down_int_res = kkk2;
      arma::mat gz_up_sigma_2_int_res = kkk2;
      for(int v = 0;v<length2;++v){
        gz_up_sigma_int_res.row(v) = gz_up_sigma_int_res.row(v) * (bbb_2(v));
        gz_down_int_res.row(v) = gz_down_int_res.row(v);
        gz_up_sigma_2_int_res.row(v) = gz_up_sigma_2_int_res.row(v) * (((-1) * bbb_3(v) + bbb_2(v) * bbb_2(v)));
      }
      //
      gz_up_sigma_int.subcube(0,u,0, (sample_size_use - 1),u,(length2 - 1)) = gz_up_sigma_int_res.t();
      gz_down_int.subcube(0,u,0, (sample_size_use - 1),u,(length2 - 1)) = gz_down_int_res.t();
      gz_up_sigma_2_int.subcube(0,u,0, (sample_size_use - 1),u,(length2 - 1)) = gz_up_sigma_2_int_res.t();

    }
    arma::vec gradiant_up_sigma_int(sample_size_use);
    arma::vec gradiant_down_int(sample_size_use);
    arma::vec Hessian_up_sigma_int(sample_size_use);
    arma::vec gradiant_int0(sample_size_use);
    arma::vec hessian_int0(sample_size_use);
    //
    for(int sample_index = 0;sample_index < sample_size_use; ++sample_index){
      arma::mat mat_incube1 = gz_up_sigma_int.subcube(sample_index,0,0, sample_index,(length1 - 1),(length2 - 1));
      arma::mat mat_incube2 = gz_down_int.subcube(sample_index,0,0, sample_index,(length1 - 1),(length2 - 1));
      arma::mat mat_incube3 = gz_up_sigma_2_int.subcube(sample_index,0,0, sample_index,(length1 - 1),(length2 - 1));
      //
      arma::mat mat_temp1 = (m1_root).t() * mat_incube1 * m2_root;
      double value1 = (double)(mat_temp1(0,0));
      arma::mat mat_temp2 = (m1_root).t() * mat_incube2 * m2_root;
      double value2 = (double)(mat_temp2(0,0));
      arma::mat mat_temp3 = (m1_root).t() * mat_incube3 * m2_root;
      double value3 = (double)(mat_temp3(0,0));

      gradiant_up_sigma_int(sample_index) = value1;
      //
      gradiant_down_int(sample_index) = value2;
      //
      Hessian_up_sigma_int(sample_index) = value3;

      gradiant_int0(sample_index) = (value1/value2-sigma_inv0(0,1))*cond_Z(sample_index);
      hessian_int0(sample_index) = (value3/value2-(value1/value2)*(value1/value2)+(sigma_inv0(0,0) * sigma_inv0(1,1) + sigma_inv0(0,1) * sigma_inv0(0,1)))*cond_Z(sample_index);
    }
    //
    double weight_temp= sum(cond_Z)/(double(sample_size_use));
    //gradiant_int(i,j) = sum((gradiant_up_sigma_int/gradiant_down_int)%cond_Z)/(double(sample_size_use)) - weight_temp*sigma_inv0(0,1);
    //
    //Hessian_int(i,j) = sum((Hessian_up_sigma_int /gradiant_down_int)%cond_Z)/(double(sample_size_use)) -
    //  sum((pow(gradiant_up_sigma_int /gradiant_down_int,2))%cond_Z)/(double(sample_size_use)) +
    //  (sigma_inv0(0,0) * sigma_inv0(1,1) + sigma_inv0(0,1) * sigma_inv0(0,1))*weight_temp;
    //
    gradiant_int(i,j)=sum(gradiant_int0);
    Hessian_int(i,j)=sum(hessian_int0);
  }
  omp_destroy_lock(&lock_temp);
  //
  return Rcpp::List::create(Rcpp::Named("gradiant_int") = gradiant_int,
                            Rcpp::Named("Hessian_int") = Hessian_int);

}
// [[Rcpp::export]]
arma::mat basic_mmle(arma::mat data_use, arma::vec S_depth, arma::vec mlemu1, arma::mat mlesigmahat1,
                       arma::mat gradiant_int, arma::mat Hessian_int,int core_num){
  int dim_use = data_use.n_cols;
  int sample_size_use = data_use.n_rows;
  //
  arma::mat sigma0_mat(dim_use,dim_use,fill::zeros);
  //
  //
  arma::mat mlesigmahat_res(dim_use,dim_use,fill::zeros);
  sigma0_mat = mlesigmahat1;
  omp_init_lock(&lock_temp);
  omp_set_lock(&lock_temp);

  omp_set_num_threads(core_num);
#pragma omp parallel for
  for(int ii = 0;ii<(dim_use - 1);++ii){
    arma::vec data1 = data_use.col(ii);
    double mu1 = mlemu1(ii);
    for(int jj = (ii+1);jj<dim_use;++jj){
      arma::vec data2 = data_use.col(jj);
      double mu2 = mlemu1(jj);
      //
      arma::vec mu_vec(2);
      mu_vec(0) = mu1;
      mu_vec(1) = mu2;
      double sigma1 = mlesigmahat1(ii,ii);
      double sigma2 = mlesigmahat1(jj,jj);
      // arma::vec vec1 = (data1 % data2)/pow(S_depth,2);
      // double sigmaint0 = log(sum(vec1)/(double)sample_size_use) - (mu1 + sigma1*0.5) - (mu2 + sigma2*0.5);
      double sigmaint0 = sigma0_mat(ii,jj);

      arma::mat sigma_matrix0(2,2);
      sigma_matrix0(0,0) = sigma1;
      sigma_matrix0(1,1) = sigma2;
      sigma_matrix0(1,0) = sigmaint0;
      sigma_matrix0(0,1) = sigmaint0;

      double det_sigma_matrix0 = det(sigma_matrix0);
      if(det_sigma_matrix0 < 0){
        sigmaint0 = sign(sigmaint0) * (sqrt(sigma1 * sigma2) * 0.99);
        sigma_matrix0(1,0) = sigmaint0;
        sigma_matrix0(0,1) = sigmaint0;
      }

      //
      mlesigmahat_res(ii,jj) = sigmaint0 - gradiant_int(ii,jj)/Hessian_int(ii,jj);
      //
      arma::mat sigma_matrix1(2,2);
      sigma_matrix1(0,0)=sigma1;
      sigma_matrix1(1,1)=sigma2;
      sigma_matrix1(1,0)=mlesigmahat_res(ii,jj);
      sigma_matrix1(0,1)=mlesigmahat_res(ii,jj);

      double det_sigma_matrix1 = det (sigma_matrix1);
      if ( det_sigma_matrix1<0 ){
        mlesigmahat_res(ii,jj)=sigmaint0;
      }
    }
  }
  omp_destroy_lock(&lock_temp);
  //
  return(mlesigmahat_res);
}

arma::mat EigenToArma_mat(Eigen::MatrixXd eigen_A) {
  arma::mat arma_B = arma::mat(eigen_A.data(), eigen_A.rows(), eigen_A.cols());
  return arma_B;
}

arma::vec EigenToArma_vec(Eigen::VectorXd eigenVec) {
  arma::vec armaVec(eigenVec.data(), eigenVec.size());
  return armaVec;
}

// [[Rcpp::export]]
Rcpp::List mmle_newton_fun(Eigen::MatrixXd data_use,Eigen::VectorXd cond_Z_input, Eigen::VectorXd S_depth, Eigen::VectorXd res_init_mu,
                              Eigen::MatrixXd res_init_sigma,Eigen::VectorXd t_root_vec, Eigen::VectorXd omega_root_vec,int k_max, int core_num,arma::mat position_input,int cut_num){
  //covert to arma
  arma::mat data_use_arma = EigenToArma_mat(data_use);
  arma::vec cond_Z_input_arma = EigenToArma_vec(cond_Z_input);
  arma::vec S_depth_arma = EigenToArma_vec(S_depth);
  arma::vec res_init_mu_arma = EigenToArma_vec(res_init_mu);
  arma::mat res_init_sigma_arma = EigenToArma_mat(res_init_sigma);
  arma::vec t_root_vec_arma = EigenToArma_vec(t_root_vec);
  arma::vec omega_root_vec_arma = EigenToArma_vec(omega_root_vec);
  Rcpp::List temp=mmle_weighted_step1(data_use_arma,cond_Z_input_arma,S_depth_arma,res_init_mu_arma,res_init_sigma_arma,t_root_vec_arma,omega_root_vec_arma,k_max,core_num,cut_num);
  Eigen::VectorXd convergence_vec=Rcpp::as<VectorXd>(temp["if_convergence_feature"]);
  arma::vec convergence_vec_arma=EigenToArma_vec(convergence_vec);
  Eigen::MatrixXd update_iter_mat = Rcpp::as<MatrixXd>(temp["update_iter_mat"]);
  arma::mat update_iter_arma=EigenToArma_mat(update_iter_mat);
  Rcpp::List mu_sigma=mmle_weighted_step2(data_use_arma,cond_Z_input_arma,S_depth_arma,convergence_vec_arma,update_iter_arma,
                                          res_init_mu_arma,res_init_sigma_arma,t_root_vec_arma,omega_root_vec_arma,k_max,core_num,cut_num);
  arma::vec mlemu1 = Rcpp::as<arma::vec>(mu_sigma["mlemu"]);
  arma::mat mlesigmahat1 = res_init_sigma_arma;

  arma::vec mlesigmahat1_diag = Rcpp::as<arma::vec>(mu_sigma["mlesigmahat"]);

  mlesigmahat1.diag() = mlesigmahat1_diag;

  bool allDiagonalPositive = arma::all(mlesigmahat1_diag > 0);
  if (!allDiagonalPositive) {
    cout << "the diagonal of sigmahat is not all positive" << endl;
  }
  Rcpp::List integrate_list=Update_sigma12_fun(data_use_arma,cond_Z_input_arma, S_depth_arma,mlemu1, mlesigmahat1,t_root_vec_arma, omega_root_vec_arma,core_num,position_input);

  arma::mat gradiant_int=Rcpp::as<arma::mat>(integrate_list["gradiant_int"]);
  arma::mat Hessian_int=Rcpp::as<arma::mat>(integrate_list["Hessian_int"]);
  gradiant_int.elem(find_nonfinite(gradiant_int)).zeros();
  Hessian_int.elem(find_nonfinite(Hessian_int)).ones();
  arma::mat mlesigmahat_res=basic_mmle(data_use_arma,  S_depth_arma, mlemu1,mlesigmahat1,gradiant_int,Hessian_int,core_num);

  mlesigmahat_res = (mlesigmahat_res + mlesigmahat_res.t());
  // Set diagonal of mlesigmahat1 to mlesigmahat1_diag
  mlesigmahat_res.diag() = mlesigmahat1_diag;

  return Rcpp::List::create(Rcpp::Named("mlemu") = mlemu1,
                            Rcpp::Named("mlesigmahat") = mlesigmahat_res);
}

