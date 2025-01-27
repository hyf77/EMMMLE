#include <RcppArmadillo.h>
#include <omp.h>
#include <Eigen/Sparse>
#include <math.h>
#include <RcppEigen.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace std;
using namespace arma;
using namespace Eigen;

static omp_lock_t lock_temp_emmle;


// [[Rcpp::export]]
arma::mat log_cul_mat_fun2(arma::mat mat_input, int core_num) {
  int num_row = mat_input.n_rows;
  int num_col = mat_input.n_cols;
  arma::mat mat_res(num_row, num_col);
  for (int i = 0; i < num_row; ++i) {
    for (int j = 0; j < num_col; ++j) {
      double value = mat_input(i, j);
      if (value > 0) {
        arma::vec vec_cur(value);
        for (int k = 0; k < value; ++k) {
          vec_cur[k] = k + 1;
        }
        mat_res(i, j) = arma::sum(arma::log(vec_cur));
      }
    }
  }
  return mat_res;
}

// [[Rcpp::export]]
arma::cube likelihood_weight(arma::mat data_use, arma::vec S_depth, arma::vec mlemu1, arma::mat mlesigmahat1,
                       arma::vec t_root_vec, arma::vec omega_root_vec, int core_num){
  //def
  arma::mat a_mat("0 ; 1");
  arma::mat b_mat("1 ; 0");
  int dim_use = data_use.n_cols;
  int sample_size_use = data_use.n_rows;
  arma::vec t_root_vec1 = t_root_vec * sqrt(2);
  arma::vec oemga_troot_vec = omega_root_vec;
  int length1 = t_root_vec1.n_elem;
  int length2 = oemga_troot_vec.n_elem;
  //
  //arma::mat gradiant_int(dim_use,dim_use,arma::fill::zeros);
  //arma::mat Hessian_int(dim_use,dim_use,arma::fill::zeros);
  arma::cube kkk1_array(sample_size_use,length1,length2);
  //
  for (int i=0;i<(dim_use-1);++i) {
    arma::vec data_1 = data_use.col(i);
    double mu1 = mlemu1(i);
    double sigma1 = mlesigmahat1(i,i);
    for (int j=(i+1);j<dim_use;++j) {
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
      //
    }
  }
  //
  return kkk1_array;

}

// [[Rcpp::export]]
double log_likelihood_fun(arma::mat data_input,  arma::mat mu_mat, arma::cube sigma_list,arma::vec pi_vec,
                            arma::vec S_depth,arma::vec t_root_vec,arma::vec omega_root_vec) {
  int dim_num=mu_mat.n_cols;
  int group_num = mu_mat.n_rows;
  int sample_num = data_input.n_rows;
  arma::mat cond_pro_X_mat(group_num, sample_num);
  int length_vec= t_root_vec.n_elem;
  arma::mat log_Y_prod = log_cul_mat_fun2(data_input,10);
  arma::vec log_Y_prod_ij_row = arma::sum(log_Y_prod, 1);
  arma::mat likelihood_mat(group_num, sample_num);
  arma::vec row_datasum=arma::sum(data_input, 1);
  for (int g = 0; g < group_num; ++g) {
    arma::vec mu_group = (mu_mat.row(g)).t();
    arma::mat sigma_group = sigma_list.subcube(0,0,g,(dim_num - 1),(dim_num - 1),g);
    arma::vec part_log_sum = row_datasum  % arma::log(S_depth) - log_Y_prod_ij_row;
    arma::cube gz_down_temp = likelihood_weight(data_input, S_depth, mu_group, sigma_group, t_root_vec, omega_root_vec, 10);
    arma::vec temp_sum(sample_num);
    for (int i = 0; i < sample_num; ++i) {
      arma::mat temp_w=gz_down_temp.subcube(i,0,0, i,(length_vec - 1),(length_vec - 1));
      double max_const = (temp_w).max();
      arma::mat temp_q=gz_down_temp.subcube(i,0,0, i,(length_vec - 1),(length_vec - 1))-max_const;
      arma::mat temp_sum0 = (omega_root_vec).t() * exp(temp_q) * omega_root_vec;
      double value2 = (double)(temp_sum0(0,0));
      temp_sum(i) = std::log(value2 / M_PI) + max_const + part_log_sum(i);
    }

    likelihood_mat.row(g) = (exp(temp_sum) * pi_vec(g)).t();
  }
  double log_likelihood_res = arma::sum(arma::log(arma::sum(likelihood_mat, 0)));

  return log_likelihood_res;
}

// [[Rcpp::export]]
arma::mat compute_weight_Z(arma::mat data_input,arma::vec weight_init,arma::mat mlemu_init, arma::cube mlesigmahat_init,arma::vec S_depth, arma::vec t_root_vec,arma::vec omega_root_vec) {
  int length_vec= t_root_vec.n_elem;
  int dim_num=mlemu_init.n_cols;
  int group_num = mlemu_init.n_rows;
  int sample_size = data_input.n_rows;
  arma::mat cond_pro_X_mat(sample_size, group_num);
  arma::cube gz_down_list(group_num*sample_size, length_vec,length_vec);
  for (int g = 0; g < group_num; ++g) {
    arma::vec res_init_mu = (mlemu_init.row(g)).t();
    arma::mat res_init_sigma = mlesigmahat_init.subcube(0,0,g,(dim_num - 1),(dim_num - 1),g);
    arma::cube weighted_update_res = likelihood_weight(data_input, S_depth, res_init_mu, res_init_sigma, t_root_vec, omega_root_vec, 10);
    gz_down_list.subcube(g*sample_size,0,0,(g*sample_size+sample_size-1),(length_vec - 1),(length_vec-1))= weighted_update_res;
  }

  for (int i = 0; i < sample_size; ++i) {
    arma::cube sample_group(group_num,length_vec,length_vec);
    for (int g = 0; g < group_num; ++g) {
      sample_group.subcube(g,0,0,g,(length_vec - 1),(length_vec-1))=gz_down_list.subcube((g*sample_size+i),0,0,(g*sample_size+i),(length_vec - 1),(length_vec-1));
    }
    double max_value = sample_group.max();
    arma::cube temp_array = sample_group- max_value;
    for (int g = 0; g < group_num; ++g) {
      arma::mat temo_mat=temp_array.subcube(g,0,0,g,(length_vec - 1),(length_vec-1));
      arma::mat temp_sum0 = (omega_root_vec).t() * exp(temo_mat) * omega_root_vec;
      double value2 = (double)(temp_sum0(0,0));
      cond_pro_X_mat(i,g) = value2 * weight_init(g);
    }
  }

  arma::vec sum_down = arma::sum(cond_pro_X_mat, 1);
  arma::mat cond_pro_Z(sample_size, group_num);

  for (int g = 0; g < group_num; ++g) {
    cond_pro_Z.col(g) = cond_pro_X_mat.col(g) / sum_down;
  }

  return cond_pro_Z;
}


// [[Rcpp::export]]
arma::mat EMMMLE_weighted_step1(arma::mat data_use, arma::vec cond_Z_input, arma::vec S_depth, arma::vec mu_init, arma::mat sigma_init,
                                arma::vec t_root_vec, arma::vec omega_root_vec,int k_max, int core_num, int cut_num){
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
  //
  arma::mat res_mat(k_max+1,dim_use);
  res_mat.submat(0, 0, (k_max-1), (dim_use-1))=update_iter_mat;
  res_mat.submat(k_max,0,k_max,(dim_use-1))=if_convergence_feature.t();
  return res_mat;

}

// [[Rcpp::export]]
arma::mat EMMMLE_weighted_step2(arma::mat data_use, arma::vec cond_Z_input, arma::vec S_depth,arma::mat res_mmle1, arma::vec mu_init, arma::mat sigma_init,
                                arma::vec t_root_vec, arma::vec omega_root_vec,int k_max, int core_num, int cut_num){
  int dim_use = data_use.n_cols;
  int sample_size = data_use.n_rows;
  int length = t_root_vec.n_elem;
  arma::mat mlemu(1,dim_use);
  arma::mat mlesigmahat(1,dim_use);
  arma::vec mlemu1(dim_use);
  arma::mat mlesigmahat1(dim_use,dim_use);
  arma::mat update_iter_mat=res_mmle1.submat(0, 0, (k_max-1), (dim_use-1));
  arma::vec if_convergence_feature=(res_mmle1.row(k_max)).t();
  //arma::mat update_iter_mat=update_mat;
  arma::vec t_root_vec1 = t_root_vec * sqrt(2);
  arma::vec oemga_troot_vec = omega_root_vec % exp(pow(t_root_vec,2)) * sqrt(2) ;
  //arma::vec if_convergence_feature=convergence_vec;
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
  arma::mat res_end2(k_max+2,dim_use);
  res_end2.row(0)=mlemu1.t();
  res_end2.row(1)=sigmadiag_end.t();
  res_end2.submat(2, 0, (k_max+1), (dim_use-1))=update_iter_mat;

  return res_end2;
  // return Rcpp::List::create(Rcpp::Named("mlemu") = mlemu1,
  //                           Rcpp::Named("update_iter_mat")=update_iter_mat,
  //                           Rcpp::Named("mlesigmahat") = mlesigmahat1);
}

// [[Rcpp::export]]
arma::mat EMMMLE_sigma12_fun(arma::mat data_use, arma::vec cond_Z_input, arma::vec S_depth, arma::vec mlemu1, arma::mat mlesigmahat1,
                               arma::vec t_root_vec, arma::vec omega_root_vec,int core_num){
  //def
  arma::mat a_mat("0 ; 1");
  arma::mat b_mat("1 ; 0");
  arma::vec cond_Z=cond_Z_input;
  int dim_use = data_use.n_cols;
  int sample_size_use = data_use.n_rows;
  arma::vec t_root_vec1 = t_root_vec * sqrt(2);
  arma::vec oemga_troot_vec = omega_root_vec;
  int length1 = t_root_vec1.n_elem;
  int length2 = oemga_troot_vec.n_elem;
  //
  arma::mat gradiant_int(dim_use,dim_use,arma::fill::zeros);
  arma::mat Hessian_int(dim_use,dim_use,arma::fill::zeros);
  //
  for (int i=0;i<(dim_use-1);++i) {
    arma::vec data_1 = data_use.col(i);
    double mu1 = mlemu1(i);
    double sigma1 = mlesigmahat1(i,i);
    for (int j=(i+1);j<dim_use;++j) {
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
  }
  //
  //arma::mat gradiant_int(dim_use,dim_use,arma::fill::zeros);
  //arma::mat Hessian_int(dim_use,dim_use,arma::fill::zeros);
  arma::mat res_end(2*dim_use,dim_use,arma::fill::zeros);
  res_end.submat(0,0,(dim_use-1),(dim_use-1))=gradiant_int;
  res_end.submat(dim_use,0,(2*dim_use-1),(dim_use-1))=Hessian_int;
  return res_end;
}

arma::mat basic_EMMMLE(arma::mat data_use, arma::vec S_depth, arma::vec mlemu1, arma::mat mlesigmahat1,
                       arma::mat gradiant_int, arma::mat Hessian_int,int core_num){
  int dim_use = data_use.n_cols;
  int sample_size_use = data_use.n_rows;
  //
  arma::mat sigma0_mat(dim_use,dim_use,fill::zeros);
  //
  //
  arma::mat mlesigmahat_res(dim_use,dim_use,fill::zeros);
  sigma0_mat = mlesigmahat1;

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
  //
  return(mlesigmahat_res);
}

// [[Rcpp::export]]
arma::mat EMMMLE_newton_fun(arma::mat data_use_arma,arma::vec cond_Z_input_arma, arma::vec S_depth_arma, arma::vec res_init_mu_arma,
                            arma::mat res_init_sigma_arma,arma::vec t_root_vec_arma, arma::vec omega_root_vec_arma,int k_max, int core_num,int cut_num){
  //arma::mat res_mat(k_max+1,dim_use);
  //res_mat.submat(0, 0, (k_max-1), (dim_use-1))=update_iter_mat;
  //res_mat.submat(k_max,0,k_max,(dim_use-1))=if_convergence_feature.t();
  int dim_use = data_use_arma.n_cols;
  arma::mat res_newton(dim_use+1,dim_use);
  arma::mat temp=EMMMLE_weighted_step1(data_use_arma,cond_Z_input_arma,S_depth_arma,res_init_mu_arma,res_init_sigma_arma,t_root_vec_arma,omega_root_vec_arma,k_max,core_num,cut_num);
  arma::mat mu_sigma=EMMMLE_weighted_step2(data_use_arma,cond_Z_input_arma,S_depth_arma,temp,
                                           res_init_mu_arma,res_init_sigma_arma,t_root_vec_arma,omega_root_vec_arma,k_max,core_num,cut_num);
  //arma::mat res_end2(k_max+2,dim_use);
  //res_end2.row(0)=mlemu;
  //res_end2.row(1)=mlesigmahat;
  //res_end2.submat(2, 0, (k_max+1), (dim_use-1))=update_iter_mat;
  res_newton.row(0)=mu_sigma.row(0);
  arma::vec mlemu1 = (mu_sigma.row(0)).t();
  arma::mat mlesigmahat1 = res_init_sigma_arma;

  arma::vec mlesigmahat1_diag = (mu_sigma.row(1)).t();

  mlesigmahat1.diag() = mlesigmahat1_diag;

  bool allDiagonalPositive = arma::all(mlesigmahat1_diag > 0);
  if (!allDiagonalPositive) {
    cout << "the diagonal of sigmahat is not all positive" << endl;
  }
  arma::mat integrate_list=EMMMLE_sigma12_fun(data_use_arma,cond_Z_input_arma, S_depth_arma,mlemu1, mlesigmahat1,t_root_vec_arma, omega_root_vec_arma,core_num);
  //res_end.submat(0,0,(dim_use-1),(dim_use-1))=gradiant_int;
  //res_end.submat(dim_use,0,(2*dim_use-1),(dim_use-1))=Hessian_int;

  arma::mat gradiant_int=integrate_list.submat(0,0,(dim_use-1),(dim_use-1));
  arma::mat Hessian_int=integrate_list.submat(dim_use,0,(2*dim_use-1),(dim_use-1));
  gradiant_int.elem(find_nonfinite(gradiant_int)).zeros();
  Hessian_int.elem(find_nonfinite(Hessian_int)).ones();
  arma::mat mlesigmahat_res=basic_EMMMLE(data_use_arma,  S_depth_arma, mlemu1,mlesigmahat1,gradiant_int,Hessian_int,core_num);

  mlesigmahat_res = (mlesigmahat_res + mlesigmahat_res.t());
  // Set diagonal of mlesigmahat1 to mlesigmahat1_diag
  mlesigmahat_res.diag() = mlesigmahat1_diag;
  res_newton.submat(1, 0, (dim_use), (dim_use-1))=mlesigmahat_res;

  return res_newton;
}

// [[Rcpp::export]]
arma::cube EMMMLE_main(arma::mat data_input,arma::vec pi_init,arma::mat mlemu_init, arma::cube mlesigmahat_init,arma::vec S_depth, arma::vec t_root_vec,arma::vec omega_root_vec,int iter_max,
                   double omega,arma::mat p_p_post_dist_Z,int k_max,int cut_num) {
  int group_num = mlemu_init.n_rows;
  int sample_size = data_input.n_rows;
  int dim_use = data_input.n_cols;
  arma::vec log_likelihood_vec(iter_max+1); 
  double log_likelihood_temp = log_likelihood_fun(data_input, mlemu_init, mlesigmahat_init, pi_init, S_depth,t_root_vec,omega_root_vec);
  log_likelihood_vec(0)=log_likelihood_temp;
  int iter = 1;
  arma::mat updated_mu(group_num, dim_use);
  arma::cube updated_sigma(dim_use,dim_use,group_num);
  arma::vec pi_update(group_num);
  while (iter <= iter_max) {
    arma::mat cond_pro_Z = compute_weight_Z(data_input, pi_init, mlemu_init, mlesigmahat_init, S_depth,t_root_vec,omega_root_vec);
    pi_update = (arma::mean(cond_pro_Z, 0)).t();
    //std::cout << arma::mean(cond_pro_Z, 0) << std::endl;
    for (int g = 0; g < group_num; ++g) {
      arma::vec post_dist_Z_iter = cond_pro_Z.col(g);
      arma::vec cond_Z_input = (1 - omega) * post_dist_Z_iter + omega * p_p_post_dist_Z.col(g);
      arma::vec res_init_mu = (mlemu_init.row(g)).t();
      arma::mat res_init_sigma = mlesigmahat_init.subcube(0,0, g,(dim_use - 1),(dim_use - 1),g);
      arma::mat mmle_newton_res = EMMMLE_newton_fun(data_input, cond_Z_input, S_depth, res_init_mu, res_init_sigma, t_root_vec, omega_root_vec, k_max, 20,cut_num);
      updated_mu.row(g) = mmle_newton_res.row(0);
      updated_sigma.subcube(0,0,g,(dim_use - 1),(dim_use - 1),g)= mmle_newton_res.submat(1, 0, (dim_use), (dim_use-1));
      //std::cout << g << std::endl;
      arma::mat cov_Y1 = mmle_newton_res.submat(1, 0, (dim_use), (dim_use-1));
      if (arma::det(cov_Y1) < 0) {
        while (arma::det(cov_Y1) < 0) {
          cov_Y1.diag() += 1e-5;
        }
      }
      updated_sigma.subcube(0,0, g,(dim_use - 1),(dim_use - 1),g)= cov_Y1;
    }

    double log_likelihood_temp0 = log_likelihood_fun(data_input, updated_mu, updated_sigma, pi_update, S_depth,t_root_vec,omega_root_vec);
    //double log_likelihood_temp0 = log_likelihood_fun(data_input, updated_mu, updated_sigma, pi_init, S_depth,t_root_vec,omega_root_vec);;

    //std::cout << log_likelihood_temp0 << std::endl;

    if (iter > 1 && log_likelihood_temp0 < log_likelihood_vec(iter-1)) {
      updated_mu = mlemu_init;
      updated_sigma = mlesigmahat_init;
      break;
    }

    log_likelihood_vec(iter)=log_likelihood_temp0;
    mlemu_init = updated_mu;
    mlesigmahat_init = updated_sigma;
    pi_init = pi_update;
    iter++;
  }
  arma::cube res_all(dim_use+1,dim_use,group_num);
  for (int g = 0; g < group_num; ++g) {
    res_all.subcube(0,0,g,0,(dim_use - 1),g)=updated_mu.row(g);
    res_all.subcube(1,0,g,dim_use,(dim_use - 1),g)=updated_sigma.subcube(0,0,g,(dim_use - 1),(dim_use - 1),g);
  }
  return res_all;
}


// [[Rcpp::export]]
Rcpp::List EMMMLE(arma::mat data_obs,arma::vec pi_init,arma::mat mu_initial_all,arma::cube sigma_inital_all,arma::vec S_depth,arma::vec t_root_vec,
                    arma::vec omega_root_vec,arma::mat p_p_post_dist_Z,int k_max,arma::mat position_input,int core_num,double omega,int iter_max,int cut_num) {
  int group_num = mu_initial_all.n_rows;
  int sample_size = data_obs.n_rows;
  int dim_use = data_obs.n_cols;
  int position_input_rows=position_input.n_rows;
  arma::cube res_mu_all(group_num, 2, position_input_rows);
  arma::mat res_sigma11_all(group_num, position_input_rows);
  arma::mat res_sigma22_all(group_num, position_input_rows);
  arma::mat res_sigma12_all(group_num, position_input_rows);
  omp_init_lock(&lock_temp_emmle);
  omp_set_lock(&lock_temp_emmle);
  //
  //
  omp_set_num_threads(core_num);
#pragma omp parallel for
  for (int element = 0; element < position_input_rows; ++element) {
    int i = (position_input(element, 1))-1;
    int j = (position_input(element, 0))-1;
    std::cout << "Iter for " << i << "," << j << std::endl;

    arma::mat data_input(sample_size,2);
    data_input.col(0)=data_obs.col(i);
    data_input.col(1)=data_obs.col(j);
    int dim_low = data_input.n_cols;
    arma::mat mu_old(group_num,dim_low);
    mu_old.col(0)=mu_initial_all.col(i);
    mu_old.col(1)=mu_initial_all.col(j);
    arma::cube sigma_old(dim_low, dim_low, group_num);
    for (int g = 0; g < group_num; ++g) {
      sigma_old(0,0,g) = sigma_inital_all(i,i,g);
      sigma_old(1,1,g) = sigma_inital_all(j,j,g);
      sigma_old(0,1,g) = sigma_inital_all(i,j,g);
      sigma_old(1,0,g) = sigma_inital_all(j,i,g);
    }

    arma::mat mlemu_init = mu_old;
    arma::cube mlesigmahat_init1 = sigma_old;
    arma::cube result_up = EMMMLE_main(data_input, pi_init, mlemu_init, mlesigmahat_init1, S_depth, t_root_vec, omega_root_vec, iter_max, omega, p_p_post_dist_Z, k_max,cut_num);

    for (int g = 0; g < group_num; ++g) {
      res_mu_all.subcube(g,0,element,g,1,element) = result_up.subcube(0,0,g,0,1,g);
      res_sigma11_all(g, element) = result_up(1,0,g);
      res_sigma22_all(g, element) = result_up(2,1,g);
      res_sigma12_all(g, element) = result_up(1,1,g);
    }

  }
  omp_destroy_lock(&lock_temp_emmle);
  return Rcpp::List::create(Rcpp::Named("res_mu_all") = res_mu_all,
                            Rcpp::Named("res_sigma11_all") = res_sigma11_all,
                            Rcpp::Named("res_sigma22_all") = res_sigma22_all,
                            Rcpp::Named("res_sigma12_all") = res_sigma12_all);
}
