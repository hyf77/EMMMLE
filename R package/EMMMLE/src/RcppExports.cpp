// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// log_cul_mat_fun2
arma::mat log_cul_mat_fun2(arma::mat mat_input, int core_num);
RcppExport SEXP _EMMMLE_log_cul_mat_fun2(SEXP mat_inputSEXP, SEXP core_numSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type mat_input(mat_inputSEXP);
    Rcpp::traits::input_parameter< int >::type core_num(core_numSEXP);
    rcpp_result_gen = Rcpp::wrap(log_cul_mat_fun2(mat_input, core_num));
    return rcpp_result_gen;
END_RCPP
}
// likelihood_weight
arma::cube likelihood_weight(arma::mat data_use, arma::vec S_depth, arma::vec mlemu1, arma::mat mlesigmahat1, arma::vec t_root_vec, arma::vec omega_root_vec, int core_num);
RcppExport SEXP _EMMMLE_likelihood_weight(SEXP data_useSEXP, SEXP S_depthSEXP, SEXP mlemu1SEXP, SEXP mlesigmahat1SEXP, SEXP t_root_vecSEXP, SEXP omega_root_vecSEXP, SEXP core_numSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data_use(data_useSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type S_depth(S_depthSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mlemu1(mlemu1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mlesigmahat1(mlesigmahat1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_root_vec(t_root_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type omega_root_vec(omega_root_vecSEXP);
    Rcpp::traits::input_parameter< int >::type core_num(core_numSEXP);
    rcpp_result_gen = Rcpp::wrap(likelihood_weight(data_use, S_depth, mlemu1, mlesigmahat1, t_root_vec, omega_root_vec, core_num));
    return rcpp_result_gen;
END_RCPP
}
// log_likelihood_fun
double log_likelihood_fun(arma::mat data_input, arma::mat mu_mat, arma::cube sigma_list, arma::vec pi_vec, arma::vec S_depth, arma::vec t_root_vec, arma::vec omega_root_vec);
RcppExport SEXP _EMMMLE_log_likelihood_fun(SEXP data_inputSEXP, SEXP mu_matSEXP, SEXP sigma_listSEXP, SEXP pi_vecSEXP, SEXP S_depthSEXP, SEXP t_root_vecSEXP, SEXP omega_root_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data_input(data_inputSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mu_mat(mu_matSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type sigma_list(sigma_listSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pi_vec(pi_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type S_depth(S_depthSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_root_vec(t_root_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type omega_root_vec(omega_root_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(log_likelihood_fun(data_input, mu_mat, sigma_list, pi_vec, S_depth, t_root_vec, omega_root_vec));
    return rcpp_result_gen;
END_RCPP
}
// compute_weight_Z
arma::mat compute_weight_Z(arma::mat data_input, arma::vec weight_init, arma::mat mlemu_init, arma::cube mlesigmahat_init, arma::vec S_depth, arma::vec t_root_vec, arma::vec omega_root_vec);
RcppExport SEXP _EMMMLE_compute_weight_Z(SEXP data_inputSEXP, SEXP weight_initSEXP, SEXP mlemu_initSEXP, SEXP mlesigmahat_initSEXP, SEXP S_depthSEXP, SEXP t_root_vecSEXP, SEXP omega_root_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data_input(data_inputSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weight_init(weight_initSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mlemu_init(mlemu_initSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type mlesigmahat_init(mlesigmahat_initSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type S_depth(S_depthSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_root_vec(t_root_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type omega_root_vec(omega_root_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_weight_Z(data_input, weight_init, mlemu_init, mlesigmahat_init, S_depth, t_root_vec, omega_root_vec));
    return rcpp_result_gen;
END_RCPP
}
// EMMMLE_weighted_step1
arma::mat EMMMLE_weighted_step1(arma::mat data_use, arma::vec cond_Z_input, arma::vec S_depth, arma::vec mu_init, arma::mat sigma_init, arma::vec t_root_vec, arma::vec omega_root_vec, int k_max, int core_num, int cut_num);
RcppExport SEXP _EMMMLE_EMMMLE_weighted_step1(SEXP data_useSEXP, SEXP cond_Z_inputSEXP, SEXP S_depthSEXP, SEXP mu_initSEXP, SEXP sigma_initSEXP, SEXP t_root_vecSEXP, SEXP omega_root_vecSEXP, SEXP k_maxSEXP, SEXP core_numSEXP, SEXP cut_numSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data_use(data_useSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type cond_Z_input(cond_Z_inputSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type S_depth(S_depthSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu_init(mu_initSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma_init(sigma_initSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_root_vec(t_root_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type omega_root_vec(omega_root_vecSEXP);
    Rcpp::traits::input_parameter< int >::type k_max(k_maxSEXP);
    Rcpp::traits::input_parameter< int >::type core_num(core_numSEXP);
    Rcpp::traits::input_parameter< int >::type cut_num(cut_numSEXP);
    rcpp_result_gen = Rcpp::wrap(EMMMLE_weighted_step1(data_use, cond_Z_input, S_depth, mu_init, sigma_init, t_root_vec, omega_root_vec, k_max, core_num, cut_num));
    return rcpp_result_gen;
END_RCPP
}
// EMMMLE_weighted_step2
arma::mat EMMMLE_weighted_step2(arma::mat data_use, arma::vec cond_Z_input, arma::vec S_depth, arma::mat res_mmle1, arma::vec mu_init, arma::mat sigma_init, arma::vec t_root_vec, arma::vec omega_root_vec, int k_max, int core_num, int cut_num);
RcppExport SEXP _EMMMLE_EMMMLE_weighted_step2(SEXP data_useSEXP, SEXP cond_Z_inputSEXP, SEXP S_depthSEXP, SEXP res_mmle1SEXP, SEXP mu_initSEXP, SEXP sigma_initSEXP, SEXP t_root_vecSEXP, SEXP omega_root_vecSEXP, SEXP k_maxSEXP, SEXP core_numSEXP, SEXP cut_numSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data_use(data_useSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type cond_Z_input(cond_Z_inputSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type S_depth(S_depthSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type res_mmle1(res_mmle1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu_init(mu_initSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma_init(sigma_initSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_root_vec(t_root_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type omega_root_vec(omega_root_vecSEXP);
    Rcpp::traits::input_parameter< int >::type k_max(k_maxSEXP);
    Rcpp::traits::input_parameter< int >::type core_num(core_numSEXP);
    Rcpp::traits::input_parameter< int >::type cut_num(cut_numSEXP);
    rcpp_result_gen = Rcpp::wrap(EMMMLE_weighted_step2(data_use, cond_Z_input, S_depth, res_mmle1, mu_init, sigma_init, t_root_vec, omega_root_vec, k_max, core_num, cut_num));
    return rcpp_result_gen;
END_RCPP
}
// EMMMLE_sigma12_fun
arma::mat EMMMLE_sigma12_fun(arma::mat data_use, arma::vec cond_Z_input, arma::vec S_depth, arma::vec mlemu1, arma::mat mlesigmahat1, arma::vec t_root_vec, arma::vec omega_root_vec, int core_num);
RcppExport SEXP _EMMMLE_EMMMLE_sigma12_fun(SEXP data_useSEXP, SEXP cond_Z_inputSEXP, SEXP S_depthSEXP, SEXP mlemu1SEXP, SEXP mlesigmahat1SEXP, SEXP t_root_vecSEXP, SEXP omega_root_vecSEXP, SEXP core_numSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data_use(data_useSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type cond_Z_input(cond_Z_inputSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type S_depth(S_depthSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mlemu1(mlemu1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mlesigmahat1(mlesigmahat1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_root_vec(t_root_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type omega_root_vec(omega_root_vecSEXP);
    Rcpp::traits::input_parameter< int >::type core_num(core_numSEXP);
    rcpp_result_gen = Rcpp::wrap(EMMMLE_sigma12_fun(data_use, cond_Z_input, S_depth, mlemu1, mlesigmahat1, t_root_vec, omega_root_vec, core_num));
    return rcpp_result_gen;
END_RCPP
}
// EMMMLE_newton_fun
arma::mat EMMMLE_newton_fun(arma::mat data_use_arma, arma::vec cond_Z_input_arma, arma::vec S_depth_arma, arma::vec res_init_mu_arma, arma::mat res_init_sigma_arma, arma::vec t_root_vec_arma, arma::vec omega_root_vec_arma, int k_max, int core_num, int cut_num);
RcppExport SEXP _EMMMLE_EMMMLE_newton_fun(SEXP data_use_armaSEXP, SEXP cond_Z_input_armaSEXP, SEXP S_depth_armaSEXP, SEXP res_init_mu_armaSEXP, SEXP res_init_sigma_armaSEXP, SEXP t_root_vec_armaSEXP, SEXP omega_root_vec_armaSEXP, SEXP k_maxSEXP, SEXP core_numSEXP, SEXP cut_numSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data_use_arma(data_use_armaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type cond_Z_input_arma(cond_Z_input_armaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type S_depth_arma(S_depth_armaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type res_init_mu_arma(res_init_mu_armaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type res_init_sigma_arma(res_init_sigma_armaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_root_vec_arma(t_root_vec_armaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type omega_root_vec_arma(omega_root_vec_armaSEXP);
    Rcpp::traits::input_parameter< int >::type k_max(k_maxSEXP);
    Rcpp::traits::input_parameter< int >::type core_num(core_numSEXP);
    Rcpp::traits::input_parameter< int >::type cut_num(cut_numSEXP);
    rcpp_result_gen = Rcpp::wrap(EMMMLE_newton_fun(data_use_arma, cond_Z_input_arma, S_depth_arma, res_init_mu_arma, res_init_sigma_arma, t_root_vec_arma, omega_root_vec_arma, k_max, core_num, cut_num));
    return rcpp_result_gen;
END_RCPP
}
// EMMMLE_main
arma::cube EMMMLE_main(arma::mat data_input, arma::vec pi_init, arma::mat mlemu_init, arma::cube mlesigmahat_init, arma::vec S_depth, arma::vec t_root_vec, arma::vec omega_root_vec, int iter_max, double omega, arma::mat p_p_post_dist_Z, int k_max, int cut_num);
RcppExport SEXP _EMMMLE_EMMMLE_main(SEXP data_inputSEXP, SEXP pi_initSEXP, SEXP mlemu_initSEXP, SEXP mlesigmahat_initSEXP, SEXP S_depthSEXP, SEXP t_root_vecSEXP, SEXP omega_root_vecSEXP, SEXP iter_maxSEXP, SEXP omegaSEXP, SEXP p_p_post_dist_ZSEXP, SEXP k_maxSEXP, SEXP cut_numSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data_input(data_inputSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pi_init(pi_initSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mlemu_init(mlemu_initSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type mlesigmahat_init(mlesigmahat_initSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type S_depth(S_depthSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_root_vec(t_root_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type omega_root_vec(omega_root_vecSEXP);
    Rcpp::traits::input_parameter< int >::type iter_max(iter_maxSEXP);
    Rcpp::traits::input_parameter< double >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type p_p_post_dist_Z(p_p_post_dist_ZSEXP);
    Rcpp::traits::input_parameter< int >::type k_max(k_maxSEXP);
    Rcpp::traits::input_parameter< int >::type cut_num(cut_numSEXP);
    rcpp_result_gen = Rcpp::wrap(EMMMLE_main(data_input, pi_init, mlemu_init, mlesigmahat_init, S_depth, t_root_vec, omega_root_vec, iter_max, omega, p_p_post_dist_Z, k_max, cut_num));
    return rcpp_result_gen;
END_RCPP
}
// EMMMLE
Rcpp::List EMMMLE(arma::mat data_obs, arma::vec pi_init, arma::mat mu_initial_all, arma::cube sigma_inital_all, arma::vec S_depth, arma::vec t_root_vec, arma::vec omega_root_vec, arma::mat p_p_post_dist_Z, int k_max, arma::mat position_input, int core_num, double omega, int iter_max, int cut_num);
RcppExport SEXP _EMMMLE_EMMMLE(SEXP data_obsSEXP, SEXP pi_initSEXP, SEXP mu_initial_allSEXP, SEXP sigma_inital_allSEXP, SEXP S_depthSEXP, SEXP t_root_vecSEXP, SEXP omega_root_vecSEXP, SEXP p_p_post_dist_ZSEXP, SEXP k_maxSEXP, SEXP position_inputSEXP, SEXP core_numSEXP, SEXP omegaSEXP, SEXP iter_maxSEXP, SEXP cut_numSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data_obs(data_obsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pi_init(pi_initSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mu_initial_all(mu_initial_allSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type sigma_inital_all(sigma_inital_allSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type S_depth(S_depthSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_root_vec(t_root_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type omega_root_vec(omega_root_vecSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type p_p_post_dist_Z(p_p_post_dist_ZSEXP);
    Rcpp::traits::input_parameter< int >::type k_max(k_maxSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type position_input(position_inputSEXP);
    Rcpp::traits::input_parameter< int >::type core_num(core_numSEXP);
    Rcpp::traits::input_parameter< double >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< int >::type iter_max(iter_maxSEXP);
    Rcpp::traits::input_parameter< int >::type cut_num(cut_numSEXP);
    rcpp_result_gen = Rcpp::wrap(EMMMLE(data_obs, pi_init, mu_initial_all, sigma_inital_all, S_depth, t_root_vec, omega_root_vec, p_p_post_dist_Z, k_max, position_input, core_num, omega, iter_max, cut_num));
    return rcpp_result_gen;
END_RCPP
}
// Update_log_cond_prob_X
Eigen::VectorXd Update_log_cond_prob_X(Eigen::VectorXd mu, Eigen::MatrixXd precision_mat, Eigen::MatrixXd Xdata, Eigen::MatrixXd scale, double eps, int max_iter, int core_num);
RcppExport SEXP _EMMMLE_Update_log_cond_prob_X(SEXP muSEXP, SEXP precision_matSEXP, SEXP XdataSEXP, SEXP scaleSEXP, SEXP epsSEXP, SEXP max_iterSEXP, SEXP core_numSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type mu(muSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type precision_mat(precision_matSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Xdata(XdataSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< int >::type core_num(core_numSEXP);
    rcpp_result_gen = Rcpp::wrap(Update_log_cond_prob_X(mu, precision_mat, Xdata, scale, eps, max_iter, core_num));
    return rcpp_result_gen;
END_RCPP
}
// soft
arma::mat soft(arma::mat A, arma::mat a_mat, int diag);
RcppExport SEXP _EMMMLE_soft(SEXP ASEXP, SEXP a_matSEXP, SEXP diagSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type a_mat(a_matSEXP);
    Rcpp::traits::input_parameter< int >::type diag(diagSEXP);
    rcpp_result_gen = Rcpp::wrap(soft(A, a_mat, diag));
    return rcpp_result_gen;
END_RCPP
}
// equal1
Rcpp::List equal1(arma::mat X, arma::vec lambda, arma::mat weight_mat, arma::mat zero_mat, double err, int maxIter, double rho, int diag, int core_num);
RcppExport SEXP _EMMMLE_equal1(SEXP XSEXP, SEXP lambdaSEXP, SEXP weight_matSEXP, SEXP zero_matSEXP, SEXP errSEXP, SEXP maxIterSEXP, SEXP rhoSEXP, SEXP diagSEXP, SEXP core_numSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type weight_mat(weight_matSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type zero_mat(zero_matSEXP);
    Rcpp::traits::input_parameter< double >::type err(errSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< int >::type diag(diagSEXP);
    Rcpp::traits::input_parameter< int >::type core_num(core_numSEXP);
    rcpp_result_gen = Rcpp::wrap(equal1(X, lambda, weight_mat, zero_mat, err, maxIter, rho, diag, core_num));
    return rcpp_result_gen;
END_RCPP
}
// equal2
Rcpp::List equal2(arma::mat X, arma::vec lambda, arma::Mat<double> weight_mat, arma::mat zero_mat, double err, int maxIter, double rho, int diag, int core_num);
RcppExport SEXP _EMMMLE_equal2(SEXP XSEXP, SEXP lambdaSEXP, SEXP weight_matSEXP, SEXP zero_matSEXP, SEXP errSEXP, SEXP maxIterSEXP, SEXP rhoSEXP, SEXP diagSEXP, SEXP core_numSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::Mat<double> >::type weight_mat(weight_matSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type zero_mat(zero_matSEXP);
    Rcpp::traits::input_parameter< double >::type err(errSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< int >::type diag(diagSEXP);
    Rcpp::traits::input_parameter< int >::type core_num(core_numSEXP);
    rcpp_result_gen = Rcpp::wrap(equal2(X, lambda, weight_mat, zero_mat, err, maxIter, rho, diag, core_num));
    return rcpp_result_gen;
END_RCPP
}
// log_cul_mat_fun1
Eigen::SparseMatrix<double> log_cul_mat_fun1(Eigen::SparseMatrix<int> mat);
RcppExport SEXP _EMMMLE_log_cul_mat_fun1(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<int> >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(log_cul_mat_fun1(mat));
    return rcpp_result_gen;
END_RCPP
}
// mmle_weighted_step1
Rcpp::List mmle_weighted_step1(arma::mat data_use, arma::vec cond_Z_input, arma::vec S_depth, arma::vec mu_init, arma::mat sigma_init, arma::vec t_root_vec, arma::vec omega_root_vec, int k_max, int core_num, int cut_num);
RcppExport SEXP _EMMMLE_mmle_weighted_step1(SEXP data_useSEXP, SEXP cond_Z_inputSEXP, SEXP S_depthSEXP, SEXP mu_initSEXP, SEXP sigma_initSEXP, SEXP t_root_vecSEXP, SEXP omega_root_vecSEXP, SEXP k_maxSEXP, SEXP core_numSEXP, SEXP cut_numSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data_use(data_useSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type cond_Z_input(cond_Z_inputSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type S_depth(S_depthSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu_init(mu_initSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma_init(sigma_initSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_root_vec(t_root_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type omega_root_vec(omega_root_vecSEXP);
    Rcpp::traits::input_parameter< int >::type k_max(k_maxSEXP);
    Rcpp::traits::input_parameter< int >::type core_num(core_numSEXP);
    Rcpp::traits::input_parameter< int >::type cut_num(cut_numSEXP);
    rcpp_result_gen = Rcpp::wrap(mmle_weighted_step1(data_use, cond_Z_input, S_depth, mu_init, sigma_init, t_root_vec, omega_root_vec, k_max, core_num, cut_num));
    return rcpp_result_gen;
END_RCPP
}
// mmle_weighted_step2
Rcpp::List mmle_weighted_step2(arma::mat data_use, arma::vec cond_Z_input, arma::vec S_depth, arma::vec convergence_vec, arma::mat update_mat, arma::vec mu_init, arma::mat sigma_init, arma::vec t_root_vec, arma::vec omega_root_vec, int k_max, int core_num, int cut_num);
RcppExport SEXP _EMMMLE_mmle_weighted_step2(SEXP data_useSEXP, SEXP cond_Z_inputSEXP, SEXP S_depthSEXP, SEXP convergence_vecSEXP, SEXP update_matSEXP, SEXP mu_initSEXP, SEXP sigma_initSEXP, SEXP t_root_vecSEXP, SEXP omega_root_vecSEXP, SEXP k_maxSEXP, SEXP core_numSEXP, SEXP cut_numSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data_use(data_useSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type cond_Z_input(cond_Z_inputSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type S_depth(S_depthSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type convergence_vec(convergence_vecSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type update_mat(update_matSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu_init(mu_initSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma_init(sigma_initSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_root_vec(t_root_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type omega_root_vec(omega_root_vecSEXP);
    Rcpp::traits::input_parameter< int >::type k_max(k_maxSEXP);
    Rcpp::traits::input_parameter< int >::type core_num(core_numSEXP);
    Rcpp::traits::input_parameter< int >::type cut_num(cut_numSEXP);
    rcpp_result_gen = Rcpp::wrap(mmle_weighted_step2(data_use, cond_Z_input, S_depth, convergence_vec, update_mat, mu_init, sigma_init, t_root_vec, omega_root_vec, k_max, core_num, cut_num));
    return rcpp_result_gen;
END_RCPP
}
// Update_sigma12_fun
Rcpp::List Update_sigma12_fun(arma::mat data_use, arma::vec cond_Z_input, arma::vec S_depth, arma::vec mlemu1, arma::mat mlesigmahat1, arma::vec t_root_vec, arma::vec omega_root_vec, int core_num, arma::mat position_input);
RcppExport SEXP _EMMMLE_Update_sigma12_fun(SEXP data_useSEXP, SEXP cond_Z_inputSEXP, SEXP S_depthSEXP, SEXP mlemu1SEXP, SEXP mlesigmahat1SEXP, SEXP t_root_vecSEXP, SEXP omega_root_vecSEXP, SEXP core_numSEXP, SEXP position_inputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data_use(data_useSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type cond_Z_input(cond_Z_inputSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type S_depth(S_depthSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mlemu1(mlemu1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mlesigmahat1(mlesigmahat1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_root_vec(t_root_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type omega_root_vec(omega_root_vecSEXP);
    Rcpp::traits::input_parameter< int >::type core_num(core_numSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type position_input(position_inputSEXP);
    rcpp_result_gen = Rcpp::wrap(Update_sigma12_fun(data_use, cond_Z_input, S_depth, mlemu1, mlesigmahat1, t_root_vec, omega_root_vec, core_num, position_input));
    return rcpp_result_gen;
END_RCPP
}
// basic_mmle
arma::mat basic_mmle(arma::mat data_use, arma::vec S_depth, arma::vec mlemu1, arma::mat mlesigmahat1, arma::mat gradiant_int, arma::mat Hessian_int, int core_num);
RcppExport SEXP _EMMMLE_basic_mmle(SEXP data_useSEXP, SEXP S_depthSEXP, SEXP mlemu1SEXP, SEXP mlesigmahat1SEXP, SEXP gradiant_intSEXP, SEXP Hessian_intSEXP, SEXP core_numSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data_use(data_useSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type S_depth(S_depthSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mlemu1(mlemu1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mlesigmahat1(mlesigmahat1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type gradiant_int(gradiant_intSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Hessian_int(Hessian_intSEXP);
    Rcpp::traits::input_parameter< int >::type core_num(core_numSEXP);
    rcpp_result_gen = Rcpp::wrap(basic_mmle(data_use, S_depth, mlemu1, mlesigmahat1, gradiant_int, Hessian_int, core_num));
    return rcpp_result_gen;
END_RCPP
}
// mmle_newton_fun
Rcpp::List mmle_newton_fun(Eigen::MatrixXd data_use, Eigen::VectorXd cond_Z_input, Eigen::VectorXd S_depth, Eigen::VectorXd res_init_mu, Eigen::MatrixXd res_init_sigma, Eigen::VectorXd t_root_vec, Eigen::VectorXd omega_root_vec, int k_max, int core_num, arma::mat position_input, int cut_num);
RcppExport SEXP _EMMMLE_mmle_newton_fun(SEXP data_useSEXP, SEXP cond_Z_inputSEXP, SEXP S_depthSEXP, SEXP res_init_muSEXP, SEXP res_init_sigmaSEXP, SEXP t_root_vecSEXP, SEXP omega_root_vecSEXP, SEXP k_maxSEXP, SEXP core_numSEXP, SEXP position_inputSEXP, SEXP cut_numSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type data_use(data_useSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type cond_Z_input(cond_Z_inputSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type S_depth(S_depthSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type res_init_mu(res_init_muSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type res_init_sigma(res_init_sigmaSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type t_root_vec(t_root_vecSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type omega_root_vec(omega_root_vecSEXP);
    Rcpp::traits::input_parameter< int >::type k_max(k_maxSEXP);
    Rcpp::traits::input_parameter< int >::type core_num(core_numSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type position_input(position_inputSEXP);
    Rcpp::traits::input_parameter< int >::type cut_num(cut_numSEXP);
    rcpp_result_gen = Rcpp::wrap(mmle_newton_fun(data_use, cond_Z_input, S_depth, res_init_mu, res_init_sigma, t_root_vec, omega_root_vec, k_max, core_num, position_input, cut_num));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_EMMMLE_log_cul_mat_fun2", (DL_FUNC) &_EMMMLE_log_cul_mat_fun2, 2},
    {"_EMMMLE_likelihood_weight", (DL_FUNC) &_EMMMLE_likelihood_weight, 7},
    {"_EMMMLE_log_likelihood_fun", (DL_FUNC) &_EMMMLE_log_likelihood_fun, 7},
    {"_EMMMLE_compute_weight_Z", (DL_FUNC) &_EMMMLE_compute_weight_Z, 7},
    {"_EMMMLE_EMMMLE_weighted_step1", (DL_FUNC) &_EMMMLE_EMMMLE_weighted_step1, 10},
    {"_EMMMLE_EMMMLE_weighted_step2", (DL_FUNC) &_EMMMLE_EMMMLE_weighted_step2, 11},
    {"_EMMMLE_EMMMLE_sigma12_fun", (DL_FUNC) &_EMMMLE_EMMMLE_sigma12_fun, 8},
    {"_EMMMLE_EMMMLE_newton_fun", (DL_FUNC) &_EMMMLE_EMMMLE_newton_fun, 10},
    {"_EMMMLE_EMMMLE_main", (DL_FUNC) &_EMMMLE_EMMMLE_main, 12},
    {"_EMMMLE_EMMMLE", (DL_FUNC) &_EMMMLE_EMMMLE, 14},
    {"_EMMMLE_Update_log_cond_prob_X", (DL_FUNC) &_EMMMLE_Update_log_cond_prob_X, 7},
    {"_EMMMLE_soft", (DL_FUNC) &_EMMMLE_soft, 3},
    {"_EMMMLE_equal1", (DL_FUNC) &_EMMMLE_equal1, 9},
    {"_EMMMLE_equal2", (DL_FUNC) &_EMMMLE_equal2, 9},
    {"_EMMMLE_log_cul_mat_fun1", (DL_FUNC) &_EMMMLE_log_cul_mat_fun1, 1},
    {"_EMMMLE_mmle_weighted_step1", (DL_FUNC) &_EMMMLE_mmle_weighted_step1, 10},
    {"_EMMMLE_mmle_weighted_step2", (DL_FUNC) &_EMMMLE_mmle_weighted_step2, 12},
    {"_EMMMLE_Update_sigma12_fun", (DL_FUNC) &_EMMMLE_Update_sigma12_fun, 9},
    {"_EMMMLE_basic_mmle", (DL_FUNC) &_EMMMLE_basic_mmle, 7},
    {"_EMMMLE_mmle_newton_fun", (DL_FUNC) &_EMMMLE_mmle_newton_fun, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_EMMMLE(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
