#include <Rcpp.h>

RcppExport SEXP model(SEXP T, SEXP I, SEXP K, SEXP M, SEXP ttt, SEXP SS, SEXP alpha_u, SEXP alpha_s, SEXP mu_u, SEXP mu_s, SEXP alpha, 
                      SEXP beta, SEXP gamma, SEXP n_s, SEXP n_u, SEXP varp_u, SEXP lambda_u,SEXP indi, SEXP d, SEXP ybar_s, SEXP ybar_u,
                      SEXP ys2_s, SEXP ys2_u, SEXP a, SEXP b, SEXP lambda, SEXP mk, SEXP Istar, SEXP mKstar, SEXP pp, SEXP pb1,
                      SEXP pb2, SEXP lambda_s,SEXP var_1, SEXP var_2, SEXP p_var, SEXP p_vars, SEXP var_1s, SEXP var_2s, SEXP m_s,SEXP Sigma_s,
                      SEXP p_varu, SEXP var_1u,SEXP var_2u, SEXP m_u, SEXP Sigma_u, SEXP p_vara, SEXP var_1a, SEXP var_2a, SEXP sig_alpha1, SEXP alpha1,
                      SEXP p_varb, SEXP var_1b, SEXP var_2b, SEXP sig_beta1, SEXP beta1, SEXP A_alphau,  SEXP A_alphas, SEXP A_gm, SEXP A_mus, SEXP A_muu, 
                      SEXP A_alpha, SEXP A_beta, SEXP Tune, SEXP pgamma);