#define RCPP_ARMADILLO_RETURN_ANYVEC_AS_VECTOR
#include <RcppArmadillo.h>
#include <RcppClock.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppClock)]]
using namespace Rcpp;


void m2pl_estep(const arma::mat &y,    // item responses.
                const arma::mat &x,    // quadrature nodes with a column 1.
                const arma::vec &wx,   // weights or prior of x.
                const arma::mat &beta, // item parameters including A and b.
                arma::vec &w0,         // f_{g}^{t}.
                arma::mat &w           // r_{gj}^{t} for all j.
){
  
  int N = y.n_rows; // number of subjects
  int J = y.n_cols; // number of items
  int G = x.n_rows; // number of quadrature nodes
  
  arma::mat log_Fj_xg(G,J);
  arma::mat log_Hj_xg(G,J);
  arma::vec tp(G);
  
  log_Fj_xg = x * beta;
  log_Fj_xg = 1 - 1/(1+exp(log_Fj_xg));
  
  log_Hj_xg = 1 - log_Fj_xg;
  log_Fj_xg = log(log_Fj_xg);
  log_Hj_xg = log(log_Hj_xg);
  
  w0.zeros();
  w.zeros();
  
  int i, j;
  for(i=0;i<N;i++){
    
    tp.zeros();
    
    for(j=0;j<J;j++){
      if(y(i,j)==1){
        tp += log_Fj_xg.col(j);
      }
      else{
        tp += log_Hj_xg.col(j);
      }
    }
    
    tp = arma::exp(tp) % wx;
    tp /= sum(tp);
    
    w0 += tp;
    
    for(j=0;j<J;j++){
      if(y(i,j)==1){
        w.col(j) += tp;
      }
    }
  } 
  
}


arma::mat calcu_sigma_hat(const arma::mat &x, const arma::vec &w0, const int N){
  arma::mat sigma = x.t() * (x.each_col() % w0) / N;
  return(sigma);
}


double obj_func_cpp(arma::mat sigma, arma::mat sigma_hat){
  arma::mat sigma_inv = arma::inv(sigma);
  return arma::accu( sigma_inv % sigma_hat ) + log(arma::det(sigma));
}


arma::mat calcu_sigma_cmle_cpp(arma::mat sigma_hat, arma::mat sigma0, double tol){
  
  arma::mat sigma1 = sigma0;
  arma::mat tmp = sigma0;
  double eps = 1;
  double step = 1;
  while(eps > tol){
    step = 1;
    tmp = arma::inv(sigma0);
    sigma1 = sigma0 - step * ( - tmp * sigma_hat * tmp + tmp );
    sigma1.diag().ones();
    sigma1 = arma::symmatu(sigma1);
    while(obj_func_cpp(sigma0, sigma_hat) < obj_func_cpp(sigma1, sigma_hat) ||
          min(arma::eig_sym(sigma1)) < 0){
      step *= 0.5;
      sigma1 = sigma0 - step * ( - tmp * sigma_hat * tmp + tmp );
      sigma1.diag().ones();
      sigma1 = arma::symmatu(sigma1);
    }
    eps = obj_func_cpp(sigma0, sigma_hat) - obj_func_cpp(sigma1, sigma_hat);
    sigma0 = sigma1;
  }
  return sigma0;
}


arma::vec update_beta(arma::vec beta,     // item parameter vector including bj and aj.
                      const arma::mat &x, // quadrature nodes with a column 1.
                      const arma::vec &w0,
                      const arma::vec &wj,
                      const int maxit,
                      const double tol = 1e-4
){
  
  int G = x.n_rows;
  int K = x.n_cols;
  
  arma::vec Fj(G);
  arma::vec u(G);
  arma::vec d1(K);
  arma::mat d2(K,K);
  
  int it = 0;
  double eps = 1.0;
  
  while(it < maxit && eps > tol){
    
    it += 1;
    Fj = 1 - 1/(1 + exp(x*beta));
    d1 = x.t() * (wj - w0 % Fj);
    
    eps = sqrt(sum(d1 % d1));
    
    u  = Fj % (1 - Fj) % w0;
    d2 = - x.t() * (x.each_col() % u);
    
    beta = beta - inv(d2)*d1;
    
  }
  
  return(beta);
}


double calcu_exp_ebic(const arma::vec &beta, // item parameter vector including bj and aj.
                      const arma::mat &x,    // quadrature nodes with a column 1.
                      const arma::vec &w0,   // f_{g}^{t}.
                      const arma::vec &wj,   // r_{gj}^{t}
                      const int df,          // number of non-zero elements in beta.
                      const int N,           // number of observations.
                      const int K,           // number of latent traits.
                      const double gamma     // parameter for ebic
){
  double obj_fcn = -2*sum((x*beta)%wj - log(1+exp(x*beta))%w0) + df*log(N) + 2*gamma*df*log(K);
  return(obj_fcn);
}


void submod_selection(arma::Col<int> &mod, arma::vec &beta, double &qfcn, // output
                      arma::mat &param_mat,       // output
                      arma::Mat<int> sub_mod_mat, // sub-model space, each col a sub-model.
                      const arma::mat xs,
                      const arma::vec w0,
                      const arma::vec wj,
                      const int N,
                      const double gamma,
                      const int maxiter_newton
){
  
  int S = sub_mod_mat.n_cols;
  int K = sub_mod_mat.n_rows - 1;
  
  arma::vec exbic_vec(S);
  arma::vec param_s(K+1);
  
  
  int s, df;
  for(s=0; s<S; s++){
    
    arma::uvec ind = find(sub_mod_mat.col(s) != 0);
    df = ind.n_elem;
    
    param_s = param_mat.col(s);
    param_s(ind) = update_beta(param_s(ind), xs.cols(ind), w0, wj, maxiter_newton);
    param_mat.col(s) = param_s;
    exbic_vec(s) = calcu_exp_ebic(param_s, xs, w0, wj, df, N, K, gamma);
    
  }
  
  int opt = arma::index_min(exbic_vec);
  mod  = sub_mod_mat.col(opt);
  beta = param_mat.col(opt);
  qfcn = exbic_vec(opt);
  
}


int submod_select_stepwise( arma::Col<int> &mod, arma::vec &beta, double &qfcn, // output
                            arma::mat &param_mat, // output
                            arma::Col<int> sub_mod_c,   // intial sub-model
                            arma::Mat<int> sub_mod_mat, // sub-model space
                            const arma::mat xs,
                            const arma::vec w0,
                            const arma::vec wj,
                            const int N,
                            const double gamma,
                            const int maxiter_newton,
                            unsigned int max_search = 0

){
  
  if(sub_mod_mat.n_cols == 0){
    stop("The sub-model space is empty set.");
  }
  if(max_search == 0){
    max_search = sub_mod_mat.n_cols;
  }
  
  // ---- added for Gems (init val) Oct 06, 2022 ----
  int k, K = sub_mod_mat.n_rows - 1;
  unsigned int s;

  arma::Col<int> pow2 = arma::zeros<arma::Col<int>>(K+1);
  arma::Col<int> sub_mod_no(pow(2,K));
  for(k=1; k<K+1; k++){
    pow2(k) = pow(2,k-1);
  }
  sub_mod_no.fill(-1);
  for(s=0; s<sub_mod_mat.n_cols; s++){
    sub_mod_no(dot(sub_mod_mat.col(s),pow2)) = s;
  }
  // -------------------------------------------------

  int in__ = 0;
  for(s=0; s<sub_mod_mat.n_cols; s++){
    if((sub_mod_c - sub_mod_mat.col(s)).is_zero()){
      in__ = 1;
      sub_mod_mat.shed_col(s);
      break;
    }
  }

  if(in__ == 0){
    stop("The initial sub-model is not in sub-model space.");
  }


  int df;
  arma::vec param_c(K+1);
  double    exbic_c;

  // -- param estimation & expected bic for current sub-model --
  arma::uvec ind = find(sub_mod_c != 0);
  df = ind.n_elem;

  param_c = param_mat.col(sub_mod_no(dot(sub_mod_c,pow2)));
  param_c(ind) = update_beta(param_c(ind), xs.cols(ind), w0, wj, maxiter_newton);
  exbic_c = calcu_exp_ebic(param_c, xs, w0, wj, df, N, K, gamma);
  param_mat.col(sub_mod_no(dot(sub_mod_c,pow2))) = param_c;

  // ---- step-wise ----
  unsigned int n_used = 1;
  unsigned int i;
  unsigned int stop__, add__, del__, min__;
  double exbic_min;
  arma::Mat<int> sub_mod_add(K+1,K);
  arma::mat      param_add(K+1,K);
  arma::vec      exbic_add(K+1);
  arma::Mat<int> sub_mod_del(K+1,K);
  arma::mat      param_del(K+1,K);
  arma::vec      exbic_del(K+1);
  arma::Col<int> sub_mod_s(K+1);
  arma::vec      param_s(K+1);

  while(n_used < max_search && sub_mod_mat.n_cols > 0){

    stop__ = 0;

    // --- forward stage ---
    if(sum(sub_mod_c) == K+1){
      stop__++;  // if full model, can not add any variable
    }
    else{
      // -- find all candidate sub-models (in sub_mod_add) --
      arma::uvec ind = find(sub_mod_c == 0);
      sub_mod_add.zeros();
      add__ = 0;

      for(i=0; i<ind.n_elem; i++){
        sub_mod_s = sub_mod_c;
        sub_mod_s(ind(i)) = 1;

        if(sub_mod_mat.n_cols > 0){
          for(s=0; s<sub_mod_mat.n_cols; s++){
            if((sub_mod_s - sub_mod_mat.col(s)).is_zero()){
              sub_mod_add.col(add__++) = sub_mod_s;
              sub_mod_mat.shed_col(s);
              break;
            }
          }
        }
      }

      if(add__ == 0){
        stop__++;  // if sub_mod_add is empty, can not add any variable
      }
      else{

        // -- compute expected bic for all sub-model in sub_mod_add.
        for(s=0; s<add__; s++){
          sub_mod_s = sub_mod_add.col(s);
          arma::uvec ind = find(sub_mod_s != 0);
          df = ind.n_elem;

          param_s = param_mat.col(sub_mod_no(dot(sub_mod_s,pow2)));
          param_s(ind)  = update_beta(param_s(ind), xs.cols(ind), w0, wj, maxiter_newton);
          param_add.col(s) = param_s;
          exbic_add(s) = calcu_exp_ebic(param_s, xs, w0, wj, df, N, K, gamma);
          param_mat.col(sub_mod_no(dot(sub_mod_s,pow2))) = param_s;
        }

        min__ = arma::index_min(exbic_add.head(add__));
        exbic_min = exbic_add(min__);

        if(exbic_min < exbic_c){
          sub_mod_c = sub_mod_add.col(min__);
          param_c   = param_add.col(min__);
          exbic_c   = exbic_add(min__);
        }
        else{
          stop__++;
        }

        n_used += add__;
      }
    } // end for forward stage

    // --- backward stage ---
    if(sum(sub_mod_c) == 1){ // note: contain an intercept.
      stop__++;  // if null model, can not remove any variable.
    }
    else{
      // -- find all candidate sub-models (in sub_mod_del) --
      arma::uvec ind = find(sub_mod_c.tail(K) == 1) + 1;
      sub_mod_del.zeros();
      del__ = 0;

      for(i=0; i<ind.n_elem; i++){
        sub_mod_s = sub_mod_c;
        sub_mod_s(ind(i)) = 0;

        if(sub_mod_mat.n_cols > 0){
          for(s=0; s<sub_mod_mat.n_cols; s++){
            if((sub_mod_s - sub_mod_mat.col(s)).is_zero()){
              sub_mod_del.col(del__++) = sub_mod_s;
              sub_mod_mat.shed_col(s);
              break;
            }
          }
        }
      }

      if(del__ == 0){
        stop__++;  // if sub_mod_del is empty, can not remove any variable.
      }
      else{
        // -- compute expected bic for all sub-model in sub_mod_del.
        for(s=0; s<del__; s++){
          sub_mod_s = sub_mod_del.col(s);
          arma::uvec ind = find(sub_mod_s != 0);
          df = ind.n_elem;

          param_s = param_mat.col(sub_mod_no(dot(sub_mod_s,pow2)));
          param_s(ind)  = update_beta(param_s(ind), xs.cols(ind), w0, wj, maxiter_newton);
          param_del.col(s) = param_s;
          exbic_del(s) = calcu_exp_ebic(param_s, xs, w0, wj, df, N, K, gamma);
          param_mat.col(sub_mod_no(dot(sub_mod_s,pow2))) = param_s;
        }

        min__ = arma::index_min(exbic_del.head(del__));
        exbic_min = exbic_del(min__);

        if(exbic_min < exbic_c){
          sub_mod_c = sub_mod_del.col(min__);
          param_c   = param_del.col(min__);
          exbic_c   = exbic_del(min__);
        }
        else{
          stop__++;
        }

        n_used += del__;
      }
    } // end for backward stage

    if(stop__ == 2){
      break;
    }

  } // end while

  mod  = sub_mod_c;
  beta = param_c;
  qfcn = exbic_c;
  // sub_mod_mat has updated in the iteration.
  
  int used = n_used;
  return(used);
}


// [[Rcpp::export]]
Rcpp::List m2pl_gems(
    const arma::mat &y,  // item response with dimension N*K.
    const arma::mat &x,  // quadrature nodes with dimension G*K.
    const arma::mat &wx, // quadrature weights with length G.
    
    arma::mat  beta,     // item parameters with dimension K*J.
    arma::mat  sigma,    // correlation of latent traits with dimension K*K.
    arma::Mat<int> mod,  // all possible sub-models for each item. 
    
    const arma::field<arma::Mat<int>> &sub_mods_list, // candidate sub-models for each item.
    arma::field<arma::mat>       beta_list,           // parameters for each sub-model.
    arma::field<arma::vec>       qfcn_list, // not used
    
    const double gamma = 0.0,  // ebic parameter, 
    
    const int is_sigmaknown = 0, // not used, whether sigma is known.
    const int minimal_test  = 1, // whether test that Q-function achieve the minimum.
    
    const int maxiter_gems   = 50,
    const int maxiter_newton = 1,
    const double tol_newton = 1e-4,
    const double tol_param  = 1e-3,
    const double tol_qfcn   = 1e-4

){
  
  Rcpp::Clock clock;
  
  int N = y.n_rows;
  int J = y.n_cols;
  int G = x.n_rows;
  int K = x.n_cols;
  
  arma::mat      new_beta  = beta;
  arma::mat      new_sigma = sigma;
  arma::Mat<int> new_mod   = mod;
  double         new_qfcn =0.0, qfcn = 0.0; // Q(\psi^{t+1}|\psi^{t}) and Q(\psi^{t}|\psi^{t})
  
  // -- for e-step --
  arma::mat xs = arma::ones(G, K+1);
  arma::vec w0(G);
  arma::mat w(G,J);
  
  // -- for ms-step --
  int j;
  arma::mat sigma_hat(K,K);
  double q0 = 0;
  arma::Col<int> mod_j(K+1);
  arma::vec      beta_j(K+1);
  double         qfcn_j;
  arma::vec      qfcn_vec(J);
  
  arma::Mat<int> used_mat(maxiter_gems,J); // record the number of the computed sub-models for each j in MS-step.
  used_mat.fill(0);
  // int used__j;
  
  // -- for stop criterion --
  double err_beta = 1.0, err_qfcn = 1.0; // err_sigma = 1.0
  
  
  int iter_gems = 0;
  arma::vec qfcn_seq(maxiter_gems);
  
  while(iter_gems < maxiter_gems){
    
    iter_gems += 1;
    Rprintf("it:%03d,  e-step.      \r", iter_gems);
    
    // ---- E-step: update xs, w0 & w ----
    clock.tick("e-step");
    sigma = 0.5 * sigma + 0.5 * sigma.t();
    xs.cols(1,K) = x * chol(sigma);
    m2pl_estep(y, xs, wx, beta, w0, w);
    clock.tock("e-step");
    
    // ---- MS-step: ----
    // ---- 1. update sigma -----
    clock.tick("ms-update-sig");
    sigma_hat = calcu_sigma_hat(xs.cols(1,K), w0, N);
    new_sigma = calcu_sigma_cmle_cpp(sigma_hat, sigma, 1e-4);
    new_sigma = arma::symmatu(new_sigma);
    q0 = N * obj_func_cpp(new_sigma, sigma_hat);
    clock.tock("ms-update-sig");
    
    // ---- 2. update beta for each item ----
    clock.tick("ms-update-beta");
    
    for(j=0;j<J;j++){
      Rprintf("it:%03d, ms-step j:%03d.\r", iter_gems, j+1);
      arma::mat beta_mat = beta_list(j);
      submod_selection(mod_j, beta_j, qfcn_j, beta_mat, sub_mods_list(j), xs, w0, w.col(j), N, gamma, maxiter_newton);
      used_mat(iter_gems-1,j) = beta_mat.n_cols;
      new_mod.col(j)  = mod_j;
      new_beta.col(j) = beta_j;
      qfcn_vec(j)     = qfcn_j;
      beta_list(j)    = beta_mat;
    }
    
    new_qfcn = q0 + sum(qfcn_vec);
    qfcn_seq(iter_gems-1) = new_qfcn;
    clock.tock("ms-update-beta");
    
    // ---- stop criterion ----
    if((new_mod - mod).is_zero()){
      
      clock.tick("calcu_qtt");
      qfcn = 0;
      qfcn += N * obj_func_cpp(sigma, sigma_hat);
      for(j=0;j<J;j++){
        qfcn += calcu_exp_ebic(beta.col(j), xs, w0, w.col(j), sum(mod.col(j)), N, K, gamma);
      }
      clock.tock("calcu_qtt");
      
      err_qfcn  = abs((new_qfcn-qfcn)/qfcn);
      
      if( err_qfcn < tol_qfcn ){
        
        if(minimal_test == 0){
          break;
        }
        else{
          
          Rprintf("it:%03d, minimal test. \r", iter_gems);
          clock.tick("ms-origin");
          
          beta = new_beta;
          mod  = new_mod;
          qfcn = new_qfcn;
          
          for(j=0;j<J;j++){
            arma::mat beta_mat = beta_list(j);
            submod_selection(mod_j, beta_j, qfcn_j, beta_mat, sub_mods_list(j), xs, w0, w.col(j), N, gamma, 50);
            used_mat(iter_gems-1,j) = beta_mat.n_cols;
            new_mod.col(j)  = mod_j;
            new_beta.col(j) = beta_j;
            qfcn_vec(j)     = qfcn_j;
            beta_list(j)    = beta_mat;
          }
          
          new_qfcn = q0 + sum(qfcn_vec);
          qfcn_seq(iter_gems-1) = new_qfcn;
          
          clock.tock("ms-origin");
          
          err_beta = norm(new_beta - beta, "fro") / norm(beta, "fro");
          
          if( (new_mod - mod).is_zero() && err_beta < tol_param ){
            break;
          }
          
        } // end minimal.test
      } 
    } // end stop criterion
    
    beta  = new_beta;
    sigma = new_sigma;
    mod   = new_mod;
    
  } // end while
  
  Rcpp::List beta_Rcpp_list = Rcpp::wrap( Rcpp::RcppArmadillo::FieldImporter< arma::Mat<double> >( beta_list ) );
  Rcpp::List qfcn_Rcpp_list = Rcpp::wrap( Rcpp::RcppArmadillo::FieldImporter< arma::Col<double> >( qfcn_list ) );
  
  clock.stop("cpu.time");
  
  List output = List::create(Rcpp::Named("beta_opt")  = new_beta,
                             Rcpp::Named("sigma_opt") = new_sigma,
                             Rcpp::Named("qfcn")      = qfcn,
                             Rcpp::Named("qfcn_seq")  = qfcn_seq.subvec(0,iter_gems-1),
                             Rcpp::Named("used_mat")  = used_mat.rows(0,iter_gems-1),
                             Rcpp::Named("iter_gems") = iter_gems,
                             
                             // test term
                             Rcpp::Named("xs") = xs,
                             Rcpp::Named("w0") = w0,
                             Rcpp::Named("w")  = w,
                             Rcpp::Named("beta_list") = beta_Rcpp_list,
                             Rcpp::Named("qfcn_list") = qfcn_Rcpp_list,
                             Rcpp::Named("qfcn_vec")  = qfcn_vec,
                             Rcpp::Named("qfcn")  = qfcn,
                             Rcpp::Named("new_mod")  = new_mod
                               
  );
  return output;
}


// [[Rcpp::export]]
Rcpp::List m2pl_gems_stepwise(
  const arma::mat &y,
  const arma::mat &x,
  const arma::mat &wx,
  
  arma::mat      beta,
  arma::mat      sigma,
  arma::Mat<int> mod,
  
  const arma::field<arma::Mat<int>> &sub_mods_list,
        arma::field<arma::mat>       beta_list,
        arma::field<arma::vec>       qfcn_list, // not used
  
  const double gamma,
  
  const int is_sigmaknown = 0,
  const int minimal_test  = 1,
  
  const int maxiter_gems   = 50,
  const int maxiter_newton = 1,
  const double tol_newton  = 1e-4,
  const double tol_param   = 1e-3,
  const double tol_qfcn    = 1e-4
                    
){
  
  Rcpp::Clock clock;
  
  int N = y.n_rows;
  int J = y.n_cols;
  int G = x.n_rows;
  int K = x.n_cols;
  
  arma::mat      new_beta  = beta;
  arma::mat      new_sigma = sigma;
  arma::Mat<int> new_mod   = mod;
  double         new_qfcn =0.0, qfcn = 0.0; // Q(\psi^{t+1}|\psi^{t}) and Q(\psi^{t}|\psi^{t})
  
  // -- for e-step --
  arma::mat xs = arma::ones(G, K+1);
  arma::vec w0(G);
  arma::mat w(G,J);
  
  // -- for ms-step --
  int j;
  arma::mat sigma_hat(K,K);
  double q0 = 0;
  arma::Col<int> mod_j(K+1);
  arma::vec      beta_j(K+1);
  double         qfcn_j;
  arma::vec      qfcn_vec(J);
  
  arma::Mat<int> used_mat(maxiter_gems,J); // record the number of the computed sub-models for each j in MS-step.
  used_mat.fill(0);
  int used__j;
  
  // -- for stop criterion --
  double err_beta  = 1.0, err_qfcn  = 1.0; // err_sigma = 1.0
  
  
  int iter_gems = 0;
  arma::vec qfcn_seq(maxiter_gems);

  while(iter_gems < maxiter_gems){
    
    iter_gems += 1;
    Rprintf("it:%03d,  e-step.      \r", iter_gems);
    
    clock.tick("e-step");
    // ---- E-step: update xs, w0 & w ----
    sigma = 0.5 * sigma + 0.5 * sigma.t();
    xs.cols(1,K) = x * chol(sigma);
    m2pl_estep(y, xs, wx, beta, w0, w);
    clock.tock("e-step");
    
    // ---- MS-step: ----
    // ---- 1. update sigma -----
    clock.tick("ms-update-sig");
    sigma_hat = calcu_sigma_hat(xs.cols(1,K), w0, N);
    new_sigma = calcu_sigma_cmle_cpp(sigma_hat, sigma, 1e-4);
    new_sigma = arma::symmatu(new_sigma);
    q0 = N * obj_func_cpp(new_sigma, sigma_hat);
    clock.tock("ms-update-sig");
    
    // ---- 2. update beta for each item ----
    clock.tick("ms-update-beta");
    
    if(iter_gems == 1){  // if iter_gems == 1, update all submodels
      for(j=0;j<J;j++){
        Rprintf("it:%03d, ms-step j:%03d.\r", iter_gems, j+1);
        arma::mat beta_mat = beta_list(j);
        submod_selection(mod_j, beta_j, qfcn_j, beta_mat, sub_mods_list(j), xs, w0, w.col(j), N, gamma, maxiter_newton);
        used_mat(iter_gems-1,j) = beta_mat.n_cols;
        new_mod.col(j)  = mod_j;
        new_beta.col(j) = beta_j;
        qfcn_vec(j)     = qfcn_j;
        beta_list(j)    = beta_mat;
      }
    }
    else{  // if iter_gems > 1, update by step wise
      for(j=0;j<J;j++){
        Rprintf("it:%03d, ms-step j:%03d.\r", iter_gems, j+1);
        arma::mat beta_mat = beta_list(j);
        used__j = submod_select_stepwise(mod_j, beta_j, qfcn_j, beta_mat, mod.col(j), sub_mods_list(j), xs, w0, w.col(j), N, gamma, maxiter_newton);
        used_mat(iter_gems-1,j) = used__j;
        new_mod.col(j)  = mod_j;
        new_beta.col(j) = beta_j;
        qfcn_vec(j)     = qfcn_j;
        beta_list(j)    = beta_mat;
      }
    }
    
    new_qfcn = q0 + sum(qfcn_vec);
    qfcn_seq(iter_gems-1) = new_qfcn;
    clock.tock("ms-update-beta");
    
    // ---- stop criterion ----
    if((new_mod - mod).is_zero()){
    
      clock.tick("calcu_qtt");
      qfcn = 0;
      qfcn += N * obj_func_cpp(sigma, sigma_hat);
      for(j=0;j<J;j++){
        qfcn += calcu_exp_ebic(beta.col(j), xs, w0, w.col(j), sum(mod.col(j)), N, K, gamma);
      }
      clock.tock("calcu_qtt");
  
      err_qfcn  = abs((new_qfcn-qfcn)/qfcn);
  
      if( err_qfcn < tol_qfcn ){
      
        if(minimal_test == 0){
          break;
        }
        else{
      
          Rprintf("it:%03d, minimal test. \r", iter_gems);
          
          clock.tick("ms-origin");
    
          beta = new_beta;
          mod  = new_mod;
          qfcn = new_qfcn;
      
          for(j=0;j<J;j++){
            arma::mat beta_mat = beta_list(j);
            submod_selection(mod_j, beta_j, qfcn_j, beta_mat, sub_mods_list(j), xs, w0, w.col(j), N, gamma, 50);
            used_mat(iter_gems-1,j) = beta_mat.n_cols;
            new_mod.col(j)  = mod_j;
            new_beta.col(j) = beta_j;
            qfcn_vec(j)     = qfcn_j;
            beta_list(j)    = beta_mat;
          }
      
          new_qfcn = q0 + sum(qfcn_vec);
          qfcn_seq(iter_gems-1) = new_qfcn;
      
          clock.tock("ms-origin");
      
          err_beta = norm(new_beta - beta, "fro") / norm(beta, "fro");

          if( (new_mod - mod).is_zero() && err_beta < tol_param ){
            break;
          }
      
        } // end minimal.test
      } 
    } // end stop criterion
    
    beta  = new_beta;
    sigma = new_sigma;
    mod   = new_mod;
    
  } // end while
  
  Rcpp::List beta_Rcpp_list = Rcpp::wrap( Rcpp::RcppArmadillo::FieldImporter< arma::Mat<double> >( beta_list ) );
  Rcpp::List qfcn_Rcpp_list = Rcpp::wrap( Rcpp::RcppArmadillo::FieldImporter< arma::Col<double> >( qfcn_list ) );
  
  clock.stop("cpu.time");
  
  List output = List::create(Rcpp::Named("beta_opt")  = new_beta,
                             Rcpp::Named("sigma_opt") = new_sigma,
                             Rcpp::Named("qfcn")      = qfcn,
                             Rcpp::Named("qfcn_seq")  = qfcn_seq.subvec(0,iter_gems-1),
                             Rcpp::Named("used_mat")  = used_mat.rows(0,iter_gems-1),
                             Rcpp::Named("iter_gems") = iter_gems,

                             // test term
                             Rcpp::Named("xs") = xs,
                             Rcpp::Named("w0") = w0,
                             Rcpp::Named("w")  = w,
                             Rcpp::Named("beta_list") = beta_Rcpp_list,
                             Rcpp::Named("qfcn_list") = qfcn_Rcpp_list,
                             Rcpp::Named("qfcn_vec")  = qfcn_vec,
                             Rcpp::Named("qfcn")  = qfcn,
                             Rcpp::Named("new_mod")  = new_mod
                             
                             );
  return output;
}


// [[Rcpp::export]]
double calcu_obs_ebic(
  const arma::mat &y,
  const arma::mat &x,
  const arma::mat &wx,
  arma::mat &beta,
  arma::mat &sigma,
  double gamma = 0.0
){
  
  int N = y.n_rows;
  int J = y.n_cols;
  int G = x.n_rows;
  int K = x.n_cols;
  int df = arma::accu(beta!=0) - J;
  
  arma::mat xs = arma::ones(G, K+1);
  xs.cols(1,K) = x * chol(sigma);
  
  arma::mat log_Pj_xg(G,J);
  arma::mat log_Qj_xg(G,J);
  arma::vec tp(G);
  
  log_Pj_xg = xs * beta;
  log_Pj_xg = 1 - 1/(1+exp(log_Pj_xg));
  log_Qj_xg = 1 - log_Pj_xg;
  
  log_Pj_xg = log(log_Pj_xg);
  log_Qj_xg = log(log_Qj_xg);
  
  arma::vec py(G);
  double obs_ebic = 0.0;
  
  int i, j;
  for(i=0;i<N;i++){
    
    py.zeros();
    
    for(j=0;j<J;j++){
      if(y(i,j)==1){
        py += log_Pj_xg.col(j);
      }
      else{
        py += log_Qj_xg.col(j);
      }
    }
    
    py = exp(py) % wx;
    obs_ebic += log(sum(py));
  }
  
  obs_ebic = -2 * obs_ebic + df*log(N) + 2*gamma*df*log(K);
  
  return(obs_ebic);
}
