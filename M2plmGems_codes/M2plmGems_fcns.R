#' -----------------------------------------------------------------------------
#' title:  "GEMS & StepGEMS for latent variable selection in M2PL model"
#' author: "Laixu3107"
#' date:   "2022.10.08"
#' -----------------------------------------------------------------------------

if(sys.nframe() == 0L){rm(list=ls()); gc()}

# ---- requaired packages ------------------------------------------------------
library(mvQuad)     # GH quadrature points and weights
Rcpp::sourceCpp("./M2plmGems_algorithm.cpp")

# ---- some useful function ----------------------------------------------------
list_all_candidate_submod <- function(K){
  # List all candidate sub-models.
  # Laixu3107, Aug 9, 2022
  # K:      number of latent traits.
  # Output: (K+1) * (2^K) mat, each col denotes a sub-model, include intercept.
  
  models <- matrix(data=0, nrow=K, ncol=2^K)
  s <- 0
  for(k in 0:K){
    com_k <- combn(K, k)
    for(l in 1:ncol(com_k)){
      s <- s + 1
      models[com_k[,l],s] <- 1
    }
  }
  return(rbind(1,models))
}


calcu_CR <- function(A_t, A_opt, fixed, col_swap=F){
  # Calculate CR for column swapping case.
  # 2022.07.02, Laixu3017.
  if(!is.null(fixed)){ # with constraints on A
    A_t   <- A_t[,-fixed]
    A_opt <- A_opt[,-fixed]
  }
  J <- ncol(A_t)
  K <- nrow(A_t)
  
  CR    <- 0
  permu <- 1:K
  opt_permu <- 0
  
  if(!col_swap){
    CR <- sum((A_opt!=0)==(A_t!=0))/K/J
  }
  else{ # without constraints on A
    
    all_permu <- gtools::permutations(K,K,1:K)
    
    for(i in 1:nrow(all_permu)){
      
      CR_permu <- sum((A_opt[all_permu[i,],]!=0)==(A_t!=0))/K/J
      
      if(CR_permu > CR){
        CR <- CR_permu
        opt_permu <- i
      }
      
    }
    
    permu <- all_permu[opt_permu,]
  }
  
  return(list(CR=CR, permu=permu))
}


# ---- GEMS for M2pl model -----------------------------------------------------
M2pl_GEMS <- function(
    y,                # N*J mat, observed item responses.
    A_init,           # K*J mat, initial value of A.
    b_init,           # J*1 vec, initial value of b.
    Sigma_init,       # K*K mat, initial value of Sigma.
    beta_list = NULL, # J  list, initial value of beta for all j and s.
    fixed,            # fixed item for identification.
    n_nodes = 5,      # number of quadrature points per dimension
    gamma = 0,        # param for EBIC, 0 <= gamma <= 1. if 0, BIC.  
    
    is.Sigmaknown = 0,
    Minimal.Test  = 1,
    
    MaxIter.GEMS   = 50,
    MaxIter.Newton = 1,
    Tol.para   = 1e-3,
    Tol.qfcn   = 1e-4,
    Tol.newton = 1e-4
    
){
  
  J <- ncol(A_init)  # number of items
  K <- nrow(A_init)  # number of latent traits
  
  # ---- GH nodes and weights ----
  GH <- createNIGrid(dim=K, type="GHN", level=n_nodes)
  x  <- GH$nodes
  wx <- GH$weights
  
  # ---- initial model & beta ----
  beta_init <- matrix(0, K+1, J)
  beta_init[ 1,] <- b_init
  beta_init[-1,] <- A_init
  
  mod_init <- matrix(0, K+1, J)
  mod_init[ 1,] <- 1
  mod_init[-1,] <- (A_init!=0)*1
  
  # ---- all sub-models ----
  sub_mods <- list_all_candidate_submod(K)
  
  sub_mods_list <- vector(mode="list", length=J)
  for(j in 1:J){
    if(j %in% fixed){
      sub_mods_list[[j]] <- mod_init[,j,drop=F]
    }
    else{
      sub_mods_list[[j]] <- sub_mods
    }
  }
  
  # ---- beta list ----
  if(is.null(beta_list)){
    beta_list <- vector(mode="list", length=J)
    for(j in 1:J){
      if(j %in% fixed){
        beta_list[[j]] <- beta_init[,j,drop=F]
      }
      else{
        beta_list[[j]] <- sub_mods*matrix(beta_init[,j], K+1, ncol(sub_mods))
      }
    }
  }
  
  # ---- Q-function list ----
  qfcn_list <- vector(mode="list", length=J)
  for(j in 1:J){
    if(j %in% fixed){
      qfcn_list[[j]] <- 0
    }
    else{
      qfcn_list[[j]] <- rep(0, ncol(sub_mods))
    }
  }
  
  # ---- call m2pl_gems ----
  time_gems <- proc.time()
  output <- m2pl_gems(y  = y,
                      x  = x,
                      wx = wx,
                      beta  = beta_init,
                      sigma = Sigma_init,
                      mod   = mod_init,
                      sub_mods_list = sub_mods_list,
                      beta_list     = beta_list,
                      qfcn_list     = qfcn_list,
                      gamma = gamma,
                      
                      is_sigmaknown = is.Sigmaknown,
                      minimal_test  = Minimal.Test,
                      
                      maxiter_gems   = MaxIter.GEMS,
                      maxiter_newton = MaxIter.Newton,
                      tol_newton     = Tol.newton,
                      tol_param      = Tol.para,
                      tol_qfcn       = Tol.qfcn
  )
  time_gems <- proc.time() - time_gems
  time_gems <- as.numeric(time_gems[3])
  
  return(list(A_opt     = output$beta_opt[-1,],
              b_opt     = output$beta_opt[1 ,],
              Sigma_opt = output$sigma_opt,
              qfcn_seq  = output$qfcn_seq,
              iter_gems = output$iter_gems,
              time_gems = time_gems,
              all_cpu_time = summary(cpu.time)
  ))
  
}

# ---- Step GEMS for M2pl model ------------------------------------------------
M2pl_StepGEMS <- function(
    y,                # N*J mat, observed item responses.
    A_init,           # K*J mat, initial value of A.
    b_init,           # J*1 vec, initial value of b.
    Sigma_init,       # K*K mat, initial value of Sigma.
    beta_list = NULL, # J  list, initial value of beta for all j and s.
    fixed,            # fixed item for identification.
    n_nodes = 5,      # number of quadrature points per dimension
    gamma = 0,        # param for EBIC, 0 <= gamma <= 1. if 0, BIC.  
    
    is.Sigmaknown = 0,
    Minimal.Test  = 1,
    
    MaxIter.GEMS   = 50,
    MaxIter.Newton = 1,
    Tol.para   = 1e-3,
    Tol.qfcn   = 1e-4,
    Tol.newton = 1e-4
    
){
  
  J <- ncol(A_init)  # number of items
  K <- nrow(A_init)  # number of latent traits
  
  # ---- GH nodes and weights ----
  GH <- createNIGrid(dim=K, type="GHN", level=n_nodes)
  x  <- GH$nodes
  wx <- GH$weights
  
  # ---- initial model & beta ----
  beta_init <- matrix(0, K+1, J)
  beta_init[ 1,] <- b_init
  beta_init[-1,] <- A_init
  
  mod_init <- matrix(0, K+1, J)
  mod_init[ 1,] <- 1
  mod_init[-1,] <- (A_init!=0)*1
  
  # ---- all sub-models ----
  sub_mods <- list_all_candidate_submod(K)
  
  sub_mods_list <- vector(mode="list", length=J)
  for(j in 1:J){
    if(j %in% fixed){
      sub_mods_list[[j]] <- mod_init[,j,drop=F]
    }
    else{
      sub_mods_list[[j]] <- sub_mods
    }
  }
  
  # ---- beta list ----
  if(is.null(beta_list)){
    beta_list <- vector(mode="list", length=J)
    for(j in 1:J){
      if(j %in% fixed){
        beta_list[[j]] <- beta_init[,j,drop=F]
      }
      else{
        beta_list[[j]] <- sub_mods*matrix(beta_init[,j], K+1, ncol(sub_mods))
      }
    }
  }
  
  # ---- Q-function list ----
  qfcn_list <- vector(mode="list", length=J)
  for(j in 1:J){
    if(j %in% fixed){
      qfcn_list[[j]] <- 0
    }
    else{
      qfcn_list[[j]] <- rep(0, ncol(sub_mods))
    }
  }
  
  # ---- call mepl_gems ----
  time_gems <- proc.time()
  output <- m2pl_gems_stepwise(
                      y  = y,
                      x  = x,
                      wx = wx,
                      beta  = beta_init,
                      sigma = Sigma_init,
                      mod   = mod_init,
                      sub_mods_list = sub_mods_list,
                      beta_list     = beta_list,
                      qfcn_list     = qfcn_list,
                      gamma = gamma,
                      
                      is_sigmaknown = is.Sigmaknown,
                      minimal_test  = Minimal.Test,
                      
                      maxiter_gems   = MaxIter.GEMS,
                      maxiter_newton = MaxIter.Newton,
                      tol_newton     = Tol.newton,
                      tol_param      = Tol.para,
                      tol_qfcn       = Tol.qfcn
  )
  time_gems <- proc.time() - time_gems
  time_gems <- as.numeric(time_gems[3])
  
  return(list(A_opt     = output$beta_opt[-1,],
              b_opt     = output$beta_opt[1 ,],
              Sigma_opt = output$sigma_opt,
              qfcn_seq  = output$qfcn_seq,
              iter_gems = output$iter_gems,
              used_mat  = output$used_mat,
              time_gems = time_gems,
              all_cpu_time = summary(cpu.time)
  ))
  
}
