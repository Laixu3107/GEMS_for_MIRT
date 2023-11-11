#' -----------------------------------------------------------------------------
#' title:  "Examples for GEMS & StepGEMS for M2PL models"
#' author: "Laixu3107"
#' date:   "2023.11.11"
#' -----------------------------------------------------------------------------


# ---- install the required packages ----
if(!require("magrittr", quietly = T)){
  install.packages("magrittr")
}

if(!require("MASS", quietly = T)){
  install.packages("MASS")
}

if(!require("mvQuad", quietly = T)){
  install.packages("mvQuad")
}

if(!require("Rcpp", quietly = T)){
  install.packages("Rcpp")
}

if(!require("RcppArmadillo", quietly = T)){
  install.packages("RcppArmadillo")
}

if(!require("RcppClock", quietly = T)){
  install.packages("RcppClock")
}


rm(list = ls())
# ---- required packages & source code ----
source("./M2plmGems_fcns.R")
library(magrittr)
library(MASS)

# ---- true model & parameters ----
A_t <- c(2.0,1.5,1.0,0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
         0.0,0.0,0.0,0.0,2.0,1.5,1.0,0.5,0.0,0.0,0.0,0.0,
         0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,2.0,1.5,1.0,0.5)
A_t <- matrix(A_t, nrow=3, ncol=12, byrow=T) # transpose of the true loading matrix
b_t <- rnorm(12,0,1)
Sigma_t <- matrix(.1,3,3); diag(Sigma_t) <- 1

# ---- generate random observed item responses ----
set.seed(2)
N <- 2000  # number of subjects
J <- ncol(A_t)
K <- nrow(A_t)
theta <- mvrnorm(n=N, mu=rep(0,K), Sigma=Sigma_t) # latent traits, here just only used to generate item responses
y <- theta %>%
  `%*%` (A_t) %>%
  `+` (matrix(data=b_t,nrow=N,ncol=J,byrow=T)) %>%
  plogis(q=.) %>%
  rbinom(n=N*J, size=1, prob=.) %>%
  matrix(data=., nrow=N, ncol=J, byrow=F)

# ---- GEMS initial parameter ----
fixed  <- c(1,5,9)  # for identification
A_init <- matrix(data=1/J, nrow=K, ncol=J, byrow=TRUE)
A_init[,fixed] <- diag(1,K) # note: the zero elements in the fixed items must be set to 0. 
b_init <- rep(0,J)
Sigma_init <- diag(1,K)

# ---- GEMS settings ----
n_nodes <-  5
is.Sigmaknown  <- 0
MaxIter.GEMS   <- 50
MaxIter.Newton <- 1
Minimal.Test   <- 1
Tol.para   <- 1e-3
Tol.qfcn   <- 1e-4
Tol.newton <- 1e-4

# ---- GEMS ----
output <- M2pl_GEMS(y = y,                  # N*J mat, observed item responses.
                    A_init = A_init,        # K*J mat, initial value of A.
                    b_init = b_init,        # K*1 vec, initial value of b.
                    Sigma_init = Sigma_init,# K*K mat, initial value of Sigma.
                    beta_list = NULL,
                    fixed = fixed,          # fixed item for identifiability.
                    n_nodes = n_nodes,      # number of quadrature points per dimension
                    gamma = 0,              # param for EBIC, 0 <= gamma <= 1. if 0, BIC.  
                    is.Sigmaknown = is.Sigmaknown,
                    Minimal.Test  = Minimal.Test,
                    MaxIter.GEMS   = MaxIter.GEMS,
                    MaxIter.Newton = MaxIter.Newton,
                    Tol.para   = Tol.para,
                    Tol.qfcn   = Tol.qfcn,
                    Tol.newton = Tol.newton
)

CR <- calcu_CR(output$A_opt, A_t, fixed, col_swap=F)$CR

cat("A_opt:\n");  print(output$A_opt);     cat("\n")
cat("b_opt:\n");  print(output$b_opt);     cat("\n")
cat("S_opt:\n");  print(output$Sigma_opt); cat("\n")

cat("bic_seq:\n");  print(output$qfcn_seq); cat("\n")

cat("iter:",        output$iter_gems,    "\n")
cat("time.total:",  output$time_gems,    "\n")
cat("CR:", CR, "\n")

print(output$all_cpu_time)


# ---- StepGEMS ----
output <- M2pl_StepGEMS(
  y = y,
  A_init = A_init,
  b_init = b_init,
  Sigma_init = Sigma_init,
  beta_list = NULL,
  fixed = fixed,     
  n_nodes = n_nodes,
  gamma = 0,
  
  is.Sigmaknown = is.Sigmaknown,
  Minimal.Test = Minimal.Test,
  
  MaxIter.GEMS   = MaxIter.GEMS,
  MaxIter.Newton = MaxIter.Newton,
  Tol.para   = Tol.para,
  Tol.qfcn   = Tol.qfcn,
  Tol.newton = Tol.newton
)

CR <- calcu_CR(output$A_opt, A_t, fixed, col_swap=F)$CR

cat("A_opt:\n");  print(output$A_opt);     cat("\n")
cat("b_opt:\n");  print(output$b_opt);     cat("\n")
cat("S_opt:\n");  print(output$Sigma_opt); cat("\n")

cat("bic_seq:\n");  print(output$qfcn_seq); cat("\n")

cat("iter:",        output$iter_gems,    "\n")
cat("time.total:",  output$time_gems,    "\n")
cat("CR:", CR, "\n")

print(output$all_cpu_time)
