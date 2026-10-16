# =============================================================================
# Covariance structures for the simulated design matrices (Appendix B.1)
#   1. compound symmetry          generate_compound_sigma()
#   2. AR(1)                      generate_ar_sigma()
#   3. sparse random precision    random_matrix_sigma_sparse()
#   4. scale-free (Barabasi-Albert) precision  generate_ba_hub_sigma()
#   5. block-diagonal precision   generate_block_sigma()
# Sourced by R/data_generation.R.
# =============================================================================

# Compound ----------------------------------------------------------------

# function of Sigma
#Compound symmetry
generate_compound_sigma<-function(p,rho){
  sigma<-matrix(rep(0,p^2),ncol=p)
  for (i in 1:p) {
    for (j in 1:p) {
      if(i==j){
        sigma[i,i]<-1
      }else{
        sigma[i,j]<-rho
      }
    }
  }
  return(sigma)
}


# AR1 ---------------------------------------------------------------------


generate_ar_sigma<-function(p,rho){
  sigma<-matrix(rep(0,p^2),ncol=p)
  for (i in 1:p) {
    for (j in 1:p) {
      if(i==j){
        sigma[i,i]<-1
      }else if(i>j){
        sigma[i,j]<-rho^(i-j)
        sigma[j,i]<-rho^(i-j)
      }
    }
  }
  return(sigma)
}




# random --------------------------------------------------------------------
sp_runif<-function(n,seed=1){
  #Uniform distribution: -0.75~-0.25,0.25~0.75 
  set.seed(seed)
  a<-runif(n,-0.5,0.5)
  a_modified <- ifelse(a < 0, a - 0.25, a + 0.25)
  return(a_modified)
}



random_matrix_sigma_sparse <- function(p, sprate = 0.5, seed = 1) {
  set.seed(seed)
  E <- matrix(sp_runif(p^2, seed), p, p)
  
  non_diag_indices <- which(row(E) != col(E))
  

  num_zeros <- round(length(non_diag_indices) * sprate)
  
  if (num_zeros > 0) {
    zero_indices <- sample(non_diag_indices, num_zeros, replace = FALSE)
    E[zero_indices] <- 0
  }
  
  E_thilde <- (E + t(E)) / 2
  lmin <- min(eigen(E_thilde)$values)
  omega_thilde <- E_thilde + (0.1 - lmin) * diag(p)
  
  L_sq <- matrix(0, nrow = p, ncol = p)
  diag(L_sq) <- sqrt(diag(solve(omega_thilde)))
  
  omega <- L_sq %*% omega_thilde %*% L_sq
  
  
  sigma <- solve(omega)
  return(sigma)
}

# hub ---------------------------------------------------------------------

library(igraph)

generate_ba_hub_sigma <- function(p, m = 2, seed = 1) {
  set.seed(seed)
  
  
  g <- sample_pa(n = p, m = m, directed = FALSE)
  
  
  adj <- as_adjacency_matrix(g, sparse = FALSE)
  E <- matrix(0, p, p)
  
  upper_idx <- which(upper.tri(adj) & adj == 1)
  if(length(upper_idx) > 0) {
    weights <- sp_runif(length(upper_idx), seed)
    E[upper_idx] <- weights
    E <- E + t(E)
  }
  
  
  lmin <- min(eigen(E)$values)
  omega_tilde <- E + (abs(lmin) + 0.1) * diag(p)
  
  
  L_sq <- diag(sqrt(diag(solve(omega_tilde))))
  omega <- L_sq %*% omega_tilde %*% L_sq
  
  return(solve(omega))
}


# block -------------------------------------------------------------------

library(igraph)

generate_block_sigma <- function(p, n_blocks = 5, within_prob = 0.3, seed = 1) {
  set.seed(seed)
  
  
  block_sizes <- as.numeric(table(cut(1:p, breaks = n_blocks, labels = FALSE)))
  
  
  pref_matrix <- diag(within_prob, n_blocks)
  g <- sample_sbm(n = p, pref.matrix = pref_matrix, block.sizes = block_sizes, directed = FALSE)
  
  
  adj <- as_adjacency_matrix(g, sparse = FALSE)
  E <- matrix(0, p, p)
  upper_idx <- which(upper.tri(adj) & adj == 1)
  
  if(length(upper_idx) > 0) {
    weights <- sp_runif(length(upper_idx), seed)
    E[upper_idx] <- weights
    E <- E + t(E)
  }
  
  
  lmin <- min(eigen(E)$values)
  omega_tilde <- E + (abs(lmin) + 0.1) * diag(p)
  
  
  L_sq <- diag(sqrt(diag(solve(omega_tilde))))
  omega <- L_sq %*% omega_tilde %*% L_sq
  
  return(solve(omega))
}
