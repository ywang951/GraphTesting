source("./functions.R")
#' @param n: number of vertices, selected from n_list
#' @param eps: epsilon, selected from eps_list
#' @param alpha: significance level
#' @param rep: number of replications
#' @return power: a matrix recording power in 3rd column
#' @return pval_dist: a matrix that each column records a list of p-values
n_list <- c(50, 150, 200)
eps_list <- seq(0, 0.1, 0.01)
n_len <- length(n_list)
eps_len <- length(eps_list)
alpha <- 0.05
rep <- 100
power <- matrix(NA, nrow = n_len * eps_len, ncol = 3)
pval_dist <- matrix(NA, nrow = rep + 1, ncol = n_len * eps_len)
idx <- 0
set.seed(110)
for (n in n_list) {
  for (eps in eps_list) {
    idx <- idx + 1
    p_val_list <- c()
    for (N in 1:rep) {
      prob <- GenerateProb(n, eps)
      Alist <- list()
      Blist <- list()
      for (j in 1:20) {
        adj <- GenerateAdj(prob$P, prob$Q)
        Alist[[j]] <- adj$A
        Blist[[j]] <- adj$B
      }
      p_val <- PermutationTest(Alist, Blist, d = 3)
      p_val_list <- c(p_val_list, p_val)
    }
    pow <- sum(p_val_list <= alpha) / length(p_val_list)
    power[idx, ] <- c(n, eps, pow)
    pval_dist[1, idx] <- paste0("n", n, "_eps", eps)
    pval_dist[2:101, idx] <- p_val_list
  }
}
