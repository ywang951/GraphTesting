# path <- "Your path"
# setwd(path)
source("./functions.R")
#' @param n: number of vertices, selected from n_list
#' @param eps: epsilon, selected from eps_list
#' @param alpha: significance level
#' @param rep: number of replications
#' @return power: a matrix recording power in 3rd column
#' @return pval_dist: a matrix that each column records a list of p-values
n_list <- c(50, 100, 150)
eps_list <- c(0, 0.02, 0.04, 0.06, 0.08, 0.1)
n_len <- length(n_list)
eps_len <- length(eps_list)
alpha <- 0.05
rep <- 100
power <- matrix(NA, nrow = n_len * eps_len, ncol = 3)
power[, 1] <- rep(n_list, each = eps_len)
power[, 2] <- rep(eps_list, n_len)
pval_dist <- matrix(NA, nrow = rep + 1, ncol = n_len * eps_len)
idx <- 0
set.seed(110)
for (n in n_list) {
  for (eps in eps_list) {
    print(c(n, eps))
    idx <- idx + 1
    p_val_list <- c()
    for (N in 1:rep) {
      print(N)
      prob <- GenerateProb(n, eps)
      adj <- GenerateAdj(prob$P, prob$Q)
      p_val <- Bootstrap(adj$A, adj$B, b_num = 1000, d = 3)
      p_val_list <- c(p_val_list, p_val)
    }
    pow <- sum(p_val_list <= alpha) / length(p_val_list)
    power[idx, 3] <- pow
    pval_dist[1, idx] <- paste0("n", n, "_eps", eps)
    pval_dist[2:101, idx] <- p_val_list
  }
}
