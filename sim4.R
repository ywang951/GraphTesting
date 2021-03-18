source("./functions.R")
#' @param n: number of vertices, selected from n_list
#' @param eps: epsilon, selected from eps_list
#' @param alpha: significance level
#' @param rep: number of replications
#' @param rho: pararmeter that controls the sparsity
#' @param d: low rank used in USVT
#' @return power
rep <- 100
d <- 3
rho <- 3
alpha <- 0.5
set.seed(110)
for (n in c(50, 100, 200, 500, 1000)){
  s <- rho * log(n) / n
  prob0 <- GenerateProb(n, eps = 0, s)
  for (eps in c(0.02, 0.1, 0.2, 0.5, 1)){
    rho0.hat <- c()
    for (k in 1:rep){
      adj0 <- GenerateAdj(prob0$P, prob0$Q)
      rho0.hat <- c(rho0.hat, Spearman(adj0$A, adj0$B, num = d))
    }
    qt <- quantile(rho0.hat, alpha)
    prob <- GenerateProb(n, eps, s)
    rho.hat <- c()
    for (k in 1:rep){
      adj <- GenerateAdj(prob$P, prob$Q)
      rho.hat <- c(rho.hat, Spearman(adj$A, adj$B, num = d))
    }
    pow <- sum(rho.hat < qt) / length(rho.hat)
    print(c(n, s * n / log(n), eps, pow))
  }
}
