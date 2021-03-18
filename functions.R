library(irlba)
library(vegan)
library(ecodist)
GenerateProb <- function(n, eps) {
  #' @param n: number of vertices
  #' @param eps: epsilon
  #' @return output: vectors of upper trianguler of true edge-probability matrices
  X <- matrix(rnorm(2 * n), nrow = n, ncol = 2)
  Z <- matrix(rnorm(2 * n, mean = 0, sd = sqrt(eps)), nrow = n, ncol = 2)
  Y <- X + Z
  D1 <- as.matrix(dist(X)^2)
  D2 <- as.matrix(dist(Y)^2)
  P <- exp(-D1^2)
  Q <- exp(-D2^2 / 4)
  output <- list()
  output$P <- P
  output$Q <- Q
  return(output)
}

GenerateAdj <- function(P, Q) {
  #' @param P,Q: edge-probability matrices
  #' @return output: adjcency matrices
  n <- dim(P)[1]
  upper.tri.ind <- upper.tri(P)
  p.upper <- P[upper.tri.ind]
  q.upper <- Q[upper.tri.ind]
  A.upper <- rbinom(n * (n - 1) / 2, 1, p.upper)
  B.upper <- rbinom(n * (n - 1) / 2, 1, q.upper)
  A <- matrix(0, ncol = n, nrow = n)
  B <- matrix(0, ncol = n, nrow = n)
  A[upper.tri.ind] <- A.upper
  B[upper.tri.ind] <- B.upper
  A <- A + t(A)
  B <- B + t(B)
  output <- list()
  output$A <- A
  output$B <- B
  return(output)
}

USVT <- function(M, num, res = TRUE) {
  # required packages: irlba
  #' @param M: (average) adjacency matrix
  #' @param num: low rank used in truncated SVD
  #' @param res: if restrict the USVT estimator within [0, 1]
  #' @return output: estimated edge probability matrix
  rst <- irlba(as.matrix(M), nv = num, nu = num, maxit = 1000)
  n <- dim(M)[1]
  if (num > 1) {
    M.svt <- rst$u %*% diag(rst$d) %*% t(rst$v)
    if (!res) {
      return(M.svt)
    } else {
      M.fin <- matrix(0, n, n)
      upper.tri.ind <- upper.tri(M)
      m.upper <- M.svt[upper.tri.ind]
      m.upper[which(m.upper > 1)] <- 1
      m.upper[which(m.upper < 0)] <- 0
      M.fin[upper.tri.ind] <- m.upper
      M.fin <- M.fin + t(M.fin)
      return(M.fin)
    }
  }
}

Spearman <- function(A, B, num, res = TRUE) {
  #' @param A, B: two adjacency matrix
  #' @param d: low rank used in USVT
  #' @param res: if restrict the USVT estimator within [0, 1]
  #' @return output: spearman rank correlation coefficient
  Phat <- USVT(A, num, res)
  Qhat <- USVT(B, num, res)
  phat.upper <- Phat[upper.tri(Phat)]
  qhat.upper <- Qhat[upper.tri(Qhat)]
  output <- cor(phat.upper, qhat.upper, method = "spearman")
  return(output)
}

Bootstrap <- function(A, B, b_num = 1000, d = 3) {
  #' @param A: (average) adjacency matrix from the first latent position
  #' @param B: (average) adjacency matrix from the second latent position
  #' @param b_num: number of bootstrapping
  #' @param d: low rank in USVT
  #' @return p_val: p-value
  n <- dim(A)[1]
  P.hat <- USVT(A, num = d)
  Q.hat <- USVT(B, num = d)
  rho0 <- Spearman(A, B, num = d)
  rho1 <- c()
  rho2 <- c()
  for (k in 1:b_num) {
    A_bts <- GenerateAdj(P.hat, P.hat)
    B_bts <- GenerateAdj(Q.hat, Q.hat)
    s1 <- Spearman(A_bts$A, A_bts$B, num = d)
    s2 <- Spearman(B_bts$A, B_bts$B, num = d)
    rho1 <- c(rho1, s1)
    rho2 <- c(rho2, s2)
    }
  p_val <- max(sum(rho0 > rho1) / b_num, sum(rho0 > rho2) / b_num)
  return(p_val)
}

PermutationTest <- function(list1, list2, p_num = 1000, d = 3) {
  #' @param list1, list2: lists of adjacency matrices
  #' @param p_num: number of permutations
  #' @param d: low rank in USVT
  #' @return pvalue: p-value
  `%notin%` <- Negate(`%in%`)
  nlist <- length(list1)
  n <- dim(list1[[1]])[1]
  A_bar_1 <- matrix(0, n, n)
  A_bar_2 <- matrix(0, n, n)
  for (i in 1:nlist) {
    A_bar_1 <- A_bar_1 + list1[[i]]
    A_bar_2 <- A_bar_2 + list2[[i]]
  }
  A_bar_1 <- A_bar_1 / nlist
  A_bar_2 <- A_bar_2 / nlist
  rho <- Spearman(A_bar_1, A_bar_2, num = d, res = FALSE)
  rho.hat <- c()
  combine.list <- c(list1, list2)
  for (i in 1:p_num) {
    relabel <- sample(1:(nlist * 2), nlist, replace = FALSE)
    A.perm <- matrix(0, n, n)
    B.perm <- matrix(0, n, n)
    for (j in 1:(nlist * 2)) {
      if (j %in% relabel) {
        A.perm <- A.perm + combine.list[[j]]
      }else{
        B.perm <- B.perm + combine.list[[j]]
      }
    }
    A.perm_bar <- A.perm / nlist
    B.perm_bar <- B.perm / nlist
    rst <- Spearman(A.perm_bar, B.perm_bar, num = d, res = FALSE)
    rho.hat <- c(rho.hat, rst)
  }
  pvalue <- sum(rho.hat < rho) / length(rho.hat)
  return(pvalue)
}

Nonmetric_MDS <- function(A, B) {
  #' @param A, B: average adjacency matrices
  #' @return outout: NMDS test statistic
  # required packages: vegan, ecodist
  A <- as.matrix(A)
  B <- as.matrix(B)
  A_mds <- 1 - A
  B_mds <- 1 - B
  Xhat <- nmds(A_mds[lower.tri(A_mds)], mindim = 2, maxdim = 2)
  Yhat <- nmds(B_mds[lower.tri(B_mds)], mindim = 2, maxdim = 2)
  yconf <- Yhat$conf[[length(Yhat$conf)]]
  xconf <- Xhat$conf[[length(Xhat$conf)]]
  Tx <- procrustes(xconf, yconf)
  Ty <- procrustes(yconf, xconf)
  output <- min(sqrt(Tx$ss), sqrt(Ty$ss))
  return(output)
}
