library(irlba)
library(vegan)
library(ecodist)
library(MCMCpack)
library(pracma)
library(igraph)
library(doBy)
library(gtools)
permutation <- function(network, namesequence) {
  # This function is from the R package "snatm", which is not compatible with my system.

  if (dim(network)[1] != length(namesequence)) {
    stop("Dimension of network and length of namesequence must be equal.")
  }

  if (length(unique(namesequence)) < dim(network)[1]) {
    stop("Each name in namesequence must be different.")
  }

  newnetwork <- network
  while (!all(rownames(newnetwork) == namesequence)) {

    for (i in 1:dim(network)[1]) {

      if (rownames(network)[i] != namesequence[i]) {
        columnposition <- which(rownames(network) == namesequence[i])
        newnetwork[, i] <- network[, columnposition]
        newnetwork[, columnposition] <- network[, i]
        tempcolumn <- newnetwork[i, ]
        newnetwork[i, ] <- newnetwork[columnposition, ]
        newnetwork[columnposition, ] <- tempcolumn
        temprowname <- colnames(newnetwork)[i]
        colnames(newnetwork)[i] <- rownames(network)[columnposition]
        rownames(newnetwork)[i] <- rownames(network)[columnposition]
        colnames(newnetwork)[columnposition] <- temprowname
        rownames(newnetwork)[columnposition] <- temprowname
        network <- newnetwork
      }
    }
  }

  newnetwork
}

Generate.A <- function(P) {
  #' @param P: edge-probability matrices
  #' @return output: adjcency matrices

  n <- dim(P)[1]
  upper.tri.ind <- upper.tri(P)
  p.upper <- P[upper.tri.ind]
  A.upper <- rbinom(n * (n - 1) / 2, 1, p.upper)
  A <- matrix(0, ncol = n, nrow = n)
  A[upper.tri.ind] <- A.upper
  A <- A + t(A)

  return(A)
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

Bootstrap <- function(A, B, b_num = 1000, d, n0) {
  #' @param A: list of adjacency matrices from the first latent position
  #' @param B: list of adjacency matrices from the second latent position
  #' @param b_num: number of bootstrapping
  #' @param d: low rank in USVT
  #' @return p_val: p-value

  m <- length(A)
  n <- dim(A[[1]])[1]
  A.avg <- B.avg <- matrix(0, n, n)

  for (k in 1:m) {
    A.avg <- A.avg + A[[k]]
    B.avg <- B.avg + B[[k]]
  }

  A.avg <- A.avg / m
  B.avg <- B.avg / m
  rho0 <- Spearman(A.avg, B.avg, num = d, res = FALSE)
  rhovec1 <- rhovec2 <- numeric(b_num)

  for (i in 1:b_num) {
    A.perm <- B.perm <- matrix(0, n, n)

    for (k in 1:m) {
      A.temp <- A[[k]]
      B.temp <- B[[k]]
      rownames(A.temp) <- colnames(A.temp) <- seq(n)
      rownames(B.temp) <- colnames(B.temp) <- seq(n)
      a <- seq(n)
      s <- sample(seq(n), n0)
      a[sort(s)] <- s
      A.perm <- A.perm + permutation(A.temp, a)
      B.perm <- B.perm + permutation(B.temp, a)
    }

    rhovec1[i] <- Spearman(A.avg, A.perm / m, num = d, res = FALSE)
    rhovec2[i] <- Spearman(B.avg, B.perm / m, num = d, res = FALSE)
  }

  pval1 <- sum(rho0 > rhovec1) / b_num
  pval2 <- sum(rho0 > rhovec2) / b_num

  return(max(pval1, pval2))
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

Simulation1 <- function(n = 100, m = 1, eps = 0, d = 2, rep = 100, n0 = 8, b_num = 1000) {
  #' @param n: number of vertices
  #' @param m: number of graphs sampled by one population
  #' @param eps: variance of the distribution difference
  #' @param d: low rank in USVT
  #' @param rep: number of Monte Carlo replicates
  #' @param n0: number of permutations in the bootstrap algorithm
  #' @param b_num: number of bootstrap samplings
  #' @return pval_list: a list of p-values

  set.seed(110)
  pval_list <- c()

  for (N in 1:rep) {
    X <- matrix(rnorm(2 * n), nrow = n, ncol = 2)
    D1 <- as.matrix(dist(X)^2)
    Z <- matrix(rnorm(2 * n, sd = sqrt(eps)), nrow = n, ncol = 2)
    D2 <- as.matrix(dist(X + Z)^2)
    P <- exp(-D1^2)
    Q <- exp(-D2^2 / 2)
    A <- B <- list()

    for (k in 1:m) {
      A[[k]] <- Generate.A(P)
      B[[k]] <- Generate.A(Q)
    }

    pval <- Bootstrap(A, B, b_num, d, n0)
    pval_list <- c(pval_list, pval)
  }

  return(pval_list)
}

Simulation2.Bootstrap <- function(eps, n = 100, d = 2, rep = 100, n0 = 8, alpha = 0.05, b_num = 1000) {
  #' @param n: number of vertices
  #' @param d: low rank in USVT
  #' @param rep: number of Monte Carlo replicates
  #' @param n0: number of permutations in the bootstrap algorithm
  #' @param alpha: significance level
  #' @param b_num: number of bootstrap samplings
  #' @param eps: number of graphs sampled by one population
  #' @return power: power of the testing with significance level alpha

  set.seed(110)
  pval_list <- c()

  for (N in 1:rep) {
    X <- matrix(rnorm(2 * n), nrow = n, ncol = 2)
    D1 <- as.matrix(dist(X)^2)
    Z <- matrix(rnorm(2 * n, sd = sqrt(eps)), nrow = n, ncol = 2)
    D2 <- as.matrix(dist(X + Z)^2)
    P <- exp(-D1^2)
    Q <- exp(-D2^2 / 2)
    A[[1]] <- Generate.A(P)
    B[[1]] <- Generate.A(Q)
    pval <- Bootstrap(A, B, b_num, d, n0)
    pval_list <- c(pval_list, pval)
  }
  power <- sum(pval_list <= alpha) / length(pval_list)

  return(power)
}

Simulation2.NMDS <- function(n_list, eps_list, alpha = 0.05, rep = 100) {
  #' @param n_list: a list of n, number of vertices
  #' @param eps_list: a list of epsilon, variance of the distribution difference
  #' @param alpha: significance level
  #' @param rep: number of Monte Carlo replicates
  #' @return power: power of the testing with significance level alpha

  n_len <- length(n_list)
  eps_len <- length(eps_list)

  if (!(0 %in% eps_list)) {
    stop("eps_list must contain 0.")
  }

  power <- matrix(NA, nrow = n_len * (eps_len - 1) + 1, ncol = 3)
  power[1, ] <- c("n", "epsilon", "Non-metric embedding")
  idx <- 1

  for (n in n_list) {
    X <- matrix(rnorm(2 * n), nrow = n, ncol = 2)
    D1 <- as.matrix(dist(X)^2)
    P <- exp(-D1^2)

    for (eps in eps_list) {
      Z <- matrix(rnorm(2 * n, sd = sqrt(eps)), nrow = n, ncol = 2)
      D2 <- as.matrix(dist(X + Z)^2)
      Q <- exp(-D2^2 / 2)
      nmds.hat <- c()
      for (N in 1:rep) {
        A_bar <- matrix(0, ncol = n, nrow = n)
        B_bar <- matrix(0, ncol = n, nrow = n)
        for (b in 1:20) {
          A_bar <- A_bar + Generate.A(P)
          B_bar <- B_bar + Generate.A(Q)
        }
        nmds.hat <- c(nmds.hat, Nonmetric_MDS(A_bar / 20, B_bar / 20))
      }

      if (eps == 0) {
        qt_nmds <- quantile(nmds.hat, 1 - alpha)
        power[idx, ] <- c(n, eps, alpha)
      } else {
        idx <- idx + 1
        pow_nmds <- sum(nmds.hat > qt_nmds) / length(nmds.hat)
        power[idx, ] <- c(n, eps, pow_nmds)
      }

    }
  }

  power <- as.data.frame(power)
  names(power) <- power[1, ]
  power <- power[-1, ]

  return(power)
}

## Additional Simulations in the Supplementary

Simulation3 <- function(n_list, eps_list, rho, alpha = 0.05, rep = 100, d = 3) {
  #' @param n_list: a list of n
  #' @param eps_list: a list of epsilon
  #' @param rho: parameter controling sparsity
  #' @param alpha: significant level
  #' @param rep: number of replications
  #' @param d: low rank used in USVT
  #' @return power

  n_len <- length(n_list)
  eps_len <- length(eps_list)

  if (!(0 %in% eps_list)) {
    stop("eps_list must contain 0.")
  }

  power <- matrix(NA, nrow = n_len * (eps_len - 1) + 1, ncol = 4)
  power[1, ] <- c("n", "epsilon", "sparsity", "power")
  idx <- 1

  for (n in n_list) {
    s <- rho * log(n) / n
    X <- matrix(rnorm(2 * n), nrow = n, ncol = 2)
    D1 <- as.matrix(dist(X)^2)
    P <- s * exp(-D1^2)

    for (eps in eps_list) {
      Z <- matrix(rnorm(2 * n, sd = sqrt(eps)), nrow = n, ncol = 2)
      D2 <- as.matrix(dist(X + Z)^2)
      Q <- s * exp(-D2^2 / 4)
      rho.hat <- c()

      for (N in 1:rep) {
        rho0 <- Spearman(Generate.A(P), Generate.A(Q), num = d)
        rho.hat <- c(rho.hat, rho0)
      }

      if (eps == 0) {
        qt <- quantile(rho.hat, alpha)
      } else {
        idx <- idx + 1
        pow_spm <- sum(rho.hat < qt) / length(rho.hat)
        power[idx, ] <- c(n, eps, s * n / log(n), pow_spm)
      }

    }
  }

  return(power)
}

Simulation4 <- function(nlist, slist, d = 3, tau = 0.7, rep = 100) {
  #' @param nlist: a list of n
  #' @param slist: a list of s
  #' @param d: low rank used in USVT
  #' @param tau: probability in quantile
  #' @param rep: number of replicates

  df <- list()
  df$n <- c()
  df$prop <- c()
  df$mean <- c()
  df$meanrank <- c()
  df$meanrankall <- c()
  df$sd <- c()
  df$sdrank <- c()
  df$sdrankall <- c()

  for (n in nlist) {
    X <- matrix(rnorm(2 * n), nrow = n, ncol = 2)
    Z <- matrix(rnorm(2 * n), nrow = n, ncol = 2)
    Y <- X + Z
    Y_trans <- MCMCpack::procrustes(Y, X)$X.new
    diff_X_Ytrans <- apply(Y_trans - X, 1, Norm)
    true_impt_node <- which(diff_X_Ytrans > quantile(diff_X_Ytrans, 0.7))
    D1 <- as.matrix(dist(X)^2)
    D2 <- as.matrix(dist(Y)^2)
    P <- exp(-D1^2)
    Q <- exp(-D2^2 / 4)

    for (s in slist) {

      if (s <= n) {
        corr_count <- c()
        rank_impt <- c()
        rank_impt_all <- c()

        for (N in 1:rep) {
          P.hat <- USVT(Generate.A(P), num = d)
          Q.hat <- USVT(Generate.A(Q), num = d)
          phat.upper <- P.hat[upper.tri(P.hat)]
          qhat.upper <- Q.hat[upper.tri(Q.hat)]
          phat.rank <- rank(phat.upper, ties.method = "average")
          qhat.rank <- rank(qhat.upper, ties.method = "average")
          diff.rank <- abs(phat.rank - qhat.rank)
          A.star <- matrix(0, ncol = n, nrow = n)
          A.star[upper.tri(A.star)] <- (diff.rank > quantile(diff.rank, tau))
          A.star <- A.star + t(A.star)
          esti_impt_node <- which.maxn(colSums(A.star), n = s)
          corr_count <- c(corr_count, sum(esti_impt_node %in% true_impt_node))
          rank_impt <- c(rank_impt, median(rank(diff_X_Ytrans[true_impt_node])[match(esti_impt_node[esti_impt_node %in% true_impt_node], true_impt_node)] / length(true_impt_node)))
          rank_impt_all <- c(rank_impt_all, median(rank(diff_X_Ytrans)[esti_impt_node] / n))
        }

        df$n <- c(df$n, n)
        df$prop <- c(df$prop, s)
        df$mean <- c(df$mean, round(mean(corr_count), 2))
        df$sd <- c(df$sd, round(sd(corr_count), 4))
        df$meanrank <- c(df$meanrank, round(mean(rank_impt), 4))
        df$sdrank <- c(df$sdrank, round(sd(rank_impt), 4))
        df$meanrankall <- c(df$meanrankall, round(mean(rank_impt_all), 4))
        df$sdrankall <- c(df$sdrankall, round(sd(rank_impt_all), 4))
      }
    }
  }

  return(as.data.frame(df))
}

Simulation5 <- function(n_list, phi, link,
                        alpha = 0.05, rep = 100, d = 3) {
  #' @param n_list: list of n
  #' @param phi: bandwidth in gaussian kernel
  #' @param link: link function with two options, logistic or gaussian
  #' @param alpha: significant level
  #' @param rep: number of replications
  #' @param d: low rank used in USVT
  #' @return power

  n_len <- length(n_list)

  if (!(0 %in% eps_list)) {
    stop("eps_list must contain 0.")
  }

  power <- matrix(NA, nrow = n_len + 1, ncol = 3)
  power[1, ] <- c("n", "link", "power")
  idx <- 1

  for (n in n_list) {
    set.seed(110)
    X <- matrix(rnorm(2 * n), nrow = n, ncol = 2)
    D <- as.matrix(dist(X)^2)
    P <- exp(-D^2)

    if (link == "logistic") {
      Q <- exp(-D) / (1 + exp(-D))
    } else {
      Q <- exp(-D^2 / phi)
    }

    rho.hat <- c()

    for (N in 1:rep) {
      P.hat <- USVT(Generate.A(P), num = d)
      Q.hat <- USVT(Generate.A(P), num = d)
      rho.hat <- c(rho.hat, norm(P.hat - Q.hat))
    }

    qt <- quantile(rho.hat, 1 - alpha)
    rho.hat <- c()

    for (N in 1:rep) {
      P.hat <- USVT(Generate.A(P), num = d)
      Q.hat <- USVT(Generate.A(Q), num = d)
      rho.hat <- c(rho.hat, norm(P.hat - Q.hat))
    }

    idx <- idx + 1
    pow <- sum(rho.hat > qt) / length(rho.hat)

    if (link == "logistic") {
      power[idx, ] <- c(n, link, pow)
    } else {
      power[idx, ] <- c(n, paste0(link, "_", phi), pow)
    }

  }

  return(power)
}

Simulation6 <- function(n = 100, rep = 100, b_num = 1000, d = 2) {
  # required package "gtools"
  #' @param n: number of vertices
  #' @param d: low rank in USVT
  #' @param rep: number of Monte Carlo replicates
  #' @param b_num: number of bootstrap samplings
  #' @return pval_list: a list of p-values

  set.seed(110)
  pval_list <- c()

  for (N in 1:rep) {
    X <- matrix(rnorm(2 * n), nrow = n, ncol = 2)
    Z1 <- matrix(rnorm(2 * n, mean = 0, sd = sqrt(0)), nrow = n, ncol = 2)
    Y <- X + Z1
    Z2 <- matrix(rnorm(2 * n, mean = 0, sd = sqrt(0)), nrow = n, ncol = 2)
    W <- X + Z2
    D1 <- as.matrix(dist(X)^2)
    D2 <- as.matrix(dist(Y)^2)
    D3 <- as.matrix(dist(W)^2)
    P <- exp(-D1^2)
    Q <- exp(-D3^2 / 4)
    R <- exp(-D2^2 / 2)
    rho0 <- Spearman(Generate.A(P), Generate.A(R), d)
    rhovec1 <- numeric(b_num)
    rhovec2 <- numeric(b_num)

    for (i in 1:b_num) {
      rhovec1[i] <- Spearman(Generate.A(P), Generate.A(Q), d)
      rhovec2[i] <- Spearman(Generate.A(R), Generate.A(Q), d)
    }

    pval1 <- sum(rho0 > rhovec1) / b_num
    pval2 <- sum(rho0 > rhovec2) / b_num
    pval_list <- c(pval_list, max(pval1, pval2))
  }

  return(pval_list)
}

Simulation7 <- function(gamma, n = 100, eps = 0, d = 2, rep = 100, n0 = 8, alpha = 0.05, b_num = 1000) {
  #' @param n: number of vertices
  #' @param eps: variance of the distribution difference
  #' @param d: low rank in USVT
  #' @param rep: number of Monte Carlo replicates
  #' @param n0: number of permutations in the bootstrap algorithm
  #' @param alpha: significance level
  #' @param b_num: number of bootstrap samplings
  #' @param gamma: parameter to control the misspecification, ranging from 0 to 1
  #' @return power: power of the testing with significance level alpha

  pval_list <- c()
  A <- list()
  B <- list()

  for (N in 1:rep) {
    X <- matrix(rnorm(2 * n), nrow = n, ncol = 2)
    D1 <- as.matrix(dist(X)^2)
    Z <- matrix(rnorm(2 * n, sd = sqrt(eps)), nrow = n, ncol = 2)
    Y <- X + Z
    D2 <- as.matrix(dist(Y)^2)
    P <- gamma * exp(-D1^2) + (1 - gamma) * exp(X %*% t(X)) / (exp(X %*% t(X)) + 1)
    Q <- gamma * exp(-D2^2 / 2) + (1 - gamma) * exp(Y %*% t(Y)) / (exp(Y %*% t(Y)) + 1)
    A[[1]] <- Generate.A(P)
    B[[1]] <- Generate.A(Q)
    pval <- Bootstrap(A, B, b_num, d, n0)
    pval_list <- c(pval_list, pval)
  }

  power <- sum(pval_list <= alpha) / length(pval_list)

  return(power)
}
