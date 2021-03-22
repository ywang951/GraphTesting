library(irlba)
library(vegan)
library(ecodist)
library(MCMCpack)
library(pracma)
library(igraph)
library(doBy)

GenerateProb <- function(n, eps, s = 1) {
  #' @param n: number of vertices
  #' @param eps: epsilon
  #' @return output: vectors of upper trianguler of true edge-probability matrices
  X <- matrix(rnorm(2 * n), nrow = n, ncol = 2)
  Z <- matrix(rnorm(2 * n, mean = 0, sd = sqrt(eps)), nrow = n, ncol = 2)
  Y <- X + Z
  D1 <- as.matrix(dist(X)^2)
  D2 <- as.matrix(dist(Y)^2)
  P <- s * exp(-D1^2)
  Q <- s * exp(-D2^2 / 4)
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

Simulation_1 <- function(n_list, eps_list, d = 3, alpha = 0.05, rep = 100) {
  #' @param n_list: list of n
  #' @param eps_list: list of epsilon
  #' @param d: low rank used in USVT
  #' @param alpha: significant level
  #' @param rep: number of replications
  #' @return power
  n_len <- length(n_list)
  eps_len <- length(eps_list)
  power <- matrix(NA, nrow = n_len * eps_len * 2 + 1, ncol = 4)
  power[1, ] <- c("n", "epsilon", "power", "methods")
  pval_dist <- matrix(NA, nrow = rep + 1, ncol = n_len * 2)
  idx <- 1
  idx2 <- 0
  for (n in n_list) {
    for (eps in eps_list) {
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
        p_val <- PermutationTest(Alist, Blist, d)
        p_val_list <- c(p_val_list, p_val)
      }
      pow <- sum(p_val_list <= alpha) / length(p_val_list)
      idx <- idx + 1
      power[idx, ] <- c(n, eps, pow, "Permutation Test")
      if (eps == 0) {
        idx2 <- idx2 + 1
        pval_dist[1, idx2] <- paste0("n", n, "_perm")
        pval_dist[2:(rep + 1), idx2] <- p_val_list
      }
      p_val_list <- c()
      for (N in 1:rep) {
        prob <- GenerateProb(n, eps)
        adj <- GenerateAdj(prob$P, prob$Q)
        p_val <- Bootstrap(adj$A, adj$B, b_num = 1000, d)
        p_val_list <- c(p_val_list, p_val)
      }
    pow <- sum(p_val_list <= alpha) / length(p_val_list)
    idx <- idx + 1
    power[idx, ] <- c(n, eps, pow, "Bootstrap")
    if (eps == 0) {
        idx2 <- idx2 + 1
        pval_dist[1, idx2] <- paste0("n", n, "_boots")
        pval_dist[2:(rep + 1), idx2] <- p_val_list
      }
    }
  }
  output <- list()
  power <- as.data.frame(power)
  names(power) <- power[1, ]
  power <- power[-1, ]
  pval_dist <- as.data.frame(pval_dist)
  names(pval_dist) <- pval_dist[1, ]
  pval_dist <- pval_dist[-1, ]
  output$pow <- power
  output$pval <- pval_dist
  return(output)
}

Simulation_2 <- function (n_list, eps_list, d = 3, alpha = 0.05, rep = 100) {
  #' @param n_list: list of n
  #' @param eps_list: list of epsilon
  #' @param d: low rank used in USVT
  #' @param alpha: significant level
  #' @param rep: number of replications
  #' @return power
  n_len <- length(n_list)
  eps_len <- length(eps_list)
  if (!(0 %in% eps_list)) {
      stop("eps_list must contain 0.")
  }
  power <- matrix(NA, nrow = n_len * (eps_len - 1)  + 1, ncol = 4)
  power[1, ] <- c("n", "epsilon", "Proposed test", "Non-metric embedding")
  idx <- 1
  for (n in n_list) {
    for (eps in eps_list) {
      prob <- GenerateProb(n, eps)
      nmds.hat <- c()
      rho.hat <- c()
      for (N in 1:rep) {
        A_bar <- matrix(0, ncol = n, nrow = n)
        B_bar <- matrix(0, ncol = n, nrow = n)
        for (b in 1:20) {
          adj <- GenerateAdj(prob$P, prob$Q)
          A_bar <- A_bar + adj$A
          B_bar <- B_bar + adj$B
        }
        A_bar <- A_bar / 20
        B_bar <- B_bar / 20
        rho <- Spearman(A_bar, B_bar, num = d, res = FALSE)
        rho.hat <- c(rho.hat, rho)
        nmds.hat <- c(nmds.hat, Nonmetric_MDS(A_bar, B_bar))
      }
      if (eps == 0) {
        qt_spm <- quantile(rho.hat, alpha)
        qt_nmds <- quantile(nmds.hat, 1 - alpha)
      } else{
        idx <- idx + 1
        pow_spm <- sum(rho.hat < qt_spm) / length(rho.hat)
        pow_nmds <- sum(nmds.hat > qt_nmds) / length(nmds.hat)
        power[idx, ] <- c(n, eps, pow_spm, pow_nmds)
      }
    }
  }
  power <- as.data.frame(power)
  names(power) <- power[1, ]
  power <- power[-1, ]
  return(power)
}

Simulation_3 <- function(n_list, eps_list, rho, alpha = 0.05, rep = 100, d = 3) {
  #' @param n_list: list of n
  #' @param eps_list: list of epsilon
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
  power <- matrix(NA, nrow = n_len * (eps_len - 1)  + 1, ncol = 4)
  power[1, ] <- c("n", "epsilon", "sparsity", "power")
  idx <- 1
  for (n in n_list) {
    s <- rho * log(n) / n
    for (eps in eps_list) {
      set.seed(110)
      prob <- GenerateProb(n, eps, s)
      rho.hat <- c()
      for (N in 1:rep) {
        adj <- GenerateAdj(prob$P, prob$Q)
        rho.hat <- c(rho.hat, Spearman(adj$A, adj$B, num = d))
      }
       if (eps == 0) {
        qt <- quantile(rho.hat, alpha)
      } else{
        idx <- idx + 1
        pow_spm <- sum(rho.hat < qt) / length(rho.hat)
        power[idx, ] <- c(n, eps, s * n / log(n), pow_spm)
      }
    }
  }
  return(power)
}

Simulation_4 <- function(nlist, slist, d = 3, tau = 0.7, rep = 100) {
  #' @param nlist: list of n
  #' @param slist: list of s
  #' @param d: low rank used in USVT
  #' @param tau: probability in quantile
  #' @param rep: number of replications
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
          adj <- GenerateAdj(P, Q)
          P.hat <- USVT(adj$A, num = d)
          Q.hat <- USVT(adj$B, num = d)
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
          rank_impt <- c(rank_impt, median(rank(diff_X_Ytrans[true_impt_node])[match(esti_impt_node[esti_impt_node%in%true_impt_node],true_impt_node)]/length(true_impt_node)))
          rank_impt_all <- c(rank_impt_all, median(rank(diff_X_Ytrans)[esti_impt_node]/n))
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

Simulation_5 <- function(n_list, phi, link,
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
    }else {
      Q <- exp(-D^2 / phi)
    }
    rho.hat <- c()
    for (N in 1:rep) {
      adj <- GenerateAdj(P, P)
      P.hat <- USVT(adj$A, num = d)
      Q.hat <- USVT(adj$B, num = d)
      rho.hat <- c(rho.hat, norm(P.hat - Q.hat))
    }
    qt <- quantile(rho.hat, 1 - alpha)
    rho.hat <- c()
    for (N in 1:rep) {
      adj <- GenerateAdj(P, Q)
      P.hat <- USVT(adj$A, num = d)
      Q.hat <- USVT(adj$B, num = d)
      rho.hat <- c(rho.hat, norm(P.hat - Q.hat))
    }
    idx <- idx + 1
    pow <- sum(rho.hat > qt) / length(rho.hat)
    if (link == "logistic") {
      power[idx, ] <- c(n, link, pow)
    }else {
      power[idx, ] <- c(n, paste0(link, "_", phi), pow)
    }
  }
  return(power)
}