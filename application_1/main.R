# prepocessing
source("../functions.R")
pt <- read.table("./participants.csv", sep = ",",
                 header = TRUE, stringsAsFactors = FALSE)
n <- dim(pt)[1]
you_id <- c()
mid_id <- c()
old_id <- c()
for (i in 1:n) {
  if (pt$age[i] >= 18 & pt$age[i] <= 35){
    you_id <- c(you_id, pt$participant_id[i])
  }
  if (pt$age[i] > 35 & pt$age[i] <= 56){
    mid_id <- c(mid_id, pt$participant_id[i])
  }
  if (pt$age[i] > 56 & pt$age[i] <= 85){
    old_id <- c(old_id, pt$participant_id[i])
  }
}

library(stringr)
path <- "./counts"
setwd(path)
files <- list.files(path = getwd(), pattern = "*.csv")
dt <- list()
old_list <- list()
mid_list <- list()
you_list <- list()

for (file in files) {
  dt[[file]] <- read.csv(file, header = FALSE)
  for (i in old_id) {
    if (str_detect(file, pattern = i) == "TRUE") {
      old_list[[i]] <- file
      break
    }
  }
  for (j in mid_id) {
    if (str_detect(file, pattern = j) == "TRUE") {
      mid_list[[j]] <- file
      break
    }
  }
  for (k in you_id) {
    if (str_detect(file, pattern = k) == "TRUE") {
      you_list[[k]] <- file
      break
    }
  }
}

# The function used to choose n0
choose_n0 <- function(file_list, N, n0, d) {
  #' @param file_list: a list of files
  #' @return rst: including estimated test statistics and power
  num <- length(file_list)
  Alist <- list()
  for (k in 1:num) {
    file <- file_list[[k]]
    A_temp <- read.csv(file, header = FALSE)
    A_temp <- as.matrix(A_temp)[-c(1, 2, 60), -c(1, 2, 60)]
    diag(A_temp) <- 0
    Alist[[k]] <- A_temp
  }
  n <- dim(Alist[[1]])[1]
  relabel1 <- sample(seq(num), num, replace = TRUE)
  relabel2 <- sample(seq(num), num, replace = TRUE)
  A <- B <- list()
  A.bar <- B.bar <- matrix(0, n, n)
  for (i in 1:num) {
    idx1 <- relabel1[i]
    A[[i]] <- Alist[[idx1]]
    A.bar <- A.bar + Alist[[idx1]]
    idx2 <- relabel2[i]
    B[[i]] <- Alist[[idx2]]
    B.bar <- B.bar + Alist[[idx2]]
  }
  A.bar <- A.bar / num
  B.bar <- B.bar / num

  rho <- Spearman(A.bar, B.bar, num = d, res = FALSE)
  rhovec1 <- rhovec2 <- numeric(N)
  for (b in 1:N) {
    A.perm <- B.perm <- matrix(0, n, n)
    for (k in 1:num) {
      A.temp <- A[[k]]
      rownames(A.temp) <- colnames(A.temp) <- seq(n)
      a <- seq(n)
      s <- sample(seq(n), n0)
      a[sort(s)] <- s
      A.perm <- A.perm + permutation(A.temp, a)
      B.temp <- B[[k]]
      rownames(B.temp) <- colnames(B.temp) <- seq(n)
      B.perm <- B.perm + permutation(B.temp, a)
    }
    A.perm <- A.perm / num
    B.perm <- B.perm / num
    rhovec1[b] <- Spearman(A.bar, A.perm, num = d, res = FALSE)
    rhovec2[b] <- Spearman(B.bar, B.perm, num = d, res = FALSE)
  }
  pval1 <- sum(rhovec1 < rho) / length(rhovec1)
  pval2 <- sum(rhovec2 < rho) / length(rhovec2)
  rst <- data.frame("rho.hat1" = rhovec1,
                    "rho.hat2" = rhovec2,
                    "rho" = rho,
                    "p.value" = max(pval1, pval2))
  return(rst)
}
# An example to choose n0
d <- 3
n0 <- 5
N <- 1000
repN <- 100

for (i = 1:repN) {
   rst <- choose_n0(file_list = old_list, N, n0, d)
   pval[i] <- rst$p.value[1]
   print(c(i, pval[i]))
}

# The function used in comparison
Comparison <- function(list1, list2, rep, n0, d) {
  #' @param list1, list2: two lists of adjacency matrices
  #' @return rst: including the observed test statistic, estimated test statistics and p-value
  num1 <- length(list1)
  num2 <- length(list2)
  n <- 128
  A <- list()
  B <- list()
  A1 <- A2 <- matrix(0, n, n)
  i <- 1
  for (file in list1) {
    A_temp <- read.csv(file, header = FALSE)
    A_temp <- as.matrix(A_temp)[-c(1, 2, 60), -c(1, 2, 60)]
    diag(A_temp) <- 0
    A1 <- A1 + A_temp
    A[[i]] <- A_temp
    i <- i + 1
  }
  i <- 1
  for (file in list2) {
    A_temp <- read.csv(file, header = FALSE)
    A_temp <- as.matrix(A_temp)[-c(1, 2, 60), -c(1, 2, 60)]
    diag(A_temp) <- 0
    A2 <- A2 + A_temp
    B[[i]] <- A_temp
    i <- i + 1
  }
  Abar1 <- A1 / num1
  Abar2 <- A2 / num2
  rho <- Spearman(Abar1, Abar2, num = d, res = FALSE)
  set.seed(3962)
  rhovec1 <- rhovec2 <- numeric(rep)
  for (i in 1:rep) {
    A.perm <- B.perm <- matrix(0, n, n)
    for (k in 1:max(num1, num2)) {
      if (k <= num1 & k <= num2) {
        A.temp <- A[[k]]
        rownames(A.temp) <- colnames(A.temp) <- seq(n)
        a <- seq(n)
        s <- sample(seq(n), n0)
        a[sort(s)] <- s
        A.perm <- A.perm + permutation(A.temp, a)
        B.temp <- B[[k]]
        rownames(B.temp) <- colnames(B.temp) <- seq(n)
        B.perm <- B.perm + permutation(B.temp, a)
      }
      if (k > num1 & k <= num2) {
        B.temp <- B[[k]]
        rownames(B.temp) <- colnames(B.temp) <- seq(n)
        a <- seq(n)
        s <- sample(seq(n), n0)
        a[sort(s)] <- s
        B.perm <- B.perm + permutation(B.temp, a)
      }
      if (k <= num1 & k > num2) {
        A.temp <- A[[k]]
        rownames(A.temp) <- colnames(A.temp) <- seq(n)
        a <- seq(n)
        s <- sample(seq(n), n0)
        a[sort(s)] <- s
        A.perm <- A.perm + permutation(A.temp, a)
      }
    }
    rhovec1[i] <- Spearman(Abar1, A.perm/num1, num = d, res = FALSE)
    rhovec2[i] <- Spearman(Abar2, B.perm/num2, num = d, res = FALSE)
  }
  pval1 <- sum(rhovec1 < rho) / length(rhovec1)
  pval2 <- sum(rhovec2 < rho) / length(rhovec2)
  rst <- data.frame("rho.hat1" = rhovec1,
                    "rho.hat2" = rhovec2,
                    "rho" = rho,
                    "p.value" = max(pval1, pval2))
  return(rst)
}


# Figure 4
rep <- 1000
rst_old_you <- Comparison(old_list, you_list, rep, n0 = 4, d = 3)
rst_old_mid <- Comparison(old_list, mid_list, rep, n0 = 4, d = 3)
rst_you_mid <- Comparison(you_list, mid_list, rep, n0 = 4, d = 3)

vec1 <- rep(NA, rep)
idx1 <- which(rst_old_you$rho.hat1 < rst_old_you$rho.hat2)
idx2 <- which(rst_old_you$rho.hat1 >= rst_old_you$rho.hat2)
vec1[idx1] <- rst_old_you$rho.hat1[idx1]
vec1[idx2] <- rst_old_you$rho.hat2[idx2]
rst_old_you$rho.hat <- vec1
rst_old_you$rho[1]

vec2 <- rep(NA, rep)
idx1 <- which(rst_old_mid$rho.hat1 < rst_old_mid$rho.hat2)
idx2 <- which(rst_old_mid$rho.hat1 >= rst_old_mid$rho.hat2)
vec2[idx1] <- rst_old_mid$rho.hat1[idx1]
vec2[idx2] <- rst_old_mid$rho.hat2[idx2]
rst_old_mid$rho.hat <- vec2
rst_old_mid$rho[1]

vec3 <- rep(NA, rep)
idx1 <- which(rst_you_mid$rho.hat1 < rst_you_mid$rho.hat2)
idx2 <- which(rst_you_mid$rho.hat1 >= rst_you_mid$rho.hat2)
vec3[idx1] <- rst_you_mid$rho.hat1[idx1]
vec3[idx2] <- rst_you_mid$rho.hat2[idx2]
rst_you_mid$rho.hat <- vec3
rst_you_mid$rho[1]

rho <- c(rst_old_you$rho.hat, rst_old_mid$rho.hat, rst_you_mid$rho.hat)
method <- c(rep("oldyoung", 100), rep("oldmid", 100), rep("midyoung", 100))
dt <- data.frame(rho, method)
dt$facets <- factor(dt$method, labels = c("H[0]:~X[middle]==X[young]",
                      "H[0]:~X[old]==X[middle]", "H[0]:~X[old]==X[young]"))

dt_text <- data.frame(
  label = c("observed~T[n]==0.9932", "observed~T[n]==0.6871",
            "observed~T[n]==0.6751"),
  facets = c("H[0]:~X[middle]==X[young]",
                      "H[0]:~X[old]==X[middle]", "H[0]:~X[old]==X[young]")
)
dt_text_2 <- data.frame(
  label = c("p-value==0.11", "p-value==0", "p-value==0"),
  facets = c("H[0]:~X[middle]==X[young]",
                      "H[0]:~X[old]==X[middle]", "H[0]:~X[old]==X[young]")
)

library(ggplot2)
fig <- ggplot(dt, aes(x = rho, ..scaled..)) +
  geom_histogram(aes(y = ..count..), colour = "black", fill = "white") +
  xlab(expression(paste("Estimated ", T[paste(n)], " under the null"))) +
  ylab("Frequency")

tiff("lifespan.tiff", units = "in", width = 15/2.3, height = 5/2.3, res = 300)
fig + facet_grid(~facets, labeller = label_parsed) + 
      theme(strip.text.x = element_text(size = 10, family = "sans")) +
      geom_text(data = dt_text,
                mapping = aes(x = -Inf, y = Inf, label = label),
                hjust = -0.06, vjust = 1.5,
                family = "mono", parse = TRUE, size = 3.5, color = "red"
      ) +
      geom_text(data = dt_text_2,
                mapping = aes(x = -Inf, y = Inf, label = label),
                hjust = -0.1, vjust = 3,
                family = "mono", parse = TRUE, size = 3.5, color = "red"
      )
dev.off()

# Table C.1
oldvsyou_list <- c(old_list, you_list)
n <- length(oldvsyou_list)

num_old <- length(old_list)
num_you <- length(you_list)

A_old <- matrix(0, nrow = 131, ncol = 131)
A_you <- matrix(0, nrow = 131, ncol = 131)

for (file in old_list) {
  A_temp <- read.csv(file, header = FALSE)
  A_temp <- as.matrix(A_temp)
  diag(A_temp) <- 0
  A_old <- A_old + A_temp
}

for (file in you_list) {
  A_temp <- read.csv(file, header = FALSE)
  A_temp <- as.matrix(A_temp)
  diag(A_temp) <- 0
  A_you <- A_you + A_temp
}

A_bar_old <- A_old[-c(1, 2, 60), -c(1, 2, 60)] / num_old
A_bar_you <- A_you[-c(1, 2, 60), -c(1, 2, 60)] / num_you

P <- USVT(A_bar_old, num = 3, res = FALSE)
Q <- USVT(A_bar_you, num = 3, res = FALSE)
p.upper <- P[upper.tri(P)]
q.upper <- Q[upper.tri(Q)]

library(doBy)
phat.rank <- rank(p.upper, ties.method = "average")
qhat.rank <- rank(q.upper, ties.method = "average")
diff.rank <- abs(phat.rank - qhat.rank)
A.star <- matrix(0, ncol = 128, nrow = 128)
A.star[upper.tri(A.star)] <- (diff.rank > quantile(diff.rank, 0.7))
A.star <- A.star + t(A.star)
esti_impt_node <- which.maxn(colSums(A.star), n = 20)
