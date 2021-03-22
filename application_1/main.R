# prepocessing
#setwd("")
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

# function used in comparison
Comparison <- function(list1, list2, rep = 100) {
  #' @param list1, list2: two lists of adjacency matrices
  #' @return rst: including estimated test statistics and power
  combine_list <- c(list1, list2)
  n <- length(combine_list)
  num1 <- length(list1)
  num2 <- length(list2)

  dt <- list()
  A1 <- matrix(0, nrow = 131, ncol = 131)
  A2 <- matrix(0, nrow = 131, ncol = 131)
  for (file in list1) {
    A_temp <- read.csv(file, header = FALSE)
    A_temp <- as.matrix(A_temp)
    diag(A_temp) <- 0
    A1 <- A1 + A_temp
  }
  for (file in list2) {
    A_temp <- read.csv(file, header = FALSE)
    A_temp <- as.matrix(A_temp)
    diag(A_temp) <- 0
    A2 <- A2 + A_temp
  }
  Abar1 <- A1[-c(1, 2, 60), -c(1, 2, 60)] / num1
  Abar2 <- A2[-c(1, 2, 60), -c(1, 2, 60)] / num2
  rho <- Spearman(Abar1, Abar2, num = 3, res = FALSE)
  rho.hat <- c()
  rep <- 100
  set.seed(110)
  for (i in 1:rep) {
    relabel <- sample(1:n, num1, replace = FALSE)
    A.perm <- matrix(0, nrow = 131, ncol = 131)
    B.perm <- matrix(0, nrow = 131, ncol = 131)
    for (file in combine_list) {
      if (file %in% combine_list[relabel]) {
        A_temp <- list()
        A_temp <- read.csv(file, header = FALSE)
        A_temp <- as.matrix(A_temp)
        diag(A_temp) <- 0
        A.perm <- A.perm + A_temp
      }else{
        B_temp <- list()
        B_temp <- read.csv(file, header = FALSE)
        B_temp <- as.matrix(B_temp)
        diag(B_temp) <- 0
        B.perm <- B.perm + B_temp
      }
    }
    A.perm_bar <- A.perm[-c(1, 2, 60), -c(1, 2, 60)] / num1
    B.perm_bar <- B.perm[-c(1, 2, 60), -c(1, 2, 60)] / num2
    rst <- Spearman(A.perm_bar, B.perm_bar, num = 3, res = FALSE)
    rho.hat <- c(rho.hat, rst)
  }
  pval <- sum(rho.hat < rho) / length(rho.hat)
  rst <- data.frame("rho.hat" = rho.hat,
                    "rho" = rho,
                    "p-value" = pval)
  return(rst)
}

# Figure 4
rst_old_you <- Comparison(old_list, you_list)
rst_old_mid <- Comparison(old_list, mid_list)
rst_you_mid <- Comparison(you_list, mid_list)

library(ggplot2)
p1 <- ggplot(rst_old_you, aes(x = rho.hat, ..scaled..)) +
  geom_histogram(aes(y = ..count..), colour = "black", fill = "white") +
  xlab(expression(paste("Estimated ", T[paste(n)], " under the null"))) +
  ylab("Frequencies") +
  geom_point(aes(x = rho, y = 0), colour = "blue")
p1
p2 <- ggplot(rst_old_mid, aes(x = rho.hat)) +
  geom_histogram(aes(y = ..count..), colour = "black", fill = "white") +
  xlab(expression(paste("Estimated ", T[paste(n)], " under the null"))) +
  ylab("Frequencies") +
  geom_point(aes(x = rho, y = 0), colour = "blue")
p2
p3 <- ggplot(rst_you_mid, aes(x = rho.hat, ..scaled..)) +
   geom_histogram(aes(y = ..count..), colour = "black", fill = "white") +
  xlab(expression(paste("Estimated ", T[paste(n)], " under the null"))) +
  ylab("Frequencies") +
  geom_point(aes(x = rho, y = 0), colour = "blue")
p3

# Table C.1
oldvsyou_list <- c(old_list, you_list)
n <- length(oldvsyou_list)

num_old <- length(old_list)
num_you <- length(you_list)

A_old <- matrix(0, nrow = 131, ncol = 131)
A_you <- matrix(0, nrow = 131, ncol = 131)

for (file in old_list){
  A_temp <- read.csv(file, header = FALSE)
  A_temp <- as.matrix(A_temp)
  diag(A_temp) <- 0
  A_old <- A_old + A_temp
}

for (file in you_list){
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
