library(ggplot2)
library(ggthemes)
library(ggsci)
library(gap)
### Simulation 1
pval_list <- Simulation1(n = 100, m = 5, eps = 0, d = 2, n0 = 5)
qqunif(pval_list, type = "unif", logscale = FALSE)

### Simulation 2 - Bootstrap
for (eps in seq(0.01, 0.05, by = 0.005)) {
  power <- Simulation2.Bootstrap(n = 100, eps, d = 2, n0 = 8, alpha = 0.05)
  print(c(eps, power))
}

### Simulation 2 - NMDS
n_list <- c(50, 100, 150)
eps_list <- seq(0, 0.1, 0.01)
rst <- Simulation2.NMDS(n_list, eps_list, alpha = 0.05)
print(rst)

### Simulation 3
Simulation3(n_list = c(100, 200, 500, 1000),
             eps_list = c(0, 0.02, 0.1, 0.2, 0.5, 1), rho = 5)


### Simulation 4
dt <- Simulation4(nlist = seq(100, 1000, 1000),
                   slist = c(20, 40, 60), d = 3, tau = 0.7)

# plots for simulation 5, please change the parameters correspondingly
dt$meanrank_v2 <- dt$meanrank - (1 - min(dt$prop, (0.3 * dt$n)) / (0.3 * dt$n))
dt$meanrank_all <- dt$meanrankall - (1 - min(dt$prop, (dt$n)) / (dt$n))
dt$prop <- as.character(dt$prop)
dt$s <- dt$prop
plt <- ggplot(dt, aes(x = n, y = mean)) +
        geom_line(aes(color = s), size = 1) +
        geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                      color = "darkorange", width = 35, size = 0.7) +
        theme_bw() +
        geom_point(aes(color = s), size = 1.5) +
        labs(x = "Number of vertices, n",
              y = "Averaged number of important vertices in estimation") +
        theme(legend.position = c(0.885, 0.35),
              legend.background = element_rect(colour = "grey")) +
        scale_color_manual(values = c("brown3", "royalblue3", "darkolivegreen"))
plt

dt_20 <- subset(dt, prop == "20", select = n:meanrank_all)
dt_40 <- subset(dt, prop == "40", select = n:meanrank_all)
dt_60 <- subset(dt, prop == "60", select = n:meanrank_all)

plt2_all <- ggplot(dt_20, aes(x = n, y = meanrankall))
plt2_all +
  geom_errorbar(aes(min = meanrankall - sdrankall,
                    ymax = meanrankall + sdrankall),
                    width = 35, size = 0.7, color = "darkorange") +
  geom_point(size = 1.5, color = "indianred4") +
  geom_line(size = 1, color = "indianred4") +
  stat_function(fun = function(x) fun_all(x, s = 20),
                size = 1, color = "darkgrey") +
  ylim(c(0.6, 1)) + theme_bw() +
  labs(x = "Number of vertices, n",
        y = "Averaged median rank of Procrustes distance")

plt3_all <- ggplot(dt_40, aes(x = n, y = meanrankall))
plt3_all +
  geom_errorbar(aes(min = meanrankall - sdrankall,
                    ymax = meanrankall + sdrankall),
                    width = 35, size = 0.7, color = "darkorange") +
  geom_point(size = 1.5, color = "deepskyblue4") +
  geom_line(size = 1, color = "deepskyblue4") +
  stat_function(fun = function(x) fun_all(x, s = 40),
                size = 1, color = "darkgrey") +
  ylim(c(0.6, 1)) + theme_bw() +
  labs(x = "Number of vertices, n",
        y = "Averaged median rank of Procrustes distance")
plt4_all <- ggplot(dt_60, aes(x = n, y = meanrankall))
plt4_all +
  geom_errorbar(aes(min = meanrankall - sdrankall,
                    ymax = meanrankall + sdrankall),
                    width = 35, size = 0.7, color = "darkorange") +
  geom_point(size = 1.5, color = "darkolivegreen") +
  geom_line(size = 1, color = "darkolivegreen") +
  stat_function(fun = function(x) fun_all(x, s = 60),
                  size = 1, color = "darkgrey") +
  ylim(c(0.6, 1)) + theme_bw() +
  labs(x = "Number of vertices, n",
        y = "Averaged median rank of Procrustes distance")

### Simulation 5
for(phi in c(0.8, 1.5, 2, 4)) {
  rst <- Simulation5(n_list = c(50, 100, 200, 500, 1000),
                       phi = phi, link = "gaussian")
  print(rst)
}
Simulation5(n_list = c(50, 100, 200, 500, 1000),
              link = "logistic")
