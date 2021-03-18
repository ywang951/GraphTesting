library(ggplot2)

dt <- read.csv("sim2_pval.csv")

ggplot(dt, aes(x = boot_n150)) +
 geom_histogram(aes(y = ..density..), colour = "black", fill = "white") +
 geom_density(alpha = .2, fill = "#FF6666") + xlab("p-values") + ylab("Density")

ggplot(dt, aes(x = perm_n150)) +
 geom_histogram(aes(y = ..density..), colour = "black", fill = "white") +
 geom_density(alpha = .2, fill="#FF6666") + ylim(c(0, 8)) +
 xlab("p-values") + ylab("Density")