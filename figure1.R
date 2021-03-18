library(ggplot2)

dt <- read.csv("sim2_power.csv")
dt$n <- as.factor(dt$n)

p <- ggplot(dt, aes(x = epsilon, y = Power)) +
  geom_line(aes(color = n), size = 1) +
  geom_point(size = 2.2, aes(color = n), shape = 21, fill = "white") +
  scale_color_manual(values = c("olivedrab", "navy", "orange")) +
  labs(x = "epsilon", y = "Power") +
  scale_x_continuous(breaks = seq(0, 0.1, 0.02)) + theme_bw()
p + facet_grid(cols = vars(Method))
