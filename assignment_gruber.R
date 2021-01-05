#data import

source('preparation.R')

data <- read.table("data.txt")
x <- data$V1

hist(data$V1, breaks = 30)


plot(density(data$V1))
lines(x, 1*f_x_giv_delta(x, 2.8), col = "red")


x = seq(0,8, by = 0.1)
plot(x, 0.7*f_x_giv_delta(x, 2.8) + 0.3*f_x_giv_delta(x, 5), type = "l")
lines(x, 1*f_x_giv_delta(x, 2.8), col = "red")
lines(x, 0.3*f_x_giv_delta(x, 5), col = "red")

