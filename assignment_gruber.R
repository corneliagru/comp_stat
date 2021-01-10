
# import ------------------------------------------------------------------


data <- read.table("data.txt")
x <- data$V1

source("preparation.R")


# description -------------------------------------------------------------

#overview data
round(summary(x), digits = 2 )




# find best parameter settings. apparently there is only one group, all deltas
# are basically the same 2.37
param_k5 <- EM(x, k = 5)
param_k5$delta

# only one class
param <- EM(x, k = 1)
delta <- param$delta
p <- param$p




# plot data and density ---------------------------------------------------

xx = seq(0,8, by = 0.1)
yy <- sapply(xx, dens_plot, delta = delta, p = p)

plot(density(x), ylim = c(0,0.4), main = "")
lines(xx, yy, col = "red")
legend("topright", legend = c("observed data", "delta = 2.37"),
       col = c("black", "red"),lty = 1,  cex = 1.1)






# Metropolis hastings -----------------------------------------------------

#find settings for metropolis hastings algorithm

res <- MCMC(theta_0 = 1, delta = delta, p = p, n = 7000,  sd = 4)

#acceptance around 30% wenn sd = 4
res$acceptance

acf(res$mc, lag.max = 50)
plot(ts(res$mc))
# -> n = 7000, burn in 1000, lag 20


samp <- MCMC_iid(res$mc, n = 300, lag = 20, burn_in  = 1000)

# looks very similar
plot(density(samp))
lines(density(x), col = "red")




# confidence intervals ----------------------------------------------------

# check if non parallel works
res_bs <- bs_param(B = 10, delta = delta, p = p, parallel = FALSE)
apply(res_bs, 2, quantile, probs = c(0.025, 0.975))




#prepare
param <- EM(x, k = 1)
delta <- param$delta
p <- param$p

#parallelize 
numCores <- detectCores()
registerDoParallel(numCores)
source("preparation.R")
start <- Sys.time()
res_bs <- bs_param(B = 1000, delta = delta, p = p, theta_0 = 1, 
                   n = 300*lag + burn_in, sd = 4, lag = 20,
                   burn_in = 1000,  parallel = TRUE)
end <- Sys.time()
print(end - start)
stopImplicitCluster()





ci <- apply(res_bs, 2, quantile, probs = c(0.025, 0.975))
ci
d_lower <- ci[1,1]
d_upper <- ci[2,1]




# plot CI ----------------------------------------------------------------


xx = seq(0,8, by = 0.1)
y_estim <- sapply(xx, dens_plot, delta = delta, p = p)
yy <- sapply(xx, dens_plot, delta = d_lower, p = p)
yyy <- sapply(xx, dens_plot, delta = d_upper, p = p)


plot(density(x), ylim = c(0,0.45), main = "", lwd = 2)
lines(xx, y_estim, col = "red", lty = "solid", lwd = 2)
lines(xx, yy, col = "gray20", lty = "dotted", lwd = 2)
lines(xx, yyy, col = "gray20", lty = "longdash", lwd = 2)
legend("topright", legend =
         c(paste0("delta = ", round(d_lower, digits = 2)),
           paste0("delta = ", round(delta, digits = 2)),
           paste0("delta = ", round(d_upper, digits = 2)),
           "true density"),
       col = c("gray20","red", "gray20", "black"), 
       lty = c("dotted","solid", "longdash", "solid"), cex = 1.1, lwd = 2)
