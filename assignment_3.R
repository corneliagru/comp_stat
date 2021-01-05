
## General important stuff -------------------------------------------------


# Load data and packages --------------------------------------------------

set.seed(123)
data <- read.table("~/R/Course_Comp_Stat/Assignment 3/data.txt", 
                   quote = "\"", comment.char = "")
x <- data$V1

library(doParallel)
source("functions.R") 
# had to do this for parallelisation to work

cores <- detectCores()
cl <- makeCluster(cores[1]-2) #not to overload computer
registerDoParallel(cl)

##  Number 1 ---------------------------------------------------------------

# initial values
k <- 1   #number of components
delta <- rep(2, k)  ## initial values for deltas
p <- rep(1/k,k)   #initial values for mixing proportions

em <-  EM_fun(x, k, delta, p)
delta_estim <- em$delta
delta_estim
p_estim <- em$p
p_estim

curve(dens(x, delta_estim, p_estim), 0, 20, n = 1000, las = 1)
# can't really be trapped somewhere!

# no matter how much i play around, there seems to be no groups and delta is
# around 2.3739

##  Number 2 ---------------------------------------------------------------

mcmc <- MCMC_fun(1, delta_estim, p_estim, n = length(x))
# plot(mcmc$all_values)
# plot(mcmc$keep_values)
# hist(mcmc$all_values, 50)
# hist(mcmc$keep_values, 50)
# # mit Visualisierung muss ich mich noch auseinander setzen

# Bootstrap

Bootstrap_fun <- function(x, delta, p, B = 1000) {
  k <- length(delta)
  
  foreach(i = 1:B, .combine = 'rbind') %dopar% {
    source("~/R/Course_Comp_Stat/Assignment 3/functions.R")
    sample <- MCMC_fun(1, delta, p, n = length(x)) #get sample
    em_b <- EM_fun(x, k, delta, p) # initial values are the estim ones
    cbind(t(em_b$delta), t(em_b$p))
  }
}

bootstrap <- Bootstrap_fun(x, delta_estim, p_estim, B = 10)
# und es kommen immer die gleichen values raus


stopCluster(cl)
