# white noise and random walk

# rm(list=ls())
# library(simts)
# library(gmwm)
# 
# # ---------------------- white noise
# frequency = 100
# q = .1 # PSD intensity of the white noise
# sqrt(q)
# delta_t = 1/frequency
# sigma2 = q/delta_t
# sigma2
# n = 100000
# x1 = rnorm(n, mean = 0, sd = sqrt(sigma2))

# 
# q <- 2
# sqrt(q)
# 
# delta_t <- 1 / frequency
# # increment variance
# sigma2_inc <- q * delta_t
# sigma2_inc
# eta <- rnorm(n, mean = 0, sd = sqrt(sigma2_inc))
# # build the random walk
# x2 <- cumsum(eta)
# 
# y = x1 + x2
# # print to csv
write.csv(x1, "R/data/test_data.csv", row.names = FALSE)













# GM
# 
# 
# # ---------------------- Gauss-Markov (FOGM)
# 
# frequency <- 100
# delta_t <- 1 / frequency
# 
# beta <- 10
# 1/beta
# q <- 2
# sqrt(q)
# 
# 
# # AR(1) parameters implied by continuous GM at delta_t
# phi <- exp(-beta * delta_t)
# sigma2 <- (q / (2 * beta)) * (1 - exp(-2 * beta * delta_t))
# n=50000
# x <- simts::gen_gts(AR1(phi = phi, sigma2 = sigma2), n = n)
# write.csv(x, "R/data/test_data.csv", row.names = FALSE)
