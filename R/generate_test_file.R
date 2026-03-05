# # white noise and random walk
# 
# rm(list=ls())
# library(simts)
# library(gmwm)
# 
# # ---------------------- white noise
# frequency = 400
# q = .1 # PSD intensity of the white noise
# sqrt(q)
# delta_t = 1/frequency
# sigma2 = q/delta_t
# sigma2
# n = 100000
# x1 = rnorm(n, mean = 0, sd = sqrt(sigma2))
# 
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
# write.csv(y, "R/data/test_data.csv", row.names = FALSE)