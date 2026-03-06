rm(list=ls())
library(simts)
library(gmwm)

# ---------------------- white noise
frequency = 100 
q = 2 # PSD intensity of the white noise
delta_t = 1/frequency
sigma2 = q/delta_t
B = 1000
n = 1000
mat_res = matrix(NA, nrow = B, ncol = 2)
for(b in seq(B)){

  set.seed(123 + b)
  # simulate white noise
  x = rnorm(n, mean = 0, sd = sqrt(sigma2))

  # compute wv
  wv_obj = wv::wvar(x)

  fit = gmwm::gmwm(model = WN(sigma2 = sigma2), input = wv_obj)
  
  q_hat = fit$estimate[1]/frequency 
  mat_res[b,] = c( fit$estimate[1], q_hat)
}
boxplot(mat_res[,1])
abline(h = sigma2)
boxplot(mat_res[,2])
abline(h = q )



# ----------------------------------random walk
frequency <- 100
q <- 2
delta_t <- 1 / frequency
# increment variance
sigma2_inc <- q * delta_t
B =1000
n=10000
mat_res = matrix(NA, nrow = B, ncol = 2)
for(b in seq(B)){
  set.seed(123 + b)
  # simulate increments
  eta <- rnorm(n, mean = 0, sd = sqrt(sigma2_inc))
  # build the random walk
  x <- cumsum(eta)
  # compute wavelet variance
  wv_obj <- wv::wvar(x)
  # fit random walk with GMWM
  fit <- gmwm::gmwm(model = RW(), input = wv_obj)


  # recover q
  q_hat <- fit$estimate[1] *frequency
  mat_res[b, ] = c(fit$estimate[1], q_hat)
  
}

boxplot(mat_res[,1])
abline(h = sigma2_inc)
boxplot(mat_res[,2])
abline(h = q)


# 
# 
# # verify AR1 correcttly recover AR1 parameters
# 
# B=100
# mat_res = matrix(NA, nrow = B, ncol = 2)
# phi = .9
# sigma2 = 2
# for(b in seq(B)){
#   # simulate process
#   n <- 10000
#   set.seed(123 + b)
#   x <- simts::gen_gts(AR1(phi =phi, sigma2 =sigma2), n = n)
#   
#   wv_obj <- wv::wvar(x)
#   fit <- gmwm::gmwm(model = AR1(phi = phi, sigma2 = sigma2), input = wv_obj)
#   
#   # AR1 estimates
#   phi_hat <- fit$estimate[1]
#   sigma2_hat <- fit$estimate[2]
#   
#   # # Convert back to continuous GM params
#   # beta_hat <- -log(phi_hat) / delta_t
#   # q_hat <- (2 * beta_hat * sigma2_hat) / (1 - exp(-2 * beta_hat * delta_t))
#   
#   # save in matrix
#   mat_res[b,] = c(phi_hat, sigma2_hat)
# }
# 
# boxplot(mat_res[,1])
# abline(h=phi)
# boxplot(mat_res[,2])
# abline(h=sigma2)




# ---------------------- Gauss-Markov (FOGM)

frequency <- 100
delta_t <- 1 / frequency

beta <- 10
q <- 2

# AR(1) parameters implied by continuous GM at delta_t
phi <- exp(-beta * delta_t)
sigma2 <- (q / (2 * beta)) * (1 - exp(-2 * beta * delta_t))

B = 1000
mat_res = matrix(NA, nrow = B, ncol = 4)
n=10000
for(b in seq(B)){
  # b=1
  set.seed(123 + b)
  x <- simts::gen_gts(AR1(phi = phi, sigma2 = sigma2), n = n)
  wv_obj <- wv::wvar(x)
  # plot(wv_obj)
  fit <- gmwm::gmwm(model = AR1(phi = phi, sigma2 = sigma2), input = wv_obj)
  # fit
  # lines(x=fit$wv$scales, y = fit$decomp.theo)
  
  # AR1 estimates
  phi_hat <- fit$estimate[1]
  sigma2_hat <- fit$estimate[2]
  
  # # Compare AR1 params
  # print(c(phi_true = phi, phi_hat = phi_hat))
  # print(c(sigma2_true = sigma2, sigma2_hat = sigma2_hat))
  # 
  # Convert back to continuous GM params
  beta_hat <- -log(phi_hat) / delta_t
  q_hat <- (2 * beta_hat * sigma2_hat) / (1 - exp(-2 * beta_hat * delta_t))
  
  mat_res[b,] = c(beta_hat, q_hat, phi_hat, sigma2_hat)  
}


boxplot(mat_res[,1])
abline(h=beta)
boxplot(mat_res[,2])
abline(h=q)
boxplot(mat_res[,3])
abline(h=phi)
boxplot(mat_res[,4])
abline(h=sigma2)



#---------------------------------------------- drift

frequency <- 100
delta_t <- 1/frequency

omega <- 2     # continuous drift rate (diamond / s)
n <- 10000
B <- 1000

mu = omega * delta_t  # drift per step (diamond)
x = mu * 1:n
mat_res <- matrix(NA, nrow=B, ncol=2)
colnames(mat_res) <- c("omega_hat", "mu_hat")

for (b in 1:B) {
  # no need seed, process is determinstic
  x = mu * 1:n
  wv_obj <- wv::wvar(x)
  # plot(wv_obj)
  # Try DR() first:
  fit <- gmwm::gmwm(model = DR(), input = wv_obj)
  
  mu_hat <- fit$estimate[1]        # drift per step (diamond)
  omega_hat <- mu_hat * frequency  # back to diamond/s
  
  mat_res[b,] <- c(omega_hat, mu_hat)
}

boxplot(mat_res[, "omega_hat"], main="omega_hat"); abline(h=omega, col="red", lwd=2)
boxplot(mat_res[, "mu_hat"], main="mu_hat"); abline(h=mu, col="red", lwd=2)


mod = DR(2)
n=10000
x=gen_gts(mod, n=n)


# QN
mod = QN(q2 = 2)
n=10000
x = gen_gts(mod, n=n)
fit = gmwm::gmwm(model = QN(), input = wv::wvar(x))
fit$estimate














# test boot strap
n = 10000
sigma2 = 1
x = rnorm(n = n,mean = 0, sd = sqrt(sigma2))
wv_emp = wv::wvar(x)
plot(wv_emp)
fit =gmwm::gmwm(model = WN(), input = wv_emp) 
summary(fit, inference = T)
gmwm::summary.gmwm
simts:::summary.gmwm(fit, inference = T)


