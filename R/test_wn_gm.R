

# ---------------------- Gauss-Markov (FOGM)
rm(list=ls())
source("R/my_plot_gmwm.R")




frequency <- 100
delta_t <- 1 / frequency

beta_gm <- 2
q_gm <- 1
q_wn <- 0.001
sigma2_wn <- q_wn / delta_t
sigma2_wn 

# AR(1) parameters implied by continuous GM at delta_t
phi_ar1 <- exp(-beta_gm * delta_t)
sigma2_ar1 <- (q_gm / (2 * beta_gm)) * (1 - exp(-2 * beta_gm * delta_t))
sigma2_ar1
n=10000
x <- simts::gen_gts(AR1(phi = phi_ar1, sigma2 = sigma2_ar1) + WN(sigma2 = sigma2_wn), n = n) 
wv_obj <- wv::wvar(x)
plot(wv_obj)




# test if fit at true value
# note that gmwm assume a freq of 1, so the parameters are transfomed accoridng to that

# transform ar1 param to gm with freq of 1
beta_gm_freq_1 <- -log(phi_ar1) / 1
q_gm_freq_1 <- (2 * beta_gm_freq_1 * sigma2_ar1 ) / (1 - exp(-2 * beta_gm_freq_1 * 1))

fit <- gmwm::gmwm(model = WN(sigma2 = sigma2_wn)+ GM(beta = beta_gm_freq_1, sigma2_gm = q_gm_freq_1), input = wv_obj)
fit2 <- gmwm::gmwm(model = WN(sigma2 = sigma2_wn)+ AR1(phi_ar1, sigma2 = sigma2_ar1), input = wv_obj)
fit3 = gmwm::gmwm(model = WN()+ GM(), input = wv_obj)
fit4 = gmwm::gmwm(model = WN()+ AR1(), input = wv_obj)

my_plot_gmwm(fit)
my_plot_gmwm(fit2)
my_plot_gmwm(fit3)
my_plot_gmwm(fit4)

# # try with smart starting values
# fit5 <- gmwm::gmwm(model = WN(sigma2 = .)+ AR1(phi = phi_ar1, sigma2 = sigma2_ar1), input = wv_obj)
# my_plot_gmwm(fit5)
# 
# 
# model = WN(1)+ GM(1,2)

trans_from_real_to_minus_1_and_1 <- function(x) {
  eps = 1e-6
  (1 - eps) * tanh(x)
}

inv_trans_from_real_to_minus_1_and_1 <- function(x) {
  eps = 1e-6
  # safety clamp in case user gives boundary value
  x <- pmin(pmax(x, -(1 - eps)), (1 - eps))
  atanh(x / (1 - eps))
}

# inv_trans_from_real_to_minus_1_and_1(trans_from_real_to_minus_1_and_1(-10))

# manual implementation

loss_fn_gmwm_wn_ar1 = function(theta, wv_obj, omega =NULL){
  # theta = sigma2 wn, phi ar1 , sigma2 wn
  sigma2_wn = exp(theta[1])
  phi_ar1 = trans_from_real_to_minus_1_and_1(theta[2])
  sigma2_ar1 = exp(theta[3])
  
  # obtain theo wv
  theo_wv = wv::ar1_to_wv(phi = phi_ar1, sigma2 = sigma2_ar1, tau = wv_obj$scales) + wv::wn_to_wv(sigma2 = sigma2_wn, tau = wv_obj$scales)
  
  # define omega
  nu_hat <- wv_obj$variance
  if (is.null(omega)) {
    omega <- diag(1 / (wv_obj$ci_low - wv_obj$ci_high)^2)
  }
  
  difference <- nu_hat - theo_wv
  objective <- as.numeric(t(difference) %*% omega %*% difference)
  
  
  return(objective)
}


# check if we provide true value
phi_ar1 = .99
sigma2_wn = 5
sigma2_ar1 = 2
n = 20000
x <- simts::gen_gts(AR1(phi = phi_ar1, sigma2 = sigma2_ar1) + WN(sigma2 = sigma2_wn), n = n) 
wv_obj <- wv::wvar(x)
plot(wv_obj)
init_param = c(log(sigma2_wn), inv_trans_from_real_to_minus_1_and_1(phi_ar1), log(sigma2_ar1))

B =100
mat_res = matrix(NA, nrow=B, ncol=3)
for(b in seq(B)){
  x <- simts::gen_gts(AR1(phi = phi_ar1, sigma2 = sigma2_ar1) + WN(sigma2 = sigma2_wn), n = n) 
  wv_obj <- wv::wvar(x)
  init_param = c(log(sigma2_wn), inv_trans_from_real_to_minus_1_and_1(phi_ar1), log(sigma2_ar1))
  res = optim(fn = loss_fn_gmwm_wn_ar1, par =init_param, wv_obj = wv_obj)
  mat_res[b,] =  c(exp(res$par[1]), trans_from_real_to_minus_1_and_1(res$par[2]), exp(res$par[3]))
}

boxplot(mat_res[,1])
abline(h=sigma2_wn)
boxplot(mat_res[,2])
abline(h=phi_ar1)
boxplot(mat_res[,3])
abline(h=sigma2_ar1)


B =500
mat_res2 = matrix(NA, nrow=B, ncol=3)
mat_res3 = matrix(NA, nrow=B, ncol=3)
# now simu with bad starting value
for(b in seq(B)){
  # b = 1
  set.seed(123 + b)
  x <- simts::gen_gts(AR1(phi = phi_ar1, sigma2 = sigma2_ar1) + WN(sigma2 = sigma2_wn), n = n) 
  wv_obj <- wv::wvar(x)
  init_param = c(log(var(x)), runif(1,.1, .99), log(var(x)))
  res = optim(fn = loss_fn_gmwm_wn_ar1, par =init_param, wv_obj = wv_obj)
  mat_res2[b,] =  c(exp(res$par[1]), trans_from_real_to_minus_1_and_1(res$par[2]), exp(res$par[3])) 
  init_param2 =  c(var(x), init_param[2], var(x))
  res2 = gmwm::gmwm(model = WN(sigma2 = init_param2[1]) + AR1(phi = init_param2[2], sigma2 = init_param2[3]), input = wv_obj)
  mat_res3[b,] = res2$estimate
}

# check own implementation
boxplot(mat_res2[,1])
abline(h = sigma2_wn)
boxplot(mat_res2[,2])
abline(h = phi_ar1)
boxplot(mat_res2[,3])
abline(h = sigma2_ar1)



boxplot(mat_res2[,1], mat_res3[,1])
abline(h=sigma2_wn)
boxplot(mat_res2[,2], mat_res3[,2])
abline(h=phi_ar1)
boxplot(mat_res2[,3], mat_res3[,3])
abline(h=sigma2_ar1)

res = optim(fn = loss_fn_gmwm_wn_ar1, par =init_param, wv_obj = wv_obj)
res
# trasnform back
c(exp(res$par[1]), trans_from_real_to_minus_1_and_1(res$par[2]), exp(res$par[3]))
# see sensibility if we provide bad starting value

res2 = optim(fn = loss_fn_gmwm_wn_ar1, par = c(log(wv_obj$variance)), wv_obj = wv_obj)
res





# generate starting value
sigma2_start_wn = var(x)
phi_start_ar1 = runif(1, 0.05, .95) 
sigma2_start_ar1 = var(x)

# convert these to gm parameter with frequency of 1 (assumed by gmwm)
# Convert back to continuous GM params
beta_start_gm <- -log(phi_start_ar1) / frequency
q_start_gm <- (2 * beta_start_gm * sigma2_ar1 ) / (1 - exp(-2 * beta_start_gm * 1))


# try with these starting value
fit <- gmwm::gmwm(model = WN(sigma2 = sigma2_start_wn)+ GM(beta = beta_start_gm, sigma2_gm = q_start_gm*10), input = wv_obj)
fit2 = gmwm::gmwm(WN(var(x))+AR1(.99, var(x)), wv_obj)
my_plot_gmwm(fit)
my_plot_gmwm(fit2)  
  




# # Compare AR1 params
  # print(c(phi_true = phi, phi_hat = phi_hat))
  # print(c(sigma2_true = sigma2, sigma2_hat = sigma2_hat))
  # 
  # Convert back to continuous GM params
  beta_hat <- -log(phi_hat) / delta_t
  q_hat <- (2 * beta_hat * sigma2_hat) / (1 - exp(-2 * beta_hat * delta_t))
  
  mat_res[b,] = c(beta_hat, q_hat, phi_hat, sigma2_hat)  
}

