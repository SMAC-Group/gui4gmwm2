
rm(list=ls())


# Custom WN + GM (via WN + AR1 equivalence) prototype
# Goal: return a gmwm-compatible object so existing app code keeps working.

# Map R -> (-1, 1) for stable AR1 phi optimization.
# trans_from_real_to_minus_1_and_1 <- function(x) {
#   eps <- 1e-6
#   (1 - eps) * tanh(x)
# }

# Inverse map (-1, 1) -> R.
# inv_trans_from_real_to_minus_1_and_1 <- function(x) {
#   eps <- 1e-6
#   x <- pmin(pmax(x, -(1 - eps)), (1 - eps))
#   atanh(x / (1 - eps))
# }




# ---------------- okay lets do a test now where we generate data from a specific frequency

n=10000
frequency <- 100
delta_t <- 1 / frequency

beta <- 2
1/beta
q <- 2
sqrt(q)

# AR(1) parameters implied by continuous GM at delta_t
phi <- exp(-beta * delta_t)
phi
sigma2 <- (q / (2 * beta)) * (1 - exp(-2 * beta * delta_t))
sigma2

# this cause a problem when you boostrap currently with GM!!!
# simts:::gm_to_ar1(c(beta, q), freq=100)

# wn continuous param
q_wn = .005 # PSD intensity of the white noise
sqrt(q_wn)
sigma2_wn = q_wn/delta_t
sigma2_wn
x <- simts::gen_gts(AR1(phi = phi, sigma2 = sigma2) + WN(sigma2_wn), n = n)
# export in csv as test data
write.csv(x, "R/data/test_data_wn_gm.csv", row.names = FALSE)

# emp_wv = wv::wvar(x)
# plot(emp_wv)
# res = fit_gmwm_wn_gm(emp_wv)
# source("R/my_plot_gmwm.R")
# my_plot_gmwm(res)
# source("R/transform_parameters.R")
# transform_parameters(res, frequency = frequency)
# source("R/parametric_bootstrap_for_ci.R")
# compute_bootstrap_ci_from_fit(res, B=100, frequency = frequency)

