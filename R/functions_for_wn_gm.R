
real_to_unit_eps <- function(theta, eps = 1e-6) {
  eps + (1 - 2 * eps) / (1 + exp(-theta))
}

unit_eps_to_real <- function(phi, eps = 1e-6) {
  log((phi - eps) / (1 - eps - phi))
}







# GMWM objective for WN + AR1 parameterization. Here phi 
loss_fn_gmwm_wn_ar1 <- function(theta, wv_obj, omega = NULL) {
  sigma2_wn <- exp(theta[1])
  phi_ar1 <- real_to_unit_eps(theta[2])
  sigma2_ar1 <- exp(theta[3])
  
  theo_wv <-
    wv::wn_to_wv(sigma2 = sigma2_wn, tau = wv_obj$scales) +
    wv::ar1_to_wv(phi = phi_ar1, sigma2 = sigma2_ar1, tau = wv_obj$scales)
  
  nu_hat <- wv_obj$variance
  if (is.null(omega)) {
    omega <- diag(1 / (wv_obj$ci_low - wv_obj$ci_high)^2)
  }
  
  difference <- nu_hat - theo_wv
  as.numeric(t(difference) %*% omega %*% difference)
}

# Check if a ts.model corresponds to exactly one WN and one GM component.
is_exact_wn_gm_model <- function(model) {
  desc <- as.character(model$desc)
  length(desc) == 2 && sum(desc == "WN") == 1 && sum(desc == "GM") == 1
}



# define a gmwm like function for that specific noise
fit_gmwm_wn_gm = function(wv_obj){
  # Initial parameter guess: log(sigma2_wn), atanh(phi_ar1), log(sigma2_ar1)
  init_theta <- c(log(wv_obj$variance[1] / 2), unit_eps_to_real(runif(1,.1, .99)), log(wv_obj$variance[1] / 2))
  
  # Optimize the GMWM objective function for WN + AR1.
  opt_result <- optim(
    par = init_theta,
    fn = loss_fn_gmwm_wn_ar1,
    wv_obj = wv_obj,
    method = "L-BFGS-B"
  )
  
  # Extract optimized parameters and transform back to original scale.
  sigma2_wn_est <- exp(opt_result$par[1])
  phi_ar1_est <-real_to_unit_eps(opt_result$par[2])
  sigma2_ar1_est <- exp(opt_result$par[3])
  
  
  # create elements similar to a gmwm object so the plot and bootstrap works
  decomp.theo <- matrix(NA, nrow = length(wv_obj$scales), ncol = 2)
  decomp.theo[,1] <- wv::wn_to_wv(sigma2_wn_est,tau = wv_obj$scales)
  decomp.theo[,2] <- wv::ar1_to_wv(phi =phi_ar1_est, sigma2 = sigma2_ar1_est, tau = wv_obj$scales)
  theo <- decomp.theo[,1] + decomp.theo[,2]
  
  # transform AR1 parameters to GM parameters assuming a frequency of 1, this is then retreated by transform parameters later with a user specific provided frequency
  beta_est <- -log(phi_ar1_est)
  q_est <- (2 * beta_est * sigma2_ar1_est) / (1 - exp(-2 * beta_est))
  
  
  # define model
  
  
  
  
  # Create a gmwm-like object with the estimated parameters.
  out = list(
    estimate = c(sigma2_wn_est, beta_est, q_est),
    estimate_ar1_parameters = c(sigma2_wn_est, phi_ar1_est, sigma2_ar1_est),
    # estimate2 = c(sigma2_wn_est, phi_ar1_est, sigma2_ar1_est), 
    model = WN(sigma2 = sigma2_wn_est) + GM(beta = beta_est, sigma2_gm = q_est),
    N = wv_obj$N,
    wv = wv_obj,
    starting=FALSE,
    decomp.theo = decomp.theo,
    theo = theo
  )
  class(out) <- "gmwm"
  return(out)
}
