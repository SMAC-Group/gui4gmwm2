gmwm_fit_to_ts_model <- function(object) {
  # This function takes a fitted gmwm object and constructs the corresponding time series model

  #-----------------------------------------------------------------------
  # x = gen_gts(WN(sigma2 = 3)+RW(gamma2 = 2), n = 100)
  # object = gmwm(model = WN() + RW()+GM(), input = wv::wvar(x))
  # object
  
  #-----------------------------------------------------------------------
  
  
  
  
  if (!inherits(object, "gmwm")) {
    stop("object must be of class gmwm")
  }

  model_desc <- as.character(object$model$desc)
  est <- object$estimate
  if (is.matrix(est) || is.data.frame(est)) {
    if ("Estimate" %in% colnames(est)) {
      est <- est[, "Estimate"]
    }
  }
  est <- as.numeric(est)

  th_param <- 1
  model <- NULL

  for (i in seq_along(model_desc)) {
    process_name <- model_desc[[i]]

    component <- if (process_name == "WN") {
      sigma2 <- est[[th_param]]
      th_param <- th_param + 1
      WN(sigma2 = sigma2)
    } else if (process_name == "RW") {
      gamma2 <- est[[th_param]]
      th_param <- th_param + 1
      RW(gamma2 = gamma2)
    } else if (process_name == "GM") {
      beta <- est[[th_param]]
      sigma2_gm <- est[[th_param + 1]]
      th_param <- th_param + 2
      GM(beta = beta, sigma2_gm = sigma2_gm)
    } else if (process_name == "QN") {
      q2 <- est[[th_param]]
      th_param <- th_param + 1
      QN(q2 = q2)
    } else if (process_name == "DR") {
      omega <- est[[th_param]]
      th_param <- th_param + 1
      DR(omega = omega)
    } else {
      stop(paste("Unsupported process in bootstrap generation:", process_name))
    }

    if (is.null(model)) {
      model <- component
    } else {
      model <- model + component
    }
  }

  model
}


compute_bootstrap_ci_from_fit <- function(object, B = 100, frequency = 1, progress_callback = NULL) {
  
  # 
  # object = res
  # B=100
  # frequency=1
  # progress_callback = NULL
  # 

  
  
  if (!inherits(object, "gmwm")) {
    stop("object must be of class gmwm")
  }
  B <- max(1, as.integer(B))
  
  bootstrap_model <- gmwm_fit_to_ts_model(object)
  n_param <- length(object$estimate)
  mat_res <- matrix(NA_real_, nrow = B, ncol = n_param)
  
  for (b in seq_len(B)) {
    if (is.function(progress_callback)) {
      progress_callback(b, B)
    }
    
    set.seed(123 + b)
    
    if(is_exact_wn_gm_model(bootstrap_model)) {
      # if the model is WN + GM we use the custom gen function, otherwise we use the simts gen function
      mod_to_generate_from = WN(sigma2 = object$estimate_ar1_parameters[1]) + AR1(phi = object$estimate_ar1_parameters[2], sigma2 = object$estimate_ar1_parameters[3])
      generated_ts =  simts::gen_gts(
        n = object$N,
        model = mod_to_generate_from
      )
      
      
    } else {
    
    generated_ts <- simts::gen_gts(
      n = object$N,
      model = bootstrap_model
    )
    }
    wv_emp_generated_ts <- wv::wvar(as.numeric(generated_ts))
    
    # execute differently if the model is WN + GM or not, if it is WN + GM we use the custom fit function, otherwise we use the gmwm fit function
    if(is_exact_wn_gm_model(bootstrap_model)) {
      fit <- fit_gmwm_wn_gm(wv_emp_generated_ts)
    } else {
      fit <- gmwm::gmwm(model = object$model, input = wv_emp_generated_ts)
    }
    
    df_param_trans <- transform_parameters(fit, frequency = frequency)
    mat_res[b, ] <- df_param_trans$`Estimated transformed parameters`
  }
  
  apply(mat_res, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
}




