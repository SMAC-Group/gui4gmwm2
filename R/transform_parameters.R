

transform_parameters = function(gmwm_fit, frequency){
  
  
  
  
  
  # frequency given in Hertz
    #---------------------------------------------
  # load("~/github_repo/gui4gmwm2/R/data/imudata.RData")
  # library(simts)
  # gmwm_fit = gmwm::gmwm(model = WN()+ RW(), input = data[[1]][[1]] )
  # frequency = 100
  # #--------------------------
  
  
  
  
  # return a table with the name of the process and transformed parameters
  if (!inherits(gmwm_fit, "gmwm")) {
    stop("model must be of class gmwm")
  }

  model_desc = as.character(gmwm_fit$model$desc)
  number_models = length(model_desc)
  gm_total = sum(model_desc == "GM")
  gm_idx = 0
  rows = list()
  th_param = 1

  for(i in seq(number_models)){
    model_name = model_desc[i]
    model_label = model_name
    if (model_name == "GM") {
      gm_idx = gm_idx + 1
      if (gm_total > 1) {
        model_label = paste0("GM", gm_idx)
      }
    }

    if(model_name == "WN"){
      sigma2_trans = gmwm_fit$estimate[th_param] / frequency
      rows[[length(rows) + 1]] = data.frame(
        "Model" = model_label,
        "Parameter" = "\\(\\sqrt{q}\\)",
        "Estimated transformed parameters" = sqrt(sigma2_trans),
        # MathJax/LaTeX unit example for White Noise
        "Units" = "\\(\\frac{\\diamond}{\\sqrt{\\mathrm{Hz}}} \\)",
        stringsAsFactors = FALSE
      )
      th_param = th_param + 1
    }else if(model_name == "RW"){
      gamma_2_trans = gmwm_fit$estimate[th_param] * frequency
      rows[[length(rows) + 1]] = data.frame(
        "Model" = model_label,
        "Parameter" ="\\(\\sqrt{q}\\)",
        "Estimated transformed parameters" = sqrt(gamma_2_trans),
        "Units" = "\\(\\frac{\\diamond}{s \\sqrt{\\mathrm{Hz}}} \\)",
        stringsAsFactors = FALSE
      )
      th_param = th_param + 1
    }else if(model_name == "GM"){
      # transform back to AR1 parameters with frequency of 1
      
      
      
      # double check with simulation that the transformation is correct
      # frequency <- 100
      # delta_t <- 1 / frequency
      # 
      # beta <- 10
      # q <- 2
      # 
      # # AR(1) parameters implied by continuous GM at delta_t
      # phi <- exp(-beta * delta_t)
      # sigma2 <- (q / (2 * beta)) * (1 - exp(-2 * beta * delta_t))
      # 
      # n=10000
      # set.seed(123 )
      # x <- simts::gen_gts(AR1(phi = phi, sigma2 = sigma2), n = n)
      # wv_obj <- wv::wvar(x)
      # fit <- gmwm::gmwm(model = AR1(phi = phi, sigma2 = sigma2), input = wv_obj)
      # phihat = fit$estimate[1]
      # sigma2hat = fit$estimate[2]
      # 
      # 
      # 
      # betahat <- -log(phihat) * frequency
      # psdhad <- 2 * betahat * sigma2hat / (1 - exp(-2 * betahat / frequency))
      # 
      # 
      # 
      
      
      
      
      gm_beta <- gmwm_fit$estimate[th_param]
      gm_psd <- gmwm_fit$estimate[th_param + 1]

      phi <- exp(-gm_beta / 1)
      sigma2 <- gm_psd / gm_beta / 2 * (1 - exp(-2 * gm_beta / 1))

      # transform now phi and sigma2 of AR1 to GM parameters with frequency
      beta <- -log(phi) * frequency
      psd <- 2 * beta * sigma2 / (1 - exp(-2 * beta / frequency))
      


      rows[[length(rows) + 1]] = data.frame(
        "Model" = model_label,
        "Parameter" = "\\(\\frac{1}{\\beta} \\)",
        "Estimated transformed parameters" =1/beta ,
        "Units" = "s",
        stringsAsFactors = FALSE
      )
      rows[[length(rows) + 1]] = data.frame(
        "Model" = model_label,
        "Parameter" = "\\(\\sqrt{q} \\)",
        "Estimated transformed parameters" = sqrt(psd),
        "Units" = "\\(\\frac{\\diamond}{s \\sqrt{\\mathrm{Hz}}} \\)",
        stringsAsFactors = FALSE
      )
      th_param = th_param + 2
    }else if(model_name == "QN"){
      q2_hat = gmwm_fit$estimate[th_param] 
      rows[[length(rows) + 1]] = data.frame(
        "Model" = model_label,
        "Parameter" = "\\(Q\\)",
        "Estimated transformed parameters" = sqrt(q2_hat),
        "Units" = "\\(\\diamond\\)",
        stringsAsFactors = FALSE
      )
      th_param = th_param + 1
    }else if(model_name == "DR"){
      mu_hat <- gmwm_fit$estimate[th_param]
      omega_hat = mu_hat * frequency
      
      rows[[length(rows) + 1]] = data.frame(
        "Model" = model_label,
        "Parameter" = "\\(\\omega\\)",
        "Estimated transformed parameters" = omega_hat,
        "Units" =  "\\(\\frac{\\diamond}{s} \\)",
        stringsAsFactors = FALSE
      )
      th_param = th_param + 1
    }
  }

  df_transformed_parameters = do.call(rbind, rows)
  names(df_transformed_parameters) <- gsub("\\.", " ", names(df_transformed_parameters))
  return(df_transformed_parameters)
}






# library(simts)
# mod = WN(sigma2 = 1)+GM(beta = .01, sigma2_gm = 1)
# x = simts::gen_gts(mod, n=10000)
# emp_wv = wvar(x)
# plot(emp_wv)
# fit = gmwm(model = RW()+WN(), emp_wv)
# fit
# 
# # function get frequency
# get_frequency = function(){
#   
# }
# 
# str(fit$model)
# transform_parameters(fit, frequency = 100)
