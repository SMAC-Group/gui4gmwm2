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

transform_parameters = function(gmwm_fit, frequency){
  # return a table with the name of the process and the transformed parameters
  
  
  
  
  # frequency given in Hertz
  #---------------------------------------------
  # gmwm_fit = fit
  # frequency = 100
  #-----------------------------------------------
  
  
  
  if (!inherits(gmwm_fit, "gmwm")) {
    stop("model must be of class gmwm")
  }
  
  # take a gmwm fit, transform parameters for use in extended kalman filter
  number_models = length(gmwm_fit$model$desc)
  number_of_total_parameters = length(gmwm_fit$estimate)
  number_of_parameters_for_each_model = sapply(gmwm_fit$model$obj.desc, length)
  
  # create table
  df_transformed_parameters = data.frame(
    model = NA,
    parameter = NA,
    transformed_parameter =NA 
  )
  
  th_param = 1
  for(i in seq(number_models)){

    # get model name
    model_name = gmwm_fit$model$desc[i]
    
    
    if(model_name == "WN"){
      sigma2_trans = gmwm_fit$estimate[th_param] / frequency
      # append in df
      df_transformed_parameters[th_param,] = c(model_name, "SIGMA2", sigma2_trans)
      th_param = th_param + 1
    }else if(model_name == "RW"){
      gamma_2_trans = gmwm_fit$estimate[th_param] * frequency
      # append
      df_transformed_parameters[th_param,] = c(model_name, "GAMMA2", gamma_2_trans)
      th_param = th_param + 1
    }else if(model_name == "GM"){
      # transform back to AR1 parameters with frequency of 1
      gm_beta <- gmwm_fit$estimate[th_param]
      gm_psd <- gmwm_fit$estimate[th_param + 1]
      
      phi <- exp(-gm_beta / 1)
      sigma2 <- gm_psd / gm_beta / 2 * (1 - exp(-2 * gm_beta / 1))
      
      # transform now phi and sigma2 of AR1 to GM parameters with frequency
      beta <- -log(phi) * frequency
      psd <- 2 * beta * sigma2 / (1 - exp(-2 * beta / frequency))
      
      # append
      df_transformed_parameters[th_param,] = c(model_name, "BETA", beta)
      df_transformed_parameters[th_param + 1,] = c(model_name, "SIGMA2_GM", psd)
      th_param = th_param + 2
    }else if(model_name == "QN"){
      df_transformed_parameters[th_param,] = c(model_name, "Q2", 0)
      th_param = th_param + 1
    }else if(model_name == "DR"){
      df_transformed_parameters[th_param,] = c(model_name, "OMEGA", 0)
      th_param = th_param + 1
    }

  }
  return(df_transformed_parameters)
}





