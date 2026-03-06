# 
# 
# # coming from summary.gmwm
# 
# 
# compute_bootstrap_ci_from_fit = function (object, B=100, frequency=1){
#   
#   # #----------------------------------
#   # DO NOT ERASE
#   # # test boot strap
#   # n = 10000
#   # sigma2 = 2
#   # x = rnorm(n = n,mean = 0, sd = sqrt(sigma2))
#   # wv_emp = wv::wvar(x)
#   # plot(wv_emp)
#   # object =gmwm::gmwm(model = WN()+GM(), input = wv_emp) 
#   # B = 100
#   # #---------------------------------------------------
#   
#   
#   
#   
#   B=10
#   # create matrix to store estimated param in the bootstrap procedure
#   mat_res = matrix(NA, nrow = B, ncol = length(object$estimate))
#   for(b in seq(B)){
#     # b=1
#     # set seed
#     set.seed(123 + b)
#     # generate data from signal
#     generated_ts = gen_model(fit$N, theta = fit$estimate, desc =   fit$model$desc, objdesc = fit$model$obj.desc);
#     # get wv on generated ts
#     wv_emp_generated_ts = wv::wvar(generated_ts)
#     # fit gmwm on generated ts
#     fit = gmwm::gmwm(model = fit$model, input = wv_emp_generated_ts)
#     # save estimated parameters
#     mat_res[b,] = fit$estimate
#   }
#   
#   
#   # compute wv
#   # wv_emp = wv::wvar(x)
#   
#   
#   
#   # summary(fit, inference = T)
#   # gmwm::summary.gmwm(fit, inference = T)
#   # simts:::summary.gmwm(fit, inference = T)
#   # inference = TRUE
#   # #---------------------------
#   
#   
#   
#   set.seed(object$seed + 5)
#   out = object$estimate
#   N = object$N
#   auto = if (N > 10000) FALSE  else TRUE
#   if (is.null(inference)) {
#     inference = auto
#   }
#   if (is.null(bs.gof)) {
#     bs.gof = if (inference) 
#       auto
#     else F
#   }
#   if (is.null(bs.gof.p.ci)) {
#     bs.gof.p.ci = if (inference) 
#       auto
#     else F
#   }
#   if (is.null(bs.theta.est)) {
#     bs.theta.est = if (inference) 
#       auto
#     else F
#   }
#   if (is.null(bs.ci)) {
#     bs.ci = if (inference) 
#       auto
#     else F
#   }
#   if ("ARMA" %in% object$model$desc) {
#     if (bs.ci == FALSE) {
#       warning(paste0("The numerical derivative of ARMA(p,q), where p > 1 and q > 1, may be inaccurate leading to inappropriate CIs.\n", 
#                      "Consider using the bs.ci = T option on the summary function."))
#     }
#   }
#   if (inference) {
#     if (any(object$model$desc == "GM")) {
#       object$estimate[, 1] = conv.gm.to.ar1(object$estimate[, 
#                                                             1], object$model$process.desc, object$freq)
#     }
#     mm = .Call("_gmwm_get_summary", PACKAGE = "gmwm", object$estimate, 
#                object$model$desc, object$model$obj.desc, object$model.type, 
#                object$wv.empir, object$theo, object$scales, object$V, 
#                solve(object$orgV), object$obj.fun, N, object$alpha, 
#                object$robust, object$eff, inference, F, bs.gof, 
#                bs.gof.p.ci, bs.theta.est, bs.ci, B)
#   }
#   else {
#     mm = vector("list", 3)
#     mm[1:3] = NA
#   }
#   if (inference) {
#     out.coln = colnames(out)
#     out = cbind(out, mm[[1]])
#     colnames(out) = c(out.coln, "CI Low", "CI High", "SE")
#     idx_gm = (object$model$desc == "GM")
#     if (any(idx_gm)) {
#       out[, 2:3] = apply(out[, 2:3], 2, FUN = conv.ar1.to.gm, 
#                          process.desc = object$model$process.desc, freq = object$freq)
#     }
#   }
#   x = structure(list(estimate = out, testinfo = mm[[2]], inference = inference, 
#                      bs.gof = bs.gof, bs.gof.p.ci = bs.gof.p.ci, bs.theta.est = bs.theta.est, 
#                      bs.ci = bs.ci, starting = object$starting, seed = object$seed, 
#                      obj.fun = object$obj.fun, N = N, freq = object$freq), 
#                 class = "summary.gmwm")
#   x
# }