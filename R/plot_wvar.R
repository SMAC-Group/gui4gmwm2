# there are currently to plotting methods for plot a wvar object, either in simts or in gmwm
# also we use gmwm version in gmwm2
# here I develop a plot method for wvar so we dont rely on these packages


my_plot_wvar = function(wv_obj, legend_position = "bottomright", show_legend = TRUE){
  # wv_obj = data[[1]][[1]]
  # str(wv_obj)

  
  
  if(class(wv_obj) != "wvar"){
    stop("wv_object must be of class wvar")
  }
  
  # get ymin and amax for y axis
  ymin = min(c(wv_obj$variance, wv_obj$ci_low))
  ymax = max(c(wv_obj$variance, wv_obj$ci_high)) 
  k_min <- floor(log10(ymin))
  k_max <- ceiling(log10(ymax))
  


  plot(x =wv_obj$scales, y=wv_obj$variance, log="xy", ylab="",
       ylim=c(ymin, ymax),
       xlab="", main="", type="l",xaxt="n", yaxt ="n")
  mtext("Wavelet Variance", side = 2, line = 3)
  mtext("Scales", side = 1, line = 2.5)
  

  
  yticks <- 10^(k_min:k_max)
  axis(2, at = yticks,    labels = parse(text = paste0("10^", k_min:k_max)), las = 1)
  
  # for x axis we can use the scales
  # get the scale in power of 2
  xmin <- min(wv_obj$scales)
  xmax <- max(wv_obj$scales)
  
  j_min <- floor(log2(xmin))
  j_max <- ceiling(log2(xmax))
  
  xticks <- 2^(j_min:j_max)
  axis(1,
       at = xticks,
       labels = parse(text = paste0("2^", j_min:j_max)))

  # add polygon
  polygon(c(wv_obj$scales, rev(wv_obj$scales)),
          c(wv_obj$ci_low, rev(wv_obj$ci_high)),
          col = "#ccf1f8", border = NA)
  
  grid( col="grey80", lty=2)

  lines(x = wv_obj$scales, y = wv_obj$variance, col="#00008b")
  points(x = wv_obj$scales, y = wv_obj$variance, col="#00008b", pch=16, cex=1.3)
  
  if (isTRUE(show_legend)) {
    legend(
      legend_position,
      legend = c("Wavelet Variance", "95% CI"),
      col = c("#00008b", NA),
      lty = c(1, NA),
      pch = c(16, 15),
      pt.cex = c(1.3, 2.2),
      pt.bg = c(NA, "#ccf1f8"),
      bty = "n"
    )
  }

  mtext(side=3, "Empirical Haar Wavelet Variance", line=1)
  box()
}

# library(imudata)
#
# data(imu6)
# library(wv)
# MTiG = list(
#   Gyro.X = wvar(imu6[,1]),
#   Gyro.Y = wvar(imu6[,2]),
#   Gyro.Z = wvar(imu6[,3]),
#   Accel.X = wvar(imu6[,4]),
#   Accel.Y = wvar(imu6[,5]),
#   Accel.Z = wvar(imu6[,6]),
#   freq = 100)
#
#
#
# wv:::plot.wvar(MTiG$Gyro.X)
# library(simts)
# fit = gmwm::gmwm(model = WN(), input = MTiG[[1]][[1]])
# fit
# plot(fit)
# plot_wvar = function(wv_object, CI=TRUE){
#   wv_obj = MTiG$Gyro.X
#
#   if(class(wv_obj) != "wvar"){
#     stop("wv_object must be of class wvar")
#   }
#
#   # plot
#   plot(wv_obj$variance, log="xy", ylab="Wavelet Variance",xlab="Scales", main="", type="l")
#
#
#
# }
