my_colors <- c(
  "#8DA0CB",  # muted lavender-blue
  "#FC8D62",  # soft coral
  "#66C2A5",  # pastel teal
  "#E78AC3",  # dusty pink
  "#A6D854",  # gentle green
  "#FFD92F"   # soft mustard yellow
)


get_plot_colors <- function(n) {
  if (n <= 0L) return(character())
  if (n <= length(my_colors)) {
    return(my_colors[seq_len(n)])
  }
  grDevices::colorRampPalette(my_colors)(n)
}

my_plot_gmwm = function(gmwm_fit, legend_position = "bottomleft", show_legend = TRUE){
  
  
  
  #-----------
  # library(simts)
  # library(gmwm)
 
  # fit1 = gmwm::gmwm(model = WN() ,  data[[1]][[1]])
  # 
  # str(data[[1]][[1]])
  # str( wvar(rnorm(1000)))
  # fit1 = gmwm::gmwm(model = WN() , )

  # gmwm:::gmwm(model=WN(), input = data[[1]][[1]])
  
  # gmwm::gmwm(model = WN())
  # fit1$estimate
  # fit2 = gmwm::gmwm(model = WN()+GM(), data[[1]][[1]], freq=1)
  # fit2$estimate
  # gmwm_fit = fit1
  #-------------
  
  
  
  
  # first plot wavelet variance
  if(class(gmwm_fit) != "gmwm"){
    stop("wv_object must be of class gmwm")
  }
  
  # get ymin and amax for y axis
  ymin = min(c(gmwm_fit$wv$variance, gmwm_fit$wv$ci_low))
  ymax = max(c(gmwm_fit$wv$variance, gmwm_fit$wv$ci_high)) 
  k_min <- floor(log10(ymin))
  k_max <- ceiling(log10(ymax))
  
  
  
  plot(x =gmwm_fit$wv$scales, y=gmwm_fit$wv$variance, log="xy", ylab="",
       ylim=c(ymin, ymax),
       xlab="", main="", type="l",xaxt="n", yaxt ="n")
  mtext("Wavelet Variance", side = 2, line = 3)
  mtext("Scales", side = 1, line = 2.5)
  
  
  
  yticks <- 10^(k_min:k_max)
  axis(2, at = yticks,    labels = parse(text = paste0("10^", k_min:k_max)), las = 1)
  
  # for x axis we can use the scales
  # get the scale in power of 2
  xmin <- min(gmwm_fit$wv$scales)
  xmax <- max(gmwm_fit$wv$scales)
  
  j_min <- floor(log2(xmin))
  j_max <- ceiling(log2(xmax))
  
  xticks <- 2^(j_min:j_max)
  axis(1,
       at = xticks,
       labels = parse(text = paste0("2^", j_min:j_max)))
  
  # add polygon
  polygon(c(gmwm_fit$wv$scales, rev(gmwm_fit$wv$scales)),
          c(gmwm_fit$wv$ci_low, rev(gmwm_fit$wv$ci_high)),
          col = "#ccf1f8", border = NA)
  
  grid( col="grey80", lty=2)
  
  lines(x = gmwm_fit$wv$scales, y = gmwm_fit$wv$variance, col="#00008b")
  points(x = gmwm_fit$wv$scales, y = gmwm_fit$wv$variance, col="#00008b", pch=16, cex=1.3)
  

  
  mtext(side=3, "Haar Wavelet Variance Representation", line=1)
  box()
  
  
  
  
  
  
  # add theoretical wv as well as each component
  # check number of models
  n_models = ncol(gmwm_fit$decomp.theo)

    
    # get color
    my_col_plot = get_plot_colors(n_models)
    
    for(i in seq(n_models)){
      lines(x = gmwm_fit$wv$scales, y = gmwm_fit$decomp.theo[,i], col=my_col_plot[i])
    }
 

  
  # add theoretical WV
  
  lines(gmwm_fit$wv$scales, y=gmwm_fit$theo, type = "b", col="darkorange", cex=1.5)
  
  
  if (isTRUE(show_legend)) {
    
    
    comp_names <- NULL
    if (!is.null(gmwm_fit$model$desc)) {
      comp_names <- as.character(gmwm_fit$model$desc)
    } else if (!is.null(gmwm_fit$model$label)) {
      comp_names <- as.character(gmwm_fit$model$label)
    }
    if (is.null(comp_names) || length(comp_names) != n_models) {
      comp_names <- if (n_models > 0) paste("Component", seq_len(n_models)) else character()
    }

    legend_labels <- c("Empirical WV", "95% CI", "Theoretical WV", comp_names)
    legend_cols <- c("#00008b", "#ccf1f8", "darkorange", if (n_models > 0) my_col_plot else character())
    legend_lty <- c(1, NA, 1, if (n_models > 0) rep(1, n_models) else integer())
    legend_pch <- c(16, 15, 1, if (n_models > 0) rep(NA, n_models) else integer())
    legend_ptcex <- c(1.3, 2.2, 1.5, if (n_models > 0) rep(NA, n_models) else numeric())
    legend_ptbg <- c(NA, "#ccf1f8", NA, if (n_models > 0) rep(NA, n_models) else character())

    legend(
      legend_position,
      legend = legend_labels,
      col = legend_cols,
      lty = legend_lty,
      pch = legend_pch,
      pt.cex = legend_ptcex,
      pt.bg = legend_ptbg,
      bty = "n"
    )
  }
}
