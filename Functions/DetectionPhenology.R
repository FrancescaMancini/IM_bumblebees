plot_DetectionPhenology_modified <- function (model, spname = NULL, bins = 12, density_function = TRUE) 
{
  require(boot)
  require(plyr)
  require(sparta)
  
  data <- model$BUGSoutput$sims.list
  if (!"beta1" %in% names(data)) 
    stop("no phenological effect was modelled!")
  pDet1 <- rowMeans(data$dtype1.p)
  jul_dates <- seq(from = 1, to = 365, length.out = bins)
  if (density_function == TRUE) {
    if (!"beta3" %in% names(data)) 
      stop("beta3 not found. Please check that the density function method was used")
    pDet <- melt(sapply(jul_dates, function(jd) {
      pDet1 + data$beta3[, 1] * (1/((2 * pi)^0.5 * data$beta2[, 
                                                              1]) * exp(-((jd - data$beta1[, 1])^2/(2 * data$beta2[, 
                                                                                                                   1]^2))))
    }))
  }
  else {
    cjd <- jul_dates - 182
    pDet <- melt(sapply(cjd, function(jd) {
      pDet1 + jd * data$beta1[, 1] + jd^2 * data$beta2[, 
                                                       1]
    }))
  }
  names(pDet) <- c("it", "bin", "lgt_pDet")
  pDet$pDet <- inv.logit(pDet$lgt_pDet)
  pDet_summary <- ddply(pDet, .(bin), summarise, mean_pDet = mean(pDet), 
                        lower95CI = quantile(pDet, 0.025), upper95CI = quantile(pDet, 
                                                                                0.975))
  if (density_function == FALSE) {
    pDet_summary$cJulDate <- cjd[pDet_summary$bin]
  }
  pDet_summary$JulianDay <- jul_dates[pDet_summary$bin]
  gp <- ggplot(data = pDet_summary, x = JulianDay, y = mean_pDet) + 
    geom_line(aes(x = JulianDay, y = mean_pDet)) + geom_ribbon(aes(x = JulianDay, 
                                                                   ymin = lower95CI, ymax = upper95CI), alpha = 0.2) + ylab("Detection probability") + 
    ggtitle(spname) + theme_bw()
  gp
}
