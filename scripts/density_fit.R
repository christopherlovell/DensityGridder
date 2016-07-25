
library(ggplot2)

library(latex2exp)

lbox = 2168.64
grid = 300

pmass = 5.43 # M_{\odot} / h * 10^{10}

#'
#' Load density data at two epochs
#'

dir <- "../data/"

snap_1 <- "weights_002_"  # z = 4.687 
snap_2 <- "weights_022_"  # z = 0

data_1 <- scan(paste0(dir,snap_1,grid,'.bin'),what=integer())
data_2 <- scan(paste0(dir,snap_2,grid,'.bin'),what=integer())

# convert to M_{\odot} h^-1 Mpc^-3 * 10^{10}
data_1 <- data_1 * (lbox/grid)**-3 * pmass
data_2 <- data_2 * (lbox/grid)**-3 * pmass


hist_1 <- hist(data_1, breaks=100, plot=F)
hist_2 <- hist(data_2, breaks=1400, plot=F)


library(MASS)

x <- seq(0,50,length=200)

fit_params_1 <- fitdistr(data_1,"lognormal")
fit_1 <- dlnorm(x, fit_params_1$estimate['meanlog'], fit_params_1$estimate['sdlog'])

fit_params_2 <- fitdistr(data_2,"lognormal")
fit_2 <- dlnorm(x, fit_params_2$estimate['meanlog'], fit_params_2$estimate['sdlog'])

png("z4p687_density_hist.png", width = 2000, height = 2000, res = 200)
plot(x, fit_1, type="l", ylab="Density", 
     xlab=TeX("$M_{\\degree} h^{-1} Mpc^{-3} \\times 10^{10}$"), ylim=c(0,max(hist_1$density)))
title(main = "z = 4.687 density histogram with lognormal fit")
lines(hist_1$mid, hist_1$density, type="l", col="red")
legend(40,0.15,legend=c("Fit","Data"),lty=c(1,1),col=c("black","red"))
dev.off()

png("z0p0_density_hist.png", width = 2000, height = 2000, res = 200)
plot(x, fit_2, type="l", ylab="Density", 
     xlab=TeX("$M_{\\degree} h^{-1} Mpc^{-3} \\times 10^{10}$"), ylim=c(0,max(hist_2$density)))
title(main = "z = 0.0 density histogram with lognormal fit")
lines(hist_2$mid, hist_2$density, type="l", col="red")
legend(40,0.15,legend=c("Fit","Data"),lty=c(1,1),col=c("black","red"))
dev.off()


# QQ plot

quants <-seq(0,1,length=81)[2:80]

png("z4p687_QQ.png", width = 2000, height = 2000, res = 200)
plot(qlnorm(quants,fit_params_1$estimate['meanlog'], fit_params_1$estimate['sdlog']),quantile(data_1,quants), xlab="Theoretical Quantiles", ylab="Sample Quantiles")
title(main = "z = 4.687 Q-Q plot of lognormal fit against data")
abline(0,1)
dev.off()

png("z0p0_QQ.png", width = 2000, height = 2000, res = 200)
plot(qlnorm(quants,fit_params_2$estimate['meanlog'], fit_params_2$estimate['sdlog']),quantile(data_2,quants), xlab="Theoretical Quantiles", ylab="Sample Quantiles")
title(main = "z = 0.0 Q-Q plot of lognormal fit against data")
abline(0,1)
dev.off()
















