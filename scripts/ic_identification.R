if(!require('dcemriS4')){install.packages('dcemriS4')}

library('dcemriS4')

library(ggplot2)

lbox = 2168.64
grid = 600

#'
#' Load density data at two epochs
#'

dir <- "../data/"
snap <- "weights_002_"  # z = 4.687 

data <- scan(paste0(dir,snap,grid,'.bin'),what=integer())

# convert to number density Mpc^-3
data <- data * (lbox/grid)**-3

arr <- array(data,c(grid,grid,grid))

# import convolution and kernel functions
source('convolve.R')


radius = 42

conv_a <- conv_it(arr,grid,lbox=lbox,side=radius)

radius = 84

conv_b <- conv_it(arr,grid,lbox=lbox,side=radius)

radius = 168

conv_c <- conv_it(arr,grid,lbox=lbox,side=radius)



# Find rms overdensities

conv_a_rms = sqrt(mean(conv_a^2))
conv_a_sigma = sqrt(mean(abs(conv_a-conv_a_rms)^2))

conv_b_rms = sqrt(mean(conv_b^2))
conv_b_sigma = sqrt(mean(abs(conv_b-conv_b_rms)^2))

conv_c_rms = sqrt(mean(conv_c^2))
conv_c_sigma = sqrt(mean(abs(conv_c-conv_c_rms)^2))


# number of overdensities
sum((conv_a < conv_a_rms+(0.002*conv_a_sigma)) & (conv_a > conv_a_rms-(0.002*conv_a_sigma)))

sum((conv_a < conv_a_rms+conv_a_sigma+(0.002*conv_a_sigma)) & (conv_a > conv_a_rms+conv_a_sigma-(0.002*conv_a_sigma)))

sum((conv_a < conv_a_rms+(2*conv_a_sigma)+(0.002*conv_a_sigma)) & (conv_a > conv_a_rms+(2*conv_a_sigma)-(0.002*conv_a_sigma)))



zero_sigma_a = (conv_a < conv_a_rms+(0.002*conv_a_sigma)) & (conv_a > conv_a_rms-(0.002*conv_a_sigma))

zero_sigma_b = (conv_b < conv_b_rms+(0.002*conv_b_sigma)) & (conv_b > conv_b_rms-(0.002*conv_b_sigma))

zero_sigma_c = (conv_c < conv_c_rms+(0.002*conv_c_sigma)) & (conv_c > conv_c_rms-(0.002*conv_c_sigma))

sum(zero_sigma_a & zero_sigma_b & zero_sigma_c)


#sum(conv > (conv_rms+(2*conv_sigma)))


# find a boring patch of the universe
zero_sigma = (conv_a < conv_a_rms+(0.002*conv_a_sigma)) & (conv_a > conv_a_rms-(0.002*conv_a_sigma)) & (conv_b < conv_b_rms+(0.002*conv_b_sigma)) & (conv_b > conv_b_rms-(0.002*conv_b_sigma))

# how many boring regions?
sum(zero_sigma)

# grab the first one, return coordinates in simulation volume
which(zero_sigma,arr.ind=T)[1,] * (lbox/grid)


# more interesting patch of the universe
two_sigma = ((conv_a < conv_a_rms+(2*conv_a_sigma)+(0.002*conv_a_sigma)) & (conv_a > conv_a_rms+(2*conv_a_sigma)-(0.002*conv_a_sigma)))

which(two_sigma,arr.ind=T)[100,] * (lbox/grid)

# exceedingly boring region
mtwo_sigma = ((conv_a < conv_a_rms-(2*conv_a_sigma)+(0.002*conv_a_sigma)) & (conv_a > conv_a_rms-(2*conv_a_sigma)-(0.002*conv_a_sigma)))


which(mtwo_sigma,arr.ind=T)[100,] * (lbox/grid)

# overdensity test with smaller window

radius = 10

conv_t <- conv_it(arr,grid,lbox=lbox,side=radius)

which(conv_t == max(conv_t),arr.ind=T) * (lbox/grid)



