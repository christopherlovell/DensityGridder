if(!require('dcemriS4')){install.packages('dcemriS4')}

library('dcemriS4')

library(ggplot2)

library(latex2exp)

lbox = 2168.64
grid = 600

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

arr_1 <- array(data_1,c(grid,grid,grid))
arr_2 <- array(data_2,c(grid,grid,grid))



#'
#' Find overdensities
#'


# Sorting the array is *very* slow due to its size
# To avoid this, just use a fraction of the maximum array value to find the N largest values 
# If N is greater than the n largest you want, simply subset by n, otherwise reduce the fraction

top_n <- function(arr,n=10000){
    N <- 0
    frac <- 1
    while(N < n){
        frac <- frac + 1
        print(frac)
        N_arr <-  which(arr > max(arr)/frac, arr.ind=T)
        N <- nrow(N_arr)
    }

    return(N_arr[sort(arr[N_arr],decreasing=T,index.return=T)$ix[1:n],])
}


n = 1000000
n_arr_1 <- top_n(arr_1,n) # coordinates of n highest density regions at z = 4.687

#arr_1[n_arr_1] # value of highest density regions, z = 4.687
#arr_1[rbind(n_arr_1[3,])]  # subset
#arr_2[n_arr_1]  # density at z = 0 of high density regions at z = 4.687

# limits for plot
lim_max <- max(max(arr_2[n_arr_1]),max(arr_1[n_arr_1]))
lim_min <- min(min(arr_2[n_arr_1]),min(arr_1[n_arr_1]))


if(lim_min == 0){
    lim_min <- lim_min + 0.1
}


png('dens_1.png')
plot(arr_1[n_arr_1],arr_2[n_arr_1],pch=16, col=rgb(0, 0, 1, 0.1), xlim=c(lim_min,lim_max),ylim=c(lim_min,lim_max),log="xy",xlab=TeX('$M_{\\degree} h^{-1} Mpc^{-3} (z = 4.687)$'),ylab=TeX('$M_{\\degree} h^{-1} Mpc^{-3} (z = 0)$'), main=paste0('Top ',n,' density regions at z=4.687\n(no convolution, grid size = ',lbox/grid,' Mpc)'))
dev.off()

n_arr_2 <- top_n(arr_2)

# limits for plot
lim_max <- max(max(arr_2[n_arr_2]),max(arr_1[n_arr_2]))
lim_min <- min(min(arr_2[n_arr_2]),min(arr_1[n_arr_2]))

png('dens_2.png')
plot(arr_2[n_arr_2],arr_1[n_arr_2], pch=16, col=rgb(0, 0, 1, 0.1),xlim=c(lim_min,lim_max),ylim=c(lim_min,lim_max),log="xy",xlab=TeX('$M_{\\degree} h^{-1} Mpc^{-3} (z = 0)$'),ylab=TeX('$M_{\\degree} h^{-1} Mpc^{-3} (z = 4.687)$'),main=paste0('Top ',n,' density regions at z=0\n(no convolution, grid size = ',lbox/grid,' Mpc)'))
dev.off()


#'
#' Alternative plotting
#'
 
#smoothScatter(arr_1[n_arr_1],arr_2[n_arr_1], xlim=c(lim_min,lim_max),ylim=c(lim_min,lim_max),log="xy",xlab=expression(rho(4.687)),ylab=expression(rho(0)),main=paste0('Top ',n,' density regions at z=4.687\n(no convolution, grid size = ',lbox/grid,' Mpc)'))
#
#
#library(ggplot2)
#
#dat <- data.frame(arr_1[n_arr_1],arr_2[n_arr_1])
#colnames(dat) <- c('arr_1','arr_2')
#
#ggplot(dat,aes(x=arr_1,y=arr_2))+
#  scale_x_log10(limits=c(lim_min,lim_max))+scale_y_log10(limits=c(lim_min,lim_max))+
#  stat_density2d(aes(alpha=..level..), geom="polygon") +
#  scale_alpha_continuous()+
##  geom_point(colour="red",alpha=0.02)+
#  theme_bw()
#


#'
#' Plot Convolved function
#'

source('convolve.R')  # import conv_it function

sidelength = 2168.64/600 * 3

conv_1 <- conv_it(arr_1,grid=600,lbox=lbox,side=sidelength)
conv_2 <- conv_it(arr_2,grid=600,lbox=lbox,side=sidelength)

n = 1000000
n_arr_1 <- top_n(conv_1,n)

lim_max <- max(max(conv_2[n_arr_1]),max(conv_1[n_arr_1]))
lim_min <- min(min(conv_2[n_arr_1]),min(conv_1[n_arr_1]))

png('dens_conv_1.png')
plot(conv_1[n_arr_1],conv_2[n_arr_1],pch=16, col=rgb(0, 0, 1, 0.1), xlim=c(lim_min,lim_max),ylim=c(lim_min,lim_max),log="xy",xlab=TeX('$M_{\\degree} h^{-1} Mpc^{-3} (z = 4.687)$'),ylab=TeX('$M_{\\degree} h^{-1} Mpc^{-3} (z = 0)$'),main=paste0('Top ',n,' density regions at z=4.687\n(top hat convolution, radius ',sidelength,' Mpc, grid size = ',lbox/grid,' Mpc)'))
dev.off()

n_arr_2 <- top_n(conv_2,n)

lim_max <- max(max(conv_2[n_arr_2]),max(conv_1[n_arr_2]))
lim_min <- min(min(conv_2[n_arr_2]),min(conv_1[n_arr_2]))

png('dens_conv_2.png')
plot(conv_2[n_arr_2],conv_1[n_arr_2], pch=16, col=rgb(0, 0, 1, 0.1),xlim=c(lim_min,lim_max),ylim=c(lim_min,lim_max),log="xy",xlab=TeX('$M_{\\degree} h^{-1} Mpc^{-3} (z = 0)$'), ylab=TeX('$M_{\\degree} h^{-1} Mpc^{-3} (z = 4.687)$'),main=paste0('Top ',n,' density regions at z=0\n(top hat convolution, radius ',sidelength,' Mpc, grid size = ',lbox/grid,' Mpc)'))
dev.off()

#'
#' Find high density regions at both epochs
#'



n = 10000
n_arr_1 <- top_n(arr_1,n)
n_arr_2 <- top_n(arr_2,n)


#matches <- apply(n_arr_2,1,function(x) sum(rowSums(x == n_arr_1) == 3))


dat <- data.frame(arr_1[n_arr_1],arr_2[n_arr_1],arr_2[n_arr_2],arr_1[n_arr_2])
colnames(dat) <- c('tz4_z4','tz4_z0','tz0_z0','tz0_z4')


lim_max <- max(dat)
lim_min <- min(dat)

png('dens_compare.png', width = 2000, height = 2000, res = 200)

p <- ggplot(dat)
p <- p + ggtitle(paste0('Top ',n,' density regions at two redshifts\n(no convolution, grid size = ',lbox/grid,' Mpc)'))
p <- p + scale_x_log10(limits=c(lim_min,lim_max))+scale_y_log10(limits=c(lim_min,lim_max))
p <- p + xlab(TeX('$M_{\\degree} h^{-1} Mpc^{-3} (z = 4.687)$')) + ylab(TeX('$M_{\\degree} h^{-1} Mpc^{-3} (z = 0)$'))

p <- p + geom_point(aes(x=tz0_z4,y=tz0_z0,colour='blue'),alpha=0.1)
p <- p + geom_point(aes(x=tz4_z4,y=tz4_z0,colour='red'),alpha=0.1)
p <- p + scale_colour_manual(values =c('blue'='blue','red'='red'), labels = c(paste0('Top ',n,' density at z=0'),paste0('Top ',n,' density at z=4.687')))
p

dev.off()


# same for convoluted arrays

n = 10000
n_arr_1 <- top_n(conv_1,n)
n_arr_2 <- top_n(conv_2,n)

dat <- data.frame(conv_1[n_arr_1],conv_2[n_arr_1],conv_2[n_arr_2],conv_1[n_arr_2])
colnames(dat) <- c('tz4_z4','tz4_z0','tz0_z0','tz0_z4')

lim_max <- max(dat)
lim_min <- min(dat)

png('dens_conv_compare.png', width = 2000, height = 2000, res = 200)

p <- ggplot(dat)
p <- p + ggtitle(paste0('Top ',n,' density regions at two redshifts\n(top hat convolution, radius ',sidelength,' Mpc, grid size = ',lbox/grid,' Mpc)'))
p <- p + scale_x_log10(limits=c(lim_min,lim_max))+scale_y_log10(limits=c(lim_min,lim_max))
p <- p + xlab(TeX('$M_{\\degree} h^{-1} Mpc^{-3} (z = 4.687)$')) + ylab(TeX('$M_{\\degree} h^{-1} Mpc^{-3} (z = 0)$'))

p <- p + geom_point(aes(x=tz0_z4,y=tz0_z0,colour='blue'),alpha=0.1)
p <- p + geom_point(aes(x=tz4_z4,y=tz4_z0,colour='red'),alpha=0.1)
p <- p + scale_colour_manual(values =c('blue'='blue','red'='red'), labels = c(paste0('Top ',n,' density at z=0'),paste0('Top ',n,' density at z=4.687')))
p

dev.off()

