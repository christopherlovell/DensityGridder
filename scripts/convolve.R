
#if(!require('dcemriS4')){install.packages('dcemriS4')}

#library('dcemriS4')


#file = "weights_022_600"


conv_it <- function(arr,lbox=2168.64,side=42,grid=600){

    # Get kernel template
    kernel <- top_hat_kernel(grid,lbox,radius=side)

    #' dcemris4::convFFT(A, B ,C)
    #' A: the template, i.e. the kernel
    #' B: the target array
    #' C: center of the kernel
    convolution <- dcemriS4::convFFT(kernel,arr,c(grid/2,grid/2,grid/2)) 

    # normalise (for top hat only)
    convolution <- convolution / sum(kernel)
    

    return(convolution)
}


box_kernel <- function(grid,lbox,side){
    
    # initialise kernel array
    template <- array(rep(F,grid^3),c(grid,grid,grid))

    # find kernel size
    kern <- as.integer(side/(lbox/grid))
    # force to quantize inside regular grid
    if(kern %% 2 != 0) kern <- kern + 1

    # kernel limits
    lower <- as.integer(grid/2 - kern/2 + 1)
    upper <- as.integer(grid/2 + kern/2)

    template[lower:upper,lower:upper,lower:upper] <- T

    return(template)
}


top_hat_kernel <- function(grid,lbox,radius){
    print("initialising kernel array")

    template <- array(rep(F,grid^3),c(grid,grid,grid))
   
    grid_radius = radius / (lbox/grid)

    lower <- grid/2 - ceiling(grid_radius)
    upper <- grid/2 + ceiling(grid_radius)

    for(i in lower:upper){
        for(j in lower:upper){
            for(k in lower:upper){
                distance <- ((i - grid/2)^2 + (j - grid/2)^2 + (k - grid/2)^2)^0.5

                if(distance <= grid_radius){
                    template[i,j,k] <- T
                }

            }
        }
    }

    print("kernel complete")

    return(template)
}


#convolution <- conv_it(file)

#'
#'  Find Mean density regions
#'

#median_density <- which(convolution == median(convolution),arr.ind=T)

#ratio = lbox/grid
#median_density * ratio



