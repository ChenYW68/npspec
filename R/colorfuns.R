
blue2red <- function(ncol,minval=0,maxval=0.7,expon=0.7,minred=0.1){

    rgb( blue = c(seq(1,maxval^(1/expon),length.out = ncol/2), seq(maxval^(1/expon),minval,length.out=ncol/2))^expon,
         green = c(seq(minval,1,length.out=ncol/2),seq(1,minval,length.out=ncol/2))^expon,
         red = c(seq(minred,1,length.out=ncol/2), seq(1,1,length.out =ncol/2) )^expon  )

}
