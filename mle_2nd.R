
fn <- function(alphaw,tstar,y){
  m <- alphaw[4]-alphaw[1]/(1+exp(-alphaw[2]*(tstar-alphaw[3])))
  return(-sum(dnorm(y,m,alphaw[5],log=TRUE)))
}

start1 <- list(alphag=c(0.1642,3.4081,1.4433,4.7228),
               alphaw=c(10.445,   3.230,   2.033, 101.517 ),
               alphat=c(0.2805, 4.4733, 1.8297, 4.9261 ),
               alphabs=c( 18.388,   4.602 ,  1.389, 137.256 ))

gmle <- optim(c(0.1642,3.4081,1.4433,4.7228,.5),fn,tstar=Tstar$x,y=lglu,control=list(maxit=1000)) #glu
wmle <- optim(c(10.445,   3.230,   2.033, 101.517,.5),fn,tstar=Tstar$x,y=waist,control=list(maxit=1000))
tmle <- optim(c(0.2805, 4.4733, 1.8297, 4.9261,.5),fn,tstar=Tstar$x,y=ltri,control=list(maxit=1000))
bmle <- optim(c( 18.388,   4.602 ,  1.389, 137.256,.5),fn,tstar=Tstar$x,y=bps,control=list(maxit=1000))

