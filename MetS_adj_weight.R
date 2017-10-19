#function to perform weight adjustment for MetS risk factors, method is from JASA 1996 paper
# A Semiparametric Transformation approach to estimating usual daily intake dist

MetS_adj_weight <- function(MetS,weights,ncdf=5000){
  fy <- matrix(0,ncol=ncol(MetS),nrow=ncdf)
  fy2 <- matrix(0,ncol=ncol(MetS),nrow=ncdf)
  
  xwaist <- seq(from=min(MetS[,1]),to=max(MetS[,1]),length.out = ncdf)
  xglu <- seq(from=min(MetS[,2]),to=max(MetS[,2]),length.out = ncdf)
  xtri <- seq(from=min(MetS[,3]),to=max(MetS[,3]),length.out = ncdf)
  xbps <- seq(from=min(MetS[,4]),to=max(MetS[,4]),length.out = ncdf)
  xldl <- seq(from=min(MetS[,5]),to=max(MetS[,5]),length.out = ncdf)
  xbpd <- seq(from=min(MetS[,6]),to=max(MetS[,6]),length.out = ncdf)
  xhdl <- seq(from=min(MetS[,7]),to=max(MetS[,7]),length.out = ncdf)
  
  xfull <- cbind(xwaist,xglu,xtri,xbps,xldl,xbpd,xhdl)
  for(i in 1:ncdf){
    for(j in 1:7){
      fy[i,j]=sum(weights*(MetS[,j] <= (xfull[i,j])))
      fy2[i,j]=sum((MetS[,j] <= (xfull[i,j])))
      
    }
  }
  
  MetSadj <- matrix(0,nrow=nrow(MetS),ncol=ncol(MetS))
  for(i in 1:7){
    empcdf <- fy[,i]/nrow(MetS)
    probs <- (1/(nrow(MetS))*(rank(MetS[,i])-0.5))
    MetSadj[,i] <- xfull[sapply(probs,function(x) which.min(abs(x-empcdf))),i]
  }
  return(MetSadj)
}


