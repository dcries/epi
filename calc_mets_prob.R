#this function simulates metS risk factor data from posterior predictive distribution 
#and then calculates probabiity of having MetS under certain criteria and cutoff
#values that are user supplied

# library(ggplot2)
# library(MASS)
# library(gridExtra)
# setwd("C:/Users/dcries/github/epi/")
# 
# demo <- read.csv("demographics_imp.csv")
# demo <- read.csv("demographics.csv")
# 
# nsim <- 1000

calc_mets_prob <- function(out,demo,tstar,nsim,criteria=3,cutoffs=list(wstm=102,wstf=88,
                                                            glu=110,tri=150,
                                                            bps=130,bpd=85,
                                                            ldl=160,hdl=40)){
  
  index2 <- sample(1:nrow(out$beta),nsim,T)
  
  waist <- matrix(0,ncol=nrow(demo),nrow=nsim)
  glu <- matrix(0,ncol=nrow(demo),nrow=nsim)
  tri <- matrix(0,ncol=nrow(demo),nrow=nsim)
  bps <- matrix(0,ncol=nrow(demo),nrow=nsim)
  bpd <- matrix(0,ncol=nrow(demo),nrow=nsim)
  ldl <- matrix(0,ncol=nrow(demo),nrow=nsim)
  hdl <- matrix(0,ncol=nrow(demo),nrow=nsim)
  
  K <- dim(out$lambda)[3]
  for(i in 1:nsim){
    
    zeta <- out$zeta[index2[i],] +1
    alpha <- out$beta[index2[i],]
    lambda <- as.matrix(out$lambda[index2[i],,])
    sds <- as.matrix(out$sds[index2[i],,])
    cormat <- as.matrix(out$cormat[index2[i],,])
    #Sigma <- out$Sigma[,,index2[i]]
    Sigma <- array(0,dim=c(7,7,K))
    for(j in 1:K){
      diag(Sigma[,,j] ) <- sds[,j]^2
      count = 1
      for(jj in 1:6){
        for(jjj in ((jj+1):7)){
          Sigma[jj,jjj,j] <- cormat[count,j]*sqrt(Sigma[jj,jj,j]*Sigma[jjj,jjj,j])
          Sigma[jjj,jj,j] <- Sigma[jj,jjj,j]
          count=count+1
        }
      }
    }
    
    
    random <- array(0,dim=c(nrow(demo),7,K))
    random <- matrix(0,nrow=7,ncol=nrow(demo))
    # !!!!!!!!!! this isn't correct yet
    #for(j in 1:K){
    for(j in 1:nrow(demo)){
      
      #random[,,j] <- mvrnorm(nrow(Z),rep(0,7),Sigma[,,j])
      random[,j] <- mvrnorm(1,rep(0,7),Sigma[,,zeta[j]])
      
    }
    # !!!!!!!!!!!!
    
    waist[i,] <- lambda[1,zeta]-alpha[1]/(1+exp(-alpha[2]*(tstar-alpha[3]))) + random[1,]
    glu[i,] <- lambda[2,zeta]-alpha[4]/(1+exp(-alpha[5]*(tstar-alpha[6]))) + random[2,]
    tri[i,] <- lambda[3,zeta]-alpha[7]/(1+exp(-alpha[8]*(tstar-alpha[9]))) + random[3,]
    bps[i,] <- lambda[4,zeta]-alpha[10]/(1+exp(-alpha[11]*(tstar-alpha[12]))) + random[4,]
    bpd[i,] <- lambda[6,zeta] + alpha[14]*tstar + random[6,] #+ alpha[22])*tstar^2
    ldl[i,] <- lambda[5,zeta] + alpha[13]*tstar + random[5,] #+ alpha[15])*tstar^2
    hdl[i,] <- lambda[7,zeta] + alpha[15]*tstar +random[7,]
    
    #if(i%%100==0){print(i)}
    
  }
  
  overwst <- waist > cutoffs$wstm
  overwst[,demo$sex==2] <- waist[,demo$sex==2] > cutoffs$wstf
  overglu <- glu > log(cutoffs$glu)
  overtri <- tri > log(cutoffs$tri)
  overbps <- bps > cutoffs$bps
  overbpd <- bpd > cutoffs$bpd
  overldl <- ldl > cutoffs$ldl
  overhdl <- hdl < cutoffs$hdl
  
  metsprob <- rep(0,nsim)
  waistprob <- rep(0,nsim)
  gluprob <- rep(0,nsim)
  triprob <- rep(0,nsim)
  bpsprob <- rep(0,nsim)
  ldlprob <- rep(0,nsim)
  bpdprob <- rep(0,nsim)
  hdlprob <- rep(0,nsim)
  
  for(i in 1:nsim){
    tp <- overwst[i,] + overglu[i,] + overtri[i,] + overbps[i,]+overbpd[i,]+overldl[i,]+overhdl[i,]
    metsprob[i] <- sum(tp >= criteria)/nrow(demo)
    waistprob[i] <- sum(overwst[i,])/nrow(demo)
    gluprob[i] <- sum(overglu[i,])/nrow(demo)
    triprob[i] <- sum(overtri[i,])/nrow(demo)
    bpsprob[i] <- sum(overbps[i,])/nrow(demo)
    ldlprob[i] <- sum(overldl[i,])/nrow(demo)
    bpdprob[i] <- sum(overbpd[i,])/nrow(demo)
    hdlprob[i] <- sum(overhdl[i,])/nrow(demo)
    
  }
  
  
  return(list(metsprob=metsprob,indprob=cbind(waistprob,gluprob,triprob,bpsprob,ldlprob,bpdprob,hdlprob)))
}
