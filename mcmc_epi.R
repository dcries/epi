library(MASS)
library(mvtnorm)
library(LaplacesDemon)

setwd("C:/Users/dcries/github/epi/")
imp1 <- read.csv("NHANES_accel_imp1.csv")
load("C:/Users/dcries/workspace2/stanout_imp1.RData")
rmat <- as.matrix(rs)
tstar <- rmat[,31:7903]
#nhanes <- read.csv("NHANES_complete.csv")
names(imp1) <- tolower(names(imp1))
imp1 <- subset(imp1,rep!=7)

nrep <- (imp1 %>% group_by(id) %>% summarise(n=length(id)))$n
meas7 <- subset(imp1, id %in% unique(id)[nrep==6] & (!is.na(imp1$education))) #individuals with all 7 days

#meas7$Tstar <- rep(colMeans(r[,21:7938]),each=6) # is from 
#meas7$std <- apply(r[,21:7938],2,sd)

tstar2 <- tstar[,complete.cases(meas7[!duplicated(meas7$id),])]
meas7 <- meas7[(!is.na(meas7$waist)) & (!is.na(meas7$bps)) & (!is.na(meas7$bpd)) & (!is.na(meas7$hdl)) & (!is.na(meas7$ldl)) & (!is.na(meas7$glu)) & (!is.na(meas7$tri)) & (!is.na(meas7$education)),] #remove NAs for waist

waist <- meas7$waist[!duplicated(meas7$id)]
lglu <- log(meas7$glu[!duplicated(meas7$id)])
ltri <- log(meas7$tri[!duplicated(meas7$id)])
bps <- (meas7$bps[!duplicated(meas7$id)])
ldl <- meas7$ldl[!duplicated(meas7$id)]
bpd <- meas7$bpd[!duplicated(meas7$id)]
hdl <- meas7$hdl[!duplicated(meas7$id)]
MetS <- (cbind(waist,lglu,ltri,bps,ldl,bpd,hdl))

#need to figure out what beta is
logl_b <- function(y,tstar,beta,Sigma,priormean,priorcov){
  #waist,glu,tri,bps,bpd,ldl,hdl
  n <- nrow(y)
  ll <- 0
  m <- matrix(c(beta[4]-beta[1]/(1+exp(-beta[2]*(tstar-beta[3]))),
          beta[8]-beta[5]/(1+exp(-beta[6]*(tstar-beta[7]))),
          beta[12]-beta[9]/(1+exp(-beta[10]*(tstar-beta[11]))),
          beta[18]-beta[15]/(1+exp(-beta[16]*(tstar-beta[17]))),
          beta[19] + beta[20]*tstar,
          beta[13] + beta[14]*tstar,
          beta[21] + beta[22]*tstar),nrow=7,byrow=T)
  
  for(i in 1:n){
    ll <- ll + dmvnorm(y[i,],m[,i],Sigma,log=TRUE)
  }
  return(ll+dmvnorm(beta,priormean,priorcov,log=TRUE)) #xxxxxxxxxxx

}


calc_meanb <- function(tstar,beta){
  m <- matrix(c(beta[4]-beta[1]/(1+exp(-beta[2]*(tstar-beta[3]))),
                beta[8]-beta[5]/(1+exp(-beta[6]*(tstar-beta[7]))),
                beta[12]-beta[9]/(1+exp(-beta[10]*(tstar-beta[11]))),
                beta[18]-beta[15]/(1+exp(-beta[16]*(tstar-beta[17]))),
                beta[19] + beta[20]*tstar,
                beta[13] + beta[14]*tstar,
                beta[21] + beta[22]*tstar),ncol=7,byrow=FALSE)
  return(m)
}

mcmc_epi <- function(y,tstar,start,prior,nsim=1000){
  #y is an nx7 matrix
  #tstar is an Xxn matrix
  #size of beta is hard coded !!!!!!!
  n <- nrow(y) #nx7 matrix
  m <- nrow(tstar)
  beta <- matrix(0,ncol=22,nrow=nsim)
  Sigma <- array(0,dim=c(ncol(y),ncol(y),nsim))
  tstar <- tstar[sample(1:m,size=nsim,replace=TRUE),]
  bm <- prior$bm
  bcov <- prior$bcov
  Sd <- prior$d
  SD <- prior$D
  currentbeta <- start$beta
  currentSigma <- start$Sigma
  
  propcov <- diag(length(currentbeta))*0.0001
  for(i in 1:nsim){
    #``update tstar"
    propb <- mvrnorm(1,mu=currentbeta,Sigma=propcov)
    logr <- logl_b(y,tstar[i,],propb,currentSigma,bm,bcov) - logl_b(y,tstar[i,],currentbeta,currentSigma,bm,bcov)
    if(log(runif(1))<logr){
      currentbeta <- propb
    }
    
    if((i < 5000) & (i > 20) & (i%%20==0)){
      propcov <- cov(beta[(1:i-1),])
    }
    currentmeanb <- calc_meanb(tstar[i,],currentbeta)
    
    currentSigma <- rinvwishart(n+Sd,t(y-currentmeanb)%*%(y-currentmeanb)+SD)
    
    beta[i,] <- currentbeta
    Sigma[,,i] <- currentSigma
    if(i%%1000) {print(i)}
  }
  return(list(beta=beta,Sigma=Sigma,propcov=propcov))
}

start <- list(beta=c(10.445,   3.230,   2.033, 101.517,0.1642,3.4081,1.4433,4.7228,
                     0.2805, 4.4733, 1.8297, 4.9261, 123,-4,
                     18.388,   4.602 ,  1.389, 137.256 , 60,3.7,63,-3),
              Sigma=diag(7)*c(15^2,.16^2,.54^2,18^2,36^2,14^2,16^2))
prior <- list(bm=c(rep(0,22)),bcov=diag(22)*1000,d=8,D=diag(7))

out <- mcmc_epi(MetS,tstar2,start,prior,10000)
