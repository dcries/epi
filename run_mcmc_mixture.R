library(MASS)
library(mvtnorm)
library(LaplacesDemon)
library(dplyr)

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


tstar2 <- tstar[,complete.cases(meas7[!duplicated(meas7$id),]) & (meas7$bpd[!duplicated(meas7$id)] > 0)]
meas7 <- meas7[(!is.na(meas7$waist)) & (!is.na(meas7$bps)) & (!is.na(meas7$bpd)) & (!is.na(meas7$hdl)) & (!is.na(meas7$ldl)) & (!is.na(meas7$glu)) & (!is.na(meas7$tri)) & (!is.na(meas7$education)) & (meas7$bpd >0),] #remove NAs for waist

waist <- meas7$waist[!duplicated(meas7$id)]
lglu <- log(meas7$glu[!duplicated(meas7$id)])
ltri <- log(meas7$tri[!duplicated(meas7$id)])
bps <- log(meas7$bps[!duplicated(meas7$id)])
ldl <- sqrt(meas7$ldl[!duplicated(meas7$id)])
bpd <- meas7$bpd[!duplicated(meas7$id)]
hdl <- log(meas7$hdl[!duplicated(meas7$id)])
MetS <- (cbind(waist,lglu,ltri,bps,ldl,bpd,hdl))

K=2
start <- list(currentbeta=c(10.445,   3.230,   2.033 ,0.1642,3.4081,1.4433,
                            0.2805, 4.4733, 1.8297,  #tri
                            #18.388,   4.602 ,  1.389, 137.256 , #bps
                            .138,4.247,1.387, #log bps
                            #126,-6.3, #ldl
                            -0.16, #sqrt(ldl)
                            2.3, #bpd
                            #62,-3.3 #hdl
                            -.02 #log(hdl)
),
  currentlambda=matrix(c(rep(101.517,K),rep(4.7228,K),rep(4.9261,K),rep(4.908,K),
                       rep(10.84,K),rep(63,K),rep(3.99,K)),ncol=K,byrow=T),
  Sigmadiag=matrix(c(15^2,.16^2,.54^2,36^2,18^2,14^2,16^2,
                     25^2,.50^2,2^2,49^2,26^2,21^2,24^2),ncol=2,byrow=FALSE),
  currentzeta=sample(0:(K-1),nrow(MetS),replace=TRUE,rep(1/K,K)),
  currentpi=rep(1/K,K),
  propcov=diag(15)*0.00001)

start$currentlambda[,2] <- start$currentlambda[,2] + rnorm(nrow(start$currentlambda),rep(0,nrow(start$currentlambda)),0.15*start$currentlambda[,2])
start$currentlambda[,1] <- start$currentlambda[,1]*.8
start$Sigmadiag[,1] <- start$Sigmadiag[,1]*.6

prior <- list(bm=c(7,3,2.11,.16,3.6,1.4,.12,4.88,2.11,18,3,1.3,rep(0,3)),
              bcov=diag(15)*c(8,1.5,.4,.08,.7,1,3,2,.4,5,1,1,rep(100,3))^2,d=8,D=diag(7),
              lm=c(98,4.7,4.73,130,0,0,0),
              lcov=diag(7)*c(17,.1,.6,7,100,100,100)^2,
              a=rep(1,K))


out = mcmc_epi_mixture(MetS,tstar2, start, prior, K,5000,3000)

length(unique(out$beta[,1]))/nrow(out$beta)
diag(out$propcov)

# add 1 to each
# w means.col(0) = beta[3]-beta[0]/(1.0+exp(-beta[1]*(tstar-beta[2])));
# g means.col(1) = beta[7]-beta[4]/(1.0+exp(-beta[5]*(tstar-beta[6])));
# t means.col(2) = beta[11]-beta[8]/(1.0+exp(-beta[9]*(tstar-beta[10])));
# bpsmeans.col(3) =  beta[15]-beta[12]/(1.0+exp(-beta[13]*(tstar-beta[14])));
# ldl means.col(4) = beta[16] + beta[17]*tstar;
# bpd means.col(5) = beta[18] + beta[19]*tstar;
# hdl means.col(6) = beta[20] + beta[21]*tstar;

plot(out$beta[,1],type="l") #waist decrease
plot(out$beta[,2],type="l") #waist k
plot(out$beta[,3],type="l") #waist midpoint
plot(out$beta[,4],type="l") #waist start
plot(out$beta[,5],type="l") #glu decrease
plot(out$beta[,6],type="l") #glu k
plot(out$beta[,7],type="l") #glu midpoint
plot(out$beta[,8],type="l") #glu start
plot(out$beta[,9],type="l") #tri decrease
plot(out$beta[,10],type="l") #tri k
plot(out$beta[,11],type="l") #tri mid
plot(out$beta[,12],type="l") #tri decrease
plot(out$beta[,13],type="l") #bps decrease
plot(out$beta[,14],type="l") #bps k
plot(out$beta[,15],type="l") #bps mid
plot(out$beta[,16],type="l") #bps start
plot(out$beta[,17],type="l") #ldl
plot(out$beta[,18],type="l") #ldl
plot(out$beta[,19],type="l") #bpd
plot(out$beta[,20],type="l") #bpd
plot(out$beta[,21],type="l") #hdl
plot(out$beta[,22],type="l") #hdl

plot(sqrt(out$Sigma[1,1,]),type="l") #w
plot(sqrt(out$Sigma[2,2,]),type="l") #g
plot(sqrt(out$Sigma[3,3,]),type="l") #t
plot(sqrt(out$Sigma[4,4,]),type="l") #bps
plot(sqrt(out$Sigma[5,5,]),type="l") #ldl
plot(sqrt(out$Sigma[6,6,]),type="l") #bpd
plot(sqrt(out$Sigma[7,7,]),type="l") #hdl

#20000 burn
[1] 0.8163531259 0.1766860175 0.0050873910 0.3554404013 0.0026411716 0.3227163282
[7] 0.0200550120 0.0021758909 0.0007838732 0.2248870768 0.0110636706 0.0005689429
[13] 0.2707886163 0.8292985560 0.0021894336 0.1842908175 0.1085846536 0.0732017781
[19] 0.4759835595 0.1631186421 0.1017801003 0.0477123960

> diag(out$propcov) #200000
[1] 34.716681022  0.902205519  0.048650489  2.978573937  0.017890427  1.167236899  0.097495781
[8]  0.015264485  0.074775603  3.032575776  0.112325053  0.012245692 10.254759688  1.299276466
[15]  0.004121287  8.389515596  6.665046331  1.646001318  1.291927573  0.295204194  1.540120809
[22]  0.372311965

> length(unique(out$beta[,1]))/nrow(out$beta)
[1] 0.1022289

colMeans(out$beta)
[1]  17.9758038   2.5948982   2.2895764 103.6226191   0.2940551   3.0327665   1.2663068   4.8404261
[9]   0.4951950   2.7066057   2.0396992   4.9701160  22.6739272   5.2236897   1.4881128 140.8617348
[17] 124.8884013  -5.4599525  63.8575194   2.0840736  61.5225344  -2.9667019

apply(out$beta,2,sd)
[1] 5.89348400 0.94981856 0.22087268 1.72541986 0.13371881 1.08028256 0.31216528 0.12351548 0.27371365
[10] 1.74156504 0.33579762 0.11066600 3.20153747 1.13983504 0.06417864 2.89562569 2.58126363 1.28280777
[19] 1.13639032 0.54320605 1.24134937 0.61024874