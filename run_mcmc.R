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

start <- list(currentbeta=c(10.445,   3.230,   2.033, 101.517,0.1642,3.4081,1.4433,4.7228,
                      0.2805, 4.4733, 1.8297, 4.9261,  #tri
                      18.388,   4.602 ,  1.389, 137.256 , #bps
                      126,-6.3, #ldl
                      63,2.3, #bpd
                      62,-3.3 #hdl
                     ),
                           currentSigma=diag(7)*c(15^2,.16^2,.54^2,36^2,18^2,14^2,16^2),
                           propcov=diag(22)*0.00001)
prior <- list(bm=c(7,3,2.11,98,.16,3.6,1.4,4.7,.12,4.88,2.11,4.73,18,3,1.3,130,rep(0,6)),
              bcov=diag(22)*c(8,1.5,.4,17,.08,.7,1,3,.1,2,.4,.6,5,1,1,7,rep(100,6))^2,d=8,D=diag(7))

start$currentbeta = c(10.6014768,   4.1067134,   2.0354505, 101.7950530,   0.2064588 ,  3.2948488,
                      1.5767090 ,  4.7471581  , 0.2396806  , 4.1355670  , 1.9258756   ,4.8684419,
                       19.4293800,   4.7369439 ,  1.5558325 ,137.2251982 ,125.6755433  ,-6.5900784,
                       62.4165241 ,  2.7786547  ,62.4926395  ,-3.6013021)

start$propcov = diag(22)*c( 0.0165999241, 0.0488733444, 0.0023339543, 0.0212046561, 0.0004781497, 0.0291257473,
                            0.0073251812, 0.0004953769, 0.0009280714, 0.1321316209, 0.0108031286, 0.0005699109,
                            0.1068854350, 0.0218794029, 0.0014728633, 0.0361485305, 0.0778616124, 0.1142208118,
                            0.0310911264, 0.0501745743, 0.0908247222, 0.0328641870)

out = mcmc_epi(MetS,tstar2, start, prior, 100000,50000)

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