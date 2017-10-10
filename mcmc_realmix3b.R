library(Rcpp)
library(MASS)
library(mvtnorm)
library(LaplacesDemon)
library(dplyr)
library(label.switching)


setwd("/home/dcries/epi/")
Rcpp::sourceCpp('mcmc_epi_mixture.cpp')
imp1 <- read.csv("nhanes_complete.csv")
load("/home/dcries/stanout.RData")
rmat <- as.matrix(rs)
tstar <- rmat[,31:3367]
#nhanes <- read.csv("NHANES_complete.csv")
names(imp1) <- tolower(names(imp1))
imp1 <- subset(imp1,rep!=7)

nrep <- (imp1 %>% group_by(id) %>% summarise(n=length(id)))$n
meas7 <- subset(imp1, id %in% unique(id)[nrep==6] & (!is.na(imp1$education))) #individuals with all 7 days


tstar2 <- tstar[,complete.cases(meas7[!duplicated(meas7$id),c("waist","glu","tri","bps","bpd","ldl","hdl")]) & (meas7$bpd[!duplicated(meas7$id)] > 0)]
meas7 <- meas7[(!is.na(meas7$waist)) & (!is.na(meas7$bps)) & (!is.na(meas7$bpd)) & (!is.na(meas7$hdl)) & (!is.na(meas7$ldl)) & (!is.na(meas7$glu)) & (!is.na(meas7$tri)) & (!is.na(meas7$education)) & (meas7$bpd >0),] #remove NAs for waist

waist <- meas7$waist[!duplicated(meas7$id)]
lglu <- log(meas7$glu[!duplicated(meas7$id)])
ltri <- log(meas7$tri[!duplicated(meas7$id)])
bps <- (meas7$bps[!duplicated(meas7$id)])
ldl <- (meas7$ldl[!duplicated(meas7$id)])
bpd <- meas7$bpd[!duplicated(meas7$id)]
hdl <- (meas7$hdl[!duplicated(meas7$id)])
MetS <- (cbind(waist,lglu,ltri,bps,ldl,bpd,hdl))

K=3
start <- list(currentbeta=c(10.445,   3.230,   2.033 ,0.1642,3.4081,1.4433,#0.1642,3.4081,1.4433,
                            0.2805, 4.4733, 1.8297,  #log tri
                            #65.736,   1.412,   2.121, #tri
                            18.388,   4.602 ,  1.389 , #bps
                            #.138,4.247,1.387, #log bps
                            -6.3, #ldl
                            #-0.16, #sqrt(ldl)
                            2.3, #bpd
                            -3.3 #hdl
                            #-.02 #log(hdl)
),
currentlambda=matrix(c(rep(101.517,K),rep(4.7228,K),rep(4.9261,K),rep(138.480,K),
                       rep(10.84,K),rep(63,K),rep(3.99,K)),ncol=K,byrow=T),
Sigmadiag=matrix(rep(c(15^2,.16^2,.24^2,36^2,18^2,14^2,16^2),K),ncol=K,byrow=FALSE),
currentzeta=sample(0:(K-1),nrow(MetS),replace=TRUE,rep(1/K,K)),
currentpi=rep(1/K,K),
propcov=diag(15)*0.00001)

#start$currentlambda[,2] <- start$currentlambda[,2] + rnorm(nrow(start$currentlambda),rep(0,nrow(start$currentlambda)),0.15*start$currentlambda[,2])
#start$currentlambda[,1] <- start$currentlambda[,1]*.8
#start$Sigmadiag[,1] <- start$Sigmadiag[,1]*.6

prior <- list(bm=c(7,3,2.11,.16,3,2.11,.12,3,2.11,18,3,2.11,rep(0,3)),
              bcov=diag(15)*c(8,1.5,.8,.08,1.5,.8,3,1.5,.8,5,1.5,.8,rep(100,3))^2,d=8,D=diag(7),
              lm=c(98,4.7,4.73,130,0,0,0),
              lcov=diag(7)*c(17,.1,.6,7,100,100,100)^2,
              a=rep(1,K))

#start$currentbeta <- start$currentbeta * 1.1
start$currentbeta <- c(7.15,1.41,1.78,0.11,1.49,1.1,0.17,1.32,1.60,15.11,2.33,1.19,-7.05,0.44,-2.4)
start$currentlambda <- start$currentlambda * 0.7
start$Sigmadiag <- start$Sigmadiag * 0.7

out3 = mcmc_epi_mixture(MetS,tstar2, start, prior, K,650000,150000,thin=10,.2)
out3$dic
# pmat <- array(0,dim=c(nrow(out1$beta),nrow(MetS),K))
# for(i in 1:K){
#   for(j in 1:nrow(MetS)){
#     pmat[,j,i] <- out1$pmat[j,i,]
#   }
# }
out3$pmat <- NULL

# permutations=label.switching(c("ECR-ITERATIVE-1","ECR-ITERATIVE-2","STEPHENS"),
#                              p=pmat,z=out1$zeta+1,K=K)
# out1=list(out1,permutations)

save(out3,file="/home/dcries/stanout_realmix3b.RData")

# length(unique(out$beta[,1]))/nrow(out$beta)
# diag(out$propcov)

# add 1 to each
# w means.col(0) = beta[3]-beta[0]/(1.0+exp(-beta[1]*(tstar-beta[2])));
# g means.col(1) = beta[7]-beta[4]/(1.0+exp(-beta[5]*(tstar-beta[6])));
# t means.col(2) = beta[11]-beta[8]/(1.0+exp(-beta[9]*(tstar-beta[10])));
# bpsmeans.col(3) =  beta[15]-beta[12]/(1.0+exp(-beta[13]*(tstar-beta[14])));
# ldl means.col(4) = beta[16] + beta[17]*tstar;
# bpd means.col(5) = beta[18] + beta[19]*tstar;
# hdl means.col(6) = beta[20] + beta[21]*tstar;

# plot(out$beta[,1],type="l") #waist decrease
# plot(out$beta[,2],type="l") #waist k
# plot(out$beta[,3],type="l") #waist midpoint
# plot(out$beta[,4],type="l") #glu decrease
# plot(out$beta[,5],type="l") #glu k
# plot(out$beta[,6],type="l") #glu midpoint
# plot(out$beta[,7],type="l") #tri decrease
# plot(out$beta[,8],type="l") #tri k
# plot(out$beta[,9],type="l") #tri mid
# plot(out$beta[,10],type="l") #bps decrease
# plot(out$beta[,11],type="l") #bps k
# plot(out$beta[,12],type="l") #bps mid
# plot(out$beta[,13],type="l") #ldl
# plot(out$beta[,14],type="l") #bpd
# plot(out$beta[,15],type="l") #hdl
# 
# plot(sqrt(out$Sigma[1,1,]),type="l") #w
# plot(sqrt(out$Sigma[2,2,]),type="l") #g
# plot(sqrt(out$Sigma[3,3,]),type="l") #t
# plot(sqrt(out$Sigma[4,4,]),type="l") #bps
# plot(sqrt(out$Sigma[5,5,]),type="l") #ldl
# plot(sqrt(out$Sigma[6,6,]),type="l") #bpd
# plot(sqrt(out$Sigma[7,7,]),type="l") #hdl
