library(Rcpp)
library(MASS)
library(mvtnorm)
library(LaplacesDemon)
library(dplyr)
library(label.switching)


setwd("C:/Users/dcries/github/epi")
Rcpp::sourceCpp('mcmc_epi_mixture_w.cpp')
source('C:/Users/dcries/github/epi/MetS_adj_weight.R')
imp1 <- read.csv("NHANES_accel_imp1.csv")
load("../../workspace/stanout_imp1.RData")
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
bps <- (meas7$bps[!duplicated(meas7$id)])
ldl <- (meas7$ldl[!duplicated(meas7$id)])
bpd <- meas7$bpd[!duplicated(meas7$id)]
hdl <- (meas7$hdl[!duplicated(meas7$id)])
MetS <- (cbind(waist,lglu,ltri,bps,ldl,bpd,hdl))
MetSadj <- MetS_adj_weight(MetS)

K=5
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
                       rep(10.84^2,K),rep(63,K),rep(exp(3.99),K)),ncol=K,byrow=T),
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

prior$bcov <- prior$bcov#*0.1

weights <- rep(1,nrow(MetS))#(meas7$smplwt[!duplicated(meas7$id)]/sum(meas7$smplwt[!duplicated(meas7$id)]))*length(meas7$smplwt[!duplicated(meas7$id)])
out5w = mcmc_epi_mixture_w(MetSadj,tstar2, start, prior, weights, K,10000,5000,thin=1,0.15)

save(out5w,file="../../workspace/stanout_mix5w2.RData")
