library(MASS)
library(mvtnorm)
#library(LaplacesDemon)
#library(dplyr)
#library(coda)

source('C:/Users/dcries/github/epi/calc_mets_prob.R')

load("C:/Users/dcries/workspace/stanout_realmix3.RData")
load("C:/Users/dcries/workspace/stanout_mix5.RData")
setwd("C:/Users/dcries/github/epi")
demoimp <- read.csv("demographics_imp.csv")
demo <- read.csv("demographics.csv")

out=out3
outimp = out5$out5

minMVPA <- matrix(0,nrow=61,ncol=6)
minMVPAimp <- matrix(0,nrow=61,ncol=6)
minMVPAquant <- matrix(0,ncol=2*6,nrow=61)
minMVPAquantimp <- matrix(0,ncol=2*6,nrow=61)

indprob <- matrix(0,nrow=61,ncol=7)
indprobimp <- matrix(0,nrow=61,ncol=7)
indprobquant <- matrix(0,nrow=61,ncol=7*2)
indprobquantimp <- matrix(0,nrow=61,ncol=7*2)
for(i in 0:60){
  for(j in 1:6){
    cc <- calc_mets_prob(out,demo,tstar = i^.25,nsim=200,criteria = j)
    imp <- calc_mets_prob(outimp,demoimp,tstar = i^.25,nsim=200,criteria = j)
    minMVPA[i+1,j] <- mean(cc$metsprob)
    minMVPAimp[i+1,j] <- mean(imp$metsprob)
    minMVPAquant[i+1,(2*j-1):(2*j)] <- quantile(cc$metsprob,probs=c(0.025,0.975))
    minMVPAquantimp[i+1,(2*j-1):(2*j)] <- quantile(imp$metsprob,probs=c(0.025,0.975))
  }
  indprob[i+1,] <- colMeans(cc$indprob)
  indprobimp[i+1,] <- colMeans(imp$indprob)
  indprobquant[i+1,] <- c(apply(cc$indprob,2,quantile,probs=c(0.025,0.975)))
  indprobquantimp[i+1,] <- c(apply(imp$indprob,2,quantile,probs=c(0.025,0.975)))
  
  print(i)
}

outlist <- list(minMVPA=minMVPA,minMVPAimp=minMVPAimp,minMVPAquant=minMVPAquant,
                minMVPAquantimp=minMVPAquantimp,indprob=indprob,indprobimp=indprobimp,
                indprobquant=indprobquant,indprobquantimp=indprobquantimp)
save(outlist,file="C:/Users/dcries/workspace/prob_mets.RData")