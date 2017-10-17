library(MASS)
library(mvtnorm)
#library(LaplacesDemon)
#library(dplyr)
#library(coda)

source('/home/dcries/epi/calc_mets_prob.R')

load("/home/dcries/stanout_realmix3.RData")
load("/home/dcries/stanout_mix5.RData")
setwd("/home/dcries/epi")
demoimp <- read.csv("demographics_imp.csv")
demo <- read.csv("demographics.csv")

out=out3
outimp = out5$out5

minMVPA <- matrix(0,nrow=61,ncol=6)
minMVPAimp <- matrix(0,nrow=61,ncol=6)
minMVPAquant <- matrix(0,ncol=2*6,nrow=61)
minMVPAquantimp <- matrix(0,ncol=2*6,nrow=61)
for(i in 0:60){
  for(j in 1:6){
    cc <- calc_mets_prob(out,demo,tstar = i^.25,nsim=200,criteria = j)
    imp <- calc_mets_prob(outimp,demoimp,tstar = i^.25,nsim=200,criteria = j)
    minMVPA[i+1,j] <- mean(cc)
    minMVPAimp[i+1,j] <- mean(imp)
    minMVPAquant[i+1,(2*j-1):(2*j)] <- quantile(cc,probs=c(0.025,0.975))
    minMVPAquantimp[i+1,(2*j-1):(2*j)] <- quantile(imp,probs=c(0.025,0.975))
  }
  print(i)
}

outlist <- list(minMVPA=minMVPA,minMVPAimp=minMVPAimp,minMVPAquant=minMVPAquant,minMVPAquantimp=minMVPAquantimp)
save(outlist,file="/home/dcries/prob_mets.RData")