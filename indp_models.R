library(nlme)
library(dplyr)
library(MCMCglmm)
library(reshape)
library(rjags)
library(MASS)
library(ggplot2)
library(lme4)

setwd("C:\\Users\\dcries\\github\\epi")
nhanes <- read.csv("NHANES_complete.csv")
names(nhanes) <- tolower(names(nhanes))
nhanes$weekend <- 0
nhanes$weekend[nhanes$dow %in% c(1,7)] <- 1
nhanes$first5 <- 0
nhanes$first5[nhanes$rep==6] <- 1
nhanes$first5[nhanes$rep==7] <- 2

nhanes$active <- 1
nhanes$active[nhanes$modvigmin ==0] <- 0

nrep <- (nhanes %>% group_by(id) %>% summarise(n=length(id)))$n
meas7 <- subset(nhanes, id %in% unique(id)[nrep==7]) #individuals with all 7 days

meas7binom <- meas7 %>% group_by(id) %>% summarise(active=sum(active),sex=sex[1],age=age[1],race=race[1])
meas7binom$total <- 7

qqnorm(nhanes$modvigmin[nhanes$modvigmin>0]^(1/4));qqline(nhanes$modvigmin[nhanes$modvigmin>0]^(1/4));
jit <- rnorm(sum(nhanes$modvigmin >0),0,5.5)
qqnorm((nhanes$modvigmin[nhanes$modvigmin>0]+jit)^(1/4));qqline((nhanes$modvigmin[nhanes$modvigmin>0]+jit)^(1/4));

m1 <- lme(I(modvigmin^.25)~age+as.factor(sex)+as.factor(race)+weekend+as.factor(first5),
          random=~1|id,correlation=corAR1(form=~1|id),data=subset(meas7,modvigmin>0),method="ML")

m2 <- lme(I(modvigmin^.25)~age+as.factor(sex)+as.factor(race)+weekend+as.factor(first5),
          random=~1|id,correlation=corARMA(p=2,q=0,form=~1|id),data=subset(meas7,modvigmin>0),method="ML")

m3 <- glmmPQL(active~age+as.factor(sex)+as.factor(race)+weekend+as.factor(first5),
              random=~1|id,correlation=corAR1(form=~1|id),data=meas7,
              family=binomial(link="probit"))

m4 <- glmmPQL(active~age+as.factor(sex)+as.factor(race)+weekend+as.factor(first5),
              random=~1|id,correlation=corARMA(p=2,q=0,form=~1|id),data=meas7,
              family=binomial(link="probit"))

m5 <- glm(active~age+as.factor(sex)+as.factor(race)+weekend+as.factor(first5),
                           data=meas7,
                           family=binomial)

m5b <- glm(cbind(active,total-active)~age+as.factor(sex)+as.factor(race),
          data=meas7binom,
          family=binomial)

m6 <- glmer(active~age+as.factor(sex)+as.factor(race)+weekend+as.factor(first5)+(1|id),
          data=meas7,
          family=binomial(link="probit"))

m6b <- glmer(cbind(active,total-active)~age+as.factor(sex)+as.factor(race)+(1|id),
            data=meas7binom,
            family=binomial(link="probit"))

#--------------------
#simulate data from model
final1 <- matrix(0,nrow=10,ncol=2)
final2 <- matrix(0,ncol=10,nrow=nrow(meas7))
for(j in 1:10){
  simdata <- matrix(0,ncol=length(unique(meas7$id)),nrow=7)
  personeffect <- rnorm(ncol(simdata),0,1.004685)
  # for(i in 1:ncol(simdata)){
  #   simdata[,i] <- as.numeric(arima.sim(n=7,list(ar=-0.05532052),sd=0.5601218))+personeffect[i]
  # }
  
  simdata <- (c(simdata)+predict(m6,newdata=meas7,type="response")+personeffect)
  probs <- pnorm(simdata)
  active <- rbinom(length(probs),1,probs)
  final1[j,] <- table(active)
  
  simdata2 <- matrix(0,ncol=length(unique(meas7$id)),nrow=7)
  personeffect2 <- rnorm(ncol(simdata2),0,0.3891849)
  for(i in 1:ncol(simdata2)){
    simdata2[,i] <- as.numeric(arima.sim(n=7,list(ar=0.1058238),sd=0.3796684))+personeffect2[i]
  }
  
  simdata2 <- c(simdata2)+predict(m1,newdata=meas7)
  simdata2 <- simdata2[active==1]
  final2[1:length(simdata2),j] <- simdata2
  #lines(ecdf(simdata2),col="red")
}

mfinal <- melt(data.frame(final2))

ggplot() + stat_ecdf(data=subset(meas7,modvigmin>0),aes(x=I(modvigmin^(1/4))))+
  stat_ecdf(data=subset(mfinal,value>0),aes(x=value,colour=variable),alpha=0.5)

ks.test(meas7$modvigmin[meas7$modvigmin>0]^(1/4),final2[final2[,1]>1,1])
ks.test((meas7$modvigmin[meas7$modvigmin>0]+rnorm(sum(meas7$modvigmin>0)))^(1/4),final2[final2[,1]>0,1])

adj=(meas7$modvigmin[meas7$modvigmin>0]+rnorm(sum(meas7$modvigmin>0),0,5))^(1/4)
ks.test(na.omit(adj),"pnorm",mean(na.omit(adj)),sd(na.omit(adj)))

#------------------
ym <- melt(meas7[,c("modvigmin","rep","id")],id.vars=c("id","rep"))
yc <- cast(ym,id+variable~rep)

x <- meas7[,c("age","sex","race","weekend","first5")]
X <- model.matrix(~age+as.factor(sex)+as.factor(race)+weekend+as.factor(first5),data=x)

jagsmodel='
model
{
  for (i in 1:N) {
    for(j in 1:K){
      #Y[i,j] <- mu[i,j] + W[i,j]
      mu[i,j] <- inprod(X[(i-1)*K+j,],beta[]) + b[i]
    }
  Y[i,] ~ dmnorm(mu[i,], Omega[,])
  #mu[i] <- inprod(X[i,],beta[]) + b[i]
  b[i] ~ dnorm(0,tau2b)
  #W[i,1:K] ~ dmnorm(muW[], Omega[,])
  }

  tau2e ~ dgamma(0.1,0.1)
  sigma2e <- 1/tau2e

  tau2b ~ dgamma(0.1,0.1)
  sigma2b <- 1/tau2b

  rho ~ dunif(-1,1)

  beta[1:p] ~ dmnorm(b_mu[], b_var[,])
  

  for (i in 1:(K-1)) {
    H[i,i] <- sigma2e 
      for (j in (i+1):K) {
        H[i,j] <- sigma2e*rho^(abs(i-j))
        H[j,i] <- H[i,j] 
    }
  }

  H[K,K] <- sigma2e
  
  
  Omega[1:K,1:K] <- inverse(H[,])
}
'
datajags <- list(
  N      = length(unique(meas7$id)),
  p      = ncol(X),
  K      = 7,
  Y      = as.matrix(yc[,3:9])^(1/4),
  X      = X,
  b_mu  = rep(0,ncol(X)),
  b_var = diag(ncol(X))/100,
  muW   = rep(0,7)
)

jagsModel <- jags.model(textConnection(jagsmodel), data=datajags,n.adapt=1000,n.chains=3)
samples <- coda.samples(jagsModel,c("beta","sigma2e","sigma2b","rho"),n.iter=1000)


jagsmodel2='
model
{
  for (i in 1:N) {
  for(j in 1:K){
    #Y[i,j] <- mu[i,j] + W[i,j]
    probit(mu[i,j]) <- inprod(X[(i-1)*K+j,],beta[]) + b[i]
    Y[i,j] ~ dbern(mu[i,j])
  }
  #Y[i,] ~ dmnorm(mu[i,], Omega[,])
  #mu[i] <- inprod(X[i,],beta[]) + b[i]
  b[i] ~ dnorm(0,tau2b)
  #W[i,1:K] ~ dmnorm(muW[], Omega[,])
  }
  
  #tau2e ~ dgamma(0.1,0.1)
  #sigma2e <- 1/tau2e
  
  tau2b ~ dgamma(0.1,0.1)
  sigma2b <- 1/tau2b
  
  #rho ~ dunif(-1,1)
  
  beta[1:p] ~ dmnorm(b_mu[], b_var[,])
  
  
  # for (i in 1:(K-1)) {
  # H[i,i] <- sigma2e 
  # for (j in (i+1):K) {
  # H[i,j] <- sigma2e*rho^(abs(i-j))
  # H[j,i] <- H[i,j] 
  # }
  # }
  # 
  # H[K,K] <- sigma2e
  # 
  # 
  # Omega[1:K,1:K] <- inverse(H[,])
}
'
yc2 <- yc[,3:9]
yc2[yc2>0] <- 1
datajags2 <- list(
  N      = length(unique(meas7$id)),
  p      = ncol(X),
  K      = 7,
  Y      = yc2,
  X      = X,
  b_mu  = rep(0,ncol(X)),
  b_var = diag(ncol(X))/100,
  muW   = rep(0,7)
)

jagsModel2 <- jags.model(textConnection(jagsmodel2), data=datajags2,n.adapt=5000,n.chains=3)
samples2 <- coda.samples(jagsModel2,c("beta","sigma2b"),n.iter=5000)



#-------------

model <- "model{
  
  # For the ones trick
  C <- 10000
  
  # for every observation
  for(i in 1:N){
      logit(w[i]) <- zeta[i]
      zeta[i] <- gamma0 + gamma1*age[i] + gamma2*gender[i] + b[ind[i],1]
    
      eta[i] <- beta0 + beta1*age[i] + beta2*gender[i] + b[ind[i],2]
    
      logGamma[i] <- log(dnorm(y[i], eta[i], tau2))
    
    # define the total likelihood, where the likelihood is (1 - w) if y < 0.0001 (z = 0) or
    # the likelihood is w * gammalik if y >= 0.0001 (z = 1). So if z = 1, then the first bit must be
    # 0 and the second bit 1. Use 1 - z, which is 0 if y > 0.0001 and 1 if y < 0.0001
     logLik[i] <- (1 - z[i]) * log(1 - w[i]) + z[i] * ( log(w[i]) + logGamma[i] )
    
      Lik[i] <- exp(logLik[i])
    
    # Use the ones trick
      p[i] <- Lik[i] / C
      ones[i] ~ dbern(p[i])

  }
  for(i in 1:m){
    b[i,1:2] ~ dmnorm(mu,Sigma)
  }
  beta0 ~ dnorm(0, 0.0001)
  beta1 ~ dnorm(0, 0.0001)
  beta2 ~ dnorm(0, 0.0001)
  gamma0 ~ dnorm(0, 0.0001)
  gamma1 ~ dnorm(0, 0.0001)
  gamma2 ~ dnorm(0, 0.0001)
  tau2 ~ dgamma(1,1)
  Sigma ~ dwish(D,v)
}"

dat <- list(
  N      = nrow(meas7),
  y      = (meas7$modvigmin)^(1/4),
  z      = meas7$active,
  ones = rep(1,nrow(meas7)),
  age=meas7$age,
  gender=meas7$sex,
  D=diag(2),
  v=3,
  ind=rep(1:length(unique(meas7$id)),each=7),
  mu=rep(0,2),
  m=length(unique(meas7$id))
)

jm <- jags.model(textConnection(model), data=dat,n.adapt=5000,n.chains=3)
s <- coda.samples(jm,c("beta0","beta1","beta2","gamma0","gamma1","gamma2","tau2","Sigma"),n.iter=5000)


#-------------

model2 <- "model{

# For the ones trick
C <- 10000

# for every individual -observation
for(i in 1:N){
logit(w[i]) <- zeta[i]
zeta[i] <- gamma0 + gamma1*age[i] + gamma2*gender[i] + b[ind[i],1]

eta[i] <- beta0 + beta1*age[i] + beta2*gender[i] + b[ind[i],2]

logGamma[i] <- log(dmnorm(y[i], eta[i], tau2))

# define the total likelihood, where the likelihood is (1 - w) if y < 0.0001 (z = 0) or
# the likelihood is w * gammalik if y >= 0.0001 (z = 1). So if z = 1, then the first bit must be
# 0 and the second bit 1. Use 1 - z, which is 0 if y > 0.0001 and 1 if y < 0.0001
logLik[i] <- (1 - z[i]) * log(1 - w[i]) + z[i] * ( log(w[i]) + logGamma[i] )

Lik[i] <- exp(logLik[i])

# Use the ones trick
p[i] <- Lik[i] / C
ones[i] ~ dbern(p[i])

}
for(i in 1:m){
b[i,1:2] ~ dmnorm(mu,Sigma)
}


  for (i in 1:(K-1)) {
    H[i,i] <- sigma2e 
for (j in (i+1):K) {
H[i,j] <- sigma2e*rho^(abs(i-j))
H[j,i] <- H[i,j] 
}
}

H[K,K] <- sigma2e

Omega[1:K,1:K] <- inverse(H[,])

rho ~ dunif(-1,1)
tau2e <- 1/sigma2e
tau2e ~ dgamma(1,1)

beta0 ~ dnorm(0, 0.0001)
beta1 ~ dnorm(0, 0.0001)
beta2 ~ dnorm(0, 0.0001)
gamma0 ~ dnorm(0, 0.0001)
gamma1 ~ dnorm(0, 0.0001)
gamma2 ~ dnorm(0, 0.0001)
tau2 ~ dgamma(1,1)
Sigma ~ dwish(D,v)
}"

dat <- list(
  N      = nrow(meas7),
  y      = (meas7$modvigmin)^(1/4),
  z      = meas7$active,
  ones = rep(1,nrow(meas7)),
  age=meas7$age,
  gender=meas7$sex,
  D=diag(2),
  v=3,
  ind=rep(1:length(unique(meas7$id)),each=7),
  mu=rep(0,2),
  m=length(unique(meas7$id))
)

jm <- jags.model(textConnection(model), data=dat,n.adapt=5000,n.chains=3)
s <- coda.samples(jm,c("beta0","beta1","beta2","gamma0","gamma1","gamma2","tau2","Sigma"),n.iter=5000)

