library(ggplot2)
library(MASS)
library(gridExtra)
setwd("C:/Users/dcries/github/epi/")


#for imputed data
load("../../workspace2/stanout_imp1.RData")
samples <- as.matrix(rs)
samples <- samples[,1:30]
load("../../workspace2/stanout_temp.RData")
demo <- read.csv("demographics_imp.csv")
Z <- model.matrix(~age+as.factor(sex)+as.factor(race)+as.factor(education),data=demo)
MetS <- read.csv("MetS_imp.csv")
sigma2y <- .2662 + -.00806*demo[,2]+.0002099*demo[,2]^2+-1.625e-06*demo[,2]^3

nsim <- 1000
index1 <- sample(1:nrow(samples),nsim,T)
index2 <- sample(1:nrow(out$beta),nsim,T)

waist <- matrix(0,ncol=nrow(demo),nrow=nsim)
glu <- matrix(0,ncol=nrow(demo),nrow=nsim)
tri <- matrix(0,ncol=nrow(demo),nrow=nsim)
bps <- matrix(0,ncol=nrow(demo),nrow=nsim)
bpd <- matrix(0,ncol=nrow(demo),nrow=nsim)
ldl <- matrix(0,ncol=nrow(demo),nrow=nsim)
hdl <- matrix(0,ncol=nrow(demo),nrow=nsim)

for(i in 1:nsim){
  gamma <- samples[index1[i],12:22]
  beta <- samples[index1[i],1:11]
  varb1 <- samples[index1[i],27]^2
  varb2 <- samples[index1[i],28]^2
  covb <- samples[index1[i],24]*(sqrt(varb1)*sqrt(varb2))
  covmat <- matrix(c(varb1,covb,covb,varb2),ncol=2,byrow=T)
  b <- mvrnorm(nrow(Z),c(0,0),covmat)
  
  alpha <- out$beta[index2[i],]
  Sigma <- out$Sigma[,,index2[i]]
  
  tstar <- (invlogit(Z%*%gamma+b[,1])*((Z%*%beta+b[,2])^4+6*sigma2y*(Z%*%beta+b[,2])^2))^.25
  random <- mvrnorm(1,rep(0,7),Sigma)
  waist[i,] <- alpha[4]-alpha[1]/(1+exp(-alpha[2]*(tstar-alpha[3]))) + random[1]
  glu[i,] <- alpha[8]-alpha[5]/(1+exp(-alpha[6]*(tstar-alpha[7]))) + random[2]
  tri[i,] <- alpha[12]-alpha[9]/(1+exp(-alpha[10]*(tstar-alpha[11]))) + random[3]
  bps[i,] <- alpha[16]-alpha[13]/(1+exp(-alpha[14]*(tstar-alpha[15]))) + random[4]
  bpd[i,] <- alpha[19] + alpha[20]*tstar + random[6] #+ alpha[22])*tstar^2 
  ldl[i,] <- alpha[17] + alpha[18]*tstar + random[5] #+ alpha[15])*tstar^2
  hdl[i,] <- alpha[21] + alpha[22]*tstar +random[7]
  

}


residwaist <- (MetS$waist-colMeans(waist))/apply(waist,2,sd)
residglu <- (MetS$lglu-colMeans(glu))/apply(glu,2,sd)
residtri <- (MetS$ltri-colMeans(tri))/apply(tri,2,sd)
residbps <- (MetS$bps-colMeans(bps))/apply(bps,2,sd)
residbpd <- (MetS$bpd-colMeans(bpd))/apply(bpd,2,sd)
residldl <- (MetS$ldl-colMeans(ldl))/apply(ldl,2,sd)
residhdl <- (MetS$hdl-colMeans(hdl))/apply(hdl,2,sd)

p1=qplot(y=residwaist,x=MetS$m) + geom_smooth()
p2=qplot(y=residglu,x=MetS$m) + geom_smooth()
p3=qplot(y=residtri,x=MetS$m) + geom_smooth()
p4=qplot(y=residbps,x=MetS$m) + geom_smooth()
p5=qplot(y=residbpd,x=MetS$m) + geom_smooth()
p6=qplot(y=residldl,x=MetS$m) + geom_smooth()
p7=qplot(y=residhdl,x=MetS$m) + geom_smooth()
grid.arrange(p1,p2,p3,p4,p5,p6,p7,nrow=4)

par(mfrow=c(4,2))
qqnorm(residwaist);qqline(residwaist)
qqnorm(residglu);qqline(residglu)
qqnorm(residtri);qqline(residtri)
qqnorm(residbps);qqline(residbps)
qqnorm(residbpd);qqline(residbpd)
qqnorm(residldl);qqline(residldl)
qqnorm(residhdl);qqline(residhdl)
