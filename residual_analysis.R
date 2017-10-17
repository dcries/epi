#this script calculates standardized residuals for the regression model
#according to the method in my dissertation-simulate new usual values
#according to taylor series appx then simualate new observations according
#to mean fcn
#plots of standardized residuals also given

library(ggplot2)
library(MASS)
library(gridExtra)
setwd("C:/Users/dcries/github/epi/")

#for imputed data
load("../../workspace/stanout_imp1.RData")
samples <- as.matrix(rs)
samples <- samples[,1:30]

#load("../../workspace2/stanout_temp2.RData")
#out$beta <- out$beta[50000:100000,]
#out$Sigma <- out$Sigma[,,50000:100000]
load("../../workspace/stanout_mix5.RData")


demo <- read.csv("demographics_imp.csv")
demo <- read.csv("demographics.csv")

Z <- model.matrix(~age+as.factor(sex)+as.factor(race)+as.factor(education),data=demo)

MetS <- read.csv("MetS_imp.csv")
MetS <- read.csv("MetS.csv")

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

K <- dim(out$lambda)[3]
for(i in 1:nsim){
  gamma <- samples[index1[i],12:22]
  beta <- samples[index1[i],1:11]
  varb1 <- samples[index1[i],27]^2
  varb2 <- samples[index1[i],28]^2
  covb <- samples[index1[i],24]*(sqrt(varb1)*sqrt(varb2))
  covmat <- matrix(c(varb1,covb,covb,varb2),ncol=2,byrow=T)
  b <- mvrnorm(nrow(Z),c(0,0),covmat)
  
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

  
  tstar <- (invlogit(Z%*%gamma+b[,1])*((Z%*%beta+b[,2])^4+6*sigma2y*(Z%*%beta+b[,2])^2))^.25
  random <- array(0,dim=c(nrow(Z),7,K))
  random <- matrix(0,nrow=7,ncol=nrow(Z))
  # !!!!!!!!!! this isn't correct yet
  #for(j in 1:K){
  for(j in 1:nrow(Z)){
    
    #random[,,j] <- mvrnorm(nrow(Z),rep(0,7),Sigma[,,j])
    random[,j] <- mvrnorm(1,rep(0,7),Sigma[,,zeta[j]])
    
  }
  # !!!!!!!!!!!!
  
  # waist[i,] <- alpha[4]-alpha[1]/(1+exp(-alpha[2]*(tstar-alpha[3]))) + random[1]
  # glu[i,] <- alpha[8]-alpha[5]/(1+exp(-alpha[6]*(tstar-alpha[7]))) + random[2]
  # tri[i,] <- alpha[12]-alpha[9]/(1+exp(-alpha[10]*(tstar-alpha[11]))) + random[3]
  # bps[i,] <- alpha[16]-alpha[13]/(1+exp(-alpha[14]*(tstar-alpha[15]))) + random[4]
  # bpd[i,] <- alpha[19] + alpha[20]*tstar + random[6] #+ alpha[22])*tstar^2
  # ldl[i,] <- alpha[17] + alpha[18]*tstar + random[5] #+ alpha[15])*tstar^2
  # hdl[i,] <- alpha[21] + alpha[22]*tstar +random[7]
  waist[i,] <- lambda[1,zeta]-alpha[1]/(1+exp(-alpha[2]*(tstar-alpha[3]))) + random[1,]
  glu[i,] <- lambda[2,zeta]-alpha[4]/(1+exp(-alpha[5]*(tstar-alpha[6]))) + random[2,]
  tri[i,] <- lambda[3,zeta]-alpha[7]/(1+exp(-alpha[8]*(tstar-alpha[9]))) + random[3,]
  bps[i,] <- lambda[4,zeta]-alpha[10]/(1+exp(-alpha[11]*(tstar-alpha[12]))) + random[4,]
  bpd[i,] <- lambda[6,zeta] + alpha[14]*tstar + random[6,] #+ alpha[22])*tstar^2
  ldl[i,] <- lambda[5,zeta] + alpha[13]*tstar + random[5,] #+ alpha[15])*tstar^2
  hdl[i,] <- lambda[7,zeta] + alpha[15]*tstar +random[7,]
  
  if(i%%100==0){print(i)}

}


residwaist <- (MetS$waist-colMeans(waist))/apply(waist,2,sd)
residglu <- (MetS$lglu-colMeans(glu))/apply(glu,2,sd)
residtri <- (MetS$ltri-colMeans(tri))/apply(tri,2,sd)
residbps <- (MetS$bps-colMeans(bps))/apply(bps,2,sd)
residbpd <- (MetS$bpd-colMeans(bpd))/apply(bpd,2,sd)
residldl <- (MetS$ldl-colMeans(ldl))/apply(ldl,2,sd)
residhdl <- (MetS$hdl-colMeans(hdl))/apply(hdl,2,sd)

p1=qplot(y=residwaist,x=MetS$X) + ggtitle("Waist Circumference") + geom_hline(yintercept = 0) + ylab("Standardized Residual") + xlab("Usual Minutes in MVPA")  + theme_bw()#+ geom_smooth()
p2=qplot(y=residglu,x=MetS$X) + ggtitle("log Glucose") + geom_hline(yintercept = 0) + ylab("Standardized Residual") + xlab("Usual Minutes in MVPA")  + theme_bw()#+ geom_smooth()#+ geom_smooth()
p3=qplot(y=residtri,x=MetS$X) + ggtitle("log Triglycerides") + geom_hline(yintercept = 0) + ylab("Standardized Residual") + xlab("Usual Minutes in MVPA")  + theme_bw()#+ geom_smooth()#+ geom_smooth()
p4=qplot(y=residbps,x=MetS$X) + ggtitle("Systolic Blood Pressure") + geom_hline(yintercept = 0) + ylab("Standardized Residual") + xlab("Usual Minutes in MVPA")  + theme_bw()#+ geom_smooth()#+ geom_smooth()
p5=qplot(y=residbpd,x=MetS$X) + ggtitle("Diastolic Blood Pressure") + geom_hline(yintercept = 0) + ylab("Standardized Residual") + xlab("Usual Minutes in MVPA")  + theme_bw()#+ geom_smooth()#+ geom_smooth()
p6=qplot(y=residldl,x=MetS$X) + ggtitle("LDL") + geom_hline(yintercept = 0) + ylab("Standardized Residual") + xlab("Usual Minutes in MVPA")  + theme_bw()#+ geom_smooth()#+ geom_smooth()
p7=qplot(y=residhdl,x=MetS$X) + ggtitle("HDL") + geom_hline(yintercept = 0) + ylab("Standardized Residual") + xlab("Usual Minutes in MVPA")  + theme_bw()#+ geom_smooth()#+ geom_smooth()
grid.arrange(p1,p2,p3,p4,p5,p6,p7,nrow=4)

par(mfrow=c(4,2))
qqnorm(residwaist,main="Waist Circumference");qqline(residwaist)
qqnorm(residglu,main="log Glucose");qqline(residglu)
qqnorm(residtri,main="log Triglyceride");qqline(residtri)
qqnorm(residbps,main="Systolic Blood Pressure");qqline(residbps)
qqnorm(residbpd,main="Diastolic Blood Pressure");qqline(residbpd)
qqnorm(residldl,main="LDL");qqline(residldl)
qqnorm(residhdl,main="HDL");qqline(residhdl)


#----------------
lam <- colMeans(out3$lambda)
sds <- colMeans(out3$sds)
pi <- colMeans(out3$pi)
gamma0 <- rep(0,7)
for(i in 2:K){
  gamma0 <- gamma0 + pi[i]*(lam[,1] - lam[,i])
}
gamma0 <- lam[,1]-gamma0

n <- round(3000*pi)
random <- matrix(0,ncol=7,nrow=3000)
for(i in 1:7){
  random[,i] <- c(rnorm(n[1],lam[i,1]-gamma0[i],sds[i,1]),rnorm(n[2],lam[i,2]-gamma0[i],sds[i,2]),rnorm(n[3],lam[i,3]-gamma0[i],sds[i,3]))
}

qqnorm(random[,1]);qqline(random[,1])
qqplot(residwaist,random[,1])
qqplot(residglu,random[,2])
qqplot(residtri,random[,3])
