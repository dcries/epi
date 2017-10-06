library(MASS)
library(mvtnorm)
library(LaplacesDemon)
library(dplyr)
library(coda)

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

#-------------------------------------
#check mixing
load("C:/Users/dcries/workspace/stanout_realmix3.RData")
out3a=out3
load("C:/Users/dcries/workspace/stanout_realmix3b.RData")
out3b=out3
load("C:/Users/dcries/workspace/stanout_realmix3c.RData")
out3c=out3

length(unique(out3a$beta[,1]))/nrow(out3a$beta)
length(unique(out3b$beta[,1]))/nrow(out3b$beta)
length(unique(out3c$beta[,1]))/nrow(out3c$beta)

beta <- mcmc.list(list(mcmc(out3a$beta),mcmc(out3b$beta),mcmc(out3c$beta)))
gelman.diag(beta)
pi <- mcmc.list(list(mcmc(out3a$pi[,1]),mcmc(out3b$pi[,1]),mcmc(out3c$pi[,3])))
gelman.diag(pi)
lambda <- mcmc.list(list(mcmc(out3a$lambda[,2,1]),mcmc(out3b$lambda[,2,1]),mcmc(out3c$lambda[,2,3])))
gelman.diag(lambda)
sigma <- mcmc.list(list(mcmc(out3a$sds[,7,1]),mcmc(out3b$sds[,7,1]),mcmc(out3c$sds[,7,3])))
gelman.diag(sigma)

plot(out3a$beta[,1],type="l")
lines(out3b$beta[,1],type="l",col=2)
lines(out3c$beta[,1],type="l",col=3)

length(unique(out$beta[,1]))/nrow(out$beta)
diag(out3a$propcov)


plot(out$beta[,1],type="l") #waist decrease
plot(out$beta[,2],type="l") #waist k
plot(out$beta[,3],type="l") #waist midpoint
plot(out$beta[,4],type="l") #glu decrease
plot(out$beta[,5],type="l") #glu k
plot(out$beta[,6],type="l") #glu midpoint
plot(out$beta[,7],type="l") #tri decrease
plot(out$beta[,8],type="l") #tri k
plot(out$beta[,9],type="l") #tri mid
plot(out$beta[,10],type="l") #bps decrease
plot(out$beta[,11],type="l") #bps k
plot(out$beta[,12],type="l") #bps mid
plot(out$beta[,13],type="l") #ldl
plot(out$beta[,14],type="l") #bpd
plot(out$beta[,15],type="l") #hdl

maxk=3
name=c("w","G","t","bps","ldl","bpd","hdl")
par(mfrow=c(3,2))
for(i in 1:7){
  for(j in 1:6){
    if(j>maxk){
      plot.new()
    }
    else{
      plot((out$sds[,i,j]),type="l",main=name[i]) #w
    }
 }
}

for(i in 1:7){
  for(j in 1:6){
    if(j>maxk){
      plot.new()
    }
    else{
      plot((out$lambda[,i,j]),type="l",main=name[i]) #w
    }
  }
}


#---------------
#test if groups differ by sex,race,etc
demo <- read.csv("demographics_imp.csv")
#demo <- read.csv("demographics.csv")

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
groups <- apply(out$zeta,2,Mode)
demo$groups <- groups
demo$a <- 1;demo$a[demo$age>=35] <- 2;demo$a[demo$age>=50] <- 3;demo$a[demo$age>=65] <- 4

chisq.test(table(demo$groups,demo$sex))
chisq.test(table(demo$groups,demo$race))
chisq.test(table(demo$groups,demo$education))
chisq.test(table(demo$groups,demo$a))

#------------
#distributyion of errors
prop <- colMeans(out5$pi)

xxw <- seq(from=-50,to=50,by=0.1)
wsd <- colMeans(out5$sds[,1,])
wlam <- colMeans(out5$lambda[,1,]) - mean(out5$lambda[,1,]%*%prop)
wden <- prop[1]*dnorm(xxw,wlam[1],wsd[1]) + prop[2]*dnorm(xxw,wlam[2],wsd[2]) + prop[3]*dnorm(xxw,wlam[3],wsd[3]) + prop[4]*dnorm(xxw,wlam[4],wsd[4]) + prop[5]*dnorm(xxw,wlam[5],wsd[5]) 
e1=qplot(x=xxw,y=wden,geom="line") +ggtitle("waist") + theme_bw()#+ geom_line(aes(x=xxw,y=dnorm(xxw,0,12.8)),colour="red")

xxg <- seq(from=-1,to=1,by=0.01)
glam <- colMeans(out5$lambda[,2,]) - mean(out5$lambda[,2,]%*%prop)
gsd <- colMeans(out5$sds[,2,])
gden <- prop[1]*dnorm(xxg,glam[1],gsd[1]) + prop[2]*dnorm(xxg,glam[2],gsd[2]) + prop[3]*dnorm(xxg,glam[3],gsd[3]) + prop[4]*dnorm(xxg,glam[4],gsd[4]) + prop[5]*dnorm(xxg,glam[5],gsd[5]) 
e2=qplot(x=xxg,y=gden,geom="line") +ggtitle("glu")+ theme_bw()#+ geom_line(aes(x=xxg,y=dnorm(xxg,0,.12)),colour="red")

xxt <- seq(from=-2.5,to=2.5,by=0.05)
tlam <- colMeans(out5$lambda[,3,]) - mean(out5$lambda[,3,]%*%prop)
tsd <- colMeans(out5$sds[,3,])
tden <- prop[1]*dnorm(xxt,tlam[1],tsd[1]) + prop[2]*dnorm(xxt,tlam[2],tsd[2]) + prop[3]*dnorm(xxt,tlam[3],tsd[3]) + prop[4]*dnorm(xxt,tlam[4],tsd[4]) + prop[5]*dnorm(xxt,tlam[5],tsd[5]) 
e3=qplot(x=xxt,y=tden,geom="line") +ggtitle("tri")+ theme_bw()#+ geom_line(aes(x=xxt,y=dnorm(xxt,0,.5)),colour="red")

xxbps <- seq(from=-70,to=70,by=0.5)
bpslam <- colMeans(out5$lambda[,4,]) - mean(out5$lambda[,4,]%*%prop)
bpssd <- colMeans(out5$sds[,4,])
bpsden <- prop[1]*dnorm(xxbps,bpslam[1],bpssd[1]) + prop[2]*dnorm(xxbps,bpslam[2],bpssd[2]) + prop[3]*dnorm(xxbps,bpslam[3],bpssd[3]) + prop[4]*dnorm(xxbps,bpslam[4],bpssd[4]) + prop[5]*dnorm(xxbps,bpslam[5],bpssd[5]) 
e4=qplot(x=xxbps,y=bpsden,geom="line") +ggtitle("bps")+ theme_bw()#+ geom_line(aes(x=xxbps,y=dnorm(xxbps,0,.5)),colour="red")

xxbpd <- seq(from=-50,to=50,by=0.5)
bpdlam <- colMeans(out5$lambda[,6,]) - mean(out5$lambda[,6,]%*%prop)
bpdsd <- colMeans(out5$sds[,6,])
bpdden <- prop[1]*dnorm(xxbpd,bpdlam[1],bpdsd[1]) + prop[2]*dnorm(xxbpd,bpdlam[2],bpdsd[2]) + prop[3]*dnorm(xxbpd,bpdlam[3],bpdsd[3]) + prop[4]*dnorm(xxbpd,bpdlam[4],bpdsd[4]) + prop[5]*dnorm(xxbpd,bpdlam[5],bpdsd[5]) 
e5=qplot(x=xxbpd,y=bpdden,geom="line") +ggtitle("bpd")+ theme_bw()#+ geom_line(aes(x=xxbpd,y=dnorm(xxbpd,0,.5)),colour="red")

xxl <- seq(from=-120,to=120,by=0.5)
llam <- colMeans(out5$lambda[,5,]) - mean(out5$lambda[,5,]%*%prop)
lsd <- colMeans(out5$sds[,5,])
lden <- prop[1]*dnorm(xxl,llam[1],lsd[1]) + prop[2]*dnorm(xxl,llam[2],lsd[2]) + prop[3]*dnorm(xxl,llam[3],lsd[3]) + prop[4]*dnorm(xxl,llam[4],lsd[4]) + prop[5]*dnorm(xxl,llam[5],lsd[5]) 
e6=qplot(x=xxl,y=lden,geom="line") +ggtitle("ldl")+ theme_bw()#+ geom_line(aes(x=xxl,y=dnorm(xxl,0,.5)),colour="red")

xxh <- seq(from=-60,to=60,by=0.5)
hlam <- colMeans(out5$lambda[,7,]) - mean(out5$lambda[,7,]%*%prop)
hsd <- colMeans(out5$sds[,7,])
hden <- prop[1]*dnorm(xxh,hlam[1],hsd[1]) + prop[2]*dnorm(xxh,hlam[2],hsd[2]) + prop[3]*dnorm(xxh,hlam[3],hsd[3]) + prop[4]*dnorm(xxh,hlam[4],hsd[4]) + prop[5]*dnorm(xxh,hlam[5],hsd[5]) 
e7=qplot(x=xxh,y=hden,geom="line") +ggtitle("hdl")+ theme_bw()#+ geom_line(aes(x=xxh,y=dnorm(xxh,0,.5)),colour="red")

grid.arrange(e1,e2,e3,e4,e5,e6,e7,nrow=4)


#------------------