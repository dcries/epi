library(ggplot2)
library(dplyr)
library(gridExtra)
library(reshape)


setwd("C:/Users/dcries/workspace")
nhanes <- read.csv("../github/epi/NHANES_complete.csv")
# names(nhanes) <- tolower(names(nhanes))
# nhanes$a <- 1;nhanes$a[nhanes$age>=35] <- 2;nhanes$a[nhanes$age>=50] <- 3;nhanes$a[nhanes$age>=65] <- 4
# nhanes$weekend <- 0
# nhanes$weekend[nhanes$dow %in% c(1,7)] <- 1
# nhanes$first5 <- 0
# nhanes$first5[nhanes$rep==6] <- 1
# nhanes$first5[nhanes$rep==7] <- 2
# nhanes <- subset(nhanes,rep!=7)
# m1 <- lm(modvigmin~(weekend)+first5,data=nhanes)
# wbar <- mean((nhanes$modvigmin[nhanes$rep <= 5]))
# w1 <- nhanes$modvigmin
# what <- predict(m1)
# w <- (1/what)*w1*wbar
# nhanes$w <- w^.25
# 
# load("stanout_naive.RData")
# naive <- as.matrix(rs)
# naivemeans <- colMeans(naive)
# load("stanout_mix5.RData") #!!!!!!!!!!!
# out=out5
# # out$beta <- out$beta[50000:100000,]
# # out$Sigma <- out$Sigma[,,50000:100000]
# indmeans <- nhanes %>% group_by(id) %>% summarise(m=mean(w),glu=glu[1],waist=waist[1],
#                                                   tri=tri[1],bps=bps[1],bpd=bpd[1],hdl=hdl[1],ldl=ldl[1],n=length(w))
# 
# waist <- mean(out$beta[,4])-mean(out$beta[,1])/(1+exp(-mean(out$beta[,2])*(x-mean(out$beta[,3]))))
# glu <- mean(out$beta[,8])-mean(out$beta[,5])/(1+exp(-mean(out$beta[,6])*(x-mean(out$beta[,7]))))
# tri <- mean(out$beta[,12])-mean(out$beta[,9])/(1+exp(-mean(out$beta[,10])*(x-mean(out$beta[,11]))))
# bps <- mean(out$beta[,16])-mean(out$beta[,13])/(1+exp(-mean(out$beta[,14])*(x-mean(out$beta[,15]))))
# bpd <- mean(out$beta[,19]) + mean(out$beta[,20])*x #+ mean(out$beta[,22])*x^2
# ldl <- mean(out$beta[,17]) + mean(out$beta[,18])*x #+ mean(out$beta[,15])*x^2
# hdl <- mean(out$beta[,21]) + mean(out$beta[,22])*x 
load("stanout_mix5.RData") #!!!!!!!!!!!
out=out5
MetS <- read.csv("MetS_imp.csv")
MetS <- read.csv("../github/epi/MetS.csv")
x <- seq(from=0.25,to=3.5,by=0.1)

lam <- colMeans(out$lambda)
sds <- colMeans(out$sds)
pi <- colMeans(out$pi)
gamma0 <- rep(0,7)
for(i in 2:K){
  gamma0 <- gamma0 + pi[i]*(lam[,1] - lam[,i])
}
gamma0 <- lam[,1]-gamma0

waist <- gamma0[1]-mean(out$beta[,1])/(1+exp(-mean(out$beta[,2])*(x-mean(out$beta[,3]))))
glu <- gamma0[2]-mean(out$beta[,4])/(1+exp(-mean(out$beta[,5])*(x-mean(out$beta[,6]))))
tri <- gamma0[3]-mean(out$beta[,7])/(1+exp(-mean(out$beta[,8])*(x-mean(out$beta[,9]))))
bps <- gamma0[4]-mean(out$beta[,10])/(1+exp(-mean(out$beta[,11])*(x-mean(out$beta[,12]))))
bpd <- gamma0[6] + mean(out$beta[,14])*x #+ mean(out$beta[,22])*x^2
ldl <- gamma0[5] + mean(out$beta[,13])*x #+ mean(out$beta[,15])*x^2
hdl <- gamma0[7] + mean(out$beta[,15])*x 


# naivewaist <- naivemeans[4]-naivemeans[1]/(1+exp(-naivemeans[2]*(x-naivemeans[3])))
# naiveglu <- naivemeans[8]-naivemeans[5]/(1+exp(-naivemeans[6]*(x-naivemeans[7])))
# naivetri <- naivemeans[12]-naivemeans[9]/(1+exp(-naivemeans[10]*(x-naivemeans[11])))
# naivebps <- naivemeans[19]-naivemeans[16]/(1+exp(-naivemeans[17]*(x-naivemeans[18])))
# naivebpd <- naivemeans[20] + naivemeans[21]*x + naivemeans[22]*x^2
# naiveldl <- naivemeans[13] + naivemeans[14]*x + naivemeans[15]*x^2
# naivehdl <- naivemeans[23] + naivemeans[24]*x 

p1=ggplot() + geom_point(data=MetS,aes(x=X,y=waist))  + geom_line(aes(x=x,y=waist),colour="blue",size=2) + theme_bw() + xlab("Usual Minutes in MVPA") + ylab("Waist Circumference")#+ geom_line(aes(x=x,y=naivewaist),colour="red",size=2,linetype=2)
p2=ggplot() + geom_point(data=MetS,aes(x=X,y=lglu)) + geom_line(aes(x=x,y=glu),colour="blue",size=2) + theme_bw() + xlab("Usual Minutes in MVPA") + ylab("log Glucose")#  + geom_line(aes(x=x,y=naiveglu),colour="red",size=2,linetype=2)#+ ylim(c(0,400))
p3=ggplot() + geom_point(data=MetS,aes(x=X,y=ltri)) + geom_line(aes(x=x,y=tri),colour="blue",size=2) + theme_bw() + xlab("Usual Minutes in MVPA") + ylab("log Triglycerides")#  + geom_line(aes(x=x,y=naivetri),colour="red",size=2,linetype=2)#+ ylim(c(0,1000))
p4=ggplot() + geom_point(data=MetS,aes(x=X,y=(bps)))  + geom_line(aes(x=x,y=bps),colour="blue",size=2) + theme_bw() + xlab("Usual Minutes in MVPA") + ylab("Systolic Blood Pressure")#+ geom_line(aes(x=x,y=naivebps),colour="red",size=2,linetype=2)
p5=ggplot() + geom_point(data=MetS,aes(x=X,y=bpd))  + geom_line(aes(x=x,y=bpd),colour="blue",size=2)+ theme_bw() + xlab("Usual Minutes in MVPA") + ylab("Diastolic Blood Pressure")# + geom_line(aes(x=x,y=naivebpd),colour="red",size=2,linetype=2)
p6=ggplot() + geom_point(data=MetS,aes(x=X,y=(ldl)))  + geom_line(aes(x=x,y=ldl),colour="blue",size=2)+ theme_bw() + xlab("Usual Minutes in MVPA") + ylab("LDL")# + geom_line(aes(x=x,y=naiveldl),colour="red",size=2,linetype=2)
p7=ggplot() + geom_point(data=MetS,aes(x=X,y=(hdl)))  + geom_line(aes(x=x,y=hdl),colour="blue",size=2)+ theme_bw() + xlab("Usual Minutes in MVPA") + ylab("HDL")# + geom_line(aes(x=x,y=naivehdl),colour="red",size=2,linetype=2)
grid.arrange(p1,p2,p3,p4,p5,p6,p7,nrow=4)


#-----------------------
load("C:/Users/dcries/workspace/prob_mets.RData")

probcc <- melt(outlist$minMVPA)
probimp <- melt(outlist$minMVPAimp)
#probcc$X2 <- paste(probcc$X2, " High MetS rf ")
#probimp$X2 <- paste(probimp$X2, " High MetS rf ")
cicc <- data.frame()
ciimp <- data.frame()
for(i in 1:6){
  cicc <- rbind(cicc,cbind(outlist$minMVPAquant[,(2*i-1):(2*i)],rep(i,nrow(outlist$minMVPAquant)),0:(nrow(outlist$minMVPAquant)-1)))
  ciimp <- rbind(ciimp,cbind(outlist$minMVPAquantimp[,(2*i-1):(2*i)],rep(i,nrow(outlist$minMVPAquantimp)),0:(nrow(outlist$minMVPAquant)-1)))
}

probcc <- cbind(probcc,cicc)
probimp <- cbind(probimp,ciimp)


ggplot(data=probcc) + geom_line(aes(x=X1-1,y=value,colour=as.factor(X2),group=as.factor(X2))) +
  geom_ribbon(aes(x=V4,ymin=V1,ymax=V2,group=as.factor(X2),fill=as.factor(X2)),alpha=0.2) + theme_bw() +
  guides(colour=guide_legend(title="No. High Risk Factors"),fill=FALSE) + xlab("Daily Minutes in MVPA")+ylab("Probability")#+
  # geom_line(data=probimp,aes(x=X1-1,y=value,colour=as.factor(X2),group=as.factor(X2))) +
  # geom_ribbon(data=probimp,aes(x=V4,ymin=V1,ymax=V2,group=as.factor(X2),fill=as.factor(X2)),alpha=0.2)


indprobcc <- melt(outlist$indprob)
indprobimp <- melt(outlist$indprobimp)
indprobcc$X2[indprobcc$X2==1] <- "P(High Waist Circum)";indprobcc$X2[indprobcc$X2==2] <- "P(High Glucose)";indprobcc$X2[indprobcc$X2==3] <- "P(High Triglycerides)";indprobcc$X2[indprobcc$X2==4] <- "P(High Sys Blood Press)";indprobcc$X2[indprobcc$X2==5] <- "P(High LDL)";indprobcc$X2[indprobcc$X2==6] <- "P(High Dias Blood Press)";indprobcc$X2[indprobcc$X2==7] <- "P(Low HDL)";
indprobimp$X2[indprobimp$X2==1] <- "P(High Waist Circum)";indprobimp$X2[indprobimp$X2==2] <- "P(High Glucose)";indprobimp$X2[indprobimp$X2==3] <- "P(High Triglycerides)";indprobimp$X2[indprobimp$X2==4] <- "P(High Sys Blood Press)";indprobimp$X2[indprobimp$X2==5] <- "P(High LDL)";indprobimp$X2[indprobimp$X2==6] <- "P(High Dias Blood Press)";indprobimp$X2[indprobimp$X2==7] <- "P(Low HDL)";
ciccind <- data.frame()
ciimpind <- data.frame()
for(i in 1:7){
  ciccind <- rbind(ciccind,cbind(outlist$indprobquant[,(2*i-1):(2*i)],rep(i,nrow(outlist$indprobquant)),0:(nrow(outlist$indprobquant)-1)))
  ciimpind <- rbind(ciimpind,cbind(outlist$indprobquantimp[,(2*i-1):(2*i)],rep(i,nrow(outlist$indprobquantimp)),0:(nrow(outlist$indprobquantimp)-1)))
}

indcc <- cbind(indprobcc,ciccind)
indimp <- cbind(indprobimp,ciimpind)

names(indcc)[2] <- "MetS_RF"
names(indimp)[2] <- "MetS_RF"

ggplot(data=indcc) + geom_line(aes(x=X1-1,y=value,colour=MetS_RF,group=MetS_RF,linetype=MetS_RF),size=1) +theme_bw() +
  xlab("Daily Minutes in MVPA")+ylab("Probability")  #+guides(colour=guide_legend(title=""))
  #geom_ribbon(aes(x=V4,ymin=V1,ymax=V2,group=as.factor(X2),fill=as.factor(X2)),alpha=0.2) 

ggplot(data=indimp) + geom_line(aes(x=X1-1,y=value,colour=MetS_RF,group=MetS_RF,linetype=MetS_RF),size=1) +theme_bw() +
  xlab("Daily Minutes in MVPA")+ylab("Probability")  #+guides(colour=guide_legend(title=""))

