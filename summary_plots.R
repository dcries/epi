library(ggplot2)
library(dplyr)
library(gridExtra)

setwd("C:/Users/dcries/workspace")
nhanes <- read.csv("../github/epi/NHANES_complete.csv")
names(nhanes) <- tolower(names(nhanes))
nhanes$a <- 1;nhanes$a[nhanes$age>=35] <- 2;nhanes$a[nhanes$age>=50] <- 3;nhanes$a[nhanes$age>=65] <- 4
nhanes$weekend <- 0
nhanes$weekend[nhanes$dow %in% c(1,7)] <- 1
nhanes$first5 <- 0
nhanes$first5[nhanes$rep==6] <- 1
nhanes$first5[nhanes$rep==7] <- 2
nhanes <- subset(nhanes,rep!=7)
m1 <- lm(modvigmin~(weekend)+first5,data=nhanes)
wbar <- mean((nhanes$modvigmin[nhanes$rep <= 5]))
w1 <- nhanes$modvigmin
what <- predict(m1)
w <- (1/what)*w1*wbar
nhanes$w <- w^.25

load("stanout_naive.RData")
naive <- as.matrix(rs)
naivemeans <- colMeans(naive)
load("stanout_mix5.RData") #!!!!!!!!!!!
out=out5
# out$beta <- out$beta[50000:100000,]
# out$Sigma <- out$Sigma[,,50000:100000]
indmeans <- nhanes %>% group_by(id) %>% summarise(m=mean(w),glu=glu[1],waist=waist[1],
                                                  tri=tri[1],bps=bps[1],bpd=bpd[1],hdl=hdl[1],ldl=ldl[1],n=length(w))

x <- seq(from=0,to=4,by=0.1)
# waist <- mean(out$beta[,4])-mean(out$beta[,1])/(1+exp(-mean(out$beta[,2])*(x-mean(out$beta[,3]))))
# glu <- mean(out$beta[,8])-mean(out$beta[,5])/(1+exp(-mean(out$beta[,6])*(x-mean(out$beta[,7]))))
# tri <- mean(out$beta[,12])-mean(out$beta[,9])/(1+exp(-mean(out$beta[,10])*(x-mean(out$beta[,11]))))
# bps <- mean(out$beta[,16])-mean(out$beta[,13])/(1+exp(-mean(out$beta[,14])*(x-mean(out$beta[,15]))))
# bpd <- mean(out$beta[,19]) + mean(out$beta[,20])*x #+ mean(out$beta[,22])*x^2
# ldl <- mean(out$beta[,17]) + mean(out$beta[,18])*x #+ mean(out$beta[,15])*x^2
# hdl <- mean(out$beta[,21]) + mean(out$beta[,22])*x 

waist <- mean(out$lambda[,1,1])-mean(out$beta[,1])/(1+exp(-mean(out$beta[,2])*(x-mean(out$beta[,3]))))
glu <- mean(out$lambda[,2,1])-mean(out$beta[,4])/(1+exp(-mean(out$beta[,5])*(x-mean(out$beta[,6]))))
tri <- mean(out$lambda[,3,1])-mean(out$beta[,7])/(1+exp(-mean(out$beta[,8])*(x-mean(out$beta[,9]))))
bps <- mean(out$lambda[,4,1])-mean(out$beta[,10])/(1+exp(-mean(out$beta[,11])*(x-mean(out$beta[,12]))))
bpd <- mean(out$lambda[,6,1]) + mean(out$beta[,14])*x #+ mean(out$beta[,22])*x^2
ldl <- mean(out$lambda[,5,1]) + mean(out$beta[,13])*x #+ mean(out$beta[,15])*x^2
hdl <- mean(out$lambda[,7,1]) + mean(out$beta[,15])*x 


naivewaist <- naivemeans[4]-naivemeans[1]/(1+exp(-naivemeans[2]*(x-naivemeans[3])))
naiveglu <- naivemeans[8]-naivemeans[5]/(1+exp(-naivemeans[6]*(x-naivemeans[7])))
naivetri <- naivemeans[12]-naivemeans[9]/(1+exp(-naivemeans[10]*(x-naivemeans[11])))
naivebps <- naivemeans[19]-naivemeans[16]/(1+exp(-naivemeans[17]*(x-naivemeans[18])))
naivebpd <- naivemeans[20] + naivemeans[21]*x + naivemeans[22]*x^2
naiveldl <- naivemeans[13] + naivemeans[14]*x + naivemeans[15]*x^2
naivehdl <- naivemeans[23] + naivemeans[24]*x 

p1=ggplot() + geom_point(data=indmeans,aes(x=m,y=waist))  + geom_line(aes(x=x,y=waist),colour="blue",size=2) + theme_bw()+ geom_line(aes(x=x,y=naivewaist),colour="red",size=2,linetype=2)
p2=ggplot() + geom_point(data=indmeans,aes(x=m,y=log(glu))) + geom_line(aes(x=x,y=glu),colour="blue",size=2) + theme_bw()  + geom_line(aes(x=x,y=naiveglu),colour="red",size=2,linetype=2)#+ ylim(c(0,400))
p3=ggplot() + geom_point(data=indmeans,aes(x=m,y=log(tri))) + geom_line(aes(x=x,y=tri),colour="blue",size=2) + theme_bw()  + geom_line(aes(x=x,y=naivetri),colour="red",size=2,linetype=2)#+ ylim(c(0,1000))
p4=ggplot() + geom_point(data=indmeans,aes(x=m,y=(bps)))  + geom_line(aes(x=x,y=bps),colour="blue",size=2) + theme_bw()+ geom_line(aes(x=x,y=naivebps),colour="red",size=2,linetype=2)
p5=ggplot() + geom_point(data=indmeans,aes(x=m,y=bpd))  + geom_line(aes(x=x,y=bpd),colour="blue",size=2)+ theme_bw() + geom_line(aes(x=x,y=naivebpd),colour="red",size=2,linetype=2)
p6=ggplot() + geom_point(data=indmeans,aes(x=m,y=(ldl)))  + geom_line(aes(x=x,y=ldl),colour="blue",size=2)+ theme_bw() + geom_line(aes(x=x,y=naiveldl),colour="red",size=2,linetype=2)
p7=ggplot() + geom_point(data=indmeans,aes(x=m,y=(hdl)))  + geom_line(aes(x=x,y=hdl),colour="blue",size=2)+ theme_bw() + geom_line(aes(x=x,y=naivehdl),colour="red",size=2,linetype=2)
grid.arrange(p1,p2,p3,p4,p5,p6,p7,nrow=4)

#load("stanout_second.RData") ???
# dat <- as.matrix(rs)
# postmeans <- colMeans(dat)


# indmeans <- nhanes %>% group_by(id) %>% summarise(m=mean(w),glu=glu[1],waist=waist[1],
#                         tri=tri[1],bps=bps[1],bpd=bpd[1],hdl=hdl[1],ldl=ldl[1],n=length(w))
# 
# x <- seq(from=0,to=4,by=0.1)
# waist <- postmeans[4]-postmeans[1]/(1+exp(-postmeans[2]*(x-postmeans[3])))
# glu <- postmeans[8]-postmeans[5]/(1+exp(-postmeans[6]*(x-postmeans[7])))
# tri <- postmeans[12]-postmeans[9]/(1+exp(-postmeans[10]*(x-postmeans[11])))
# bps <- postmeans[19]-postmeans[16]/(1+exp(-postmeans[17]*(x-postmeans[18])))
# bpd <- postmeans[20] + postmeans[21]*x #+ postmeans[22]*x^2
# ldl <- postmeans[13] + postmeans[14]*x #+ postmeans[15]*x^2
# hdl <- postmeans[23] + postmeans[24]*x 
# 
# naivewaist <- naivemeans[4]-naivemeans[1]/(1+exp(-naivemeans[2]*(x-naivemeans[3])))
# naiveglu <- naivemeans[8]-naivemeans[5]/(1+exp(-naivemeans[6]*(x-naivemeans[7])))
# naivetri <- naivemeans[12]-naivemeans[9]/(1+exp(-naivemeans[10]*(x-naivemeans[11])))
# naivebps <- naivemeans[19]-naivemeans[16]/(1+exp(-naivemeans[17]*(x-naivemeans[18])))
# naivebpd <- naivemeans[20] + naivemeans[21]*x + naivemeans[22]*x^2
# naiveldl <- naivemeans[13] + naivemeans[14]*x + naivemeans[15]*x^2
# naivehdl <- naivemeans[23] + naivemeans[24]*x 
# 
# ggplot() + geom_point(data=indmeans,aes(x=m,y=waist)) + geom_line(aes(x=x,y=waist),colour="red",size=2) + geom_line(aes(x=x,y=naivewaist),colour="blue",size=2) + theme_bw()
# ggplot() + geom_point(data=indmeans,aes(x=m,y=log(glu))) + geom_line(aes(x=x,y=glu),colour="red",size=2) + geom_line(aes(x=x,y=naiveglu),colour="blue",size=2) + theme_bw()
# ggplot() + geom_point(data=indmeans,aes(x=m,y=log(tri))) + geom_line(aes(x=x,y=tri),colour="red",size=2) + geom_line(aes(x=x,y=naivetri),colour="blue",size=2) + theme_bw()
# ggplot() + geom_point(data=indmeans,aes(x=m,y=bps)) + geom_line(aes(x=x,y=bps),colour="red",size=2) + geom_line(aes(x=x,y=naivebps),colour="blue",size=2) + theme_bw()
# ggplot() + geom_point(data=indmeans,aes(x=m,y=bpd)) + geom_line(aes(x=x,y=bpd),colour="red",size=2) + geom_line(aes(x=x,y=naivebpd),colour="blue",size=2) + theme_bw()
# ggplot() + geom_point(data=indmeans,aes(x=m,y=ldl)) + geom_line(aes(x=x,y=ldl),colour="red",size=2) + geom_line(aes(x=x,y=naiveldl),colour="blue",size=2) + theme_bw()
# ggplot() + geom_point(data=indmeans,aes(x=m,y=hdl)) + geom_line(aes(x=x,y=hdl),colour="red",size=2) + geom_line(aes(x=x,y=naivehdl),colour="blue",size=2) + theme_bw()


#----------
#resiudal calculations
# tstar <- read.csv("../github/epi/tstar.csv")
# 
# nhanes <- subset(nhanes,rep!=7)
# m1 <- lm(modvigmin~(weekend)+first5,data=nhanes)
# wbar <- mean((nhanes$modvigmin[nhanes$rep <= 5]))
# w1 <- nhanes$modvigmin
# what <- predict(m1)
# w <- (1/what)*w1*wbar
# nhanes$w <- w^.25
# nrep <- (nhanes %>% group_by(id) %>% summarise(n=length(id)))$n
# meas7 <- subset(nhanes, id %in% unique(id)[nrep==6]) #individuals with all 7 days
# indmeans2 <- meas7 %>% group_by(id) %>% summarise(m=mean(w),glu=glu[1],waist=waist[1],
#                                                   tri=tri[1],bps=bps[1],bpd=bpd[1],hdl=hdl[1],ldl=ldl[1],n=length(w),age=age[1])
# 
# metsrf <- indmeans2[(!is.na(indmeans2$waist)) & (!is.na(indmeans2$bps)) & (!is.na(indmeans2$bpd)) & (!is.na(indmeans2$hdl)) & (!is.na(indmeans2$ldl)) & (!is.na(indmeans2$glu)) & (!is.na(indmeans2$tri)) ,] #remove NAs for waist
# metsrf$tstar <- tstar$x
# 
# metsrf$pwaist <- mean(out$beta[,4])-mean(out$beta[,1])/(1+exp(-mean(out$beta[,2])*(metsrf$tstar-mean(out$beta[,3]))))
# metsrf$pglu <- mean(out$beta[,8])-mean(out$beta[,5])/(1+exp(-mean(out$beta[,6])*(metsrf$tstar-mean(out$beta[,7]))))
# metsrf$ptri <- mean(out$beta[,12])-mean(out$beta[,9])/(1+exp(-mean(out$beta[,10])*(metsrf$tstar-mean(out$beta[,11]))))
# metsrf$pbps <- mean(out$beta[,16])-mean(out$beta[,13])/(1+exp(-mean(out$beta[,14])*(metsrf$tstar-mean(out$beta[,15]))))
# metsrf$pldl <- mean(out$beta[,17]) + mean(out$beta[,18])*metsrf$tstar #+ mean(out$beta[,22])*metsrf$tstar^2
# metsrf$pbpd <- mean(out$beta[,19]) + mean(out$beta[,20])*metsrf$tstar #+ mean(out$beta[,15])*metsrf$tstar^2
# metsrf$phdl <- mean(out$beta[,21]) + mean(out$beta[,22])*metsrf$tstar 
# swaist <- mean(sqrt(out$Sigma[1,1,]));sglu <- mean(sqrt(out$Sigma[2,2,]));stri <- mean(sqrt(out$Sigma[3,3,]));sbps <- mean(sqrt(out$Sigma[4,4,]));
# sbpd <- mean(sqrt(out$Sigma[6,6,]));sldl <- mean(sqrt(out$Sigma[5,5,]));shdl <- mean(sqrt(out$Sigma[7,7,]));
# 
# 
# # metsrf$pwaist <- postmeans[4]-postmeans[1]/(1+exp(-postmeans[2]*(metsrf$tstar-postmeans[3])))
# # metsrf$pglu <- postmeans[8]-postmeans[5]/(1+exp(-postmeans[6]*(metsrf$tstar-postmeans[7])))
# # metsrf$ptri <- postmeans[12]-postmeans[9]/(1+exp(-postmeans[10]*(metsrf$tstar-postmeans[11])))
# # metsrf$pbps <- postmeans[19]-postmeans[16]/(1+exp(-postmeans[17]*(metsrf$tstar-postmeans[18])))
# # metsrf$pbpd <- postmeans[20] + postmeans[21]*metsrf$tstar + postmeans[22]*metsrf$tstar^2
# # metsrf$pldl <- postmeans[13] + postmeans[14]*metsrf$tstar + postmeans[15]*metsrf$tstar^2
# # metsrf$phdl <- postmeans[23] + postmeans[24]*metsrf$tstar 
# # swaist <- postmeans[74];sglu <- postmeans[75];stri <- postmeans[76];sbps <- postmeans[77];
# #sbpd <- postmeans[78];sldl <- postmeans[79];shdl <- postmeans[80];
# 
# p1=qplot(data=metsrf,y=(waist-pwaist)/swaist,x=tstar) + geom_smooth()
# p2=qplot(data=metsrf,y=(log(glu)-pglu)/sglu,x=tstar) + geom_smooth()
# p3=qplot(data=metsrf,y=(log(tri)-ptri)/stri,x=tstar) + geom_smooth()
# p4=qplot(data=metsrf,y=(bps-pbps)/sbps,x=tstar) + geom_smooth()
# p5=qplot(data=metsrf,y=(bpd-pbpd)/sbpd,x=tstar) + geom_smooth()
# p6=qplot(data=metsrf,y=(ldl-pldl)/sldl,x=tstar) + geom_smooth()
# p7=qplot(data=metsrf,y=(hdl-phdl)/shdl,x=tstar) + geom_smooth()
# grid.arrange(p1,p2,p3,p4,p5,p6,p7,nrow=4)
# 
# par(mfrow=c(4,2))
# qqnorm((metsrf$waist-metsrf$pwaist)/swaist);qqline((metsrf$waist-metsrf$pwaist)/swaist)
# qqnorm((log(metsrf$glu)-metsrf$pglu)/sglu);qqline((log(metsrf$glu)-metsrf$pglu)/sglu)
# qqnorm((log(metsrf$tri)-metsrf$ptri)/stri);qqline((log(metsrf$tri)-metsrf$ptri)/stri)
# qqnorm((metsrf$bps-metsrf$pbps)/sbps);qqline((metsrf$bps-metsrf$pbps)/sbps)
# qqnorm((metsrf$bpd-metsrf$pbpd)/sbpd);qqline((metsrf$bpd-metsrf$pbpd)/sbpd)
# qqnorm((metsrf$ldl-metsrf$pldl)/sldl);qqline((metsrf$ldl-metsrf$pldl)/sldl)
# qqnorm((metsrf$hdl-metsrf$phdl)/shdl);qqline((metsrf$hdl-metsrf$phdl)/shdl)
