library(ggplot2)
library(dplyr)

setwd("C:/Users/dcries/workspace2")
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

#load("stanout_second.RData") ???
indmeans <- nhanes %>% group_by(id) %>% summarise(m=mean(w),glu=glu[1],waist=waist[1],
                        tri=tri[1],bps=bps[1],bpd=bpd[1],hdl=hdl[1],ldl=ldl[1],n=length(w))
dat <- as.matrix(rs)
postmeans <- colMeans(dat)

x <- seq(from=0,to=4,by=0.1)
waist <- postmeans[4]-postmeans[1]/(1+exp(-postmeans[2]*(x-postmeans[3])))
glu <- postmeans[8]-postmeans[5]/(1+exp(-postmeans[6]*(x-postmeans[7])))
tri <- postmeans[12]-postmeans[9]/(1+exp(-postmeans[10]*(x-postmeans[11])))
bps <- postmeans[19]-postmeans[16]/(1+exp(-postmeans[17]*(x-postmeans[18])))
bpd <- postmeans[20] + postmeans[21]*x + postmeans[22]*x^2
ldl <- postmeans[13] + postmeans[14]*x + postmeans[15]*x^2
hdl <- postmeans[23] + postmeans[24]*x 

ggplot() + geom_point(data=indmeans,aes(x=m,y=waist)) + geom_line(aes(x=x,y=waist),colour="red",size=2) + theme_bw()
ggplot() + geom_point(data=indmeans,aes(x=m,y=log(glu))) + geom_line(aes(x=x,y=glu),colour="red",size=2) + theme_bw()
ggplot() + geom_point(data=indmeans,aes(x=m,y=log(tri))) + geom_line(aes(x=x,y=tri),colour="red",size=2) + theme_bw()
ggplot() + geom_point(data=indmeans,aes(x=m,y=bps)) + geom_line(aes(x=x,y=bps),colour="red",size=2) + theme_bw()
ggplot() + geom_point(data=indmeans,aes(x=m,y=bpd)) + geom_line(aes(x=x,y=bpd),colour="red",size=2) + theme_bw()
ggplot() + geom_point(data=indmeans,aes(x=m,y=ldl)) + geom_line(aes(x=x,y=ldl),colour="red",size=2) + theme_bw()
ggplot() + geom_point(data=indmeans,aes(x=m,y=hdl)) + geom_line(aes(x=x,y=hdl),colour="red",size=2) + theme_bw()
