library(ggplot2)
library(lme4)
library(dplyr)
library(reshape)
library(car)
library(nlme)

setwd("C:\\Users\\dcries\\github\\epi")
nhanes <- read.csv("NHANES_complete.csv")

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


#-------------
a=nhanes %>% group_by(id) %>% summarise(m=mean(w),
        g=sex[1],r=race[1],age=age[1],n=length(sex),s=(n-1)*var(w),
        glu=glu[1],tri=tri[1],waist=waist[1],bps=bps[1],bpd=bpd[1],ldl=ldl[1],hdl=hdl[1],a=a[1])

#for mean function needing to split parameters by gender, etc
qplot(data=a,x=m,y=waist,group=as.factor(g)) +geom_smooth()
qplot(data=a,x=m,y=waist,group=as.factor(r)) +geom_smooth()
qplot(data=a,x=m,y=waist,group=as.factor(a)) +geom_smooth()

qplot(data=a,x=m,y=log(tri),group=as.factor(g)) +geom_smooth()
qplot(data=a,x=m,y=log(tri),group=as.factor(r)) +geom_smooth()
qplot(data=a,x=m,y=log(tri),group=as.factor(a)) +geom_smooth()

qplot(data=a,x=m,y=log(glu),group=as.factor(g)) +geom_smooth()
qplot(data=a,x=m,y=log(glu),group=as.factor(r)) +geom_smooth()
qplot(data=a,x=m,y=log(glu),group=as.factor(a)) +geom_smooth()

qplot(data=a,x=m,y=bps,group=as.factor(g)) +geom_smooth()
qplot(data=a,x=m,y=bps,group=as.factor(r)) +geom_smooth()
qplot(data=a,x=m,y=bps,group=as.factor(a)) +geom_smooth()

qplot(data=a,x=m,y=bpd,group=as.factor(g)) +geom_smooth()
qplot(data=a,x=m,y=bpd,group=as.factor(r)) +geom_smooth()
qplot(data=a,x=m,y=bpd,group=as.factor(a)) +geom_smooth()

qplot(data=a,x=m,y=ldl,group=as.factor(g)) +geom_smooth()
qplot(data=a,x=m,y=ldl,group=as.factor(r)) +geom_smooth()
qplot(data=a,x=m,y=ldl,group=as.factor(a)) +geom_smooth()

qplot(data=a,x=m,y=hdl,group=as.factor(g)) +geom_smooth()
qplot(data=a,x=m,y=hdl,group=as.factor(r)) +geom_smooth()
qplot(data=a,x=m,y=hdl,group=as.factor(a)) +geom_smooth()

#for variance
nhanes2 <- subset(nhanes,w>0)
d = nhanes2 %>% group_by(id) %>% summarise(m=mean(w),glu=glu[1],waist=waist[1],tri=tri[1],bps=bps[1],bpd=bpd[1],hdl=hdl[1],ldl=ldl[1],n=length(w))

mfull <- lme(w~age+as.factor(race)+as.factor(sex),data=nhanes2,random=~1|id,correlation=corAR1(form=~1|id),method="REML")
rf <- NULL
rfs <- unlist(ranef(mfull))
for(i in 1:length(unique(nhanes2$id))){
  rf <- c(rf,rep(rfs[i],d$n[i]))
}
nhanes2$w2 <- nhanes2$w - predict(mfull) - rf
agevar=nhanes2 %>% group_by(age) %>% summarise(v=var(w2),weights=smplwt[1])

qplot(data=agevar,x=age,y=v)  + geom_smooth()

m1 <- lm(v~age+I(age^2)+I(age^3),data=agevar)
num <- 18:85
eval <- predict(m1)
# m1 <- nls(v~L2+L/(1+exp(-k*(age-x0))),start=list(L=.21,k=2,x0=55,L2=.24),data=agevar)
# eval <- coef(m1)["L2"]+coef(m1)["L"]/(1+exp(-coef(m1)["k"]*(num-coef(m1)["x0"])))
qplot(data=agevar,x=age,y=v)  + geom_smooth() + geom_line(aes(x=num,y=eval),colour="red")


#----------------------
#don't know about these
a1 <- a[a$n>1,]
a1$a <- 1;a1$a[a1$age>=35] <- 2;a1$a[a1$age>=50] <- 3;a1$a[a1$age>=65] <- 4
a1$b <- 1;a1$b[a1$age>=35] <- 2;a1$b[a1$age>=50] <- 3;a1$b[a1$age>=65] <- 4

a1 %>% group_by(g) %>% summarise(m=mean(s))
a1 %>% group_by(r) %>% summarise(m=mean(s))
a1 %>% group_by(a) %>% summarise(m=mean(s))
a1 %>% group_by(g,r) %>% summarise(m=mean(s))
a1 %>% group_by(g,a) %>% summarise(m=mean(s))
a1 %>% group_by(a,r) %>% summarise(m=mean(s))

a2 <- a1 %>% group_by(age) %>% summarise(m=mean(s))
a3 <- a1 %>% group_by(age,r) %>% summarise(m=mean(s))

#plot age versus avg variance, appears to be s shape curve
qplot(data=a2,x=age,y=m)  + geom_smooth()
qplot(data=subset(a3,r==1),x=age,y=m)  + geom_smooth()


# m1 <- nls(m~L2+L/(1+exp(-k*(age-x0))),start=list(L=.6,k=2,x0=52,L2=.7),data=a2)
# num <- 18:85
# eval <- coef(m1)["L2"]+coef(m1)["L"]/(1+exp(-coef(m1)["k"]*(num-coef(m1)["x0"])))
# qplot(data=a2,x=age,y=m)  + geom_smooth() + geom_line(aes(x=num,y=eval),colour="red")
# 
# m2 <- lm(m~age,data=a2)
# plot(predict(m2),resid(m2));abline(a=0,b=0)
# plot(a2$age,a2$m);lines(a2$age,predict(m2),col="red")
# m3 <- lm(m~age+I(age^2),data=a2)
# anova(m3,m2)
# plot(predict(m3),resid(m3));abline(a=0,b=0)
# plot(a2$age,a2$m);lines(a2$age,predict(m3),col="red")

#different values of rho
nhanes2$b <- 1;nhanes2$b[nhanes2$bmi>=18] <- 2;nhanes2$b[nhanes2$bmi>=25] <- 3;nhanes2$b[nhanes2$bmi>=30] <- 4;nhanes2$b[nhanes2$bmi>=35] <- 5

m1c <- lme(w~1,data=nhanes2,random=~1|id,correlation=corAR1(form=~1|id),method="ML")

mm <- lme(w~1,data=subset(nhanes2,sex==1),random=~1|id,correlation=corAR1(form=~1|id),method="ML")
mf <- lme(w~age+as.factor(race),data=subset(nhanes2,sex==2),random=~1|id,correlation=corAR1(form=~1|id),method="ML")

mr1 <- lme(w~1,data=subset(nhanes2,race==1),random=~1|id,correlation=corAR1(form=~1|id),method="ML")
mr2 <- lme(w~1,data=subset(nhanes2,race==2),random=~1|id,correlation=corAR1(form=~1|id),method="ML")
mr3 <- lme(w~1,data=subset(nhanes2,race==3),random=~1|id,correlation=corAR1(form=~1|id),method="ML")
mr4 <- lme(w~1,data=subset(nhanes2,race==4),random=~1|id,correlation=corAR1(form=~1|id),method="ML")
mr5 <- lme(w~1,data=subset(nhanes2,race==5),random=~1|id,correlation=corAR1(form=~1|id),method="ML")

ma1 <- lme(w~1,data=subset(nhanes2,a==1),random=~1|id,correlation=corAR1(form=~1|id),method="ML")
ma2 <- lme(w~1,data=subset(nhanes2,a==2),random=~1|id,correlation=corAR1(form=~1|id),method="ML")
ma3 <- lme(w~1,data=subset(nhanes2,a==3),random=~1|id,correlation=corAR1(form=~1|id),method="ML")
ma4 <- lme(w~1,data=subset(nhanes2,a==4),random=~1|id,correlation=corAR1(form=~1|id),method="ML")

mb1 <- lme(w~1,data=subset(nhanes2,b==1),random=~1|id,correlation=corAR1(form=~1|id),method="ML")
mb2 <- lme(w~1,data=subset(nhanes2,b==2),random=~1|id,correlation=corAR1(form=~1|id),method="ML")
mb3 <- lme(w~1,data=subset(nhanes2,b==3),random=~1|id,correlation=corAR1(form=~1|id),method="ML")
mb4 <- lme(w~1,data=subset(nhanes2,b==4),random=~1|id,correlation=corAR1(form=~1|id),method="ML")
mb5 <- lme(w~1,data=subset(nhanes2,b==5),random=~1|id,correlation=corAR1(form=~1|id),method="ML")

