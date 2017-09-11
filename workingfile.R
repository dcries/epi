library(ggplot2)
library(lme4)
library(dplyr)
library(reshape)
library(car)
library(nlme)

setwd("C:\\Users\\dcries\\github\\epi")
nhanes <- read.csv("NHANES_complete.csv")

names(nhanes) <- tolower(names(nhanes))
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

d = nhanes %>% group_by(id) %>% summarise(m=mean(w),glu=glu[1],waist=waist[1],tri=tri[1],bps=bps[1],bpd=bpd[1],hdl=hdl[1],ldl=ldl[1])

#-------------

a=nhanes %>% group_by(id) %>% summarise(m=mean(modvigmin^.25),
        g=sex[1],r=race[1],age=age[1],n=length(sex),s=(n-1)*var(modvigmin^.25))

a1 <- a[a$n>1,]
a1$a <- 1;a1$a[a1$age>=35] <- 2;a1$a[a1$age>=50] <- 3;a1$a[a1$age>=65] <- 4
# sum(a1$s[a1$g==1])/sum(a1$g==1);sum(a1$s[a1$g==2])/sum(a1$g==2)
# sum(a1$s[a1$r==1])/sum(a1$r==1);sum(a1$s[a1$r==2])/sum(a1$r==2);sum(a1$s[a1$r==3])/sum(a1$r==3);sum(a1$s[a1$r==4])/sum(a1$r==4);sum(a1$s[a1$r==5])/sum(a1$r==5)
# sum(a1$s[a1$age %in% c(18:34)])/sum(a1$age %in% c(18:34));sum(a1$s[a1$age %in% c(35:50)])/sum(a1$age %in% c(35:50));sum(a1$s[a1$age %in% c(51:65)])/sum(a1$age %in% c(51:65));sum(a1$s[a1$age %in% c(66:85)])/sum(a1$age %in% c(66:85));

a1 %>% group_by(g) %>% summarise(m=mean(s))
a1 %>% group_by(r) %>% summarise(m=mean(s))
a1 %>% group_by(a) %>% summarise(m=mean(s))
a1 %>% group_by(g,r) %>% summarise(m=mean(s))
a1 %>% group_by(g,a) %>% summarise(m=mean(s))
a1 %>% group_by(a,r) %>% summarise(m=mean(s))

a2 <- a1 %>% group_by(age) %>% summarise(m=mean(s))
m2 <- lm(m~age,data=a2)
plot(predict(m2),resid(m2));abline(a=0,b=0)
plot(a2$age,a2$m);lines(a2$age,predict(m2),col="red")
m3 <- lm(m~age+I(age^2),data=a2)
anova(m3,m2)
plot(predict(m3),resid(m3));abline(a=0,b=0)
plot(a2$age,a2$m);lines(a2$age,predict(m3),col="red")

#different values of rho
m1c <- lme(w~1,data=nhanes,random=~1|id,correlation=corAR1(form=~1|id),method="ML")

mm <- lme(w~1,data=subset(nhanes,sex==1),random=~1|id,correlation=corAR1(form=~1|id),method="ML")
mf <- lme(w~1,data=subset(nhanes,sex==2),random=~1|id,correlation=corAR1(form=~1|id),method="ML")

