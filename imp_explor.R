

setwd("C:\\Users\\dcries\\github\\epi")
imp1 <- read.csv("NHANES_accel_imp1.csv")
imp2 <- read.csv("NHANES_accel_imp2.csv")
imp3 <- read.csv("NHANES_accel_imp3.csv")
imp4 <- read.csv("NHANES_accel_imp4.csv")
imp5 <- read.csv("NHANES_accel_imp5.csv")
nhanes <- read.csv("nhanes_complete.csv")

#for impx
imp1$weekend <- 0
imp1$weekend[imp1$dow %in% c(1,7)] <- 1
imp1$first5 <- 0
imp1$first5[imp1$rep==6] <- 1
imp1$first5[imp1$rep==7] <- 2

imp1 <- subset(imp1,rep!=7)
m1 <- lm(modvigmin^.25~(weekend)+first5,data=imp1)
wbar <- mean((imp1$modvigmin[imp1$rep <= 5])^.25)
w1 <- imp1$modvigmin^.25
what <- predict(m1)
w <- (1/what)*w1*wbar
imp1$modvigmin2 <- w

indlevel2 <- imp1 %>% group_by(id) %>% summarise(m=mean(modvigmin^(1/4)),m2=mean(modvigmin2),s=sd(modvigmin^(1/4)),glu=glu[1],waist=waist[1],ldl=ldl[1],hdl=hdl[1],bps=bps[1],bpd=bpd[1],tri=tri[1],gender=sex[1],race=race[1],age=age[1],n=length(id))
qplot(data=indlevel2,x=m2,y=log(glu)) + geom_smooth()
qplot(data=indlevel2,x=m2,y=log(tri)) + geom_smooth()
qplot(data=indlevel2,x=m2,y=waist) + geom_smooth()
qplot(data=indlevel2,x=m2,y=ldl) + geom_smooth()
qplot(data=indlevel2,x=m2,y=hdl) + geom_smooth()
qplot(data=indlevel2,x=m2,y=bps) + geom_smooth()
qplot(data=indlevel2,x=m2,y=bpd) + geom_smooth()

mglu2 <- nls(log(glu)~L2-L/(1+exp(-k*(m-x0))),start=list(L=.5,k=3,x0=1.2,L2=log(120)),data=indlevel2)
mwaist2 <- nls(waist~L2-L/(1+exp(-k*(m-x0))),start=list(L=20,k=3,x0=1.2,L2=110),data=indlevel2)
mtri2 <- nls(log(tri)~L2-L/(1+exp(-k*(m-x0))),start=list(L=.5,k=3,x0=1.5,L2=5),data=indlevel2)
mbps2 <- nls(bps~L2-L/(1+exp(-k*(m-x0))),start=list(L=20,k=5,x0=1.2,L2=140),data=indlevel2)

df=data.frame(x=seq(from=0,to=4,by=0.1))
df$bps=coef(mbps2)[4]-coef(mbps2)[1]/(1+exp(-coef(mbps2)[2]*(df$x-coef(mbps2)[3])))
qplot(data=indlevel,x=m,y=bps) + geom_smooth() + geom_line(data=df,aes(x=x,y=bps),col="red") 
df$glu=coef(mglu)["L2"]-coef(mglu)["L"]/(1+exp(-coef(mglu)["k"]*(df$x-coef(mglu)["x0"])))
qplot(data=indlevel,x=m,y=log(glu)) + geom_smooth() + geom_line(data=df,aes(x=x,y=glu),col="red") 
df$waist=coef(mwaist)["L2"]-coef(mwaist)["L"]/(1+exp(-coef(mwaist)["k"]*(df$x-coef(mwaist)["x0"])))
qplot(data=indlevel,x=m,y=waist) + geom_smooth() + geom_line(data=df,aes(x=x,y=waist),col="red") 
df$tri=coef(mtri)["L2"]-coef(mtri)["L"]/(1+exp(-coef(mtri)["k"]*(df$x-coef(mtri)["x0"])))
qplot(data=indlevel,x=m,y=log(tri)) + geom_smooth() + geom_line(data=df,aes(x=x,y=tri),col="red") #+ ylim(c(0,750))
