setwd("C:/Users/dcries/workspace2")


#-----------------
#do for complete cases
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
nrep <- (nhanes %>% group_by(id) %>% summarise(n=length(id)))$n
meas7 <- subset(nhanes, id %in% unique(id)[nrep==6]) #individuals with all 7 days


load("stanout.RData")
r=as.matrix(rs)
tstar=colMeans(r[,31:3367])
std = apply(r[,31:3367],2,sd)
meas7$tstar <- rep(tstar,each=6)
meas7$std <- rep(std,each=6)
meas7 <- meas7[(!is.na(meas7$waist)) & (!is.na(meas7$bps)) & (!is.na(meas7$bpd)) & (!is.na(meas7$hdl)) & (!is.na(meas7$ldl)) & (!is.na(meas7$glu)) & (!is.na(meas7$tri)) & (!is.na(meas7$education)),] #remove NAs for waist

comp <- meas7[!duplicated(meas7$id),]
write.csv(comp$tstar,file="C:/Users/dcries/github/epi/tstar.csv",row.names=FALSE)
write.csv(comp$std,file="C:/Users/dcries/github/epi/sdT.csv",row.names=FALSE)



#------------
imp1 <- read.csv("../github/epi/NHANES_accel_imp1.csv")
names(imp1) <- tolower(names(imp1))
imp1 <- imp1[!is.na(imp1$education),]
imp1$weekend <- 0
imp1$weekend[imp1$dow %in% c(1,7)] <- 1
imp1$first5 <- 0
imp1$first5[imp1$rep==6] <- 1
imp1$first5[imp1$rep==7] <- 2
imp1$active <- 1
imp1$active[imp1$modvigmin ==0] <- 0
imp1 <- subset(imp1,rep!=7)
m1 <- lm(modvigmin~(weekend)+first5,data=imp1)
wbar <- mean((imp1$modvigmin[imp1$rep <= 5]))
w1 <- imp1$modvigmin
what <- predict(m1)
w <- (1/what)*w1*wbar
imp1$modvigmin2 <- w^.25
nrep <- (imp1 %>% group_by(id) %>% summarise(n=length(id)))$n
meas7 <- subset(imp1, id %in% unique(id)[nrep==6]) #individuals with all 7 days

load("stanout_imp1.RData")
r=as.matrix(rs)
tstar=colMeans(r[,31:7903])
std = apply(r[,31:7903],2,sd)
meas7$tstar <- rep(tstar,each=6)
meas7$std <- rep(std,each=6)
meas7 <- meas7[(!is.na(meas7$waist)) & (!is.na(meas7$bps)) & (!is.na(meas7$bpd)) & (!is.na(meas7$hdl)) & (!is.na(meas7$ldl)) & (!is.na(meas7$glu)) & (!is.na(meas7$tri)) & (!is.na(meas7$education)),] #remove NAs for waist

comp <- meas7[!duplicated(meas7$id),]
write.csv(comp$tstar,file="C:/Users/dcries/github/epi/tstar_imp1.csv",row.names=FALSE)
write.csv(comp$std,file="C:/Users/dcries/github/epi/sdT_imp1.csv",row.names=FALSE)


#------------
imp1 <- read.csv("../github/epi/NHANES_accel_imp2.csv")
names(imp1) <- tolower(names(imp1))
imp1 <- imp1[!is.na(imp1$education),]
imp1$weekend <- 0
imp1$weekend[imp1$dow %in% c(1,7)] <- 1
imp1$first5 <- 0
imp1$first5[imp1$rep==6] <- 1
imp1$first5[imp1$rep==7] <- 2
imp1$active <- 1
imp1$active[imp1$modvigmin ==0] <- 0
imp1 <- subset(imp1,rep!=7)
m1 <- lm(modvigmin~(weekend)+first5,data=imp1)
wbar <- mean((imp1$modvigmin[imp1$rep <= 5]))
w1 <- imp1$modvigmin
what <- predict(m1)
w <- (1/what)*w1*wbar
imp1$modvigmin2 <- w^.25
nrep <- (imp1 %>% group_by(id) %>% summarise(n=length(id)))$n
meas7 <- subset(imp1, id %in% unique(id)[nrep==6]) #individuals with all 7 days

load("stanout_imp2.RData")
r=as.matrix(rs)
tstar=colMeans(r[,31:7903])
std = apply(r[,31:7903],2,sd)
meas7$tstar <- rep(tstar,each=6)
meas7$std <- rep(std,each=6)
meas7 <- meas7[(!is.na(meas7$waist)) & (!is.na(meas7$bps)) & (!is.na(meas7$bpd)) & (!is.na(meas7$hdl)) & (!is.na(meas7$ldl)) & (!is.na(meas7$glu)) & (!is.na(meas7$tri)) & (!is.na(meas7$education)),] #remove NAs for waist

comp <- meas7[!duplicated(meas7$id),]
write.csv(comp$tstar,file="C:/Users/dcries/github/epi/tstar_imp2.csv",row.names=FALSE)
write.csv(comp$std,file="C:/Users/dcries/github/epi/sdT_imp2.csv",row.names=FALSE)

#------------
imp1 <- read.csv("../github/epi/NHANES_accel_imp3.csv")
names(imp1) <- tolower(names(imp1))
imp1 <- imp1[!is.na(imp1$education),]
imp1$weekend <- 0
imp1$weekend[imp1$dow %in% c(1,7)] <- 1
imp1$first5 <- 0
imp1$first5[imp1$rep==6] <- 1
imp1$first5[imp1$rep==7] <- 2
imp1$active <- 1
imp1$active[imp1$modvigmin ==0] <- 0
imp1 <- subset(imp1,rep!=7)
m1 <- lm(modvigmin~(weekend)+first5,data=imp1)
wbar <- mean((imp1$modvigmin[imp1$rep <= 5]))
w1 <- imp1$modvigmin
what <- predict(m1)
w <- (1/what)*w1*wbar
imp1$modvigmin2 <- w^.25
nrep <- (imp1 %>% group_by(id) %>% summarise(n=length(id)))$n
meas7 <- subset(imp1, id %in% unique(id)[nrep==6]) #individuals with all 7 days

load("stanout_imp3.RData")
r=as.matrix(rs)
tstar=colMeans(r[,31:7903])
std = apply(r[,31:7903],2,sd)
meas7$tstar <- rep(tstar,each=6)
meas7$std <- rep(std,each=6)
meas7 <- meas7[(!is.na(meas7$waist)) & (!is.na(meas7$bps)) & (!is.na(meas7$bpd)) & (!is.na(meas7$hdl)) & (!is.na(meas7$ldl)) & (!is.na(meas7$glu)) & (!is.na(meas7$tri)) & (!is.na(meas7$education)),] #remove NAs for waist

comp <- meas7[!duplicated(meas7$id),]
write.csv(comp$tstar,file="C:/Users/dcries/github/epi/tstar_imp3.csv",row.names=FALSE)
write.csv(comp$std,file="C:/Users/dcries/github/epi/sdT_imp3.csv",row.names=FALSE)


#------------
imp1 <- read.csv("../github/epi/NHANES_accel_imp4.csv")
names(imp1) <- tolower(names(imp1))
imp1 <- imp1[!is.na(imp1$education),]
imp1$weekend <- 0
imp1$weekend[imp1$dow %in% c(1,7)] <- 1
imp1$first5 <- 0
imp1$first5[imp1$rep==6] <- 1
imp1$first5[imp1$rep==7] <- 2
imp1$active <- 1
imp1$active[imp1$modvigmin ==0] <- 0
imp1 <- subset(imp1,rep!=7)
m1 <- lm(modvigmin~(weekend)+first5,data=imp1)
wbar <- mean((imp1$modvigmin[imp1$rep <= 5]))
w1 <- imp1$modvigmin
what <- predict(m1)
w <- (1/what)*w1*wbar
imp1$modvigmin2 <- w^.25
nrep <- (imp1 %>% group_by(id) %>% summarise(n=length(id)))$n
meas7 <- subset(imp1, id %in% unique(id)[nrep==6]) #individuals with all 7 days

load("stanout_imp4.RData")
r=as.matrix(rs)
tstar=colMeans(r[,31:7903])
std = apply(r[,31:7903],2,sd)
meas7$tstar <- rep(tstar,each=6)
meas7$std <- rep(std,each=6)
meas7 <- meas7[(!is.na(meas7$waist)) & (!is.na(meas7$bps)) & (!is.na(meas7$bpd)) & (!is.na(meas7$hdl)) & (!is.na(meas7$ldl)) & (!is.na(meas7$glu)) & (!is.na(meas7$tri)) & (!is.na(meas7$education)),] #remove NAs for waist

comp <- meas7[!duplicated(meas7$id),]
write.csv(comp$tstar,file="C:/Users/dcries/github/epi/tstar_imp4.csv",row.names=FALSE)
write.csv(comp$std,file="C:/Users/dcries/github/epi/sdT_imp4.csv",row.names=FALSE)


#------------
imp1 <- read.csv("../github/epi/NHANES_accel_imp5.csv")
names(imp1) <- tolower(names(imp1))
imp1 <- imp1[!is.na(imp1$education),]
imp1$weekend <- 0
imp1$weekend[imp1$dow %in% c(1,7)] <- 1
imp1$first5 <- 0
imp1$first5[imp1$rep==6] <- 1
imp1$first5[imp1$rep==7] <- 2
imp1$active <- 1
imp1$active[imp1$modvigmin ==0] <- 0
imp1 <- subset(imp1,rep!=7)
m1 <- lm(modvigmin~(weekend)+first5,data=imp1)
wbar <- mean((imp1$modvigmin[imp1$rep <= 5]))
w1 <- imp1$modvigmin
what <- predict(m1)
w <- (1/what)*w1*wbar
imp1$modvigmin2 <- w^.25
nrep <- (imp1 %>% group_by(id) %>% summarise(n=length(id)))$n
meas7 <- subset(imp1, id %in% unique(id)[nrep==6]) #individuals with all 7 days

load("stanout_imp5.RData")
r=as.matrix(rs)
tstar=colMeans(r[,31:7903])
std = apply(r[,31:7903],2,sd)
meas7$tstar <- rep(tstar,each=6)
meas7$std <- rep(std,each=6)
meas7 <- meas7[(!is.na(meas7$waist)) & (!is.na(meas7$bps)) & (!is.na(meas7$bpd)) & (!is.na(meas7$hdl)) & (!is.na(meas7$ldl)) & (!is.na(meas7$glu)) & (!is.na(meas7$tri)) & (!is.na(meas7$education)),] #remove NAs for waist

comp <- meas7[!duplicated(meas7$id),]
write.csv(comp$tstar,file="C:/Users/dcries/github/epi/tstar_imp5.csv",row.names=FALSE)
write.csv(comp$std,file="C:/Users/dcries/github/epi/sdT_imp5.csv",row.names=FALSE)
