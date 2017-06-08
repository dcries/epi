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

sum(nhanes$modvigmetsb==0)/nrow(nhanes)
sum(nhanes$modvigmin==0)/nrow(nhanes) #not necessarily in a bout, although we have that data
sum(nhanes$vigmin==0)/nrow(nhanes)
sum(nhanes$lightmin==0)/nrow(nhanes)

#no diff in modvigmin by day really except 7th
nhanes %>% group_by(rep) %>% summarise(m=mean(modvigmin),s=sd(modvigmin)/sqrt(length(modvigmin)))
#no diff in modmin by day really except 7th
nhanes %>% group_by(rep) %>% summarise(m=mean(modmin),s=sd(modmin)/sqrt(length(modmin)))
#no diff in modmin by day really except 7th and 2nd
nhanes %>% group_by(rep) %>% summarise(m=mean(lightmin),s=sd(lightmin)/sqrt(length(lightmin)))
#mondays most vig, fri-sun least, tues-thurs similar
nhanes %>% group_by(dow) %>% summarise(m=mean(vigmin),s=sd(vigmin)/sqrt(length(lightmin)))
#weekends lower
nhanes %>% group_by(dow) %>% summarise(m=mean(modvigmin),s=sd(modvigmin)/sqrt(length(modvigmin)))
#weekends lower
nhanes %>% group_by(dow) %>% summarise(m=mean(modmin),s=sd(modmin)/sqrt(length(modmin)))
#weekends and thursday lower, friday higher
nhanes %>% group_by(dow) %>% summarise(m=mean(lightmin),s=sd(lightmin)/sqrt(length(lightmin)))
#first day vig
nhanes %>% group_by(rep) %>% summarise(m=mean(vigmin),s=sd(vigmin)/sqrt(length(lightmin)))

qplot(data=nhanes,x=as.factor(rep),y=modvigmin,geom="boxplot")
anova(lm(modvigmin~as.factor(rep),data=subset(nhanes,rep<6))) #first 5 days are similar, last two

qplot(data=nhanes,x=as.factor(dow),y=modvigmin,geom="boxplot")
anova(lm(modvigmin~as.factor(dow),data=subset(nhanes,dow%in%c(2:6)))) #M-F similar

#test for weekend effects
anova(lm(modvigmin ~ as.factor(weekend),data=nhanes))
anova(lm(modmin ~ as.factor(weekend),data=nhanes))
anova(lm(lightmin ~ as.factor(weekend),data=nhanes))
anova(lm(vigmin ~ as.factor(weekend),data=nhanes))

#cor between wear time and exercise time
cor(nhanes$wearmin,nhanes$lightmin);plot(nhanes$wearmin,nhanes$lightmin)
cor(nhanes$wearmin,nhanes$modmin);plot(nhanes$wearmin,nhanes$modmin)
cor(nhanes$wearmin,nhanes$modvigmin);plot(nhanes$wearmin,nhanes$modvigmin)
cor(nhanes$wearmin,nhanes$vigmin);plot(nhanes$wearmin,nhanes$vigmin)


nrep <- (nhanes %>% group_by(id) %>% summarise(n=length(id)))$n

m1 <- lmer(modvigmin~(1|id),data=nhanes,REML=FALSE)
#m1b <- gls(modvigmin~1,data=nhanes,correlation=corAR1(form=~1|id),method="ML")
m1c <- lme(modvigmin~1,data=nhanes,random=~1|id,correlation=corAR1(form=~1|id),method="ML")

re <- rep(unlist(ranef(m1)),nrep)
nhanes$error <- nhanes$modvigmin - re - fixef(m1)

id = unique(nhanes$id)
acf(nhanes$error[nhanes$id==id[1]])

eht <- nhanes %>% group_by(id) %>% summarise(e1=error[1],e2=error[2],e3=error[3],e4=error[4],e5=error[5],e6=error[6],e7=error[7])
#eht <- nhanes %>% group_by(id) %>% summarise(e1=modmin[1],e2=modmin[2],e3=modmin[3],e4=modmin[4],e5=modmin[5],e6=modmin[6],e7=modmin[7])

ind <- which(!is.na(eht$e7))
boxtest <- rep(0,length(ind))
for(i in 1:length(ind)){
  boxtest[i] <- Box.test(unlist(eht[ind[i],-1]),lag=2,type="Ljung-Box")$p.value
}

#lag 1 difference
o1a = c(eht$e1[nrep > 1],eht$e2[nrep > 2],eht$e3[nrep > 3],eht$e4[nrep > 4],eht$e5[nrep > 5],eht$e6[nrep > 6])
o1b = c(eht$e2[nrep > 1],eht$e3[nrep > 2],eht$e4[nrep > 3],eht$e5[nrep > 4],eht$e6[nrep > 5],eht$e7[nrep > 6])
cor(o1a,o1b);cor.test(o1a,o1b)
plot(o1a,o1b)
#lag 2 diff
o2a = c(eht$e1[nrep > 2],eht$e2[nrep > 3],eht$e3[nrep > 4],eht$e4[nrep > 5],eht$e5[nrep > 6])
o2b = c(eht$e3[nrep > 2],eht$e4[nrep > 3],eht$e5[nrep > 4],eht$e6[nrep > 5],eht$e7[nrep > 6])
cor(o2a,o2b);cor.test(o2a,o2b)
plot(o2a,o2b)
#lag 3 diff
o3a = c(eht$e1[nrep > 3],eht$e2[nrep > 4],eht$e3[nrep > 5],eht$e4[nrep > 6])
o3b = c(eht$e4[nrep > 3],eht$e5[nrep > 4],eht$e6[nrep > 5],eht$e7[nrep > 6])
cor(o3a,o3b);cor.test(o3a,o3b)
plot(o3a,o3b)
#lag 4 diff
o4a = c(eht$e1[nrep > 4],eht$e2[nrep > 5],eht$e3[nrep > 6])
o4b = c(eht$e5[nrep > 4],eht$e6[nrep > 5],eht$e7[nrep > 6])
cor(o4a,o4b);cor.test(o4a,o4b)
plot(o4a,o4b)
#lag 5 diff
o5a = c(eht$e1[nrep > 5],eht$e2[nrep > 6])
o5b = c(eht$e6[nrep > 5],eht$e7[nrep > 6])
cor(o5a,o5b);cor.test(o5a,o5b)
plot(o5a,o5b)

#test if error larger on day 1 compared to others
t.test(eht$e1[nrep > 1],eht$e2[nrep > 1])
all7 <- eht[complete.cases(eht),]
mall7 <- melt(data.frame(all7[,-1]))
m2 = lm(value~as.factor(variable),data=mall7)
summary(m2)

#-------------------------------------------
#look at relationships among mets risk factors
y <- nhanes %>% filter(!duplicated(id))
pairs(y[,c("glu","waist","tri","ldl","hdl","bps","bpd","modvigmin","lightmin","vigmin","predv02","estv02","fitlev")])

#correlation plot of variables
cormat <- cor(y[y$vigmin>0,c("glu","waist","tri","ldl","hdl","bps","bpd","modvigmin","lightmin","vigmin","predv02","estv02","fitlev")],use="pairwise.complete.obs")
cormat[lower.tri(cormat)] <- 0
diag(cormat) <- 0
df <- expand.grid(x=1:ncol(cormat),y=1:ncol(cormat))
df$covar <- c(cormat)
df$x <- factor(df$x,labels=c("glu","waist","tri","ldl","hdl","bps","bpd","modvigmin","lightmin","vigmin","predv02","estv02","fitlev"))
df$y <- factor(df$y,labels=c("glu","waist","tri","ldl","hdl","bps","bpd","modvigmin","lightmin","vigmin","predv02","estv02","fitlev"))
df$corr <- as.character(c(round(cormat,3)))
df$corr[df$corr=="0"] <- ""

ggplot(data=df,aes(x=as.factor(x),y=as.factor(y))) + geom_tile(aes(fill=covar)) + 
  geom_text(aes(label=corr)) +
  #geom_text(aes(label=count)) + scale_fill_gradient(low = "white", high = "red") 
  scale_fill_gradient2(high="blue",low="red",mid="white")

#-------------------------
#multivariate regression on mean observed modvig, light
a <- nhanes %>% group_by(id) %>% summarise(modvig=mean(modvigmin),light=mean(lightmin))
y$modvig <- a$modvig
y$light <- a$light

m3 <- lm(with(y,cbind(glu,waist,tri,hdl,bps,bpd))~y$modvig)
coef(m3)

m4 <- lm(with(y,cbind(glu,waist,tri,hdl,bps,bpd))~y$light)
coef(m4)

m5 <- lm(with(y,cbind(glu,waist,tri,hdl,bps,bpd))~y$modvig+y$estv02)
coef(m5)


#inidivual distributions
qplot(log(y$glu))
qplot(log(y$tri))
qqnorm(log(y$tri));qqline(log(y$tri))
qqnorm(log(y$hdl));qqline(log(y$hdl))
qqnorm(sqrt(y$ldl));qqline(sqrt(y$ldl))
qplot(y$bpd)
qqnorm((y$bpd));qqline((y$bpd))


#---------------------------------
library(rstan)
model = "
data{
  int<lower = 0> N; //number of individuals
  int<lower = 0> K; //number of obs per individual

  vector[K] Y[N];
  vector[N] D;

  cov_matrix[K] R;
}

parameters{
  vector[K] a;
  vector[K] b;

  cov_matrix[K] Omega;
}

transformed parameters{
  vector[K] mu[N];

for(k in 1:K){
  for(i in 1:N){
    mu[i, k] <- a[k] + b[k] * D[i];
    }
  }
}

model{
  for(i in 1:N){
    Y[i] ~ multi_normal(mu[i], Omega);
  }

  Omega ~ inv_wishart(K+1, R);
}
"
model = stan_model(model_code=model)
results = sampling(model,data=dataList)
r=as.matrix(results)

