library(rstan)
library(dplyr)
setwd("/home/dcries")
imp1 <- read.csv("epi/NHANES_accel_imp1.csv")
#nhanes <- read.csv("NHANES_complete.csv")
names(imp1) <- tolower(names(imp1))
imp1$weekend <- 0
imp1$weekend[imp1$dow %in% c(1,7)] <- 1
imp1$first5 <- 0
imp1$first5[imp1$rep==6] <- 1
imp1$first5[imp1$rep==7] <- 2

imp1$active <- 1
imp1$active[imp1$modvigmin ==0] <- 0

imp1 <- subset(imp1,rep!=7)
m1 <- lm(modvigmin^.25~(weekend)+first5,data=imp1)
wbar <- mean((imp1$modvigmin[imp1$rep <= 5])^.25)
w1 <- imp1$modvigmin^.25
what <- predict(m1)
w <- (1/what)*w1*wbar
imp1$modvigmin2 <- w

nrep <- (imp1 %>% group_by(id) %>% summarise(n=length(id)))$n
meas7 <- subset(imp1, id %in% unique(id)[nrep==6]) #individuals with all 7 days



meas7 <- meas7[(!is.na(meas7$waist)) & (!is.na(meas7$bps)) & (!is.na(meas7$bpd)) & (!is.na(meas7$hdl)) & (!is.na(meas7$ldl)) & (!is.na(meas7$glu)) & (!is.na(meas7$tri)),] #remove NAs for waist

#Tstar <- ((meas7 %>% group_by(id) %>% summarise(m=mean(modvigmin)))$m)^.25
Tstar <- ((meas7 %>% group_by(id) %>% summarise(m=mean(modvigmin2)))$m)

models <- "
data{
  int<lower=0> N; //number of individuals
  vector[N] Tstar; // usual ^ 1/4
  //vector[N] sdT;//

  real waist[N];
  real lglu[N];
  real ltri[N];
  real bps[N];
  real ldl[N];
  real hdl[N];
  real bpd[N];
  //real<lower=0> age[N];
  //int<lower=0> gender[N];

}
transformed data{
  //vector[2] zeros;
  //zeros = rep_vector(0,2);
}
parameters{
  vector[4] alphaw;
  vector[4] alphag;
  vector[4] alphat;
  vector[4] alphabs;
  vector[3] alphal;
  vector[3] alphabd;
  vector[2] alphah;


  real<lower=0> sigmawaist;
  real<lower=0> sigmaglu;
  real<lower=0> sigmatri;
  real<lower=0> sigmabps;
  real<lower=0> sigmabpd;
  real<lower=0> sigmaldl;
  real<lower=0> sigmahdl;

  //vector[N] appx;
}


model{

  vector[N] muwaist;
  vector[N] mulglu;
  vector[N] multri;
  vector[N] mubps;
  vector[N] mubpd;
  vector[N] muhdl;
  vector[N] muldl;

  //vector[N] mu[7]; //7 for number of regressions

for(i in 1:N){
    muwaist[i] = alphaw[4]-alphaw[1]/(1+exp(-alphaw[2]*(Tstar[i]-alphaw[3])));
    mulglu[i] = alphag[4]-alphag[1]/(1+exp(-alphag[2]*(Tstar[i]-alphag[3])));
    multri[i] = alphat[4]-alphat[1]/(1+exp(-alphat[2]*(Tstar[i]-alphat[3])));
    mubps[i] = alphabs[4]-alphabs[1]/(1+exp(-alphabs[2]*(Tstar[i]-alphabs[3])));
    muldl[i] = alphal[1] + alphal[2]*Tstar[i] + alphal[3]*pow(Tstar[i],2);
    mubpd[i] = alphabd[1] + alphabd[2]*Tstar[i] + alphabd[3]*pow(Tstar[i],2);
    muhdl[i] = alphah[1] + alphah[2]*Tstar[i];

  waist[i] ~ normal(muwaist[i],sigmawaist);
  lglu[i] ~ normal(mulglu[i],sigmaglu);
  ltri[i] ~ normal(multri[i],sigmatri);
  bps[i] ~ normal(mubps[i],sigmabps);
  ldl[i] ~ normal(muldl[i],sigmaldl);
  hdl[i] ~ normal(muhdl[i],sigmahdl);
  bpd[i] ~ normal(mubpd[i],sigmabpd);

}

  //muwaist = alphaw[4]-alphaw[1]/(1+exp(-alphaw[2]*((Tstar+appx-alphaw[3])));
  //mulglu = alphag[4]-alphag[1]/(1+exp(-alphag[2]*(Tstar+appx-alphag[3])));
  //multri = alphat[4]-alphat[1]/(1+exp(-alphat[2]*(Tstar+appx-alphat[3])));
  //mubps = alphabs[4]-alphabs[1]/(1+exp(-alphabs[2]*(Tstar+appx-alphabs[3])));
  //muldl = alphal[1] + alphal[2]*(Tstar+appx) + alphal[3]*pow(Tstar+appx,2);
  //mubpd = alphabd[1] + alphabd[2]*(Tstar+appx) + alphabd[3]*pow(Tstar+appx,2);
  //muhdl = alphah[1] + alphah[2]*(Tstar+appx);

  //appx ~ normal(0,sdT);

  sigmawaist ~ cauchy(0,1);
  sigmaglu ~ cauchy(0,1);
  sigmatri ~ cauchy(0,1);
  sigmabps ~ cauchy(0,1);
  sigmabpd ~ cauchy(0,1);
  sigmaldl ~ cauchy(0,1);
  sigmahdl ~ cauchy(0,1);


  alphaw ~ normal(0,100);
  alphag ~ normal(0,100);
  alphat ~ normal(0,100);
  alphabs ~ normal(0,100);
  alphal ~ normal(0,100);
  alphabd ~ normal(0,100);
  alphah ~ normal(0,100);
}
"

waist <- meas7$waist[!duplicated(meas7$id)]
lglu <- log(meas7$glu[!duplicated(meas7$id)])
ltri <- log(meas7$tri[!duplicated(meas7$id)])
bps <- (meas7$bps[!duplicated(meas7$id)])
ldl <- meas7$ldl[!duplicated(meas7$id)]
hdl <- meas7$hdl[!duplicated(meas7$id)]
bpd <- meas7$bpd[!duplicated(meas7$id)]

dat=list(N      = length(unique(meas7$id)),
         #age= meas7$age[!duplicated(meas7$id)],
         #gender= meas7$sex[!duplicated(meas7$id)],nu=3,D=diag(2),
         waist=waist,lglu=lglu,ltri=ltri,bps=bps,ldl=ldl,
         hdl=hdl,bpd=bpd,Tstar=Tstar
)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

start1 <- list(alphag=c(0.1642,3.4081,1.4433,4.7228),
  alphaw=c(10.445,   3.230,   2.033, 101.517 ),
alphat=c(0.2805, 4.4733, 1.8297, 4.9261 ),
alphabs=c( 18.388,   4.602 ,  1.389, 137.256 ))
start2 <- list(alphag=c(0.1242,4.4081,1.1433,3.7228),
  alphaw=c(8.445,   2.230,   1.033, 80.517 ),
alphat=c(0.1805, 3.4733, 1.1297, 3.9261 ),
alphabs=c( 12.388,   3.602 ,  1.189, 120.256 ))
start3 <- list(alphag=c(0.3642,5.4081,1.9433,5.7228),
  alphaw=c(11.445,   4.230,   2.933, 111.517 ),
alphat=c(0.3805, 5.4733, 2.8297,5.9261 ),
alphabs=c( 19.388,   5.602 ,  2.389, 147.256 ))
start4 <- list(alphag=c(0.2642,2.4081,2.4433,5.7228),
               alphaw=c(9.445,   4.230,   1.033, 121.517 ),
               alphat=c(0.2205, 6.4733, 2.8297, 5.9261 ),
               alphabs=c( 28.388,   3.602 ,  2.389, 157.256 ))

ms <- stan_model(model_code=models)
rs <- sampling(ms,dat,c("alphaw",
                        "alphag","alphat","alphal",
                        "alphabs","alphabd",
                        "alphah",
                        "sigmawaist","sigmabps",
                        "sigmaglu","sigmatri","sigmaldl",
                        "sigmahdl","sigmabpd"
),
iter=2000,chains=4,init=list(start1,start2,start3,start4),
control=list(adapt_delta=0.99,max_treedepth=15))
#summary(rs)
save(rs,file="stanout_2nd_imp1.RData")
