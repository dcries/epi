library(rstan)
library(dplyr)
library(reshape)

#setwd("C:\\Users\\dcries\\github\\epi")
setwd("/home/dcries")
nhanes <- read.csv("epi/nhanes_complete.csv")
#nhanes <- read.csv("NHANES_complete.csv")
names(nhanes) <- tolower(names(nhanes))
nhanes$weekend <- 0
nhanes$weekend[nhanes$dow %in% c(1,7)] <- 1
nhanes$first5 <- 0
nhanes$first5[nhanes$rep==6] <- 1
nhanes$first5[nhanes$rep==7] <- 2

nhanes$active <- 1
nhanes$active[nhanes$modvigmin ==0] <- 0

nhanes <- subset(nhanes,rep!=7)
m1 <- lm(modvigmin^.25~(weekend)+first5,data=nhanes)
wbar <- mean((nhanes$modvigmin[nhanes$rep <= 5])^.25)
w1 <- nhanes$modvigmin^.25
what <- predict(m1)
w <- (1/what)*w1*wbar
nhanes$modvigmin2 <- w

nrep <- (nhanes %>% group_by(id) %>% summarise(n=length(id)))$n
meas7 <- subset(nhanes, id %in% unique(id)[nrep==6]) #individuals with all 7 days


meas7 <- meas7[(!is.na(meas7$waist)) & (!is.na(meas7$bps)) & (!is.na(meas7$bpd)) & (!is.na(meas7$hdl)) & (!is.na(meas7$ldl)) & (!is.na(meas7$glu)) & (!is.na(meas7$tri)),] #remove NAs for waist

models <- "
data{
  int<lower=0> N; //number of individuals
  int<lower=0> k; // number of obs per individual
  int<lower=0> n2; //total number of non-zero obs
  //real<lower=0> y[N];
  vector[k] y[N];
  vector[n2] y2;
  int pk; //number of covariates plus intercept
  matrix[N,pk] X;
  real waist[N];
  //real lglu[N];
  //real ltri[N];
  real bps[N];
  //real ldl[N];
  real hdl[N];
  real bpd[N];
  //real<lower=0> age[N];
  //int<lower=0> gender[N];
  int<lower=0> numnonzeros[N]; //number of nonzero minutes days for each individual
  matrix[N,k] nonzeropos; //position of nonzero minutes for each indivudal
  real nu;
  matrix[2,2] D;
}
transformed data{
  vector[2] zeros;
  zeros = rep_vector(0,2);
}
parameters{
  vector[pk] gamma;
//    real gamma1;
//real gamma2;
  vector[pk] beta;
  vector[4] alphaw;
  //vector[4] alphag;
  //vector[4] alphat;
  //vector[4] alphabs;
  //vector[3] alphal;
  //vector[3] alphabd;
  vector[2] alphah;

//real beta1;
//real beta2;
  cov_matrix[2] Sigma;
  real<lower=0> sigmae;
  real<lower=-1,upper=1> rho;
  matrix[N,2] b;
  real<lower=0> sigmawaist;
  //real<lower=0> sigmaglu;
  //real<lower=0> sigmatri;
  //real<lower=0> sigmabps;
  //real<lower=0> sigmabpd;
  //real<lower=0> sigmaldl;
  real<lower=0> sigmahdl;

}
transformed parameters{
  cov_matrix[k] ar1mat[N];

for(i in 1:N){
  for (m in 1:k){
    ar1mat[i,m,m] = pow(sigmae,2.0);
  }
}
//respecify ar1mat, nonzeropos, numnonzeros
for(i in 1:N){
  for (m in 1:(k-1)) {
    for (n in (m+1):k) {
      ar1mat[i,m,n] = pow(sigmae,2.0) * if_else(n<=numnonzeros[i],pow(rho,nonzeropos[i,n]-nonzeropos[i,m]),0);
      ar1mat[i,n,m] = ar1mat[i,m,n];
    }
  }
}

}
model{

//matrix[N,k] mu;
  vector[N] T; //usual
  vector[N] mu;
  vector[N] Tstar; // usual ^ 1/4
  vector[N] p;
  vector[N] muwaist;
  //vector[N] mulglu;
  //vector[N] multri;
  //vector[N] mubps;
  //vector[N] mubpd;
  vector[N] muhdl;
  //vector[N] muldl;


  int pos;
  pos = 1;

  for(i in 1:N){
    p[i] = Phi(X[i,]*gamma+b[i,1]);
    mu[i] = X[i,]*beta + b[i,2];
    T[i] = pow(mu[i],4.0) + 6*pow(sigmae,2.0)*pow(mu[i],2.0);
    Tstar[i] = p[i]*pow(T[i],0.25);

    muwaist[i] = alphaw[4]-alphaw[1]/(1+exp(-alphaw[2]*(Tstar[i]-alphaw[3])));
    //mulglu[i] = alphag[4]-alphag[1]/(1+exp(-alphag[2]*(Tstar[i]-alphag[3])));
    //multri[i] = alphat[4]-alphat[1]/(1+exp(-alphat[2]*(Tstar[i]-alphat[3])));
    //mubps[i] = alphabs[4]-alphabs[1]/(1+exp(-alphabs[2]*(Tstar[i]-alphabs[3])));
    //muldl[i] = alphal[1] + alphal[2]*Tstar[i] + alphal[3]*pow(Tstar[i],2);
    //mubpd[i] = alphabd[1] + alphabd[2]*Tstar[i] + alphabd[3]*pow(Tstar[i],2);
    muhdl[i] = alphah[1] + alphah[2]*Tstar[i];

    waist[i] ~ normal(muwaist[i],sigmawaist);
    //lglu[i] ~ normal(mulglu[i],sigmaglu);
    //ltri[i] ~ normal(multri[i],sigmatri);
    //bps[i] ~ normal(mubps[i],sigmabps);
    //ldl[i] ~ normal(muldl[i],sigmaldl);
    hdl[i] ~ normal(muhdl[i],sigmahdl);
    //bpd[i] ~ normal(mubpd[i],sigmabpd);


    b[i] ~ multi_normal(zeros,Sigma);

    for(j in 1:k){
      if(y[i][j] == 0){
        //increment_log_prob(bernoulli_log(0,p[i]));
        target += bernoulli_lpmf(0|p[i]);
      }
      else{
        //increment_log_prob(bernoulli_log(1,p[i]));
        target += bernoulli_lpmf(1|p[i]);// + normal_lpdf(y[i][j]|mu[i],sigmae);
      }
    }
    if(numnonzeros[i]>0){
      target += multi_normal_lpdf(segment(y2,pos,numnonzeros[i])|rep_vector(mu[i],numnonzeros[i]),block(ar1mat[i],1,1,numnonzeros[i],numnonzeros[i]));
      pos = pos + numnonzeros[i];
    }
  }




//  for(j in 1:k){
//    for(i in 1:N){
//      mu[i,j] = beta0 + beta1*age[i] + beta2*gender[i] + b[i,2];
//    }
//  }


//  rho ~ uniform(-1,1);
  sigmae ~ cauchy(0,1);
  sigmawaist ~ cauchy(0,1);
  //sigmaglu ~ cauchy(0,1);
  //sigmatri ~ cauchy(0,1);
  //sigmabps ~ cauchy(0,1);
  //sigmabpd ~ cauchy(0,1);
  //sigmaldl ~ cauchy(0,1);
  sigmahdl ~ cauchy(0,1);

 beta ~ normal(0,100); //implies independent priors for beta ?
//  beta1 ~ normal(0,100);
//  beta2 ~ normal(0,100);
  gamma ~ normal(0,100);
//  gamma1 ~ normal(0,100);
//  gamma2 ~ normal(0,100);
  alphaw ~ normal(0,100);
  //alphag ~ normal(0,100);
  //alphat ~ normal(0,100);
  //alphabs ~ normal(0,100);
  //alphal ~ normal(0,100);
  //alphabd ~ normal(0,100);
  alphah ~ normal(0,100);
  Sigma ~ inv_wishart(nu,D);
}
"

ym <- melt(meas7[,c("modvigmin2","rep","id")],id.vars=c("id","rep"))
yc <- cast(ym,id+variable~rep)
nonzeros <- (meas7 %>% group_by(id) %>% summarise(n=sum(active)))$n
a1=meas7$rep[meas7$active==1]
a2=cumsum(nonzeros)

nonzeropos <- matrix(0,ncol=length(unique(meas7$id)),nrow=6)
nonzeropos[,1] <- a1[1:nonzeros[1]]
for(i in 2:ncol(nonzeropos)){
    temp <- a1[(a2[i-1]+1):a2[i]]
    if(length(temp)<6){
      n0 <- 6-length(temp)
      temp <- c(temp,rep(0,n0))
    }
    nonzeropos[,i] <- temp
}

x <- model.matrix(~age+as.factor(sex)+as.factor(race),data=meas7[!duplicated(meas7$id),c("age","sex","race")])

waist <- meas7$waist[!duplicated(meas7$id)]
lglu <- log(meas7$glu[!duplicated(meas7$id)])
ltri <- log(meas7$tri[!duplicated(meas7$id)])
bps <- (meas7$bps[!duplicated(meas7$id)])
ldl <- meas7$ldl[!duplicated(meas7$id)]
hdl <- meas7$hdl[!duplicated(meas7$id)]
bpd <- meas7$bpd[!duplicated(meas7$id)]


dat=list(y=(yc[,3:8]),  N      = length(unique(meas7$id)),
         k      = 6, age= meas7$age[!duplicated(meas7$id)],
         gender= meas7$sex[!duplicated(meas7$id)],nu=3,D=diag(2),
         numnonzeros=nonzeros,nonzeropos=t(nonzeropos),
         y2=(meas7$modvigmin2[meas7$modvigmin2>0]),n2=sum(meas7$modvigmin>0),
         X=x,pk=ncol(x),waist=waist,lglu=lglu,ltri=ltri,bps=bps,ldl=ldl,
         hdl=hdl,bpd=bpd
)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

start1 <- list(#alphag=c(0.1642,3.4081,1.4433,4.7228),
               alphaw=c(10.445,   3.230,   2.033, 101.517 ))#,
               #alphat=c(0.2805, 4.4733, 1.8297, 4.9261 ),
               #alphabs=c( 18.388,   4.602 ,  1.389, 137.256 ))
start2 <- list(#alphag=c(0.1242,4.4081,1.1433,3.7228),
               alphaw=c(8.445,   2.230,   1.033, 80.517 ))#,
               #alphat=c(0.1805, 3.4733, 1.1297, 3.9261 ),
               #alphabs=c( 12.388,   3.602 ,  1.189, 120.256 ))
start3 <- list(#alphag=c(0.3642,5.4081,1.9433,5.7228),
               alphaw=c(11.445,   4.230,   2.933, 111.517 ))#,
               #alphat=c(0.3805, 5.4733, 2.8297,5.9261 ),
               #alphabs=c( 19.388,   5.602 ,  2.389, 147.256 ))
start4 <- list(alphag=c(0.2642,2.4081,2.4433,5.7228),
               alphaw=c(9.445,   4.230,   1.033, 121.517 ),
               alphat=c(0.2205, 6.4733, 2.8297, 5.9261 ),
               alphabs=c( 28.388,   3.602 ,  2.389, 157.256 ))

ms <- stan_model(model_code=models)
rs <- sampling(ms,dat,c("beta","gamma","sigmae","Sigma","rho","alphaw",
                        #"alphag","alphat","alphal",
                        #"alphabs","alphabd",
                        "alphah",
                        "sigmawaist",#"sigmabps",
                        #"sigmaglu","sigmatri","sigmaldl",
                        "sigmahdl"#,"sigmabpd"
                        ),
                       iter=200,chains=3,init=list(start1,start2,start3),
               control=list(adapt_delta=0.99,max_treedepth=15))
summary(rs)
save(rs,file="stanout_full.RData")

