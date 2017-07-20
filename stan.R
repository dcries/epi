library(rstan)
library(dplyr)
library(reshape)

#setwd("C:\\Users\\dcries\\github\\epi")
setwd("/home/dcries/epi")
nhanes <- read.csv("nhanes_complete.csv")
names(nhanes) <- tolower(names(nhanes))
nhanes$weekend <- 0
nhanes$weekend[nhanes$dow %in% c(1,7)] <- 1
nhanes$first5 <- 0
nhanes$first5[nhanes$rep==6] <- 1
nhanes$first5[nhanes$rep==7] <- 2

nhanes$active <- 1
nhanes$active[nhanes$modvigmin ==0] <- 0

nrep <- (nhanes %>% group_by(id) %>% summarise(n=length(id)))$n
meas7 <- subset(nhanes, id %in% unique(id)[nrep==7]) #individuals with all 7 days


models <- "
data{
  int<lower=0> N; //number of individuals
  int<lower=0> k; // number of obs per individual
  //real<lower=0> y[N];
  matrix[N,k] y;
  real<lower=0> age[N];
  int<lower=0> gender[N];
}
parameters{
//  real gamma0;
//    real gamma1;
//real gamma2;
real beta0;
real beta1;
real beta2;
//cov_matrix[2] Sigma;
real<lower=0> sigma2e;
real<lower=-1,upper=1> rho;
}
transformed parameters{
cov_matrix[k] ar1mat;


for (m in 1:k){
  ar1mat[m,m] <- sigma2e;
}

for (m in 1:(k-1)) {
  for (n in (m+1):k) {
    ar1mat[m,n] = sigma2e * pow(rho,abs(m-n));
    ar1mat[n,m] = ar1mat[m,n];
  }
}

}
model{

matrix[N,k] mu;

  for(i in 1:N){
    for(j in 1:k){
      mu[i,j] = beta0 + beta1*age[i] + beta2*gender[i];
    }
  }

  for(i in 1:N){
    y[i,] ~ multi_normal(mu[i,],ar1mat);
  }

  rho ~ uniform(-1,1);
  sigma2e ~ inv_gamma(1,1);
  beta0 ~ normal(0,100);
  beta1 ~ normal(0,100);
  beta2 ~ normal(0,100);
}
"

ym <- melt(meas7[,c("modvigmin","rep","id")],id.vars=c("id","rep"))
yc <- cast(ym,id+variable~rep)

x <- meas7[,c("age","sex","race","weekend","first5")]

dat=list(y=yc[,3:9],  N      = length(unique(meas7$id)),
         k      = 7, age= meas7$age[!duplicated(meas7$id)],
         gender= meas7$sex[!duplicated(meas7$id)]
)

ms <- stan_model(model_code=models)
rs <- sampling(ms,dat,c("beta0","beta1","beta2","sigma2e","rho"))
