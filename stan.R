library(rstan)
library(dplyr)
library(reshape)

#setwd("C:\\Users\\dcries\\github\\epi")
setwd("/home/dcries")
nhanes <- read.csv("epi/nhanes_complete.csv")
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
  int<lower=0> n2; //total number of non-zero obs
  //real<lower=0> y[N];
  vector[k] y[N];
  vector[n2] y2;
  int pk; //number of covariates plus intercept
  matrix[N,pk] X;
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
//real beta1;
//real beta2;
cov_matrix[2] Sigma;
real<lower=0> sigma2e;
real<lower=-1,upper=1> rho;
matrix[N,2] b;
//real mu2;
//real<lower=0> sigma2;
}
transformed parameters{
cov_matrix[k] ar1mat[N];


//for (m in 1:k){
//  ar1mat[m,m] = sigma2e;
//}

//for (m in 1:(k-1)) {
//  for (n in (m+1):k) {
//    ar1mat[m,n] = sigma2e * pow(rho,abs(m-n));
//    ar1mat[n,m] = ar1mat[m,n];
//  }
//}

//for(i in 1:N){
  //for (m in 1:k){
    //ar1mat[i,m,m] = sigma2e;
 // }
//}
//respecify ar1mat, nonzeropos, numnonzeros
//for(i in 1:N){
  //for (m in 1:(numnonzeros[i]-1)) {
    //for (n in (m+1):numnonzeros[i]) {
      //ar1mat[i,m,n] = 0;//sigma2e * pow(rho,nonzeropos[i,n]-nonzeropos[i,m]);
      //ar1mat[i,n,m] = 0;//ar1mat[i,m,n];
    //}
  //}
//}

for(i in 1:N){
  for (m in 1:k){
    ar1mat[i,m,m] = sigma2e;
    //for(n in (m+1):k){
    //ar1mat[i,m,n] = sigma2e*pow(rho,abs(m-n));
    //ar1mat[i,n,m] = sigma2e*pow(rho,abs(m-n));
    //}
  }
}
//respecify ar1mat, nonzeropos, numnonzeros
for(i in 1:N){
  for (m in 1:(k-1)) {
    for (n in (m+1):k) {
      ar1mat[i,m,n] = sigma2e * if_else(n<=numnonzeros[i],pow(rho,nonzeropos[i,n]-nonzeropos[i,m]),0);
      ar1mat[i,n,m] = ar1mat[i,m,n];
    }
  }
}

}
model{

//matrix[N,k] mu;
vector[N] mu;
vector[N] p;
  int pos;
  pos = 1;

  for(i in 1:N){
    p[i] = Phi(X[i,]*gamma+b[i,1]);
    mu[i] = X[i,]*beta + b[i,2];
    b[i] ~ multi_normal(zeros,Sigma);

    for(j in 1:k){
      if(y[i][j] == 0){
        //increment_log_prob(bernoulli_log(0,p[i]));
        target += bernoulli_lpmf(0|p[i]);
      }
      else{
        //increment_log_prob(bernoulli_log(1,p[i]));
        target += bernoulli_lpmf(1|p[i]);// + normal_lpdf(y[i][j]|mu[i],sigma2e);
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

// for(i in 1:N){
    //y[i] ~ multi_normal(row(mu,i),ar1mat);
   // b[i] ~ multi_normal(zeros,Sigma);
  //}

//would need to change dimension of mu
//  for(i in 1:N){
  //    mu[i] = beta0 + beta1*age[i] + beta2*gender[i];
 // }

// need to define ind, change y--so it only includes non-zeros

  //for(i in 1:N){
   // segment(y2,pos,numnonzeros[i]) ~ multi_normal(rep_vector(mu[i],numnonzeros[i]),block(ar1mat[i],1,1,numnonzeros[i],numnonzeros[i]));
    //segment(y[i],pos,numnonzeros[i]) ~ multi_normal(rep_vector(mu[i],numnonzeros[i]),block(ar1mat[i],1,1,numnonzeros[i],numnonzeros[i]));
    //pos = pos + numnonzeros[i];
  //}
//for(i in 1:n2){
//y2[i] ~ normal(mu2,sigma2);
//}

//mu2 ~ normal(0,1000);
//sigma2 ~ inv_gamma(1,1);
//  rho ~ uniform(-1,1);
  sigma2e ~ inv_gamma(1,1);
 beta ~ normal(0,100); //implies independent priors for beta ?
//  beta1 ~ normal(0,100);
//  beta2 ~ normal(0,100);
  gamma ~ normal(0,100);
//  gamma1 ~ normal(0,100);
//  gamma2 ~ normal(0,100);
  Sigma ~ inv_wishart(nu,D);
}
"

ym <- melt(meas7[,c("modvigmin","rep","id")],id.vars=c("id","rep"))
yc <- cast(ym,id+variable~rep)
nonzeros <- (meas7 %>% group_by(id) %>% summarise(n=sum(active)))$n
a1=meas7$rep[meas7$active==1]
a2=cumsum(nonzeros)

nonzeropos <- matrix(0,ncol=length(unique(meas7$id)),nrow=7)
nonzeropos[,1] <- a1[1:nonzeros[1]]
for(i in 2:ncol(nonzeropos)){
    temp <- a1[(a2[i-1]+1):a2[i]]
    if(length(temp)<7){
      n0 <- 7-length(temp)
      temp <- c(temp,rep(0,n0))
    }
    nonzeropos[,i] <- temp
}

x <- model.matrix(~age+sex+as.factor(race)+weekend+first5,data=meas7[!duplicated(meas7$id),c("age","sex","race","weekend","first5")])


dat=list(y=(yc[,3:9])^(1/4),  N      = length(unique(meas7$id)),
         k      = 7, age= meas7$age[!duplicated(meas7$id)],
         gender= meas7$sex[!duplicated(meas7$id)],nu=3,D=diag(2),
         numnonzeros=nonzeros,nonzeropos=t(nonzeropos),
         y2=(meas7$modvigmin[meas7$modvigmin>0])^(1/4),n2=sum(meas7$modvigmin>0),
         X=x,pk=ncol(x)
)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

ms <- stan_model(model_code=models)
rs <- sampling(ms,dat,c("beta","gamma","sigma2e","Sigma","rho"),
                       iter=2000)
summary(rs)
save(rs,file="stanout.RData")

