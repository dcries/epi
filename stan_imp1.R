library(rstan)
library(dplyr)
library(reshape)

#setwd("C:\\Users\\dcries\\github\\epi")
setwd("/home/dcries")
imp1 <- read.csv("epi/NHANES_accel_imp1.csv")
#nhanes <- read.csv("NHANES_complete.csv")
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
# m1 <- lm(modvigmin^.25~(weekend)+first5,data=nhanes)
# wbar <- mean((nhanes$modvigmin[nhanes$rep <= 5])^.25)
# w1 <- nhanes$modvigmin^.25
# what <- predict(m1)
# w <- (1/what)*w1*wbar
m1 <- lm(modvigmin~(weekend)+first5,data=imp1)
wbar <- mean((imp1$modvigmin[imp1$rep <= 5]))
w1 <- imp1$modvigmin
what <- predict(m1)
w <- (1/what)*w1*wbar
imp1$modvigmin2 <- w^.25
nrep <- (imp1 %>% group_by(id) %>% summarise(n=length(id)))$n
meas7 <- subset(imp1, id %in% unique(id)[nrep==6]) #individuals with all 7 days
meas7$a <- 1;meas7$a[meas7$age>=35] <- 2;meas7$a[meas7$age>=50] <- 3;meas7$a[meas7$age>=65] <- 4
meas7$corrg <- 1;meas7$corrg[meas7$a == 4] <- 2

#meas7 <- meas7[(!is.na(meas7$waist)) & (!is.na(meas7$bps)) & (!is.na(meas7$bpd)) & (!is.na(meas7$hdl)),] #remove NAs for waist

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
//real waist[N];
//real lglu[N];
//real ltri[N];
//real bps[N];
//real ldl[N];
//real hdl[N];
//real bpd[N];
real<lower=0> age[N];
//int<lower=0> gender[N];
int<lower=0> numnonzeros[N]; //number of nonzero minutes days for each individual
matrix[N,k] nonzeropos; //position of nonzero minutes for each indivudal
real nu;
matrix[2,2] D;
  vector[4] theta;
int cg[N]; //correlation group

}
transformed data{
vector[2] zeros;
vector<lower=0>[N] sigmae;
zeros = rep_vector(0,2);

  for(i in 1:N){
    //sigmae[i] = theta[4]+theta[1]/(1+exp(-theta[2]*(age[i]-theta[3])));
   sigmae[i] = theta[1] + theta[2]*age[i] + theta[3]*pow(age[i],2.0) + theta[4]*pow(age[i],3.0);
  }
}
parameters{
vector[pk] gamma;
//    real gamma1;
//real gamma2;
vector[pk] beta;
//vector[4] alphaw;
//vector[4] alphag;
//vector[4] alphat;
//vector[4] alphabs;
//vector[3] alphal;
//vector[3] alphabd;
//vector[2] alphah;

//real beta1;
//real beta2;
//cov_matrix[2] Sigma;
corr_matrix[2] L;
vector<lower=0>[2] sigmab;
//real<lower=0> sigmae;
//real<lower=-1,upper=1> rho;
vector <lower=-1,upper=1>[2] rho;

matrix[N,2] b;
//real<lower=0> sigma2waist;
//real<lower=0> sigma2glu;
//real<lower=0> sigma2tri;
//real<lower=0> sigma2bps;
//real<lower=0> sigma2bpd;
//real<lower=0> sigma2ldl;
//real<lower=0> sigma2hdl;
 //vector[4] theta;


}
transformed parameters{
cov_matrix[k] ar1mat[N];
vector[N] Tstar; // usual ^ 1/4
vector[N] T; //usual
vector[N] mu;
vector[N] p;
//vector<lower=0>[N] sigmae;

for(i in 1:N){
  //sigmae[i] = theta[4]+theta[1]/(1+exp(-theta[2]*(age[i]-theta[3])));
   //sigmae[i] = theta[1] + theta[2]*age[i] + theta[3]*pow(age[i],2.0) + theta[4]*pow(age[i],3.0);
for (m in 1:k){
//ar1mat[i,m,m] = pow(sigmae,2.0);
ar1mat[i,m,m] = sigmae[i];
}
}
//respecify ar1mat, nonzeropos, numnonzeros
for(i in 1:N){
for (m in 1:(k-1)) {
for (n in (m+1):k) {
      //ar1mat[i,m,n] = pow(sigmae[i],2.0) * if_else(n<=numnonzeros[i],pow(rho,nonzeropos[i,n]-nonzeropos[i,m]),0);
      ar1mat[i,m,n] = sigmae[i] * if_else(n<=numnonzeros[i],pow(rho[cg[i]],nonzeropos[i,n]-nonzeropos[i,m]),0);
      ar1mat[i,n,m] = ar1mat[i,m,n];
}
}
}

for(i in 1:N){
p[i] = inv_logit(X[i,]*gamma+b[i,1]);
mu[i] = X[i,]*beta + b[i,2];
    //T[i] = p[i]*(pow(mu[i],4.0) + 6*pow(sigmae,2.0)*pow(mu[i],2.0));
    T[i] = p[i]*(pow(mu[i],4.0) + 6*sigmae[i]*pow(mu[i],2.0));
Tstar[i] = pow(T[i],0.25);
}
}
model{

//matrix[N,k] mu;

//vector[N] muwaist;
//vector[N] mulglu;
//vector[N] multri;
//vector[N] mubps;
//vector[N] mubpd;
//vector[N] muhdl;
//vector[N] muldl;


int pos;
pos = 1;

for(i in 1:N){

//muwaist[i] = alphaw[1]-alphaw[2]/(1+exp(-alphaw[3]*(Tstar[i]-alphaw[4])));
//mulglu[i] = alphag[1]-alphag[2]/(1+exp(-alphag[3]*(Tstar[i]-alphag[4])));
//multri[i] = alphat[1]-alphat[2]/(1+exp(-alphat[3]*(Tstar[i]-alphat[4])));
//mubps[i] = alphabs[1]-alphabs[2]/(1+exp(-alphabs[3]*(Tstar[i]-alphabs[4])));
//muldl[i] = alphal[1] + alphal[2]*Tstar[i] + alphal[3]*pow(Tstar[i],2);
//mubpd[i] = alphabd[1] + alphabd[2]*Tstar[i] + alphabd[3]*pow(Tstar[i],2);
//muhdl[i] = alphah[1] + alphah[2]*Tstar[i];

//waist[i] ~ normal(muwaist[i],sigma2waist);
//lglu[i] ~ normal(mulglu[i],sigma2glu);
//ltri[i] ~ normal(multri[i],sigma2tri);
//bps[i] ~ normal(mubps[i],sigma2bps);
//ldl[i] ~ normal(muldl[i],sigma2ldl);
//hdl[i] ~ normal(muhdl[i],sigma2hdl);
//bpd[i] ~ normal(mubpd[i],sigma2bpd);


b[i] ~ multi_normal(zeros,diag_matrix(sigmab)*L*diag_matrix(sigmab));

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
//sigmae ~ cauchy(0,1);
//sigma2waist ~ inv_gamma(1,1);
//sigma2glu ~ inv_gamma(1,1);
//sigma2tri ~ inv_gamma(1,1);
//sigma2bps ~ inv_gamma(1,1);
//sigma2bpd ~ inv_gamma(1,1);
//sigma2ldl ~ inv_gamma(1,1);
//sigma2hdl ~ inv_gamma(1,1);

beta ~ normal(0,100); //implies independent priors for beta ?
//  beta1 ~ normal(0,100);
//  beta2 ~ normal(0,100);
gamma ~ normal(0,100);
//  gamma1 ~ normal(0,100);
//  gamma2 ~ normal(0,100);
//alphaw ~ normal(0,100);
//alphag ~ normal(0,100);
//alphat ~ normal(0,100);
//alphabs ~ normal(0,100);
//alphal ~ normal(0,100);
//alphabd ~ normal(0,100);
//alphah ~ normal(0,100);
//Sigma ~ inv_wishart(nu,D);
L ~ lkj_corr(1.0);
sigmab ~ cauchy(0,1);

  //theta[1] ~ normal(0,5);
  //theta[2] ~ normal(0,5);
  //theta[3] ~ normal(0,5);
  //theta[4] ~ normal(0,5);
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

x <- model.matrix(~age+bmi+as.factor(sex)+as.factor(race)+as.factor(education),data=meas7[!duplicated(meas7$id),c("age","bmi","sex","race","education","weekend","first5")])

waist <- meas7$waist[!duplicated(meas7$id)]
lglu <- log(meas7$glu[!duplicated(meas7$id)])
ltri <- log(meas7$tri[!duplicated(meas7$id)])
bps <- (meas7$bps[!duplicated(meas7$id)])
ldl <- meas7$ldl[!duplicated(meas7$id)]
hdl <- meas7$hdl[!duplicated(meas7$id)]
bpd <- meas7$bpd[!duplicated(meas7$id)]
corrgroup <- meas7$corrg[!duplicated(meas7$id)]


dat=list(y=(yc[,3:8]),  N      = length(unique(meas7$id)),
         k      = 6, age= meas7$age[!duplicated(meas7$id)],
         gender= meas7$sex[!duplicated(meas7$id)],nu=3,D=diag(2),
         numnonzeros=nonzeros,nonzeropos=t(nonzeropos),
         y2=(meas7$modvigmin2[meas7$modvigmin2>0]),n2=sum(meas7$modvigmin>0),
         X=x,pk=ncol(x),theta=c(.2662,-.00806,.0002099,-1.625e-06),
         hdl=hdl,bpd=bpd,cg=corrgroup
)
#theta=c(.175,.161,56.67,.209),
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

start1 <- list(theta=c(2.662e-01,   -8.060e-03,    2.099e-04,   -1.625e-06))
start2 <- list(theta=c(.62,.08,45.90,.49))
start3 <- list(theta=c(.82,.18,65.90,1.69))
start4 <- list(theta=c(.82,.10,50.90,.99))

ms <- stan_model(model_code=models)
rs <- sampling(ms,dat,c("beta","gamma","L","sigmab","rho","Tstar"#,"alphaw",
                        #"alphag","alphat","alphal",
                        #"alphabs","alphabd","alphah",
                        #"sigma2waist","sigma2bps",
                        #"sigma2glu","sigma2tri","sigma2ldl",
                        #"sigma2hdl","sigma2bpd"
),chains=8,control=list(adapt_delta=0.99),
iter=2000)
(rs)
save(rs,file="/ptmp/dcries/stanout_imp1.RData")

