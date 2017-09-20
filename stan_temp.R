library(rstan)
library(dplyr)
setwd("/home/dcries")
imp1 <- read.csv("epi/NHANES_accel_imp1.csv")
Tstar <- read.csv("epi/tstar_imp1.csv")
sdT <- read.csv("epi/sdT_imp1.csv")

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

#meas7$Tstar <- rep(colMeans(r[,21:7938]),each=6) # is from 
#meas7$std <- apply(r[,21:7938],2,sd)


meas7 <- meas7[(!is.na(meas7$waist)) & (!is.na(meas7$bps)) & (!is.na(meas7$bpd)) & (!is.na(meas7$hdl)) & (!is.na(meas7$ldl)) & (!is.na(meas7$glu)) & (!is.na(meas7$tri)) & (!is.na(meas7$education)),] #remove NAs for waist

#Tstar <- ((meas7 %>% group_by(id) %>% summarise(m=mean(modvigmin)))$m)^.25
#Tstar <- ((meas7 %>% group_by(id) %>% summarise(m=mean(modvigmin2)))$m)



models <- "
data{
int<lower=0> N; //number of individuals
vector[N] Tstar; // usual ^ 1/4
vector[N] sdT;// vector of standard deviation parameters of Tstar 
vector[7] MetS[N]; //7xN matrix of MetS data

//real waist[N];
//real lglu[N];
//real ltri[N];
//real bps[N];
//real ldl[N];
//real hdl[N];
//real bpd[N];
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


//real<lower=0> sigmawaist;
//real<lower=0> sigmaglu;
//real<lower=0> sigmatri;
//real<lower=0> sigmabps;
//real<lower=0> sigmabpd;
//real<lower=0> sigmaldl;
//real<lower=0> sigmahdl;

vector[N] appx;
vector<lower=0>[7] sigmar;
corr_matrix[7] Lr;
}


model{

//vector[N] muwaist;
//vector[N] mulglu;
//vector[N] multri;
//vector[N] mubps;
//vector[N] mubpd;
//vector[N] muhdl;
//vector[N] muldl;

vector[7] mu[N]; //7 for number of regressions

for(i in 1:N){
//muwaist[i] = alphaw[4]-alphaw[1]/(1+exp(-alphaw[2]*(Tstar[i]-alphaw[3])));
//mulglu[i] = alphag[4]-alphag[1]/(1+exp(-alphag[2]*(Tstar[i]-alphag[3])));
//multri[i] = alphat[4]-alphat[1]/(1+exp(-alphat[2]*(Tstar[i]-alphat[3])));
//mubps[i] = alphabs[4]-alphabs[1]/(1+exp(-alphabs[2]*(Tstar[i]-alphabs[3])));
//muldl[i] = alphal[1] + alphal[2]*Tstar[i] + alphal[3]*pow(Tstar[i],2);
//mubpd[i] = alphabd[1] + alphabd[2]*Tstar[i] + alphabd[3]*pow(Tstar[i],2);
//muhdl[i] = alphah[1] + alphah[2]*Tstar[i];

//waist[i] ~ normal(muwaist[i],sigmawaist);
//lglu[i] ~ normal(mulglu[i],sigmaglu);
//ltri[i] ~ normal(multri[i],sigmatri);
//bps[i] ~ normal(mubps[i],sigmabps);
//ldl[i] ~ normal(muldl[i],sigmaldl);
//hdl[i] ~ normal(muhdl[i],sigmahdl);
//bpd[i] ~ normal(mubpd[i],sigmabpd);


mu[i,1] = alphaw[4]-alphaw[1]/(1+exp(-alphaw[2]*(Tstar[i]+appx[i]-alphaw[3])));
mu[i,2] = alphag[4]-alphag[1]/(1+exp(-alphag[2]*(Tstar[i]+appx[i]-alphag[3])));
mu[i,3] = alphat[4]-alphat[1]/(1+exp(-alphat[2]*(Tstar[i]+appx[i]-alphat[3])));
mu[i,4] = alphabs[4]-alphabs[1]/(1+exp(-alphabs[2]*(Tstar[i]+appx[i]-alphabs[3])));
mu[i,5] = alphal[1] + alphal[2]*(Tstar[i]+appx[i]);// + alphal[3]*pow(Tstar[i]+appx[i],2);
mu[i,6] = alphabd[1] + alphabd[2]*(Tstar[i]+appx[i]);// + alphabd[3]*pow(Tstar[i]+appx[i],2);
mu[i,7] = alphah[1] + alphah[2]*(Tstar[i]+appx[i]);

MetS[i] ~ multi_normal(mu[i],diag_matrix(sigmar)*Lr*diag_matrix(sigmar));


}

//muwaist = alphaw[4]-alphaw[1]/(1+exp(-alphaw[2]*((Tstar+appx-alphaw[3])));
//mulglu = alphag[4]-alphag[1]/(1+exp(-alphag[2]*(Tstar+appx-alphag[3])));
//multri = alphat[4]-alphat[1]/(1+exp(-alphat[2]*(Tstar+appx-alphat[3])));
//mubps = alphabs[4]-alphabs[1]/(1+exp(-alphabs[2]*(Tstar+appx-alphabs[3])));
//muldl = alphal[1] + alphal[2]*(Tstar+appx) ;//+ alphal[3]*pow(Tstar+appx,2);
//mubpd = alphabd[1] + alphabd[2]*(Tstar+appx) ;//+ alphabd[3]*pow(Tstar+appx,2);
//muhdl = alphah[1] + alphah[2]*(Tstar+appx);

appx ~ normal(0,sdT);

//sigmawaist ~ cauchy(0,1);
//sigmaglu ~ cauchy(0,1);
//sigmatri ~ cauchy(0,1);
//sigmabps ~ cauchy(0,1);
//sigmabpd ~ cauchy(0,1);
//sigmaldl ~ cauchy(0,1);
//sigmahdl ~ cauchy(0,1);


alphaw[1] ~ normal(7,8);
alphaw[2] ~ normal(3,1.5);
alphaw[3] ~ normal(2.11,.4);
alphaw[4] ~ normal(98,17);

alphabs[1] ~ normal(5,10);
alphabs[2] ~ normal(4,2);
alphabs[3] ~ normal(2.11,.4);
alphabs[4] ~ normal(125,22);

alphag[1] ~ normal(0.05,.15);
alphag[2] ~ normal(3.6,2);
alphag[3] ~ normal(2.11,.4);
alphag[4] ~ normal(4.62,.25);

alphat[1] ~ normal(.27,.2);//normal(.12,.2); 
alphat[2] ~ normal(4.88,2);//normal(4.88,2);
alphat[3] ~ normal(1.81,1);//normal(2.11,.4);
alphat[4] ~ normal(4.92,2);//normal(4.73,.6);

//alphaw[1] ~ normal(10,8);
//alphaw[2] ~ normal(3,1.5);
//alphaw[3] ~ normal(2,1);
//alphaw[4] ~ normal(100,15);

//alphabs[1] ~ normal(18,10);
//alphabs[2] ~ normal(4,2);
//alphabs[3] ~ normal(1.3,1);
//alphabs[4] ~ normal(137,20);

//alphag[1] ~ normal(0.16,.15);
//alphag[2] ~ normal(3.6,2);
//alphag[3] ~ normal(1.4,1);
//alphag[4] ~ normal(4.7,3);

//alphat[1] ~ normal(.27,.2);
//alphat[2] ~ normal(4.88,2);
//alphat[3] ~ normal(1.81,1);
//alphat[4] ~ normal(4.92,2);

alphal ~ normal(0,100);
alphabd ~ normal(0,100);
alphah ~ normal(0,100);

sigmar ~ cauchy(0,1);
Lr ~ lkj_corr(1.0);
}
"

waist <- meas7$waist[!duplicated(meas7$id)]
lglu <- log(meas7$glu[!duplicated(meas7$id)])
ltri <- log(meas7$tri[!duplicated(meas7$id)])
bps <- (meas7$bps[!duplicated(meas7$id)])
ldl <- meas7$ldl[!duplicated(meas7$id)]
bpd <- meas7$bpd[!duplicated(meas7$id)]
hdl <- meas7$hdl[!duplicated(meas7$id)]

#Tstar <- meas7$Tstar[!duplicated(meas7$id)]
#sdT <- meas7$std[!duplicated(meas7$id)]
MetS <- (cbind(waist,lglu,ltri,bps,ldl,bpd,hdl))


dat=list(N      = length(unique(meas7$id)),
         #age= meas7$age[!duplicated(meas7$id)],
         #gender= meas7$sex[!duplicated(meas7$id)],nu=3,D=diag(2),
         Tstar=Tstar$x,sdT=sdT$x,MetS=MetS
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
                        "Lr","sigmar"
),
iter=1000,chains=3,init=list(start1,start2,start3),
control=list(adapt_delta=0.99,max_treedepth=15))
(rs)
save(rs,file="stanout_temp.RData")
