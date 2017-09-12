library(xtable)

setwd("C:/Users/dcries/workspace2")


#naive analysis
load("stanout_naive.RData")
#load("stanout_second.RData") ???
naive <- as.matrix(rs)
mem <- as.matrix(rs)

#for nonlinear function parameters
naivenonlinear <- naive[,c(1:12,16:19)]
memnonlinear <- mem[,c(1:12,16:19)]
tab <- matrix(0,nrow=32,ncol=3)
ind1 <- seq(from=1,to=31,by=2);ind2 <- seq(from=2,to=32,by=2)
ind3 <- c(1,5,9,13,2,6,10,14,3,7,11,15,4,8,12,16)
for(i in 1:(nrow(tab)/2)){
  tab[ind1[i],] <- round(c(mean(memnonlinear[,ind3[i]]),quantile(memnonlinear[,ind3[i]],probs=c(0.025,0.975))),3)
  tab[ind2[i],] <- round(c(mean(naivenonlinear[,ind3[i]]),quantile(naivenonlinear[,ind3[i]],probs=c(0.025,0.975))),3)
}
ci <- paste0("(",tab[,2],",",tab[,3],")")
metsfactor <- rep(c("Waist","","Glucose","","Triglyceride","","Sys Blood Press",""),4)
parameter <- c("L",rep("",7),"K",rep("",7),"B",rep("",7),"M",rep("",7))
df <- data.frame(cbind(parameter,metsfactor,rep(c("MEM","Naive"),8),tab[,1],ci))
names(df) <- c("Parameter","MetS RF","Model","Post Mean","95% CI")

#for linear function parameters
naivelinear <- naive[,c(13:15,20:22,23:24)]
memlinear <- mem[,c(13:15,20:22,23:24)]
tab2 <- matrix(0,nrow=16,ncol=3)
ind4 <- seq(from=1,to=15,by=2);ind5 <- seq(from=2,to=16,by=2)
ind6 <- c(1,4,7,2,5,8,3,6)
for(i in 1:(nrow(tab2)/2)){
  tab2[ind4[i],] <- round(c(mean(memlinear[,ind6[i]]),quantile(memlinear[,ind6[i]],probs=c(0.025,0.975))),3)
  tab2[ind5[i],] <- round(c(mean(naivelinear[,ind6[i]]),quantile(naivelinear[,ind6[i]],probs=c(0.025,0.975))),3)
}
ci2 <- paste0("(",tab2[,2],",",tab2[,3],")")
metsfactor2 <- c(rep(c("LDL","","Dias Blood Pressure","","HDL",""),2),"LDL","","Dias Blood Pressure","")
parameter2 <- c("b0",rep("",5),"b1",rep("",5),"b2",rep("",3))
df2 <- data.frame(cbind(parameter2,metsfactor2,rep(c("MEM","Naive"),8),tab2[,1],ci2))
names(df2) <- c("Parameter","MetS RF","Model","Post Mean","95% CI")


#correlation matrix, order w,g,t,ldl,bps,bpd,hdl
cormatnaive <- naive[,25:73]
cormatmem <- mem[,25:73]

naivemean <- colMeans(cormatnaive)
naivemem <- colMeans(cormatmem)

corn <- matrix(naivemean,ncol=7,nrow=7)
corm <- matrix(naivemem,ncol=7,nrow=7)

#variance terms
varnaive <- naive[,74:80]
varmem <- mem[,74:80]

colMeans(varnaive)
colMeans(varmem)



#--------------------
