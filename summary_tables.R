library(xtable)

setwd("C:/Users/dcries/workspace")

load("stanout_realmix3.RData")
#naive analysis
#load("stanout_naive.RData")
#load("stanout_second.RData") ???
#naive <- as.matrix(rs)
#mem <- as.matrix(rs)

#for nonlinear function parameters
#naivenonlinear <- naive[,c(1:12,16:19)]


#-----------------------------
#nonlinear parameters

gamma0 <- matrix(0,ncol=7,nrow=nrow(out3$lambda))
for(i in 1:nrow(out3$lambda)){
  for(j in 2:K){
    gamma0[i,] <- gamma0[i,] + out$pi[i,j]*(out3$lambda[i,,1] - out3$lambda[i,,j])
  }
  gamma0[i,] <- out3$lambda[i,,1]-gamma0[i,]
}


tab <- matrix(0,nrow=16,ncol=3)
for(i in 1:4){
  tab[i,] <- round(c(mean(gamma0[,i]),quantile(gamma0[,i],probs=c(0.025,0.975))),3)
}
ind1 <- c(1,4,7,10,2,5,8,11,3,6,9,12)
for(i in 1:length(ind1)){
  tab[i+4,] <- round(c(mean(out3$beta[,ind1[i]]),quantile(out3$beta[,ind1[i]],probs=c(0.025,0.975))),3)
}


ci <- paste0("(",tab[,2],",",tab[,3],")")
metsfactor <- rep(c("Waist","Glucose","Triglyceride","Sys Blood Press"),4)
parameter <- c("$\\gamma_0$",rep("",3),"L",rep("",3),"K",rep("",3),"B",rep("",3))
df <- data.frame(cbind(parameter,metsfactor,tab[,1],ci))
names(df) <- c("Parameter","MetS RF","Post Mean","95\\% CI")

print(xtable(df,align="rll|ll"),sanitize.text.function=function(x){x},include.rownames = FALSE,
      hline.after = c(-1,0,4,8,12,16))
#--------------------

#for linear function parameters
tab2 <- matrix(0,nrow=6,ncol=3)
for(i in 1:3){
  tab2[i,] <- round(c(mean(gamma0[,i+4]),quantile(gamma0[,i+4],probs=c(0.025,0.975))),3)
}

ind2 <- c(13,14,15)
for(i in 1:length(ind2)){
  tab2[i+3,] <- round(c(mean(out3$beta[,ind2[i]]),quantile(out3$beta[,ind2[i]],probs=c(0.025,0.975))),3)
}


ci2 <- paste0("(",tab2[,2],",",tab2[,3],")")
metsfactor2 <- rep(c("LDL","Dias Blood Press","HDL"),2)
parameter2 <- c("$\\gamma_0$",rep("",2),"$\\gamma_1$",rep("",2))
df2 <- data.frame(cbind(parameter2,metsfactor2,tab2[,1],ci2))
names(df2) <- c("Parameter","MetS RF","Post Mean","95\\% CI")

print(xtable(df2,align="rll|ll"),sanitize.text.function=function(x){x},include.rownames = FALSE,
      hline.after = c(-1,0,3,6))

#--------------------------
#mean mat of lambda values
lammat <- data.frame(round(cbind(colMeans(out3$pi),t(apply(out3$lambda,c(2,3),mean))),2))
names(lammat) <- c("$p$","$\\lambda_{wst}$","$\\lambda_{glu}$",
                   "$\\lambda_{tri}$","$\\lambda_{bps}$","$\\lambda_{ldl}$",
                   "$\\lambda_{bpd}$","$\\lambda_{hdl}$")

print(xtable(lammat,align="r|l|lllllll"),sanitize.text.function=function(x){x},include.rownames = TRUE)

#--------------------------
#mean mat of sds values
sigmamat <- data.frame(round(cbind(colMeans(out3$pi),t(apply(out3$sds,c(2,3),mean))),2))
names(sigmamat) <- c("$p$","$\\sigma_{wst}$","$\\sigma_{glu}$",
                   "$\\sigma_{tri}$","$\\sigma_{bps}$","$\\sigma_{ldl}$",
                   "$\\sigma_{bpd}$","$\\sigma_{hdl}$")

print(xtable(sigmamat,align="r|l|lllllll"),sanitize.text.function=function(x){x},include.rownames = TRUE)

#-----------------------
#correlation matrix
sds <- t(apply(out3$sds,c(2,3),mean))
cormat <- array(0,dim=c(7,7,dim(out3$lambda)[3]))

for(i in 1:dim(cormat)[3]){
  count = 1
  for(j in 1:(nrow(cormat)-1)){
    for(k in (j+1):ncol(cormat)){
      cormat[j,k,i] = mean(out3$cormat[,count,i])
      count = count + 1
    }
  }
}

diag(cormat[,,1]) <- sds[1,]
diag(cormat[,,2]) <- sds[2,]
diag(cormat[,,3]) <- sds[3,]

cordf1 <- data.frame(round(cormat[,,1],3))
names(cordf1) <- c("wst","glu","tri","bps","ldl","bpd","hdl")
rownames(cordf1) <- c("wst","glu","tri","bps","ldl","bpd","hdl")
cordf1[cordf1==0] <- ""
print(xtable(cordf1,align="r|lllllll"))

cordf1 <- data.frame(round(cormat[,,2],3))
names(cordf1) <- c("wst","glu","tri","bps","ldl","bpd","hdl")
rownames(cordf1) <- c("wst","glu","tri","bps","ldl","bpd","hdl")
cordf1[cordf1==0] <- ""
print(xtable(cordf1,align="r|lllllll"))

cordf1 <- data.frame(round(cormat[,,3],3))
names(cordf1) <- c("wst","glu","tri","bps","ldl","bpd","hdl")
rownames(cordf1) <- c("wst","glu","tri","bps","ldl","bpd","hdl")
cordf1[cordf1==0] <- ""
print(xtable(cordf1,align="r|lllllll"))

#--------------------
