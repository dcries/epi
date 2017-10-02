library(reshape)
library(data.table)
setwd("C:\\Users\\dcries\\data")
data <- read.csv("nhanes_casted.csv")
full <- fread("nhanes_raw.csv")


load("C:/Users/dcries/workspace/impute1.RData")
results1 <- results
load("C:/Users/dcries/workspace/impute2.RData")
results2 <- results
load("C:/Users/dcries/workspace/impute3.RData")
results3 <- results
load("C:/Users/dcries/workspace/impute4.RData")
results4 <- results


impute1 <- cbind(data[,1:3],rbind(results1[[1]][[1]],results2[[1]][[1]],results3[[1]][[1]],results4[[1]][[1]]))
impute2 <- cbind(data[,1:3],rbind(results1[[2]][[1]],results2[[2]][[1]],results3[[2]][[1]],results4[[2]][[1]]))
impute3 <- cbind(data[,1:3],rbind(results1[[3]][[1]],results2[[3]][[1]],results3[[3]][[1]],results4[[3]][[1]]))
impute4 <- cbind(data[,1:3],rbind(results1[[4]][[1]],results2[[4]][[1]],results3[[4]][[1]],results4[[4]][[1]]))
impute5 <- cbind(data[,1:3],rbind(results1[[5]][[1]],results2[[5]][[1]],results3[[5]][[1]],results4[[5]][[1]]))

rm(results1);rm(results2);rm(results3);rm(results4);rm(results)



imp1 <- data.table(melt(impute1,id.vars=c("SEQN","PAXDAY","rep")))
imp1[,minute:=(as.numeric(substring(variable,2)))]
imp1[,PAXHOUR:=floor(minute/60)]
imp1[,PAXMINUT:=minute-PAXHOUR*60]
imp1 <- imp1[order(SEQN, PAXDAY,PAXHOUR,PAXMINUT)]
imp1 <- cbind(full[,c("PAXN","PAXSTAT","PAXCAL")],imp1)
fwrite(imp1,file="impute1.csv")

rm(impute1);rm(imp1)

imp2 <- data.table(melt(impute2,id.vars=c("SEQN","PAXDAY","rep")))
imp2[,minute:=(as.numeric(substring(variable,2)))]
imp2[,PAXHOUR:=floor(minute/60)]
imp2[,PAXMINUT:=minute-PAXHOUR*60]
imp2 <- imp2[order(SEQN, PAXDAY,PAXHOUR,PAXMINUT)]
imp2 <- cbind(full[,c("PAXN","PAXSTAT","PAXCAL")],imp2)
fwrite(imp2,file="impute2.csv")

rm(impute2);rm(imp2)


imp3 <- data.table(melt(impute3,id.vars=c("SEQN","PAXDAY","rep")))
imp3[,minute:=(as.numeric(substring(variable,2)))]
imp3[,PAXHOUR:=floor(minute/60)]
imp3[,PAXMINUT:=minute-PAXHOUR*60]
imp3 <- imp3[order(SEQN, PAXDAY,PAXHOUR,PAXMINUT)]
imp3 <- cbind(full[,c("PAXN","PAXSTAT","PAXCAL")],imp3)
fwrite(imp3,file="impute3.csv")

rm(impute3);rm(imp3)

imp4 <- data.table(melt(impute4,id.vars=c("SEQN","PAXDAY","rep")))
imp4[,minute:=(as.numeric(substring(variable,2)))]
imp4[,PAXHOUR:=floor(minute/60)]
imp4[,PAXMINUT:=minute-PAXHOUR*60]
imp4 <- imp4[order(SEQN, PAXDAY,PAXHOUR,PAXMINUT)]
imp4 <- cbind(full[,c("PAXN","PAXSTAT","PAXCAL")],imp4)
fwrite(imp4,file="impute4.csv")

rm(impute4);rm(imp4)

imp5 <- data.table(melt(impute5,id.vars=c("SEQN","PAXDAY","rep")))
imp5[,minute:=(as.numeric(substring(variable,2)))]
imp5[,PAXHOUR:=floor(minute/60)]
imp5[,PAXMINUT:=minute-PAXHOUR*60]
imp5 <- imp5[order(SEQN, PAXDAY,PAXHOUR,PAXMINUT)]
imp5 <- cbind(full[,c("PAXN","PAXSTAT","PAXCAL")],imp5)
fwrite(imp5,file="impute5.csv")

rm(impute5);rm(imp5)

