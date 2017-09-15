library(ggplot2)
library(dplyr)
library(reshape)
library(foreign)
library(randomForest)

setwd("\\\\my.files.iastate.edu\\Users\\dcries\\Downloads")

demo1 <- read.xport("DEMO_C.XPT")
demo2 <- read.xport("DEMO_D.XPT")
#names(demo1)[1] <- "id";names(demo2)[1] <- "id";
vars <- c("SEQN","RIAGENDR","RIDAGEYR","DMDMARTL","DMDHHSIZ","INDFMINC","INDFMPIR","DMDHREDU")

demo <- rbind(demo1[,vars],demo2[,vars])
names(demo) <- c("id","gender","age","married","pplhshld","income","pir","education")
demo$education[demo$education > 5] <- NA
demo$education <- as.factor(demo$education)
demo$gender <- as.factor(demo$gender)
demo$married <- as.factor(demo$married)
demo$income <- as.factor(demo$income)

rf <- randomForest(education~.,data=demo[complete.cases(demo),-c(1,2,4)],importance=TRUE,ntree=100)
rf
varImpPlot(rf)

# demo2 <- subset(demo,education %in% c(1,3,5))
# demo2$education = as.factor(as.character(demo2$education))
# 
# rf <- randomForest(education~.,data=demo2[complete.cases(demo2),-c(1,2)],importance=TRUE,ntree=100)
# rf
# varImpPlot(rf)
setwd("C:\\Users\\dcries\\github\\epi")
nhanes <- read.csv("NHANES_complete.csv")

ids <- demo$id[which(demo$id %in% nhanes$id)]
newdemo <- demo[demo$id %in% ids,]

nas <- newdemo[is.na(newdemo$education) & (newdemo$id %in% ids),c("id","age","married","pplhshld","income","pir")]
preds <- predict(rf,newdata=nas[complete.cases(nas),-1]) #227 individuals had education level imputed
#nas$id[is.na(nas$pir)|is.na(nas$income)] #remove 47 individuals
nonaid <- nas[complete.cases(nas),"id"]

for(i in 1:length(preds)){
  newdemo$education[newdemo$id==nonaid[i]] <- preds[i]
}

d=left_join(nhanes[,1:37],newdemo[,c("id","education")])
d2=d[!is.na(d$education),]#remove 47 individuals
write.csv(d2,file="nhanes_complete.csv",row.names=FALSE)


imp1 <- read.csv("NHANES_accel_imp1.csv")
d3=left_join(imp1,newdemo[,c("id","education")])
write.csv(d3,file="NHANES_accel_imp1.csv",row.names=FALSE)

imp2 <- read.csv("NHANES_accel_imp2.csv")
d3=left_join(imp2,newdemo[,c("id","education")])
write.csv(d3,file="NHANES_accel_imp2.csv",row.names=FALSE)

imp3 <- read.csv("NHANES_accel_imp3.csv")
d3=left_join(imp3,newdemo[,c("id","education")])
write.csv(d3,file="NHANES_accel_imp3.csv",row.names=FALSE)

imp4 <- read.csv("NHANES_accel_imp4.csv")
d3=left_join(imp4,newdemo[,c("id","education")])
write.csv(d3,file="NHANES_accel_imp4.csv",row.names=FALSE)

imp5 <- read.csv("NHANES_accel_imp5.csv")
d3=left_join(imp5,newdemo[,c("id","education")])
write.csv(d3,file="NHANES_accel_imp5.csv",row.names=FALSE)
