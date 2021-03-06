library(accelmissing)
library(mice); library(pscl)

setwd("C:\\Users\\dcries\\data")
data <- read.csv("nhanes_casted.csv")
nhanes <- read.csv("C:\\Users\\dcries\\github\\epi\\NHANES_complete.csv")
all.equal(unique(data$SEQN),unique(nhanes$id))

flag60 = create.flag(data[,-c(1:3)], window=60)
label <- data[,1:2]
demo <- nhanes[!duplicated(nhanes$id),c("sex","age","race","bmi")]
demo$sex <- as.factor(demo$sex)
demo$race <- as.factor(demo$race)

#mr = missing.rate(label, flag60, mark.missing=0, time.range=c("09:00", "20:59"))
#mr$total #22.6 percent

# wearing proportion over time
#wear.time.plot(data[,-c(1:3)], label, flag60)

accelimp = accel.impute(PA=data[,-c(1:3)], label=label, flag=flag60, demo=demo,
 method="zipln", time.range=c("09:00","20:59"), m=5)

save(accelimp,file="impute.RData")
