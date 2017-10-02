library(dplyr)

setwd("C:\\Users\\dcries\\github\\epi")
imp1 <- read.csv("NHANES_accel_imp1.csv")
imp2 <- read.csv("NHANES_accel_imp2.csv")
imp3 <- read.csv("NHANES_accel_imp3.csv")
imp4 <- read.csv("NHANES_accel_imp4.csv")
imp5 <- read.csv("NHANES_accel_imp5.csv")
nhanes <- read.csv("nhanes_complete.csv")

names(imp1) <- tolower(names(imp1))
names(imp2) <- tolower(names(imp2))
names(imp3) <- tolower(names(imp3))
names(imp4) <- tolower(names(imp4))
names(imp5) <- tolower(names(imp5))


epi <- nhanes[,c("id","glu","waist","tri","ldl","hdl","bps","bpd")]

a=inner_join(imp1,epi)
imp1 <- a[!duplicated(a),]
a=inner_join(imp2,epi)
imp2 <- a[!duplicated(a),]
a=inner_join(imp3,epi)
imp3 <- a[!duplicated(a),]
a=inner_join(imp4,epi)
imp4 <- a[!duplicated(a),]
a=inner_join(imp5,epi)
imp5 <- a[!duplicated(a),]

write.csv(imp1,file="NHANES_accel_imp1.csv",row.names=F)
write.csv(imp2,file="NHANES_accel_imp2.csv",row.names=F)
write.csv(imp3,file="NHANES_accel_imp3.csv",row.names=F)
write.csv(imp4,file="NHANES_accel_imp4.csv",row.names=F)
write.csv(imp5,file="NHANES_accel_imp5.csv",row.names=F)
