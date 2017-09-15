library(ggplot2)
library(lme4)
library(dplyr)
library(reshape)
library(foreign)

setwd("C:\\Users\\dcries\\github\\epi")
nhanes <- read.csv("NHANES_accel.csv")
names(nhanes) <- tolower(names(nhanes))

setwd("\\\\my.files.iastate.edu\\Users\\dcries\\Downloads")

bp1 <- read.xport("BPX_C.xpt")
bp2 <- read.xport("BPX_D.xpt")

waist1 <- read.xport("BMX_C.xpt")
waist2 <- read.xport("BMX_D.xpt")

hdl1 <- read.xport("L13_C.xpt")
hdl2 <- read.xport("HDL_D.xpt")

tri1 <- read.xport("L13AM_C.xpt")
tri2 <- read.xport("TRIGLY_D.xpt")

glu1 <- read.xport("L10AM_C.xpt")
glu2 <- read.xport("glu_D.xpt")

fitness1 <- read.xport("CVX_C.XPT")

meds1 <- read.xport("BPQ_C.xpt")
meds2 <- read.xport("BPQ_D.xpt")

#blood pressure systotlic and diastolic
bp <- rbind(bp1[,c("SEQN","BPXSY1","BPXSY2","BPXSY3","BPXSY4","BPXDI1","BPXDI2","BPXDI3","BPXDI4")],
bp2[,c("SEQN","BPXSY1","BPXSY2","BPXSY3","BPXSY4","BPXDI1","BPXDI2","BPXDI3","BPXDI4")])
names(bp)[-1] <- c(paste0("bps",1:4),paste0("bpd",1:4))

#waist circumference
waist <- rbind(waist1[,c("SEQN","BMXWAIST")],
waist2[,c("SEQN","BMXWAIST")])
names(waist)[-1] <- "waist"

#tryglyceride mg, LDL mg
tri <- rbind(tri1[,c("SEQN","LBXTR","LBDLDL")],
tri2[,c("SEQN","LBXTR","LBDLDL")])
names(tri)[-1] <- c("tri","ldl")

#hdl mg, 
names(hdl1)[which(names(hdl1)=="LBXHDD")] <- "LBDHDD"
hdl <- rbind(hdl1[,c("SEQN","LBDHDD")],
hdl2[,c("SEQN","LBDHDD")])
names(hdl)[-1] <- "hdl"

#fasting glucose mg
glu <- rbind(glu1[,c("SEQN","LBXGLU")],
glu2[,c("SEQN","LBXGLU")])
names(glu)[-1] <- "glu"

#predicted v02, estimated v02, and by class
fit <- fitness1[,c("SEQN","CVDVOMAX","CVDESVO2","CVDFITLV")]
names(fit) <- c("id","predv02","estv02","fitlev")

a=left_join(bp,glu)
b=left_join(waist,tri)
d=left_join(a,b)
full=left_join(d,hdl)
full$bps <- rowMeans(full[,2:5],na.rm=TRUE)
full$bpd <- rowMeans(full[,6:9],na.rm=TRUE)
sum(complete.cases(full[,10:16]))

names(full)[1] <- "id"
complete=inner_join(nhanes,full[,c(1,10:16)])
sum(complete.cases(complete))
length(unique(complete$id[complete.cases(complete)]))


#---------------
meds <- rbind(meds1[,c("SEQN","BPQ040A")],meds2[,c("SEQN","BPQ040A")])
names(meds) <- c("id","bpmed")
complete=left_join(nhanes,meds)

write.csv(complete,file="nhanes_complete.csv",row.names=FALSE)


#----------------
setwd("\\\\my.files.iastate.edu\\Users\\dcries\\Downloads")

# demo1 <- read.xport("DEMO_C.XPT")
# demo2 <- read.xport("DEMO_D.XPT")
# names(demo1)[1] <- "id";names(demo2)[1] <- "id";
# demo <- rbind(demo1[,c("id","DMDHREDU")],demo2[,c("id","DMDHREDU")])
# d=left_join(nhanes,demo)
# names(d)[38] <- "education"
# write.csv(d,file="nhanes_complete.csv",row.names=FALSE)
