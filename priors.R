library(ggplot2)
library(lme4)
library(dplyr)
library(reshape)
library(foreign)

setwd("C:\\Users\\dcries\\github\\epi")
nhanes <- read.csv("NHANES_accel.csv")
names(nhanes) <- tolower(names(nhanes))

setwd("\\\\my.files.iastate.edu\\Users\\dcries\\Downloads")

pa <- read.xport("PAQ_B.xpt")
bp <- read.xport("BPX_B.xpt")
glu <- read.xport("L10AM_B.XPT")
tri <- read.xport("L13_2_B.XPT")
waist <- read.xport("BMX_B.XPT")

#blood pressure systotlic and diastolic
bp <- bp[,c("SEQN","BPXSY1","BPXDI1")]
names(bp)[-1] <- c(paste0("bps",1),paste0("bpd",1))

#waist circumference
waist <- waist[,c("SEQN","BMXWAIST")]
names(waist)[-1] <- "waist"

#tryglyceride mg, LDL mg
tri <- tri[,c("SEQN","LB2TR","LB2LDL","LB2HDL")]
names(tri)[-1] <- c("tri","ldl","hdl")

#fasting glucose mg
glu <- glu[,c("SEQN","LBXGLU")]
names(glu)[-1] <- "glu"

#pa, mod actvitiy y/n, walked or biked y/n, level of PA 1-4
pa <- pa[,c("SEQN","PAD320","PAD020","PAQ180")]
names(pa)[-1] <- c("modact","walked","paqlevel")

#blood pressure
pabp <- left_join(pa,bp)
pabp[complete.cases(pabp),] %>% group_by(modact) %>% summarise(ms = mean(bps1),md=mean(bpd1),n=length(bps1),s=sd(bps1))
pabp[complete.cases(pabp),] %>% group_by(walked) %>% summarise(ms = mean(bps1),md=mean(bpd1),n=length(bps1),s=sd(bps1))
pabp[complete.cases(pabp),] %>% group_by(paqlevel) %>% summarise(ms = mean(bps1),md=mean(bpd1),n=length(bps1),s=sd(bps1))


#waist
pawaist <- left_join(pa,waist)
pawaist[complete.cases(pawaist),] %>% group_by(modact) %>% summarise(ms = mean(waist),n=length(waist),s=sd(waist))
pawaist[complete.cases(pawaist),] %>% group_by(walked) %>% summarise(ms = mean(waist),n=length(waist),s=sd(waist))
pawaist[complete.cases(pawaist),] %>% group_by(paqlevel) %>% summarise(ms = mean(waist),n=length(waist),s=sd(waist))

#tri
patri <- left_join(pa,tri)
patri[complete.cases(patri),] %>% group_by(modact) %>% summarise(mt = mean(log(tri)),ml=mean(ldl),md=mean(hdl),n=length(tri),s=sd(log(tri)))
patri[complete.cases(patri),] %>% group_by(walked) %>% summarise(mt = mean(log(tri)),ml=mean(ldl),md=mean(hdl),n=length(tri),s=sd(log(tri)))
patri[complete.cases(patri),] %>% group_by(paqlevel) %>% summarise(mt = mean(log(tri)),ml=mean(ldl),md=mean(hdl),n=length(tri),s=sd(log(tri)))

#glu
paglu <- left_join(pa,glu)
paglu[complete.cases(paglu),] %>% group_by(modact) %>% summarise(mt = mean(log(glu)),n=length(glu),s=sd(log(glu)))
paglu[complete.cases(paglu),] %>% group_by(walked) %>% summarise(mt = mean(log(glu)),n=length(glu),s=sd(log(glu)))
paglu[complete.cases(paglu),] %>% group_by(paqlevel) %>% summarise(mt = mean(log(glu)),n=length(glu),s=sd(log(glu)))
