library(foreign)
library(reshape)
library(data.table)

setwd("C:\\Users\\dcries\\data")

nhanes <- read.csv("C:\\Users\\dcries\\github\\epi\\NHANES_complete.csv")
#----------------------------
# nhanes1 <- read.xport("paxraw_c.xpt")
# nhanes1dt <- data.table(nhanes1)
# nhanes2 <- read.xport("paxraw_d.xpt")
# nhanes2dt <- data.table(nhanes2)
# 
# #reduce data to only adults and those wiht 10080 accelerometry minutes, so it matches data given from nci, see sas code
# #keeps only reliable data and in calibration
# nhanes1red <- nhanes1dt[SEQN %in% unique(nhanes$id)]
# nhanes2red <- nhanes2dt[SEQN %in% unique(nhanes$id),-"PAXSTEP"]
# final <- rbindlist(list(nhanes1red,nhanes2red))
# 
# fwrite(final,file="nhanes_raw.csv")

#--------------------
full <- fread("nhanes_raw.csv")
#make sure its sorted
full <- full[order(SEQN,PAXN)]

full[,minute := 60*PAXHOUR+PAXMINUT]
#based on replicate days
full[,rep:= 1];full[PAXN>1440, rep:=2];full[PAXN>2*1440, rep:=3];full[PAXN>3*1440, rep:=4];full[PAXN>4*1440, rep:=5];full[PAXN>5*1440, rep:=6];full[PAXN>6*1440, rep:=7]

full2 <- data.frame(full)
fullm <- reshape::melt(full2[,c("SEQN","PAXDAY","PAXINTEN","minute","rep")],id.vars=c("SEQN","PAXDAY","minute","rep"))

fullc <- cast(fullm,SEQN+PAXDAY+rep~minute)

fullc1 <- cast(fullm[1:40320000,],SEQN+PAXDAY+rep~minute)
fullc2 <- cast(fullm[40320001:nrow(fullm),],SEQN+PAXDAY+rep~minute)


fwrite(data.table(rbindlist(list(fullc1,fullc2))),file="nhanes_casted.csv")
#-----------------------------


# nh <- full[SEQN %in% c(21005,21010),]
# #add minute on 1440 scale
# nh[,minute := 60*nh$PAXHOUR+nh$PAXMINUT]
# #based on replicate days
# nh[,rep:= 1];nh[PAXN>1440, rep:=2];nh[PAXN>2*1440, rep:=3];nh[PAXN>3*1440, rep:=4];nh[PAXN>4*1440, rep:=5];nh[PAXN>5*1440, rep:=6];nh[PAXN>6*1440, rep:=7]
# 
# nh2 <- data.frame(nh)
# nhm <- reshape::melt(nh2[,c("SEQN","PAXDAY","PAXINTEN","minute","rep")],id.vars=c("SEQN","PAXDAY","minute","rep"))
# 
# nhc <- cast(nhm,SEQN+PAXDAY+rep~minute)
