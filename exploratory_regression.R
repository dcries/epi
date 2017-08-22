indlevel <- nhanes %>% group_by(id) %>% summarise(m=mean(modvigmin^(1/4)),s=sd(modvigmin^(1/4)),glu=glu[1],waist=waist[1],ldl=ldl[1],hdl=hdl[1],bps=bps[1],bpd=bpd[1],tri=tri[1],gender=sex[1],race=race[1],age=age[1],n=length(id))
indlevel$agegroup <- 1;indlevel$agegroup[indlevel$age>29] <- 2;indlevel$agegroup[indlevel$age>45] <- 3;indlevel$agegroup[indlevel$age>63] <- 4;


a <- nhanes %>% group_by(id) %>% summarise(m=mean(modvigmin),s=sd(modvigmin))
q1=qplot(data=a,x=m) + xlab('Observed Mean Minutes in MVPA') + theme_bw()
q2=qplot(data=indlevel,x=m,y=s) + xlab('Observed Mean Minutes in MVPA (Transformed)') + ylab("Standard Deviation") +theme_bw()
grid.arrange(q1,q2,nrow=1)

p1 <- qplot(data=indlevel,x=m,y=log(glu))+ylab("log Glucose") + xlab("Mean Transformed MVPA Minutes") + geom_smooth() + theme_bw()
p2 <- qplot(data=indlevel,x=m,y=waist)+ylab("Waist") + xlab("Mean Transformed MVPA Minutes") + geom_smooth() + theme_bw()
p3 <- qplot(data=indlevel,x=m,y=log(tri))+ylab("log Triglycerides") + xlab("Mean Transformed MVPA Minutes") + geom_smooth() + theme_bw()
p4 <- qplot(data=indlevel,x=m,y=ldl)+ylab("LDL") + xlab("Mean Transformed MVPA Minutes") + geom_smooth() + theme_bw()
p5 <- qplot(data=indlevel,x=m,y=hdl)+ylab("HDL") + xlab("Mean Transformed MVPA Minutes") + geom_smooth() + theme_bw()
p6 <- qplot(data=indlevel,x=m,y=bps)+ylab("Blood Pressure (S)") + xlab("Mean Transformed MVPA Minutes") + geom_smooth() + theme_bw()
p7 <- qplot(data=indlevel,x=m,y=bpd)+ylab("Blood Pressure (D)") + xlab("Mean Transformed MVPA Minutes") + geom_smooth() + theme_bw()
grid.arrange(p1,p2,p3,p4,p5,p6,p7,nrow=4)

plot(a$m,a$s)


mbps <- nls(bps~L2-L/(1+exp(-k*(m-x0))),start=list(L=20,k=5,x0=1.2,L2=140),data=indlevel)

df=data.frame(x=seq(from=0,to=4,by=0.1))
df$bps=137.256-18.388/(1+exp(-4.602*(df$x-1.389)))
qplot(data=indlevel,x=m,y=bps) + geom_smooth() + geom_line(data=df,aes(x=x,y=bps),col="red") 

#mglu <- nls(glu~L2-L/(1+exp(-k*(m-x0))),start=list(L=20,k=3,x0=1.2,L2=120),data=indlevel)
mglu <- nls(log(glu)~L2-L/(1+exp(-k*(m-x0))),start=list(L=.5,k=3,x0=1.2,L2=log(120)),data=indlevel)
df$glu=coef(mglu)["L2"]-coef(mglu)["L"]/(1+exp(-coef(mglu)["k"]*(df$x-coef(mglu)["x0"])))
qplot(data=indlevel,x=m,y=log(glu)) + geom_smooth() + geom_line(data=df,aes(x=x,y=glu),col="red") 

mwaist <- nls(waist~L2-L/(1+exp(-k*(m-x0))),start=list(L=20,k=3,x0=1.2,L2=110),data=indlevel)
df$waist=coef(mwaist)["L2"]-coef(mwaist)["L"]/(1+exp(-coef(mwaist)["k"]*(df$x-coef(mwaist)["x0"])))
qplot(data=indlevel,x=m,y=waist) + geom_smooth() + geom_line(data=df,aes(x=x,y=waist),col="red") 

#mtri <- nls(tri~L2-L/(1+exp(-k*(m-x0))),start=list(L=30,k=3,x0=1.5,L2=150),data=indlevel)
mtri <- nls(log(tri)~L2-L/(1+exp(-k*(m-x0))),start=list(L=.5,k=3,x0=1.5,L2=5),data=indlevel)
df$tri=coef(mtri)["L2"]-coef(mtri)["L"]/(1+exp(-coef(mtri)["k"]*(df$x-coef(mtri)["x0"])))
qplot(data=indlevel,x=m,y=log(tri)) + geom_smooth() + geom_line(data=df,aes(x=x,y=tri),col="red") #+ ylim(c(0,750))

qplot(data=indlevel,x=m,y=ldl) + geom_smooth()
qplot(data=indlevel,x=m,y=hdl) + geom_smooth()
qplot(data=indlevel,x=m,y=bpd) + geom_smooth()

qplot(data=nhanes,x=modvigmin,y=ldl) + geom_smooth()
qplot(data=nhanes,x=modvigmin,y=hdl) + geom_smooth()
qplot(data=nhanes,x=modvigmin,y=bpd) + geom_smooth()

summary(lm(hdl~m,data=indlevel))
