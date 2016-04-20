##Does total cells of each symbiont change with presence of other symbionts 
##total A in Sym_ID - A vs. AB vs. AD

##Make subset variables
NCdata=read.table(file.choose(),header=T,sep="\t") ##Data for R cont rem
sub_A=subset(NCdata,NCdata$Sym_ID=="A"|NCdata$Sym_ID=="AB"|NCdata$Sym_ID=="AD")

##A ANOVA with Light treatment nested
modelA=aov(sub_A$Log_A~sub_A$Sym_ID*sub_A$Env+Error(sub_A$Rep/sub_A$Sym_ID/sub_A$Env))
summary(modelA)

##Make variables for histogram
A_A=subset(sub_A,sub_A$Sym_ID=="A")
A_AB=subset(sub_A,sub_A$Sym_ID=="AB")
A_AD=subset(sub_A,sub_A$Sym_ID=="AD")

##Histogram to show difference in distribution
par(mfrow=c(3,1))
hist(A_A$Log_A,breaks=4,xlim=c(0,6),main="A Only Treatment",xlab="",col="gray40")
hist(A_AB$Log_A,breaks=6,xlim=c(0,6),main="AB Treatment",xlab="",col="gray40")
hist(A_AD$Log_A,breaks=6,xlim=c(0,6),main="AD Treatment",xlab="Log Total Symbiodinium A",col="gray40")


##B ANOVA with Light treatment nested
modelB=aov(sub_B$Log_B~sub_B$Sym_ID*sub_B$Env+Error(sub_B$Rep/sub_B$Sym_ID/sub_B$Env))
summary(modelB)

##Make variables for histogram
B_B=subset(sub_B,sub_B$Sym_ID=="B")
B_AB=subset(sub_B,sub_B$Sym_ID=="AB")
B_BD=subset(sub_B,sub_B$Sym_ID=="BD")

##Histogram to show difference in distribution
par(mfrow=c(3,1))
hist(B_B$Log_B,breaks=4,xlim=c(0,6),main="B Only Treatment",xlab="")
hist(B_AB$Log_B,breaks=6,xlim=c(0,6),main="AB Treatment",xlab="")
hist(B_BD$Log_B,breaks=6,xlim=c(0,6),main="BD Treatment",xlab="Log Total Symbiodinium B")


##D ANOVA
modelD=aov(sub_D$Log_D~sub_D$Sym_ID*sub_D$Env+Error(sub_D$Rep/sub_D$Sym_ID/sub_D$Env))
summary(modelD)

D_D=subset(sub_D,sub_D$Sym_ID=="D")
D_AD=subset(sub_D,sub_D$Sym_ID=="AD")
D_BD=subset(sub_D,sub_D$Sym_ID=="BD")

par(mfrow=c(3,1))
hist(D_D$Log_D,breaks=4,xlim=c(0,6),main="D Only Treatment",xlab="")
hist(D_AD$Log_D,breaks=6,xlim=c(0,6),main="AD Treatment",xlab="")
hist(D_BD$Log_D,breaks=6,xlim=c(0,6),main="BD Treatment",xlab="Log Total Symbiodinium D")


##Time Variables
sub_A1=subset(sub_A,sub_A$Time==1)
sub_A2=subset(sub_A,sub_A$Time==2)
sub_B1=subset(sub_B,sub_B$Time==1)
sub_B2=subset(sub_B,sub_B$Time==2)
sub_D1=subset(sub_D,sub_D$Time==1)
sub_D2=subset(sub_D,sub_D$Time==2)

##Environment by Time
sub_A1L=subset(sub_A1,sub_A1$Env=="Light")
sub_A1D=subset(sub_A1,sub_A1$Env=="Dark")
sub_A2L=subset(sub_A2,sub_A2$Env=="Light")
sub_A2D=subset(sub_A2,sub_A2$Env=="Dark")

sub_B1L=subset(sub_B1,sub_B1$Env=="Light")
sub_B1D=subset(sub_B1,sub_B1$Env=="Dark")
sub_B2L=subset(sub_B2,sub_B2$Env=="Light")
sub_B2D=subset(sub_B2,sub_B2$Env=="Dark")

sub_D1L=subset(sub_D1,sub_D1$Env=="Light")
sub_D1D=subset(sub_D1,sub_D1$Env=="Dark")
sub_D2L=subset(sub_D2,sub_D2$Env=="Light")
sub_D2D=subset(sub_D2,sub_D2$Env=="Dark")

##FIGURES
##hist(~sub_A1L$Log_A|sub_A1L$Sym_ID,layout=c(1,3))
##quartz()
##hist(~sub_A1D$Log_A|sub_A1D$Sym_ID,layout=c(1,3))

#hist(~sub_A2L$Log_A|sub_A2L$Sym_ID,layout=c(1,3))
#quartz()
#hist(~sub_A2D$Log_A|sub_A2D$Sym_ID,layout=c(1,3))

#hist(~sub_B1L$Log_B|sub_B1L$Sym_ID,layout=c(1,3))
#quartz()
#hist(~sub_B1D$Log_B|sub_B1D$Sym_ID,layout=c(1,3))

#hist(~sub_B2L$Log_B|sub_B2L$Sym_ID,layout=c(1,3))
#quartz()
#hist(~sub_B2D$Log_B|sub_B2D$Sym_ID,layout=c(1,3))

#hist(~sub_D1L$Log_D|sub_D1L$Sym_ID,layout=c(1,3))
#quartz()
#hist(~sub_D1D$Log_D|sub_D1D$Sym_ID,layout=c(1,3))

#hist(~sub_D2L$Log_D|sub_D2L$Sym_ID,layout=c(1,3))
#quartz()
#hist(~sub_D2D$Log_D|sub_D2D$Sym_ID,layout=c(1,3))