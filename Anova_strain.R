##2015 June 17
##Make variables for log cell counts of individual strains, individual times, for mixed sym treatment

ANOVA_strain<-function(qPCR_all)

qPCR_all=read.table(file.choose(),header=T,sep="\t")
qPCR_all$Rep=as.factor(qPCR_all$Rep)
qPCR_T1=subset(qPCR_all,qPCR_all$Time==1)
qPCR_T2=subset(qPCR_all,qPCR_all$Time==2)

qPCR_T1_A=subset(qPCR_T1,qPCR_T1$Sym_ID=="AB"|qPCR_T1$Sym_ID=="AD")
qPCR_T2_A=subset(qPCR_T2,qPCR_T2$Sym_ID=="AB"|qPCR_T2$Sym_ID=="AD")

model_T1_A=anova(lm(qPCR_T1_A$Log_A~qPCR_T1_A$Sym_ID+qPCR_T1_A$Sym_ID/qPCR_T1_A$Env+qPCR_T1_A$Env/qPCR_T1_A$Rep))
model_T2_A=anova(lm(qPCR_T2_A$Log_A~qPCR_T2_A$Sym_ID+qPCR_T2_A$Sym_ID/qPCR_T2_A$Env+qPCR_T2_A$Env/qPCR_T2_A$Rep))

qPCR_T1_B=subset(qPCR_T1,qPCR_T1$Sym_ID=="AB"|qPCR_T1$Sym_ID=="BD")
qPCR_T2_B=subset(qPCR_T2,qPCR_T2$Sym_ID=="AB"|qPCR_T2$Sym_ID=="BD")

model_T1_B=anova(lm(qPCR_T1_B$Log_B~qPCR_T1_B$Sym_ID+qPCR_T1_B$Sym_ID/qPCR_T1_B$Env+qPCR_T1_B$Env/qPCR_T1_B$Rep))
model_T2_B=anova(lm(qPCR_T2_B$Log_B~qPCR_T2_B$Sym_ID+qPCR_T2_B$Sym_ID/qPCR_T2_B$Env+qPCR_T2_B$Env/qPCR_T2_B$Rep))

qPCR_T1_D=subset(qPCR_T1,qPCR_T1$Sym_ID=="AD"|qPCR_T1$Sym_ID=="BD")
qPCR_T2_D=subset(qPCR_T2,qPCR_T2$Sym_ID=="AD"|qPCR_T2$Sym_ID=="BD")

model_T1_D=anova(lm(qPCR_T1_D$Log_D~qPCR_T1_D$Sym_ID+qPCR_T1_D$Sym_ID/qPCR_T1_D$Env+qPCR_T1_D$Env/qPCR_T1_D$Rep))
model_T2_D=anova(lm(qPCR_T2_D$Log_D~qPCR_T2_D$Sym_ID+qPCR_T2_D$Sym_ID/qPCR_T2_D$Env+qPCR_T2_D$Env/qPCR_T2_D$Rep))

cat("\n","A at Time 1","\n")
print(model_T1_A)
cat("\n","A at Time 2","\n")
print(model_T2_A)

cat("\n","B at Time 1","\n")
print(model_T1_B)
cat("\n","B at Time 2","\n")
print(model_T2_B)

cat("\n","D at Time 1","\n")
print(model_T1_D)
cat("\n","D at Time 2","\n")
print(model_T2_D)

##PLOTS by Symbiont

A1plot<-bwplot(qPCR_T1_A$Log_A~qPCR_T1_A$Sym_ID|qPCR_T1_A$Env,ylab="Log A",main="A cells at 6 Weeks")
A2plot<-bwplot(qPCR_T2_A$Log_A~qPCR_T2_A$Sym_ID|qPCR_T2_A$Env,ylab="Log A",main="A cells at 8 Weeks")
B1plot<-bwplot(qPCR_T1_B$Log_B~qPCR_T1_B$Sym_ID|qPCR_T1_B$Env,ylab="Log B",main="B cells at 6 Weeks")
B2plot<-bwplot(qPCR_T2_B$Log_B~qPCR_T2_B$Sym_ID|qPCR_T2_B$Env,ylab="Log B",main="B cells at 8 Weeks")
D1plot<-bwplot(qPCR_T1_D$Log_D~qPCR_T1_D$Sym_ID|qPCR_T1_D$Env,ylab="Log D",main="D cells at 6 Weeks")
D2plot<-bwplot(qPCR_T2_D$Log_D~qPCR_T2_D$Sym_ID|qPCR_T2_D$Env,ylab="Log D",main="D cells at 8 Weeks")
print(A1plot,split=c(1,1,2,3),more=TRUE)
print(A2plot,split=c(2,1,2,3),more=TRUE)
print(B1plot,split=c(1,2,2,3),more=TRUE)
print(B2plot,split=c(2,2,2,3),more=TRUE)
print(D1plot,split=c(1,3,2,3),more=TRUE)
print(D2plot,split=c(2,3,2,3),more=F)

##PLOTS by Treatment
A1plot<-bwplot(qPCR_T1_A$Log_A~qPCR_T1_A$Env|qPCR_T1_A$Sym_ID,ylab="Log A",ylim=c(-0.5,6))
A2plot<-bwplot(qPCR_T2_A$Log_A~qPCR_T2_A$Env|qPCR_T2_A$Sym_ID,ylab="Log A",ylim=c(-0.5,6))
B1plot<-bwplot(qPCR_T1_B$Log_B~qPCR_T1_B$Env|qPCR_T1_B$Sym_ID,ylab="Log B",ylim=c(-0.5,6))
B2plot<-bwplot(qPCR_T2_B$Log_B~qPCR_T2_B$Env|qPCR_T2_B$Sym_ID,ylab="Log B",ylim=c(-0.5,6))
D1plot<-bwplot(qPCR_T1_D$Log_D~qPCR_T1_D$Env|qPCR_T1_D$Sym_ID,ylab="Log D",ylim=c(-0.5,6))
D2plot<-bwplot(qPCR_T2_D$Log_D~qPCR_T2_D$Env|qPCR_T2_D$Sym_ID,ylab="Log D",ylim=c(-0.5,6))
print(A1plot,split=c(1,1,2,3),more=TRUE)
print(A2plot,split=c(2,1,2,3),more=TRUE)
print(B1plot,split=c(1,2,2,3),more=TRUE)
print(B2plot,split=c(2,2,2,3),more=TRUE)
print(D1plot,split=c(1,3,2,3),more=TRUE)
print(D2plot,split=c(2,3,2,3),more=F)

##PLOT individual "only" inculations for comparison
#make variables
Aonly_T1=subset(qPCR_T1,qPCR_T1$Sym_ID=="A")
Aonly_T2=subset(qPCR_T2,qPCR_T2$Sym_ID=="A")
Bonly_T1=subset(qPCR_T1,qPCR_T1$Sym_ID=="B")
Bonly_T2=subset(qPCR_T2,qPCR_T2$Sym_ID=="B")
Donly_T1=subset(qPCR_T1,qPCR_T1$Sym_ID=="D")
Donly_T2=subset(qPCR_T2,qPCR_T2$Sym_ID=="D")
AonlyPlot_T1<-bwplot(Aonly_T1$Log_A~Aonly_T1$Env|Aonly_T1$Sym_ID,ylab="Log A",ylim=c(-0.5,6.5))
BonlyPlot_T1<-bwplot(Bonly_T1$Log_B~Bonly_T1$Env|Bonly_T1$Sym_ID,ylab="Log B",ylim=c(-0.5,6.5))
DonlyPlot_T1<-bwplot(Donly_T1$Log_D~Donly_T1$Env|Donly_T1$Sym_ID,ylab="Log D",ylim=c(-0.5,6.5))
AonlyPlot_T2<-bwplot(Aonly_T2$Log_A~Aonly_T2$Env|Aonly_T2$Sym_ID,ylab="Log A",ylim=c(-0.5,6.5))
BonlyPlot_T2<-bwplot(Bonly_T2$Log_B~Bonly_T2$Env|Bonly_T2$Sym_ID,ylab="Log B",ylim=c(-0.5,6.5))
DonlyPlot_T2<-bwplot(Donly_T2$Log_D~Donly_T2$Env|Donly_T2$Sym_ID,ylab="Log D",ylim=c(-0.5,6.5))
print(AonlyPlot_T1,split=c(1,1,2,3),more=TRUE)
print(AonlyPlot_T2,split=c(2,1,2,3),more=TRUE)
print(BonlyPlot_T1,split=c(1,2,2,3),more=TRUE)
print(BonlyPlot_T2,split=c(2,2,2,3),more=TRUE)
print(DonlyPlot_T1,split=c(1,3,2,3),more=TRUE)
print(DonlyPlot_T2,split=c(2,3,2,3),more=F)

##Make variables to plot all data together
qPCR_T1=subset(qPCR_all,qPCR_all$Time==1)
qPCR_T2=subset(qPCR_all,qPCR_all$Time==2)

plvar_T1_A=subset(qPCR_T1,qPCR_T1$Sym_ID=="AB"|qPCR_T1$Sym_ID=="AD"|qPCR_T1$Sym_ID=="A")
plvar_T2_A=subset(qPCR_T2,qPCR_T2$Sym_ID=="AB"|qPCR_T2$Sym_ID=="AD"|qPCR_T2$Sym_ID=="A")
plvar_T1_B=subset(qPCR_T1,qPCR_T1$Sym_ID=="AB"|qPCR_T1$Sym_ID=="BD"|qPCR_T1$Sym_ID=="B")
plvar_T2_B=subset(qPCR_T2,qPCR_T2$Sym_ID=="AB"|qPCR_T2$Sym_ID=="BD"|qPCR_T2$Sym_ID=="B")
plvar_T1_D=subset(qPCR_T1,qPCR_T1$Sym_ID=="AD"|qPCR_T1$Sym_ID=="BD"|qPCR_T1$Sym_ID=="D")
plvar_T2_D=subset(qPCR_T2,qPCR_T2$Sym_ID=="AD"|qPCR_T2$Sym_ID=="BD"|qPCR_T2$Sym_ID=="D")

A1plot_all<-bwplot(plvar_T1_A$Log_A~plvar_T1_A$Env|plvar_T1_A$Sym_ID,ylab="Log A",ylim=c(-0.5,6.5))
A2plot_all<-bwplot(plvar_T2_A$Log_A~plvar_T2_A$Env|plvar_T2_A$Sym_ID,ylab="Log A",ylim=c(-0.5,6.5))
B1plot_all<-bwplot(plvar_T1_B$Log_B~plvar_T1_B$Env|plvar_T1_B$Sym_ID,ylab="Log B",ylim=c(-0.5,6.5))
B2plot_all<-bwplot(plvar_T2_B$Log_B~plvar_T2_B$Env|plvar_T2_B$Sym_ID,ylab="Log B",ylim=c(-0.5,6.5))
D1plot_all<-bwplot(plvar_T1_D$Log_D~plvar_T1_D$Env|plvar_T1_D$Sym_ID,ylab="Log D",ylim=c(-0.5,6.5))
D2plot_all<-bwplot(plvar_T2_D$Log_D~plvar_T2_D$Env|plvar_T2_D$Sym_ID,ylab="Log D",ylim=c(-0.5,6.5))

print(A1plot_all,split=c(1,1,2,3),more=TRUE)
print(A2plot_all,split=c(2,1,2,3),more=TRUE)
print(B1plot_all,split=c(1,2,2,3),more=TRUE)
print(B2plot_all,split=c(2,2,2,3),more=TRUE)
print(D1plot_all,split=c(1,3,2,3),more=TRUE)
print(D2plot_all,split=c(2,3,2,3),more=F)

##PLOT ALL TOGETHER! AHH!
print(A1plot,split=c(1,1,4,3),more=TRUE)
print(AonlyPlot_T1,split=c(2,1,4,3),more=TRUE)
print(A2plot,split=c(3,1,4,3),more=TRUE)
print(AonlyPlot_T2,split=c(4,1,4,3),more=TRUE)
print(B1plot,split=c(1,2,4,3),more=TRUE)
print(BonlyPlot_T1,split=c(2,2,4,3),more=TRUE)
print(B2plot,split=c(3,2,4,3),more=TRUE)
print(BonlyPlot_T2,split=c(4,2,4,3),more=TRUE)
print(D1plot,split=c(1,3,4,3),more=TRUE)
print(DonlyPlot_T1,split=c(2,3,4,3),more=TRUE)
print(D2plot,split=c(3,3,4,3),more=TRUE)
print(DonlyPlot_T2,split=c(4,3,4,3),more=F)

##PLOT for proper figure spacing

##make individual variables
AB_T1_A=subset(qPCR_T1,qPCR_T1$Sym_ID=="AB")
AB_T2_A=subset(qPCR_T2,qPCR_T2$Sym_ID=="AB")
AD_T1_A=subset(qPCR_T1,qPCR_T1$Sym_ID=="AD")
AD_T2_A=subset(qPCR_T2,qPCR_T2$Sym_ID=="AD")
AB_T1_B=subset(qPCR_T1,qPCR_T1$Sym_ID=="AB")
AB_T2_B=subset(qPCR_T2,qPCR_T2$Sym_ID=="AB")
BD_T1_B=subset(qPCR_T1,qPCR_T1$Sym_ID=="BD")
BD_T2_B=subset(qPCR_T2,qPCR_T2$Sym_ID=="BD")
AD_T1_D=subset(qPCR_T1,qPCR_T1$Sym_ID=="AD")
AD_T2_D=subset(qPCR_T2,qPCR_T2$Sym_ID=="AD")
BD_T1_D=subset(qPCR_T1,qPCR_T1$Sym_ID=="BD")
BD_T2_D=subset(qPCR_T2,qPCR_T2$Sym_ID=="BD")
Aonly_T1=subset(qPCR_T1,qPCR_T1$Sym_ID=="A")
Aonly_T2=subset(qPCR_T2,qPCR_T2$Sym_ID=="A")
Bonly_T1=subset(qPCR_T1,qPCR_T1$Sym_ID=="B")
Bonly_T2=subset(qPCR_T2,qPCR_T2$Sym_ID=="B")
Donly_T1=subset(qPCR_T1,qPCR_T1$Sym_ID=="D")
Donly_T2=subset(qPCR_T2,qPCR_T2$Sym_ID=="D")


AonlyPlot_T1<-bwplot(Aonly_T1$Log_A~Aonly_T1$Env|Aonly_T1$Sym_ID,ylab="",ylim=c(-0.5,6.5))
ABPlot_T1_A<-bwplot(AB_T1_A$Log_A~AB_T1_A$Env|AB_T1_A$Sym_ID,ylab="",ylim=c(-0.5,6.5))
ADPlot_T1_A<-bwplot(AD_T1_A$Log_A~AD_T1_A$Env|AD_T1_A$Sym_ID,ylab="",ylim=c(-0.5,6.5))

BonlyPlot_T1<-bwplot(Bonly_T1$Log_B~Bonly_T1$Env|Bonly_T1$Sym_ID,ylab="",ylim=c(-0.5,6.5))
ABPlot_T1_B<-bwplot(AB_T1_B$Log_B~AB_T1_B$Env|AB_T1_B$Sym_ID,ylab="",ylim=c(-0.5,6.5))
BDPlot_T1_B<-bwplot(BD_T1_B$Log_B~BD_T1_B$Env|BD_T1_B$Sym_ID,ylab="",ylim=c(-0.5,6.5))

DonlyPlot_T1<-bwplot(Donly_T1$Log_D~Donly_T1$Env|Donly_T1$Sym_ID,ylab="",ylim=c(-0.5,6.5))
ADPlot_T1_D<-bwplot(AD_T1_D$Log_D~AD_T1_D$Env|AD_T1_D$Sym_ID,ylab="",ylim=c(-0.5,6.5))
BDPlot_T1_D<-bwplot(BD_T1_D$Log_D~BD_T1_D$Env|BD_T1_D$Sym_ID,ylab="",ylim=c(-0.5,6.5))

AonlyPlot_T2<-bwplot(Aonly_T2$Log_A~Aonly_T2$Env|Aonly_T2$Sym_ID,ylab="Log A",ylim=c(-0.5,6.5))
ABPlot_T2_A<-bwplot(AB_T2_A$Log_A~AB_T2_A$Env|AB_T2_A$Sym_ID,ylab="Log A",ylim=c(-0.5,6.5))
ADPlot_T2_A<-bwplot(AD_T2_A$Log_A~AD_T2_A$Env|AD_T2_A$Sym_ID,ylab="Log A",ylim=c(-0.5,6.5))

BonlyPlot_T2<-bwplot(Bonly_T2$Log_B~Bonly_T2$Env|Bonly_T2$Sym_ID,ylab="Log B",ylim=c(-0.5,6.5))
ABPlot_T2_B<-bwplot(AB_T2_B$Log_B~AB_T2_B$Env|AB_T2_B$Sym_ID,ylab="Log B",ylim=c(-0.5,6.5))
BDPlot_T2_B<-bwplot(BD_T2_B$Log_B~BD_T2_B$Env|BD_T2_B$Sym_ID,ylab="Log B",ylim=c(-0.5,6.5))

DonlyPlot_T2<-bwplot(Donly_T2$Log_D~Donly_T2$Env|Donly_T2$Sym_ID,ylab="Log D",ylim=c(-0.5,6.5))
ADPlot_T2_D<-bwplot(AD_T2_D$Log_D~AD_T2_D$Env|AD_T2_D$Sym_ID,ylab="Log D",ylim=c(-0.5,6.5))
BDPlot_T2_D<-bwplot(BD_T2_D$Log_D~BD_T2_D$Env|BD_T2_D$Sym_ID,ylab="Log D",ylim=c(-0.5,6.5))
##DATA AT 6 Weeks
print(ABPlot_T1_A,split=c(1,1,4,3),more=TRUE)
print(ADPlot_T1_A,split=c(2,1,4,3),more=TRUE)
print(AonlyPlot_T1,split=c(4,1,4,3),more=TRUE)
print(ABPlot_T1_B,split=c(1,2,4,3),more=TRUE)
print(BDPlot_T1_B,split=c(3,2,4,3),more=TRUE)
print(BonlyPlot_T1,split=c(4,2,4,3),more=TRUE)
print(ADPlot_T1_D,split=c(2,3,4,3),more=TRUE)
print(BDPlot_T1_D,split=c(3,3,4,3),more=TRUE)
print(DonlyPlot_T1,split=c(4,3,4,3),more=F)


##DATA AT 8 Weeks

print(ABPlot_T2_A,split=c(1,1,4,3),more=TRUE)
print(ADPlot_T2_A,split=c(2,1,4,3),more=TRUE)
print(AonlyPlot_T2,split=c(4,1,4,3),more=TRUE)
print(ABPlot_T2_B,split=c(1,2,4,3),more=TRUE)
print(BDPlot_T2_B,split=c(3,2,4,3),more=TRUE)
print(BonlyPlot_T2,split=c(4,2,4,3),more=TRUE)
print(ADPlot_T2_D,split=c(2,3,4,3),more=TRUE)
print(BDPlot_T2_D,split=c(3,3,4,3),more=TRUE)
print(DonlyPlot_T2,split=c(4,3,4,3),more=F)