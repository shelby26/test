
twoway1waySym<-function(qPCR_all)

qPCR_all=read.table(file.choose(),header=T,sep="\t") ## choose qPCR data for R.txt
qPCR_all$Rep=as.factor(qPCR_all$Rep)
##qPCR=subset(qPCR_all,qPCR_all$Time==1)

qPCR_Aonly=subset(qPCR,qPCR$Sym_ID=="A"&qPCR$B_Count==0&qPCR$D_Count==0)
qPCR_Bonly=subset(qPCR,qPCR$Sym_ID=="B"&qPCR$A_Count==0&qPCR$D_Count==0)
qPCR_Donly=subset(qPCR,qPCR$Sym_ID=="D"&qPCR$A_Count==0&qPCR$B_Count==0)
qPCR_ABonly=subset(qPCR,qPCR$Sym_ID=="AB"&qPCR$D_Count==0)
qPCR_BDonly=subset(qPCR,qPCR$Sym_ID=="BD"&qPCR$A_Count==0)
qPCR_ADonly=subset(qPCR,qPCR$Sym_ID=="AD"&qPCR$B_Count==0)

model2way<-aov(qPCR$Log_Total~qPCR$Sym_ID*qPCR$Env+Error(qPCR$Rep/qPCR$Sym_ID/qPCR$Env))
print(summary(model2way))

cat("\n","SymA","\n")
model1way_SymA<-aov(qPCR_Aonly$Log_Total~qPCR_Aonly$Env+Error(qPCR_Aonly$Rep/qPCR_Aonly$Env))
print(summary(model1way_SymA))

cat("\n","SymB","\n")
model1way_SymB<-aov(qPCR_Bonly$Log_Total~qPCR_Bonly$Env+Error(qPCR_Bonly$Rep/qPCR_Bonly$Env))
print(summary(model1way_SymB))

cat("\n","SymD","\n")
model1way_SymD<-aov(qPCR_Donly$Log_Total~qPCR_Donly$Env+Error(qPCR_Donly$Rep/qPCR_Donly$Env))
print(summary(model1way_SymD))

cat("\n","SymAB","\n")
model1way_SymAB<-aov(qPCR_ABonly$Log_Total~qPCR_ABonly$Env+Error(qPCR_ABonly$Rep/qPCR_ABonly$Env))
print(summary(model1way_SymAB))

cat("\n","SymAD","\n")
model1way_SymAD<-aov(qPCR_ADonly$Log_Total~qPCR_ADonly$Env+Error(qPCR_ADonly$Rep/qPCR_ADonly$Env))
print(summary(model1way_SymAD))

cat("\n","SymBD","\n")
model1way_SymBD<-aov(qPCR_BDonly$Log_Total~qPCR_BDonly$Env+Error(qPCR_BDonly$Rep/qPCR_BDonly$Env))
print(summary(model1way_SymBD))