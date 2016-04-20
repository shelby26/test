##Make variables of cells counts for cct1 and cct2 in qPCR dominance scripts pooling all times and environments into Sym_ID categories

##Separate Data for AB_only
AB_onlyT1=subset(qPCR,qPCR$Sym_ID=="AB"&qPCR$D_Count==0&qPCR$Time==1)
AD_onlyT1=subset(qPCR,qPCR$Sym_ID=="AD"&qPCR$B_Count==0&qPCR$Time==1)
BD_onlyT1=subset(qPCR,qPCR$Sym_ID=="BD"&qPCR$A_Count==0&qPCR$Time==1)
AB_onlyT2=subset(qPCR,qPCR$Sym_ID=="AB"&qPCR$D_Count==0&qPCR$Time==2)
AD_onlyT2=subset(qPCR,qPCR$Sym_ID=="AD"&qPCR$B_Count==0&qPCR$Time==2)
BD_onlyT2=subset(qPCR,qPCR$Sym_ID=="BD"&qPCR$A_Count==0&qPCR$Time==2)

ct1AB_A_Count_T1=AB_onlyT1$A_Count
ct1AB_A_Count_T1[ct1AB_A_Count_T1==0]<-1
ct2AB_B_Count_T1=AB_onlyT1$B_Count
ct2AB_B_Count_T1[ct2AB_B_Count_T1==0]<-1

ct1AD_A_Count_T1=AD_onlyT1$A_Count
ct1AD_A_Count_T1[ct1AD_A_Count_T1==0]<-1
ct2AD_D_Count_T1=AD_onlyT1$D_Count
ct2AD_D_Count_T1[ct2AD_D_Count_T1==0]<-1

ct1BD_B_Count_T1=BD_onlyT1$B_Count
ct1BD_B_Count_T1[ct1BD_B_Count_T1==0]<-1
ct2BD_D_Count_T1=BD_onlyT1$D_Count
ct2BD_D_Count_T1[ct2BD_D_Count_T1==0]<-1

ct1AB_A_Count_T2=AB_onlyT2$A_Count
ct1AB_A_Count_T2[ct1AB_A_Count_T2==0]<-1
ct2AB_B_Count_T2=AB_onlyT2$B_Count
ct2AB_B_Count_T2[ct2AB_B_Count_T2==0]<-1

ct1AD_A_Count_T2=AD_onlyT2$A_Count
ct1AD_A_Count_T2[ct1AD_A_Count_T2==0]<-1
ct2AD_D_Count_T2=AD_onlyT2$D_Count
ct2AD_D_Count_T2[ct2AD_D_Count_T2==0]<-1

ct1BD_B_Count_T2=BD_onlyT2$B_Count
ct1BD_B_Count_T2[ct1BD_B_Count_T2==0]<-1
ct2BD_D_Count_T2=BD_onlyT2$D_Count
ct2BD_D_Count_T2[ct2BD_D_Count_T2==0]<-1


##Run each subset through both qPCR dominance scripts (log based and percent based)
cat("***ALL AB ONLY T1***")
qPCR_percentdom(ct1AB_A_Count_T1,ct2AB_B_Count_T1,.95,100)
qPCR_dominance(ct1AB_A_Count_T1,ct2AB_B_Count_T1,.95,100)

cat("***ALL AB ONLY T2***")
qPCR_percentdom(ct1AB_A_Count_T2,ct2AB_B_Count_T2,.95,100)
qPCR_dominance(ct1AB_A_Count_T2,ct2AB_B_Count_T2,.95,100)

cat("***ALL AD ONLY T1***")
qPCR_percentdom(ct1AD_A_Count_T1,ct2AD_D_Count_T1,.95,100)
qPCR_dominance(ct1AD_A_Count_T1,ct2AD_D_Count_T1,.95,100)

cat("***ALL AD ONLY T2***")
qPCR_percentdom(ct1AD_A_Count_T2,ct2AD_D_Count_T2,.95,100)
qPCR_dominance(ct1AD_A_Count_T2,ct2AD_D_Count_T2,.95,100)

cat("***ALL BD ONLY***")
qPCR_percentdom(ct1BD_B_Count_T1,ct2BD_D_Count_T1,.95,100)
qPCR_dominance(ct1BD_B_Count_T1,ct2BD_D_Count_T1,.95,100)

cat("***ALL BD ONLY***")
qPCR_percentdom(ct1BD_B_Count_T2,ct2BD_D_Count_T2,.95,100)
qPCR_dominance(ct1BD_B_Count_T2,ct2BD_D_Count_T2,.95,100)

quartz()
par(mfrow=c(2,3))
hist(AB_onlyT1$Per_A,col="red",main="",xlab="Percent A in AB T1")
hist(AD_onlyT1$Per_A,col="red",main="",xlab="Percent A in AD T1")
hist(BD_onlyT1$Per_B,col="blue",main="",xlab="Percent B in BD T1")
hist(AB_onlyT1$Per_B,col="blue",main="",xlab="Percent B in AB T1")
hist(AD_onlyT1$Per_D,col="green",main="",xlab="Percent D in AD T1")
hist(BD_onlyT1$Per_D,col="green",main="",xlab="Percent D in BD T1")

quartz()
par(mfrow=c(2,3))
hist(AB_onlyT2$Per_A,col="red",main="",xlab="Percent A in AB T2")
hist(AD_onlyT2$Per_A,col="red",main="",xlab="Percent A in AD T2")
hist(BD_onlyT2$Per_B,col="blue",main="",xlab="Percent B in BD T2")
hist(AB_onlyT2$Per_B,col="blue",main="",xlab="Percent B in AB T2")
hist(AD_onlyT2$Per_D,col="green",main="",xlab="Percent D in AD T2")
hist(BD_onlyT2$Per_D,col="green",main="",xlab="Percent D in BD T2")
