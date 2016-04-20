##2 JULY 2015
##Make figure of Log Total Cells by Symbiont and Light Treatment

quartz()
par(mfrow=c(2,1),mar=c(3,4,3,2))
boxplot(qPCR_T1$Log_Total~qPCR_T1$Env*qPCR_T1$Sym_ID,col=c("grey","gold"),main="6 Weeks After Inoculation",ylab=expression(paste("Log (", italic("Symbiodinium "), "cells per polyp)")),xaxt='n')
 xlabel<-c("A","","AB","","AD","","B","","BD","","D","")
 axis(labels=xlabel,side=1,at=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5))
 abline(v=c(2.5,4.5,6.5,8.5,10.5))
 means1<-c(3.201,3.597,3.332,3.400,3.044,3.879,3.444,3.405,3.056,3.832,3.603,4.041)
 points(means1,pch=16)
##Significance
 points(x=5.5,y=4.8,pch=8,col="red")
 points(x=9.5,y=4.8,pch=8,col="red")

 boxplot(qPCR_T2$Log_Total~qPCR_T2$Env*qPCR_T2$Sym_ID,col=c("grey","gold"),main="8 Weeks After Inoculation",ylab=expression(paste("Log (", italic("Symbiodinium "), "cells per polyp)")),xaxt='n')
 xlabel<-c("A","","AB","","AD","","B","","BD","","D","")
 axis(labels=xlabel,side=1,at=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5))
 abline(v=c(2.5,4.5,6.5,8.5,10.5))
 means2<-c(4.247,4.435,3.630,4.020,4.594,4.375,4.032,4.417,4.018,4.657,4.499,4.715)
points(means2,pch=16)

##SINGLE CELL TREATEMENTS ONLY
quartz()
boxplot(qPCR_T1_1$Log_Total~qPCR_T1_1$Env*qPCR_T1_1$Sym_ID,col=c("grey","gold"),main="6 Weeks After Inoculation",ylab=expression(paste("Log (", italic("Symbiodinium "), "cells per polyp)")),xaxt='n')
 xlabel<-c("A","","B","","D","")
 axis(labels=xlabel,side=1,at=c(1.5,2.5,3.5,4.5,5.5,6.5))
 abline(v=c(2.5,4.5))
 means=c(3.200667,3.597167,3.444500,3.405429,3.602818,4.040714)
 points(means,pch=16)
 
 quartz()
boxplot(qPCR_T2_1$Log_Total~qPCR_T2_1$Env*qPCR_T2_1$Sym_ID,col=c("grey","gold"),main="8 Weeks After Inoculation",ylab=expression(paste("Log (", italic("Symbiodinium "), "cells per polyp)")),xaxt='n')
 xlabel<-c("A","","B","","D","")
 axis(labels=xlabel,side=1,at=c(1.5,2.5,3.5,4.5,5.5,6.5))
 abline(v=c(2.5,4.5))
 means=c(4.247500,4.435538,4.03225,4.41675,4.499000,4.714929)
 points(means,pch=16)