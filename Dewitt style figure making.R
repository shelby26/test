##Take subsets of data and make Dewitt Style comparisons of total cells in single vs. multiple infection treatments
##Load lattice for bwplots 


##Separate data by figure
AB_6_Dark=subset(qPCR,qPCR$FIG=="AB_6_Dark")
AB_6_Dark$FIG<-factor(AB_6_Dark$FIG)
AB_6_Dark$SYM<-factor(AB_6_Dark$SYM)
AB_6_Dark_means=tapply(AB_6_Dark$CELLS,AB_6_Dark$SYM,mean)

AB_6_Light=subset(qPCR,qPCR$FIG=="AB_6_Light")
AB_6_Light$FIG<-factor(AB_6_Light$FIG)
AB_6_Light$SYM<-factor(AB_6_Light$SYM)
AB_6_Light_means=tapply(AB_6_Light$CELLS,AB_6_Light$SYM,mean)

AD_6_Dark=subset(qPCR,qPCR$FIG=="AD_6_Dark")
AD_6_Dark$FIG<-factor(AD_6_Dark$FIG)
AD_6_Dark$SYM<-factor(AD_6_Dark$SYM)
AD_6_Dark_means=tapply(AD_6_Dark$CELLS,AD_6_Dark$SYM,mean)

AD_6_Light=subset(qPCR,qPCR$FIG=="AD_6_Light")
AD_6_Light$FIG<-factor(AD_6_Light$FIG)
AD_6_Light$SYM<-factor(AD_6_Light$SYM)
AD_6_Light_means=tapply(AD_6_Light$CELLS,AD_6_Light$SYM,mean)

BD_6_Dark=subset(qPCR,qPCR$FIG=="BD_6_Dark")
BD_6_Dark$FIG<-factor(BD_6_Dark$FIG)
BD_6_Dark$SYM<-factor(BD_6_Dark$SYM)
BD_6_Dark_means=tapply(BD_6_Dark$CELLS,BD_6_Dark$SYM,mean)

BD_6_Light=subset(qPCR,qPCR$FIG=="BD_6_Light")
BD_6_Light$FIG<-factor(BD_6_Light$FIG)
BD_6_Light$SYM<-factor(BD_6_Light$SYM)
BD_6_Light_means=tapply(BD_6_Light$CELLS,BD_6_Light$SYM,mean)

par(mfrow=c(1,2))
boxplot(AB_6_Dark$CELLS~AB_6_Dark$SYM,main="A vs. B effects in Low Light",ylab="# cells")
points(AB_6_Dark_means,col="red",pch=16)
boxplot(AB_6_Light$CELLS~AB_6_Light$SYM,main="A vs. B effects in High Light",ylab="# cells")
points(AB_6_Light_means,col="red",pch=16)