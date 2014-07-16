#Growth curve plot
par(mfrow=c(2,1))
plot(VBGF.ex[,2],VBGF.ex[,1],pch=21,bg="blue1",xlab="",ylab="")
points(VBGF.ex[,3],VBGF.ex[,1],pch=21,bg="blue4")
points(VBGF.ex[,4],VBGF.ex[,1],pch=21,bg="lightblue")
lines(c(0:13),vbgf.fits,pch=21,bg="black",lwd=3)

#Plot Age samples
barplot(rev(age.samples.dat[,2]),horiz=T, names.arg=abbreviate(rev(age.samples.dat[,1])),las=2)
barplot(rev(age.samples.dat[,3]),horiz=T,col="red",add=T,las=2)

#Left-right comps
otowt.comp.spt<-read.table('clipboard',header=T)
otowt.comp.sbl<-read.table('clipboard',header=T)
otowt.comp.ptl<-read.table('clipboard',header=T)
plot(otowt.comp.spt$Left/max(otowt.comp.spt),otowt.comp.spt$Right/max(otowt.comp.spt),xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",pch=21,col="black",bg="red")
points(otowt.comp.sbl$Left/max(otowt.comp.sbl),otowt.comp.sbl$Right/max(otowt.comp.sbl),xlab="",ylab="",pch=21,col="black",bg="yellow")
points(otowt.comp.ptl$Blind/max(otowt.comp.ptl),otowt.comp.ptl$Eyed/max(otowt.comp.ptl),xlab="",ylab="",pch=24,col="black",bg="green")
abline(a=0,b=1,col="black")

#Piecewise/Linear model filts
par(mfrow=c(2,2))
OtoAge.plot(Spp.Lms$Hake.Lms.300,2)
OtoAge.plot(Spp.Lms$Hake.Lms.300,3)
OtoAge.plot(Spp.Lms$Petrale.Lms.80,2)
OtoAge.plot(Spp.Lms$Petrale.Lms.80,3)
par(mfrow=c(2,2))
OtoAge.plot(Spp.Lms$Sablefish,2)
OtoAge.plot(Spp.Lms$Sablefish,3)
OtoAge.plot(Spp.Lms$Splitnose,2)
OtoAge.plot(Spp.Lms$Splitnose,3)

#Comparison of piecewise linear model by area 
par(mfrow=c(2,2))
OtoAge.plot(Sab.nsei_ssei.Lms,2)
points(Spp.Lms$Sablefish[[1]]$Females,pch=21,bg="red")
OtoAge.plot(Sab.nsei_ssei.Lms,3)
points(Spp.Lms$Sablefish[[1]]$Males,pch=21,bg="red")
par(mfrow=c(2,2))
LM.comp.plots(Sab.nsei_ssei.Lms,2)
LM.comp.plots(Spp.Lms$Sablefish,2,add.plot=T,col.bpts=c("red","black"))
LM.comp.plots(Sab.nsei_ssei.Lms,3)
LM.comp.plots(Spp.Lms$Sablefish,3,add.plot=T,col.bpts=c("red","black"))

#1:1 Plots
par(mfrow=c(2,2))
plot.age2age(Spp.age.pred,"Hake300",gen.type=2,age.comp.type=3,rd.type=1,col.dots="red",lab.nam=c("F","T"))
plot.age2age(Spp.age.pred,"Hake300",gen.type=3,age.comp.type=3,rd.type=1,col.dots="blue",lab.nam=c("F","F"),add.plot="T")
plot.age2age(Spp.age.pred,"Petrale80",gen.type=2,age.comp.type=3,rd.type=1,col.dots="red",lab.nam=c("F","F"))
plot.age2age(Spp.age.pred,"Petrale80",gen.type=3,age.comp.type=3,rd.type=1,col.dots="blue",lab.nam=c("F","F"),add.plot="T")
plot.age2age(Spp.age.pred,"Sablefish",gen.type=2,age.comp.type=3,rd.type=1,col.dots="red",lab.nam=c("T","T"))
plot.age2age(Spp.age.pred,"Sablefish",gen.type=3,age.comp.type=3,rd.type=1,col.dots="blue",lab.nam=c("F","F"),add.plot="T")
plot.age2age(Spp.age.pred,"Splitnose",gen.type=2,age.comp.type=3,rd.type=1,col.dots="red",lab.nam=c("T","F"))
plot.age2age(Spp.age.pred,"Splitnose",gen.type=3,age.comp.type=3,rd.type=1,col.dots="blue",lab.nam=c("F","F"),add.plot="T")

#Parameter comparisons
par(mfrow=c(2,1),mar=c(4,4,3,2))
vbgf.par.comp.plot(Hake300.vbgf,c(8:21),1,labelx="N")
vbgf.par.comp.plot(Hake300.vbgf,c(8:21),2,labelx="N")
par(mfrow=c(2,1),mar=c(4,4,3,2))
vbgf.par.comp.plot(Petrale80.vbgf,c(8:21),1,labelx="N")
vbgf.par.comp.plot(Petrale80.vbgf,c(8:21),2,labelx="N")
par(mfrow=c(2,1),mar=c(4,4,3,2))
vbgf.par.comp.plot(Sablefish.vbgf,c(8:21),1,labelx="N")
vbgf.par.comp.plot(Sablefish.vbgf,c(8:21),2,labelx="N")
par(mfrow=c(2,1),mar=c(4,4,3,2))
vbgf.par.comp.plot(Splitnose.vbgf,c(8:21),1,labelx="N")
vbgf.par.comp.plot(Splitnose.vbgf,c(8:21),2,labelx="N")

#Parameter comparisons- redux
#Linf
plotCI(c(1:8)-0.1,Spp.RE.Linf[1,],uiw=Spp.RE.Linf[2,],ylim=c(-1,1),pch=21,pt.bg="red",xlab="",ylab=expression(L[inf]),axes=F,cex=1.25)
box()
axis(1,at=c(1:dim(Spp.RE.Linf)[2]),labels=rep(c("PW","LM"),4),las=2,cex=0.75)
axis(2)
abline(h=0,lty=2,col="black")
abline(v=seq(2.5,6.5,by=2),lwd=1,col="black")
plotCI(c(1:8)+0.1,Spp.RE.Linf[3,],uiw=Spp.RE.Linf[4,],ylim=c(-0.5,0.5),pch=21,pt.bg="blue",xlab="",ylab=expression(L[inf]),cex=1.25,add=T)
#k
plotCI(c(1:8)-0.1,Spp.RE.k[1,],uiw=-Spp.RE.k[2,],ylim=c(-1,1),pch=21,pt.bg="red",xlab="",ylab="k",axes=F,cex=1.25)
box()
axis(1,at=c(1:dim(Spp.RE.k)[2]),labels=rep(c("PW","LM"),4),las=2,cex=0.75)
axis(2)
abline(h=0,lty=2,col="black")
abline(v=seq(2.5,6.5,by=2),lwd=1,col="black")
plotCI(c(1:8)+0.1,Spp.RE.k[3,],uiw=-Spp.RE.k[4,],ylim=c(-0.5,0.5),pch=21,pt.bg="blue",xlab="",ylab="k",cex=1.25,add=T)

#VBGF fits
#Hake300
par(mfrow=c(2,1),mar=c(2,4,2,3))
VBGF.fit.plot(Spp.age.pred$Hake300$F$Age,Spp.age.pred$Hake300$F$Length,ceiling(max(Spp.age.pred$Hake300$F[,c(2:4)],na.rm=T)/10)*10,ceiling(max(Spp.age.pred$Hake300$F[,5],na.rm=T)/10)*10,"red","",add.plot="F",print.vbgf="F")
VBGF.fit.plot(Spp.age.pred$Hake300$F$PW_age,Spp.age.pred$Hake300$F$Length,ceiling(max(Spp.age.pred$Hake300$F[,c(2:4)],na.rm=T)/10)*10,ceiling(max(Spp.age.pred$Hake300$F[,5],na.rm=T)/10)*10,"pink","",add.plot="T",print.vbgf="F")
VBGF.fit.plot(Spp.age.pred$Hake300$M$Age,Spp.age.pred$Hake300$M$Length,ceiling(max(Spp.age.pred$Hake300$M[,c(2:4)],na.rm=T)/10)*10,ceiling(max(Spp.age.pred$Hake300$M[,5],na.rm=T)/10)*10,"blue","",add.plot="F",print.vbgf="F")
VBGF.fit.plot(Spp.age.pred$Hake300$M$PW_age,Spp.age.pred$Hake300$M$Length,ceiling(max(Spp.age.pred$Hake300$M[,c(2:4)],na.rm=T)/10)*10,ceiling(max(Spp.age.pred$Hake300$M[,5],na.rm=T)/10)*10,"lightblue","",add.plot="T",print.vbgf="F")
#Petrale80
par(mfrow=c(2,1),mar=c(2,4,2,3))
VBGF.fit.plot(Spp.age.pred$Petrale80$F$Age,Spp.age.pred$Petrale80$F$Length,ceiling(max(Spp.age.pred$Petrale80$F[,c(2:4)],na.rm=T)/10)*10,ceiling(max(Spp.age.pred$Petrale80$F[,5],na.rm=T)/10)*10,"red","",add.plot="F",print.vbgf="F")
VBGF.fit.plot(Spp.age.pred$Petrale80$F$PW_age,Spp.age.pred$Petrale80$F$Length,ceiling(max(Spp.age.pred$Petrale80$F[,c(2:4)],na.rm=T)/10)*10,ceiling(max(Spp.age.pred$Petrale80$F[,5],na.rm=T)/10)*10,"pink","",add.plot="T",print.vbgf="F")
VBGF.fit.plot(Spp.age.pred$Petrale80$M$Age,Spp.age.pred$Petrale80$M$Length,ceiling(max(Spp.age.pred$Petrale80$M[,c(2:4)],na.rm=T)/10)*10,ceiling(max(Spp.age.pred$Petrale80$M[,5],na.rm=T)/10)*10,"blue","",add.plot="F",print.vbgf="F")
VBGF.fit.plot(Spp.age.pred$Petrale80$M$PW_age,Spp.age.pred$Petrale80$M$Length,ceiling(max(Spp.age.pred$Petrale80$M[,c(2:4)],na.rm=T)/10)*10,ceiling(max(Spp.age.pred$Petrale80$M[,5],na.rm=T)/10)*10,"lightblue","",add.plot="T",print.vbgf="F")
#Sablefish
par(mfrow=c(2,1),mar=c(2,4,2,3))
VBGF.fit.plot(Spp.age.pred$Sablefish$F$Age,Spp.age.pred$Sablefish$F$Length,ceiling(max(Spp.age.pred$Sablefish$F[,c(2:4)],na.rm=T)/10)*10,ceiling(max(Spp.age.pred$Sablefish$F[,5],na.rm=T)/10)*10,"red","",add.plot="F",print.vbgf="F")
VBGF.fit.plot(Spp.age.pred$Sablefish$F$PW_age,Spp.age.pred$Sablefish$F$Length,ceiling(max(Spp.age.pred$Sablefish$F[,c(2:4)],na.rm=T)/10)*10,ceiling(max(Spp.age.pred$Sablefish$F[,5],na.rm=T)/10)*10,"pink","",add.plot="T",print.vbgf="F")
VBGF.fit.plot(Spp.age.pred$Sablefish$M$Age,Spp.age.pred$Sablefish$M$Length,ceiling(max(Spp.age.pred$Sablefish$M[,c(2:4)],na.rm=T)/10)*10,ceiling(max(Spp.age.pred$Sablefish$M[,5],na.rm=T)/10)*10,"blue","",add.plot="F",print.vbgf="F")
VBGF.fit.plot(Spp.age.pred$Sablefish$M$PW_age,Spp.age.pred$Sablefish$M$Length,ceiling(max(Spp.age.pred$Sablefish$M[,c(2:4)],na.rm=T)/10)*10,ceiling(max(Spp.age.pred$Sablefish$M[,5],na.rm=T)/10)*10,"lightblue","",add.plot="T",print.vbgf="F")
#Splitnose
par(mfrow=c(2,1),mar=c(2,4,2,3))
VBGF.fit.plot(Spp.age.pred$Splitnose$F$Age,Spp.age.pred$Splitnose$F$Length,ceiling(max(Spp.age.pred$Splitnose$F[,c(2:4)],na.rm=T)/10)*10,ceiling(max(Spp.age.pred$Splitnose$F[,5],na.rm=T)/10)*10,"red","",add.plot="F",print.vbgf="F")
VBGF.fit.plot(Spp.age.pred$Splitnose$F$PW_age,Spp.age.pred$Splitnose$F$Length,ceiling(max(Spp.age.pred$Splitnose$F[,c(2:4)],na.rm=T)/10)*10,ceiling(max(Spp.age.pred$Splitnose$F[,5],na.rm=T)/10)*10,"pink","",add.plot="T",print.vbgf="F")
VBGF.fit.plot(Spp.age.pred$Splitnose$M$Age,Spp.age.pred$Splitnose$M$Length,ceiling(max(Spp.age.pred$Splitnose$M[,c(2:4)],na.rm=T)/10)*10,ceiling(max(Spp.age.pred$Splitnose$M[,5],na.rm=T)/10)*10,"blue","",add.plot="F",print.vbgf="F")
VBGF.fit.plot(Spp.age.pred$Splitnose$M$PW_age,Spp.age.pred$Splitnose$M$Length,ceiling(max(Spp.age.pred$Splitnose$M[,c(2:4)],na.rm=T)/10)*10,ceiling(max(Spp.age.pred$Splitnose$M[,5],na.rm=T)/10)*10,"lightblue","",add.plot="T",print.vbgf="F")

Spp.AE.name<-levels(age.error.in$Species)
par(mfrow=c(2,2),mar=c(5,2,4,2))
plot.AErr.comp(age.error.in,Spp.AE.name[1],col="black")
plot.AErr.comp(age.error.in,Spp.AE.name[2],col="blue")
par(mfrow=c(2,2),mar=c(5,2,4,2))
plot.AErr.comp(age.error.in,Spp.AE.name[3])
plot.AErr.comp(age.error.in,Spp.AE.name[4])
par(mfrow=c(2,2),mar=c(5,2,4,2))
plot.AErr.comp(age.error.in,Spp.AE.name[5])
plot.AErr.comp(age.error.in,Spp.AE.name[6])
par(mfrow=c(2,2),mar=c(5,2,4,2))
plot.AErr.comp(age.error.in,Spp.AE.name[7])
plot.AErr.comp(age.error.in,Spp.AE.name[8])