#Prep data
#Spp.dat.otos<-read.table('clipboard',header=T)
#save(Spp.dat.otos,file="C://Users//david.lin//documents-export-2014-07-02//Spp_dat_otos.DMP")
#load("C://Users//david.lin//documents-export-2014-07-02//Spp_dat_otos.DMP")
load("C://Users//DavidMarkLin//Downloads//documents-export-2014-09-07//Spp_dat_otos.DMP")
#Subset species data
Hake.dat.otos<-subset(Spp.dat.otos,Species=="Hake")
Petrale.dat.otos<-subset(Spp.dat.otos,Species=="Petrale")
Petrale.dat.otos$Sex<-as.factor(toupper(Petrale.dat.otos$Sex))
Sablefish.dat.otos<-subset(Spp.dat.otos,Species=="Sablefish")
Splitnose.dat.otos<-subset(Spp.dat.otos,Species=="Splitnose")
Splitnose.dat.otos$Sex<-as.factor(toupper(Splitnose.dat.otos$Sex))

#Ageing lab data
Age.lab.NWFSC.dat<-read.table('clipboard',header=T)
NWFSC.sexlt.dat<-read.table('clipboard',header=T)
NWFSC.Age.Sex.Lt.dat<-match.f(Age.lab.NWFSC.dat,NWFSC.sexlt.dat,"FDBarcode","OTOLITH_BARCODE",c("FISH_SEX","LENGTH_CM"))
NWFSC.Age.Sex.Lt.dat.spp<-list()
spp.names.NWFSC.dat<-unique(NWFSC.Age.Sex.Lt.dat$PacFINCode)
for(i in 1:length(spp.names.NWFSC.dat)){NWFSC.Age.Sex.Lt.dat.spp[[i]]<-subset(NWFSC.Age.Sex.Lt.dat,PacFINCode==spp.names.NWFSC.dat[i])}
names(NWFSC.Age.Sex.Lt.dat.spp)<-spp.names.NWFSC.dat
save(NWFSC.Age.Sex.Lt.dat.spp,file="C:/Users/copeja/Desktop/Current projects/otolith weight/Otolith_weight_age/oto-weight/NWFSC_Age_Sex_Lt_dat_spp.DMP")
load("C:/Users/copeja/Desktop/Current projects/otolith weight/Otolith_weight_age/oto-weight/NWFSC_Age_Sex_Lt_dat_spp.DMP")

all.spp.all.surveys<-read.table('clipboard',header=T)
hake.all<-subset(all.spp.all.surveys,PacFINCode=="PWHT")
##########################
### RUN and FIT models ###
##########################
mapply(function(x) quantile(Spp.dat.AFM.mod[,1],,na.rm=T),x=1:3)
############
### Hake ###
############
#Check summary stats of otolith weight at age
Hake.Bps.dat<-Oto.Age.Model.fits(Hake.dat.otos,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=c(10,300),steppin=10,breakpts)$Data
par(mfrow=c(2,2))
hake.summ.All<-boxplot(Hake.Bps.dat$A[,1]~Hake.Bps.dat$A[,2])
hake.summ.F<-boxplot(Hake.Bps.dat$F[,1]~Hake.Bps.dat$F[,2])
hake.summ.M<-boxplot(Hake.Bps.dat$M[,1]~Hake.Bps.dat$M[,2])
#Calculate breakpoints
Hake.Bps.1<-Oto.Age.Model.fits(Hake.dat.otos,oto.age.col=c(5,4,3),sextype="All",Bp.find="T",rngSplit=list(c(10,300),c(10,300),c(10,300)),steppin=c(1,1,1),lowbreaks=c(0,0,0),jitter=0,one.to.one="F",bps=c(0.2,0.6),quad.fit= "F",add.lengths="F")
Hake.Bps.1b<-Oto.Age.Model.fits(Hake.dat.otos,oto.age.col=c(5,4,3),sextype="All",Bp.find="T",rngSplit=list(c(10,300),c(10,300),c(10,300)),steppin=c(1,1,1),lowbreaks=c(0,0,0),jitter=0,one.to.one="F",bps=c(0.05,0.1,0.8),quad.fit= "F",add.lengths="F")
Hake.Bps.2a<-Oto.Age.Model.fits(Hake.dat.otos,oto.age.col=c(5,4,3),sextype="All",Bp.find="T",rngSplit=list(c(10,273),c(10,280),c(10,242)),steppin=c(1,1,1),lowbreaks=c(0,0,0),highbreaks=c(273,288,242),jitter=0,one.to.one="F",bps=c(0.2,0.6),quad.fit="T",add.lengths="T")
Hake.Bps.2b<-Oto.Age.Model.fits(Hake.dat.otos,oto.age.col=c(5,4,3),sextype="All",Bp.find="T",rngSplit=list(c(10,273),c(10,280),c(10,242)),steppin=c(1,1,1),lowbreaks=c(0,0,0),highbreaks=c(273,330,242),jitter=0,one.to.one="F",bps=c(0.2,0.8),quad.fit="T",add.lengths="T")
Hake.Bps.2<-Oto.Age.Model.fits(Hake.dat.otos,oto.age.col=c(5,4,3),sextype="All",Bp.find="T",rngSplit=list(c(50,273),c(90,280),c(50,242)),steppin=c(1,1,1),lowbreaks=c(0,0,0),highbreaks=c(273,288,242),jitter=0,one.to.one="F",bps=0.2,quad.fit="T",add.lengths="T")
Hake.Bps.3<-Oto.Age.Model.fits(Hake.dat.otos,oto.age.col=c(5,4,3),sextype="All",Bp.find="T",rngSplit=list(c(50,273),c(90,288),c(50,242)),steppin=c(1,1,1),lowbreaks=c(0,0,0),highbreaks=c(273,288,242),jitter=100,one.to.one="T",bps=0.2,add.lengths="F",quad.fit="F",comp.plot="F",comp.num=list(c(1,1),c(0,1),c(1,0)))
plot(Hake.Bps.3$Data_for_comps$Females[,1],Hake.Bps.3$Data_for_comps$Females[,3],pch=21,col="black",bg="red",xlab="Annuli ages",ylab="Oto Wt ages",xlim=c(0,max(Hake.Bps.3$Data_for_comps$Females)),ylim=c(0,max(Hake.Bps.3$Data_for_comps$Females)))
abline(a=0,b=1,col="black",lwd=2)
plot(Hake.Bps.3$Data_for_comps$Males[,1],Hake.Bps.3$Data_for_comps$Males[,3],pch=21,col="black",bg="blue",xlab="Annuli ages",ylab="Oto Wt ages",xlim=c(0,max(Hake.Bps.3$Data_for_comps$Females)),ylim=c(0,max(Hake.Bps.3$Data_for_comps$Females)))
abline(a=0,b=1,col="black",lwd=2)

female.hake<-subset(NWFSC.Age.Sex.Lt.dat.spp$PWHT,FISH_SEX=="f")
female.hake.mod<-female.hake$OtolithWeight[female.hake$OtolithWeight*1000<288]
femhake_pred<-predict(Hake.Bps.3$lms.out$Females,newdata=data.frame(OtoWt=female.hake.mod*1000))
plot(female.hake$FinalAge[female.hake$OtolithWeight*1000<288],femhake_pred,pch=21,col="black",bg="red",xlab="Annuli ages",ylab="Oto Wt ages",xlim=c(0,max(female.hake.mod,femhake_pred)),ylim=c(0,max(female.hake.mod,femhake_pred)))
abline(a=0,b=1,col="black",lwd=2)

femhake_pred<-predict(Hake.Bps.3$lms.out$Females,newdata=data.frame(OtoWt=hake.all$OtolithWeight*1000))
plot(hake.all$FinalAge,femhake_pred,pch=21,col="black",bg="red",xlab="Annuli ages",ylab="Oto Wt ages",xlim=c(0,max(female.hake$FinalAge,femhake_pred)),ylim=c(0,max(female.hake$FinalAge,femhake_pred)))
abline(a=0,b=1,col="black",lwd=2)


#Hake.Bps.F<-Oto.Age.Model.fits(Hake.dat.otos,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=c(91,91),steppin=1,lowbreaks=c(50,90,50),mod.select == "T")
#Hake.Bps.M<-Oto.Age.Model.fits(Hake.dat.otos,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=c(130,140),steppin=0.1,lowbreaks=c(50,90,50))
Hake.Lms<-Oto.Age.Model.fits(Hake.dat.otos,oto.age.col=c(5,4),sextype="All",Bp.find="F",rngSplit=c(0,100),steppin=10,breakpts=c(Hake.Bps.A$Bps$All$par[1],Hake.Bps.F$Bps$Females$par[1],Hake.Bps.M$Bps$Males$par[1]),lowbreaks=c(50,90,50))
#Fit linear models
par(mfrow=c(2,2))
OtoAge.plot(Hake.Lms,1)
OtoAge.plot(Hake.Lms,2)
OtoAge.plot(Hake.Lms,3)
#Hake <300mg -- use this one for now
#Check summary stats of otolith weight at age
Hake.dat.otos.300<-subset(Hake.dat.otos,OtoWt<300)
Hake.Bps.10_300.dat<-Oto.Age.Model.fits(Hake.dat.otos.300,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=c(10,300),steppin=10)$Data
par(mfrow=c(2,2))
hake300.summ.All<-boxplot(Hake.Bps.10_300.dat$A[,1]~Hake.Bps.10_300.dat$A[,2])
hake300.summ.F<-boxplot(Hake.Bps.10_300.dat$F[,1]~Hake.Bps.10_300.dat$F[,2])
hake300.summ.M<-boxplot(Hake.Bps.10_300.dat$M[,1]~Hake.Bps.10_300.dat$M[,2])
#Calculate breakpoints
Hake.Bps.A.300<-Oto.Age.Model.fits(Hake.dat.otos.300,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=c(135,140),steppin=0.01,lowbreaks=c(50,90,50))
Hake.Bps.F.300<-Oto.Age.Model.fits(Hake.dat.otos.300,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=c(120,130),steppin=0.1,lowbreaks=c(50,90,50))
Hake.Bps.M.300<-Oto.Age.Model.fits(Hake.dat.otos.300,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=c(135,140),steppin=0.01,lowbreaks=c(50,90,50))
Hake.Lms.300<-Oto.Age.Model.fits(Hake.dat.otos.300,oto.age.col=c(5,4),sextype="All",Bp.find="F",rngSplit=c(0,100),steppin=10,breakpts=c(Hake.Bps.A.300$Bps$All$par[1],Hake.Bps.F.300$Bps$Females$par[1],Hake.Bps.M.300$Bps$Males$par[1]),lowbreaks=c(50,90,50))
#Fit linear models
par(mfrow=c(2,2))
OtoAge.plot(Hake.Lms.300,1)
OtoAge.plot(Hake.Lms.300,2)
OtoAge.plot(Hake.Lms.300,3)
############

###############
### Petrale ###
###############
#Check summary stats of otolith weight at age
Petrale.Bps.dat<-Oto.Age.Model.fits(Petrale.dat.otos,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=c(0,100),steppin=10)$Data
par(mfrow=c(2,2))
petrale.summ.All<-boxplot(Petrale.Bps.dat$A[,1]~Petrale.Bps.dat$A[,2])
petrale.summ.F<-boxplot(Petrale.Bps.dat$F[,1]~Petrale.Bps.dat$F[,2])
petrale.summ.M<-boxplot(Petrale.Bps.dat$M[,1]~Petrale.Bps.dat$M[,2])
#Calculate breakpoints
Petrale.Bps.1<-Oto.Age.Model.fits(Petrale.dat.otos,oto.age.col=c(5,4,3),sextype="All",Bp.find="T",rngSplit=list(c(0,110),c(0,110),c(0,110),c(0,110)),steppin=c(0.1,0.1,0.1),lowbreaks=c(0,0,0),jitter=0,one.to.one="F",bps=c(0.2),quad.fit="T",add.lengths="T")
Petrale.Bps.2<-Oto.Age.Model.fits(Petrale.dat.otos,oto.age.col=c(5,4,3),sextype="All",Bp.find="T",rngSplit=list(c(0,110),c(0,110),c(0,110),c(0,110)),steppin=c(0.1,0.1,0.1),lowbreaks=c(0,0,0),highbreaks=c(130,130,130),jitter=1,one.to.one="T",quad.fit="T",mod.select="T",add.lengths="T")
Petrale.Bps.3<-Oto.Age.Model.fits(Petrale.dat.otos,oto.age.col=c(5,4,3),sextype="All",Bp.find="T",rngSplit=list(c(0,110),c(0,110),c(0,110),c(0,110)),steppin=c(0.1,0.1,0.1),lowbreaks=c(0,0,0),highbreaks=c(130,130,130),jitter=100,one.to.one="T",quad.fit="T",mod.select="T",add.lengths="T",comp.plot="T",comp.num=list(c(1,0),c(1,1),c(1,0)))
Petrale.Bps.F<-Oto.Age.Model.fits(Petrale.dat.otos,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=c(20,24),steppin=0.1,lowbreaks=c(0,0,0))
Petrale.Bps.M<-Oto.Age.Model.fits(Petrale.dat.otos,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=c(22,24),steppin=0.1,lowbreaks=c(0,0,0))
Petrale.Lms<-Oto.Age.Model.fits(Petrale.dat.otos,oto.age.col=c(5,4),sextype="All",Bp.find="F",rngSplit=c(0,100),steppin=10,breakpts=c(Petrale.Bps.A$Bps$All$par[1],Petrale.Bps.F$Bps$Females$par[1],Petrale.Bps.M$Bps$Males$par[1]),lowbreaks=c(0,0,0))
#Fit linear models
par(mfrow=c(2,2))
OtoAge.plot(Petrale.Lms,1)
OtoAge.plot(Petrale.Lms,2)
OtoAge.plot(Petrale.Lms,3)
#Petrale <80mg -- use this one for now
#Check summary stats of otolith weight at age
Petrale.dat.otos.80<-subset(Petrale.dat.otos,OtoWt<80)
Petrale.Bps.10_80.dat<-Oto.Age.Model.fits(Petrale.dat.otos.80,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=c(0,100),steppin=10)$Data
petrale80.summ.All<-boxplot(Petrale.Bps.10_80.dat$A[,1]~Petrale.Bps.10_80.dat$A[,2])
petrale80.summ.F<-boxplot(Petrale.Bps.10_80.dat$F[,1]~Petrale.Bps.10_80.dat$F[,2])
petrale80.summ.M<-boxplot(Petrale.Bps.10_80.dat$M[,1]~Petrale.Bps.10_80.dat$M[,2])
#Calculate breakpoints
Petrale.Bps.80.A<-Oto.Age.Model.fits(Petrale.dat.otos.80,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=c(20,26),steppin=0.1,lowbreaks=c(0,0,0))
Petrale.Bps.80.F<-Oto.Age.Model.fits(Petrale.dat.otos.80,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=c(24,30),steppin=0.1,lowbreaks=c(0,0,0))
Petrale.Bps.80.M<-Oto.Age.Model.fits(Petrale.dat.otos.80,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=c(20,26),steppin=0.1,lowbreaks=c(0,0,0))
Petrale.Lms.80<-Oto.Age.Model.fits(Petrale.dat.otos.80,oto.age.col=c(5,4),sextype="All",Bp.find="F",rngSplit=c(0,100),steppin=10,breakpts=c(Petrale.Bps.80.A$Bps$All$par[1],Petrale.Bps.80.F$Bps$Females$par[1],Petrale.Bps.80.M$Bps$Males$par[1]),lowbreaks=c(0,0,0))
#Fit linear models
par(mfrow=c(2,2))
OtoAge.plot(Petrale.Lms.80,1)
OtoAge.plot(Petrale.Lms.80,2)
OtoAge.plot(Petrale.Lms.80,3)
###############

#################
### Sablefish ###
#################
#Check summary stats of otolith weight at age
Sablefish.Bps.dat<-Oto.Age.Model.fits(Sablefish.dat.otos,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=c(0,100),steppin=10)$Data
par(mfrow=c(2,2))
sabfish.summ.All<-boxplot(Sablefish.Bps.dat$A[,1]~Sablefish.Bps.dat$A[,2])
sabfish.summ.F<-boxplot(Sablefish.Bps.dat$F[,1]~Sablefish.Bps.dat$F[,2])
sabfish.summ.M<-boxplot(Sablefish.Bps.dat$M[,1]~Sablefish.Bps.dat$M[,2])
#Calculate breakpoints
Sablefish.Bps.1<-Oto.Age.Model.fits(Sablefish.dat.otos,oto.age.col=c(5,4,3),sextype="All",Bp.find="T",rngSplit=list(c(0,80),c(0,80),c(0,80)),steppin=c(0.01,0.01,0.01),lowbreaks=c(12,12,12),jitter=0,one.to.one="F",bps=c(0.2,0.6),quad.fit="F",add.lengths="F")
Sablefish.Bps.2<-Oto.Age.Model.fits(Sablefish.dat.otos,oto.age.col=c(5,4,3),sextype="All",Bp.find="T",rngSplit=list(c(0,80),c(0,80),c(0,80)),steppin=c(0.01,0.01,0.01),lowbreaks=c(12,12,12),highbreaks=c(80,80,80),jitter=0,one.to.one="F",two.bp="T",mod.select="T",quad.fit="T",add.lengths="T")
Sablefish.Bps.3<-Oto.Age.Model.fits(Sablefish.dat.otos,oto.age.col=c(5,4,3),sextype="All",Bp.find="T",rngSplit=list(c(0,80),c(0,80),c(0,80)),steppin=c(0.01,0.01,0.01),lowbreaks=c(12,12,12),highbreaks=c(80,80,80),jitter=100,one.to.one="T",two.bp="F",mod.select="T",quad.fit="T",add.lengths="T",comp.plot="T",comp.num=list(c(1,1),c(0,1),c(0,1)))
Sablefish.Bps.A<-Oto.Age.Model.fits(Sablefish.dat.otos,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=list(c(22,24),c(22,24),c(21,23)),steppin=c(0.01,0.01,0.01),lowbreaks=c(12,12,12),jitter=100,one.to.one="T",altop="F",two.bp="T")
Sablefish.Bps.F<-Oto.Age.Model.fits(Sablefish.dat.otos,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=c(22,24),steppin=0.01,lowbreaks=c(12,12,12))
Sablefish.Bps.M<-Oto.Age.Model.fits(Sablefish.dat.otos,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=c(21,23),steppin=0.01,lowbreaks=c(12,12,12))
#Fit linear models
Sablefish.NWFSC.Lms<-Oto.Age.Model.fits(Sablefish.dat.otos,oto.age.col=c(5,4),sextype="All",Bp.find="F",rngSplit=c(15,25),steppin=0.1,breakpts=c(Sablefish.Bps.A$Bps$All$par[1],Sablefish.Bps.F$Bps$Females$par[1],Sablefish.Bps.M$Bps$Males$par[1]),lowbreaks=c(12,12,12))
par(mfrow=c(2,2))
OtoAge.plot(Sablefish.NWFSC.Lms,1)
OtoAge.plot(Sablefish.NWFSC.Lms,2)
OtoAge.plot(Sablefish.NWFSC.Lms,3)
#################

#################
### Splitnose ###
#################
#Check summary stats of otolith weight at age
Splitnose.Bps.dat<-Oto.Age.Model.fits(Splitnose.dat.otos,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=c(0,1000),steppin=50,breakpts)$Data
par(mfrow=c(2,2))
splitnose.summ.All<-boxplot(Splitnose.Bps.dat$A[,1]~Splitnose.Bps.dat$A[,2])
splitnose.summ.F<-boxplot(Splitnose.Bps.dat$F[,1]~Splitnose.Bps.dat$F[,2])
splitnose.summ.M<-boxplot(Splitnose.Bps.dat$M[,1]~Splitnose.Bps.dat$M[,2])
#Calculate breakpoints
Splitnose.Bps.1<-Oto.Age.Model.fits(Splitnose.dat.otos,oto.age.col=c(5,4,3),sextype="All",Bp.find="T",rngSplit=list(c(45,300),c(45,300),c(45,300)),steppin=c(0.01,0.01,0.01),lowbreaks=c(45,45,45),jitter=0,one.to.one="F",two.bp="T",mod.select="F",quad.fit="T",add.lengths="T")
Splitnose.Bps.2<-Oto.Age.Model.fits(Splitnose.dat.otos,oto.age.col=c(5,4,3),sextype="All",Bp.find="T",rngSplit=list(c(45,600),c(45,600),c(45,600)),steppin=c(0.01,0.01,0.01),lowbreaks=c(45,45,45),highbreaks=c(700,700,700),jitter=0,one.to.one="F",two.bp="T",mod.select="T",quad.fit="T",add.lengths="T")
Splitnose.Bps.3<-Oto.Age.Model.fits(Splitnose.dat.otos,oto.age.col=c(5,4,3),sextype="All",Bp.find="T",rngSplit=list(c(45,600),c(45,600),c(45,600)),steppin=c(0.01,0.01,0.01),lowbreaks=c(45,45,45),highbreaks=c(700,700,700),jitter=100,one.to.one="T",two.bp="F",mod.select="T",quad.fit="F",add.lengths="T",comp.plot="T")
Splitnose.Bps.A<-Oto.Age.Model.fits(Splitnose.dat.otos,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=list(c(109,111),c(127,129),c(104,106),c(15,25)),steppin=c(0.01,0.01,0.01),lowbreaks=c(45,45,45),jitter=100,one.to.one="T",altmethod = "SANN",altop="T")
Splitnose.Bps.F<-Oto.Age.Model.fits(Splitnose.dat.otos,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=c(127,129),steppin=0.01,lowbreaks=c(45,45,45))
Splitnose.Bps.M<-Oto.Age.Model.fits(Splitnose.dat.otos,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=c(104,106),steppin=0.01,lowbreaks=c(45,45,45))
Splitnose.Lms<-Oto.Age.Model.fits(Splitnose.dat.otos,oto.age.col=c(5,4),sextype="All",Bp.find="F",rngSplit=c(15,25),steppin=0.1,breakpts=c(Splitnose.Bps.A$Bps$All$par[1],Splitnose.Bps.F$Bps$Females$par[1],Splitnose.Bps.M$Bps$Males$par[1]),lowbreaks=c(45,45,45))
#Fit linear models
par(mfrow=c(2,2))
OtoAge.plot(Splitnose.Lms,1)
OtoAge.plot(Splitnose.Lms,2)
OtoAge.plot(Splitnose.Lms,3)
#################


Spp.Lms<-list(Hake.Lms,Hake.Lms.300,Petrale.Lms,Petrale.Lms.80,Sablefish.NWFSC.Lms,Splitnose.Lms)
names(Spp.Lms)<-c("Hake","Hake300","Petrale","Petrale80","Sablefish","Splitnose")
save(Spp.Lms,file="C://Users//david.lin//documents-export-2014-07-02//Spp_Lms.DMP")
#Get Numbers
Spp.Ns<-Spp.samp.size(Spp.Lms)
N.Hake300.F<-hist(Hake.Lms.300$Data$Females[,2],breaks=c(0,3,6,9,12,20))
N.Hake300.M<-hist(Hake.Lms.300$Data$Males[,2],breaks=c(0,3,6,9,12,20))
N.Petrale80.F<-hist(Petrale.Lms.80$Data$Females[,2],breaks=c(0,3,6,9,12,20))
N.Petrale80.M<-hist(Petrale.Lms.80$Data$Males[,2],breaks=c(0,3,6,9,12,20))
N.Sablefish.F<-hist(Sablefish.NWFSC.Lms$Data$Females[,2],breaks=c(0,10,20,30,50,100))
N.Sablefish.M<-hist(Sablefish.NWFSC.Lms$Data$Males[,2],breaks=c(0,10,20,30,50,100))
N.Splitnose.F<-hist(Splitnose.Lms$Data$Females[,2],breaks=c(0,5,20,40,70,100))
N.Splitnose.M<-hist(Splitnose.Lms$Data$Males[,2],breaks=c(0,5,20,40,70,105))
#Get R^2s
Spp.R2<-Spp.r2(Spp.Lms)
#Get predicted values
spp.low.ages.mat<-list(
matrix(c(50,90,1,2),nrow=2,ncol=2),matrix(c(50,1),nrow=1,ncol=2),
matrix(c(0,1),nrow=1,ncol=2),matrix(c(0,1),nrow=1,ncol=2),
matrix(c(6,12,0,1),nrow=2,ncol=2),matrix(c(6,12,0,1),nrow=2,ncol=2),
matrix(c(25,45,1,2),nrow=2,ncol=2),matrix(c(22,45,1,2),nrow=2,ncol=2))
Spp.age.pred<-Pred.ages.oto(Spp.Lms,spp.low.ages.mat)
hist(Spp.age.pred$Hake$A[,2])

#VBGF parameters
Hake.vbgf<-vbgf.pars(Spp.age.pred$Hake)
Hake300.vbgf<-vbgf.pars(Spp.age.pred$Hake300)
Petrale.vbgf<-vbgf.pars(Spp.age.pred$Petrale)
Petrale80.vbgf<-vbgf.pars(Spp.age.pred$Petrale80)
Sablefish.vbgf<-vbgf.pars(Spp.age.pred$Sablefish)
Splitnose.vbgf<-vbgf.pars(Spp.age.pred$Splitnose)

#Relative errors
Hake.RE.Linf<-Spp.RE(Hake.vbgf,1)
Hake300.RE.Linf<-Spp.RE(Hake300.vbgf,1)
Petrale.RE.Linf<-Spp.RE(Petrale.vbgf,1)
Petrale80.RE.Linf<-Spp.RE(Petrale80.vbgf,1)
Sablefish.RE.Linf<-Spp.RE(Sablefish.vbgf,1)
Splintose.RE.Linf<-Spp.RE(Splitnose.vbgf,1)
Hake.RE.k<-Spp.RE(Hake.vbgf,2)
Hake300.RE.k<-Spp.RE(Hake300.vbgf,2)
Petrale.RE.k<-Spp.RE(Petrale.vbgf,2)
Petrale80.RE.k<-Spp.RE(Petrale80.vbgf,2)
Sablefish.RE.k<-Spp.RE(Sablefish.vbgf,2)
Splintose.RE.k<-Spp.RE(Splitnose.vbgf,2)
Spp.RE.Linf<-cbind(Hake300.RE.Linf,Petrale80.RE.Linf,Sablefish.RE.Linf,Splintose.RE.Linf)
Spp.RE.k<-cbind(Hake300.RE.k,Petrale80.RE.k,Sablefish.RE.k,Splintose.RE.k)

#Plot parameter estimate comps

#VBGF fits

#1:1 plots.

#Ageing error
Hake300.F.4AE<-Punt_age_prep(Spp.age.pred$Hake300$F[,c(2,3)])
Hake300.M.4AE<-Punt_age_prep(Spp.age.pred$Hake300$M[,c(2,3)])
Petrale80.F.4AE<-Punt_age_prep(Spp.age.pred$Petrale80$F[,c(2,3)])
Petrale80.M.4AE<-Punt_age_prep(Spp.age.pred$Petrale80$M[,c(2,3)])
Sablefish.F.4AE<-Punt_age_prep(Spp.age.pred$Sablefish$F[,c(2,3)])
Sablefish.M.4AE<-Punt_age_prep(Spp.age.pred$Sablefish$M[,c(2,3)])
Splitnose.F.4AE<-Punt_age_prep(Spp.age.pred$Splitnose$F[,c(2,3)])
Splitnose.M.4AE<-Punt_age_prep(Spp.age.pred$Splitnose$M[,c(2,3)])
Spp.AgeError<-list(Hake300.F.4AE,Hake300.M.4AE,Petrale80.F.4AE,Petrale80.M.4AE,Sablefish.F.4AE,Sablefish.M.4AE,Splitnose.F.4AE,Splitnose.M.4AE)
names(Spp.AgeError)<-c("Hale_F","Hake_M","Petrale_F","Petrale_M","Sabelfish_F","Sablefish_M","Splitnose_F","Splitnose_M")
save(Spp.AgeError,file="C://Users//david.lin//documents-export-2014-07-02//Spp_AgeError.DMP")


#ASssessment comparisons
Dir<-"C:/Users/copeja/Desktop/Current projects/otolith weight/Assessment sensitivity/Sablefish/"
sablefish_annuli<-SS_output(paste(Dir,"2011 sablefish model files_annuli/",sep=""))
sablefish_otolith<-SS_output(paste(Dir,"2011 sablefish model files_otowt/",sep=""))
sablefish_getoutput<-SSgetoutput(keyvec = NULL, dirvec = c(paste(Dir,"2011 sablefish model files_annuli/",sep=""),paste(Dir,"2011 sablefish model files_otowt/",sep="")))
sablefish_summaryoutput <- SSsummarize(sablefish_getoutput)
SSplotComparisons(sablefish_summaryoutput)
(sablefish_annuli$timeseries$SpawnBio-sablefish_otolith$timeseries$SpawnBio)/sablefish_annuli$timeseries$SpawnBio


##Splitnose sole
#annuli, female and male
cbind(c(1.66,seq(2,100,1)),VBGF(30.48,0.1,-4.04,c(1.66,seq(2,100,1))))
cbind(c(1.66,seq(2,100,1)),VBGF(29.77,0.06,-9.96,c(1.66,seq(2,100,1))))
#otowt, female and male
cbind(c(1.66,seq(2,100,1)),VBGF(32.81,0.05,-9.96,c(1.66,seq(2,100,1))))
cbind(c(1.66,seq(2,100,1)),VBGF(34.25,0.03,-18.63,c(1.66,seq(2,100,1))))
Dir<-"C:/Users/copeja/Desktop/Current projects/otolith weight/Assessment sensitivity/Splintose/"
splitnose_annuli<-SS_output(paste(Dir,"Splitnose_XSSS_annuli",sep=""))
splitnose_otolith<-SS_output(paste(Dir,"Splitnose_XSSS_otowt",sep=""))
splitnose_getoutput<-SSgetoutput(keyvec = NULL, dirvec = c(paste(Dir,"Splitnose_XSSS_annuli/",sep=""),paste(Dir,"Splitnose_XSSS_otowt/",sep="")))
splitnose_summaryoutput <- SSsummarize(splitnose_getoutput)
SSplotComparisons(splitnose_summaryoutput,subplots=1:4)

##Petrale sole
#annuli, female and male
cbind(seq(0,17,1),VBGF(50.19,0.25,-0.27,seq(0,17,1)))
cbind(seq(0,17,1),VBGF(39.74,0.47,0.33,seq(0,17,1)))
#otowt, female and male
cbind(seq(0,17,1),VBGF(52.02,0.22,-0.33,seq(0,17,1)))
cbind(seq(0,17,1),VBGF(40.33,0.41,0.21,seq(0,17,1)))

Dir<-"C:/Users/copeja/Desktop/Current projects/otolith weight/Assessment sensitivity/Petrale/"
Petrale_annuli<-SS_output(paste(Dir,"Petrale2013_annuli",sep=""))
Petrale_otolith<-SS_output(paste(Dir,"Petrale2013_otowt",sep=""))
petrale_getoutput<-SSgetoutput(keyvec = NULL, dirvec = c(paste(Dir,"Petrale2013_annuli/",sep=""),paste(Dir,"Petrale2013_otowt/",sep="")))
petrale_summaryoutput <- SSsummarize(petrale_getoutput)
SSplotComparisons(petrale_summaryoutput)

##Pacific hake
#annuli, female and male
cbind(seq(0,20,1),VBGF(50.1,0.37,-0.67,seq(0,20,1)))
cbind(seq(0,20,1),VBGF(48.65,0.41,-0.58,seq(0,20,1)))
#otowt, female and male
cbind(seq(0,20,1),VBGF(52.66,0.24,-1.76,seq(0,20,1)))
cbind(seq(0,20,1),VBGF(49.58,0.33,0.9,seq(0,20,1)))

Dir<-"C:/Users/copeja/Desktop/Current projects/otolith weight/Assessment sensitivity/P_hake/"
Hake_annuli<-SS_output(paste(Dir,"Phake_annuli",sep=""))
Hake_otolith<-SS_output(paste(Dir,"Phake_otowt",sep=""))
Hake_getoutput<-SSgetoutput(keyvec = NULL, dirvec = c(paste(Dir,"Phake_annuli/",sep=""),paste(Dir,"Phake_otowt/",sep="")))
Hake_summaryoutput <- SSsummarize(Hake_getoutput)
SSplotComparisons(Hake_summaryoutput)


