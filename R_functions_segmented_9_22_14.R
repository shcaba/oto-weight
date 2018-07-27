#############################################
#############  FUNCTIONS  ###################
#############################################
library(gplots)
library(reshape2)
library(segmented)
#Fit oto weights and ages
#oto.age.col = columns that contain the otolith weight and ages, respectively
#sextype = "All" means both M and F use any unsexed individuals in the data set.
#BP.find = T- find breakpoints. BP.find = F- fit linear models.
#steppin = intervals to search over
#breakpts = breakpoints used when Bp.find = "F"
#lowbreaks = the otolith weight at which the data set is separated. The low numbers are not used in any of the linear models.
#intcpt = "Y" estimates the intercept; "N" sets it to 0.


Oto.Age.Model.fits<-function(spp.dat.in,oto.age.col=c(5,4,3),sextype="All",Bp.find="F",rngSplit=list(c(10,200),c(10,200),c(10,200)),lowbreaks=c(0,0,0),highbreaks=c(10000,10000,10000),steppin=c(10,10,10),breakpts=c(0,0,0),intcpt="Y",optmethod = "Nelder-Mead",pref.method = c(1,1,1),jitter=0,one.to.one="F",bps = 0,quad.fit="F",add.lengths="F",comp.plot="F",comp.num=list(c(1,1),c(1,1),c(1,1)),skip.check=F)
{
  #Data prep
  CVs<-dcast(spp.dat.in,Sex~Age,fun.aggregate=sd,value.var="OtoWt")[,-1]/dcast(spp.dat.in,Sex~Age,fun.aggregate=mean,value.var="OtoWt")[,-1]
  rownames(CVs)<-c("F","M")
  dat.names<-c("All","Females","Males")
  Spp.dat.AFM<-Spp.dat.Lt.AFM<-list()
    
  Spp.dat.AFM[[1]]<-na.omit(cbind(spp.dat.in[,oto.age.col[1]],spp.dat.in[,oto.age.col[2]],spp.dat.in[,oto.age.col[3]]))
  if(length(oto.age.col)==3){Spp.dat.Lt.AFM[[1]]<-na.omit(spp.dat.in[,oto.age.col])}
 
  if(sextype=="All")
  {
    Spp.dat.AFM[[2]]<-na.omit(cbind(subset(spp.dat.in,Sex!="M")[,oto.age.col[1]],subset(spp.dat.in,Sex!="M")[,oto.age.col[2]],subset(spp.dat.in,Sex!="M")[,oto.age.col[3]]))
    Spp.dat.AFM[[3]]<-na.omit(cbind(subset(spp.dat.in,Sex!="F")[,oto.age.col[1]],subset(spp.dat.in,Sex!="F")[,oto.age.col[2]],subset(spp.dat.in,Sex!="F")[,oto.age.col[3]]))
    names(Spp.dat.AFM)<-dat.names
    #Lts<-list()
    #if(length(spp.dat.in[,3])!=length(Spp.dat.AFM[[1]][,1])){Lts[[1]]<-spp.dat.in[as.numeric(-attributes(Spp.dat.AFM[[1]])$na.action),3]}
    #if(length(spp.dat.in[,3])==length(Spp.dat.AFM[[1]][,1])){Lts[[1]]<-spp.dat.in[,3]}
    #if(length(subset(spp.dat.in,Sex!="M")[,3])!=length(Spp.dat.AFM[[2]][,1])){Lts[[2]]<-subset(spp.dat.in,Sex!="M")[as.numeric(-attributes(Spp.dat.AFM[[2]])$na.action),3]}
    #if(length(subset(spp.dat.in,Sex!="M")[,3])==length(Spp.dat.AFM[[2]][,1])){Lts[[2]]<-subset(spp.dat.in,Sex!="M")[,3]}
    #if(length(subset(spp.dat.in,Sex!="F")[,3])!=length(Spp.dat.AFM[[3]][,1])){Lts[[3]]<-subset(spp.dat.in,Sex!="F")[as.numeric(-attributes(Spp.dat.AFM[[3]])$na.action),3]}
    #if(length(subset(spp.dat.in,Sex!="F")[,3])==length(Spp.dat.AFM[[3]][,1])){Lts[[3]]<-subset(spp.dat.in,Sex!="F")[,3]}
  }
  if(sextype!="All")
  {
    Spp.dat.AFM[[2]]<-na.omit(cbind(subset(spp.dat.in,Sex=="F")[,oto.age.col[1]],subset(spp.dat.in,Sex=="F")[,oto.age.col[2]]))
    Spp.dat.AFM[[3]]<-na.omit(cbind(subset(spp.dat.in,Sex=="M")[,oto.age.col[1]],subset(spp.dat.in,Sex=="M")[,oto.age.col[2]]))
    names(Spp.dat.AFM)<-dat.names
    #Lts<-list()
    #if(length(spp.dat.in[,3])!=length(Spp.dat.AFM[[1]][,1])){Lts[[1]]<-spp.dat.in[as.numeric(-attributes(Spp.dat.AFM[[1]])$na.action),3]}
    #if(length(spp.dat.in[,3])==length(Spp.dat.AFM[[1]][,1])){Lts[[1]]<-spp.dat.in[,3]}
    #if(length(subset(spp.dat.in,Sex=="F")[,3])!=length(Spp.dat.AFM[[2]][,1])){Lts[[2]]<-subset(spp.dat.in,Sex=="F")[as.numeric(-attributes(Spp.dat.AFM[[2]])$na.action),3]}
    #if(length(subset(spp.dat.in,Sex=="F")[,3])==length(Spp.dat.AFM[[2]][,1])){Lts[[2]]<-subset(spp.dat.in,Sex=="F")[,3]}
    #if(length(subset(spp.dat.in,Sex=="M")[,3])!=length(Spp.dat.AFM[[3]][,1])){Lts[[3]]<-subset(spp.dat.in,Sex=="M")[as.numeric(-attributes(Spp.dat.AFM[[3]])$na.action),3]}
    #if(length(subset(spp.dat.in,Sex=="M")[,3])==length(Spp.dat.AFM[[3]][,1])){Lts[[3]]<-subset(spp.dat.in,Sex=="M")[,3]}
  }
  #names(Lts)<-c("Lt_A","Lt_F","Lt_M")
  colnames(Spp.dat.AFM[[1]])<-colnames(Spp.dat.AFM[[2]])<-colnames(Spp.dat.AFM[[3]])<-colnames(spp.dat.in[,oto.age.col])
  Dat.Bp.Lm.out<-list(Spp.dat.AFM,lowbreaks)
  names(Dat.Bp.Lm.out)<-c("Data.in","Low_wt_breaks_used")
  
  #Find breakpoint
  if(Bp.find=="T")
  {
    Spp.bps<-Prof.out<-Low.breaks<-out.f.all<-pred.ages.out<-BP2.models<-list()
    #Create 2x2 empty panels for graphs
    
    for(i in 1:length(Spp.dat.AFM))
    {
      
      Spp.dat.AFM.mod<-(Spp.dat.AFM[[i]][(Spp.dat.AFM[[i]][,1]>lowbreaks[i])&(Spp.dat.AFM[[i]][,1]<highbreaks[i]),])
      seg.data<-as.data.frame(Spp.dat.AFM.mod)
      colnames(seg.data)<-c("OtoWt","Age")
      x.lm<-lm(Age~OtoWt,data=seg.data)
      if (pref.method[i] == 1){
        out.f<-segmented(x.lm, seg.Z=~OtoWt,psi=(sum(rngSplit[[i]])/2))
      } else if (pref.method[i] == 0){
        out.f<-x.lm
      } else if (pref.method[i] == 2){
        out.f<-segmented(x.lm, seg.Z=~OtoWt,psi=c(quantile(Spp.dat.AFM.mod[,1],0.25,na.rm=T),quantile(Spp.dat.AFM.mod[,1],0.75,na.rm=T)),jt=FALSE)
      } else {
        out.f<-NA
      }
      
      out.f.all[[i]]<-out.f
      Sort.ind<-sort(Spp.dat.AFM.mod[,1],index.return=T)
      sorted.dat<-Spp.dat.AFM.mod[Sort.ind$ix,]
      otowt.minus.NA <- as.matrix(sorted.dat[,1])
      age.minus.NA <- as.matrix(sorted.dat[,2])
      otowt.pred.ages<-predict(out.f,newdata=data.frame(OtoWt=otowt.minus.NA))
      orig.pred.val<-predict(out.f)
      pred.ages.out[[i]]<-cbind(orig.pred.val,otowt.pred.ages,age.minus.NA)
      colnames(pred.ages.out[[i]])<-c("original","org.pred","annuli_ages")
            
      par(mfrow=c(2,2))
      
      if (add.lengths=="T"){
        library(scatterplot3d)
        library(rgl)
        library(segmented)
        
        Ages.Spp<-as.matrix(Spp.dat.AFM.mod[,2])
        Weights.Spp<-as.matrix(Spp.dat.AFM.mod[,1])
        Lengths.Spp<-as.matrix(Spp.dat.AFM.mod[,3])
        out.zero<-lm(Ages.Spp~Weights.Spp+Lengths.Spp)
        #plot3d(x=Weights.Spp,y=Lengths.Spp,z=Ages.Spp,main=dat.names[i])
        #planes3d(out.zero$coefficients[[2]],out.zero$coefficients[[3]],-1,out.zero$coefficients[[1]],col="blue")
        #print(AIC(out.zero))
        
        out.one.prelim<-lm(Ages.Spp~Weights.Spp+Lengths.Spp)
        oo<-try(segmented(out.one.prelim,seg.Z=~Weights.Spp+Lengths.Spp,psi=list(Weights.Spp=median(Spp.dat.AFM.mod[,1]),Lengths.Spp=median(Spp.dat.AFM.mod[,3]))))
        
        if(class(oo)[1]=="try-error")
        {     
          iii<-0
          while(class(oo)[1]=="try-error"|iii<=10)
          {
            oo<-try(segmented(x.lm, seg.Z=~OtoWt,psi=bps.init,jt=FALSE),silent=TRUE)
            iii<-iii+1
          }
        }
        
        if(class(oo)[1]!="try-error"){
          out.one<-oo
          #plot3d(x=Weights.Spp,y=Lengths.Spp,z=Ages.Spp,main=dat.names[i])
          #planes3d(out.one$coefficients[[2]],out.one$coefficients[[3]],-1,out.one$coefficients[[1]],col="blue")
          #planes3d(out.one$coefficients[[4]],out.one$coefficients[[5]],1,out.one$coefficients[[1]]-out.one$coefficients[[4]]*out.one$psi[1,2]-out.one$coefficients[[5]]*out.one$psi[2,2],col="red")
          #print(AIC(out.one))
        } 
        else{out.one <- NA}
      }
      
      #Start jitter
      if(jitter>0){
        out.f.jitter=out.f.jitter.all=list()
        output.total= matrix(NA,nrow=jitter,ncol=length(out.f$coefficients)+length(out.f$psi[,2])+1 )
        for(x in 1:jitter){
          #out.f.jitter[[x]]<-segmented(x.lm, seg.Z=~OtoWt,psi=quantile(Spp.dat.AFM.mod[,1],0.5,na.rm=T),jt=TRUE)
          out.f.jitter[[x]]<-segmented(x.lm, seg.Z=~OtoWt,psi=(sum(rngSplit[[i]])/2),jt=TRUE)
         for(coeff.x in 1:length(out.f$coefficients))
         {
          output.total[x,coeff.x] = out.f.jitter[[x]]$coefficients[[coeff.x]] #y-intercept
         }

          for(psi.x in 1:length(out.f$psi[,2]))
          {
            output.total[x,length(out.f$coefficients)+psi.x] = out.f.jitter[[x]]$psi[[psi.x,2]] #y-intercept
          } 

          output.total[x,ncol(output.total)] = AIC(out.f.jitter[[x]]) #AIC
        }
         names(output.total)=c("Intercept","First Slope","Second Slope","Break_Pt","AIC" )
         x.spec = c(1:jitter)
         par(mfrow=c(3,2),oma=c(0,0,1,0))
         if (min(output.total[,1])<0.9*out.f$coefficients[[1]] | max(output.total[,1])>1.1*out.f$coefficients[[1]]){
           plot(x.spec,output.total[,1], xlab = "jitter number", ylab = "intercept")
         } else {
           plot(x.spec,output.total[,1], xlab = "jitter number", ylab = "intercept",ylim = c(0.9*out.f$coefficients[[1]],1.1*out.f$coefficients[[1]]))
         }
        abline(h=min(output.total[,1]), col = "red", lty = 3)
        abline(h=out.f$coefficients[[1]],col = "green")
        legend("topright", legend = c("base case","jitter minimum") , lty=c(1,2), col=c("green","red"), bty='n', cex=.75)
        if (min(output.total[,2])<0.9*out.f$coefficients[[2]] | max(output.total[,2])>1.1*out.f$coefficients[[2]]){
          plot(x.spec,output.total[,2], xlab = "jitter number", ylab = "first slope")
        } else {
          plot(x.spec,output.total[,2], xlab = "jitter number", ylab = "first slope",ylim = c(0.9*out.f$coefficients[[2]],1.1*out.f$coefficients[[2]]))
        }
        abline(h=min(output.total[,2]), col = "red", lty = 3)
        abline(h=out.f$coefficients[[2]],col = "green")
        legend("topright", legend = c("base case","jitter minimum") , lty=c(1,2), col=c("green","red"), bty='n', cex=.75)
        if (min(output.total[,3])<0.9*out.f$coefficients[[4]] | max(output.total[,3])>1.1*out.f$coefficients[[3]]){
          plot(x.spec,output.total[,3], xlab = "jitter number", ylab = "second slope")
        } else {
          plot(x.spec,output.total[,3], xlab = "jitter number", ylab = "second slope",ylim = c(0.9*out.f$coefficients[[3]],1.1*out.f$coefficients[[3]]))
        }
        abline(h=min(output.total[,3]), col = "red", lty = 3)
        abline(h=out.f$coefficients[[3]],col = "green")
        legend("topright", legend = c("base case","jitter minimum") , lty=c(1,2), col=c("green","red"), bty='n', cex=.75)
        if (min(output.total[,4])<0.9*out.f$psi[1,2] | max(output.total[,4])>1.1*out.f$psi[1,2]){
          plot(x.spec,output.total[,4], xlab = "jitter number", ylab = "breakpoint")
        } else {
          plot(x.spec,output.total[,4], xlab = "jitter number", ylab = "breakpoint",ylim = c(0.9*out.f$psi[1,2],1.1*out.f$psi[1,2]))
        }
        abline(h=min(output.total[,4]), col = "red", lty = 3)
        abline(h=out.f$psi[1,2],col = "green")
        legend("topright", legend = c("base case","jitter minimum") , lty=c(1,2), col=c("green","red"), bty='n', cex=.75)
        if (abs(min(output.total[,5])-AIC(out.f))<5){
          plot(x.spec,output.total[,5], ylim = c(AIC(out.f)-5,AIC(out.f)+5),xlab = "jitter number", ylab = "AIC value")
        } else {
          plot(x.spec,output.total[,5],xlab = "jitter number", ylab = "AIC value",ylim=c(min(output.total[,5]),max(output.total[,5])))
        }
        abline(h=AIC(out.f), col = "green")
        abline(h=min(output.total[,5]),col = "red", lty = 3)
        legend("topright", legend = c("base case","jitter minimum") , lty=c(1,2), col=c("green","red"), bty='n', cex=.75)
        points(AIC(out.f), col = "blue")
        y.min.val <- min(output.total[x,5])
        abline(h=y.min.val, col = "green")
        #Do one to one comparison
        min.val.index =(output.total[,5]-AIC(out.f))==min(output.total[,5]-AIC(out.f))  #Gets index of minimum obj fun from jitter
        #min.parm.vals=output.total[min.val.index,] #Get parameter values from min jitter run
        min.jitt.num<-c(1:jitter)[min.val.index]
        out.f.jitter.all[[i]]<-out.f.jitter[[min.jitt.num[1]]]
        if (one.to.one == "T")
        {
          
          ages.pred.min.val<-predict(out.f.jitter[[min.jitt.num[1]]])
          plot(orig.pred.val,ages.pred.min.val,xlim = c(0,max(orig.pred.val)+2),ylim = c(0,max(ages.pred.min.val)+2),xlab = "benchmark age values",ylab = "Alternative Parameter Age values")
          abline(a=0,b=1,col="green")
          
         # pred.ages.out[[i]]<-cbind(pred.ages.out,ages.pred.min.val)
        #  colnames(pred.ages.out[[i]])<-c("original","org.pred","annuli_ages", "jitter min")
          #print(pred.ages.out)
          #print(sorted.dat)
        }
        title( main = dat.names[i],outer = TRUE)  
      }
      #End jitter  
      #Spp.bps[[i]]<-out.f
      #Low.breaks[[i]]<-lowbreaks[i]
      #names(Spp.bps)[[i]]<-names(Low.breaks)[[i]]<-dat.names[i]
      #likeProf <- PRegLikeProf.fn(out.f$par,x.in=Spp.dat.AFM.mod[,1],y.in=Spp.dat.AFM.mod[,2],rngSplit=rngSplit[[i]],steppin=steppin[i])
      #likeProf$fDiff <- likeProf$f-out.f$value
      #likeProf$CI <- likeProf$fDiff-2.5
      #Prof.out[[i]]<-tmp <- likeProf
      #tmp[tmp$f==min(tmp$f),]
      #plot graphs using the modified data with piecewise curves.
      #plot(tmp$split,tmp$fDiff,type="l",xlab="Breakpoint",ylab="-LogLike Difference from minimum",lwd=2)
      #abline(h=qchisq(0.95,1)/2,col="red",lwd=1.5,lty=2)
      #title( main = dat.names[i],outer = TRUE)       
    if(quad.fit == "T"){
      x <- Spp.dat.AFM.mod[,1]
      y <- Spp.dat.AFM.mod[,2]
      par(mfrow=c(1,1))
      quad <- lm(y ~ poly(x, 2, raw= TRUE))
      #print((AIC(quad)))
      domain <- seq(rngSplit[[i]][1],rngSplit[[i]][2],0.1)
      #print(quad)
      fit = quad$coefficients[[1]]+quad$coefficients[[2]]*domain+quad$coefficients[[3]]*domain^2
      #print(c(length(domain),length(fit)))
    }
    
    lin.reg <- lm(Spp.dat.AFM.mod[,2]~Spp.dat.AFM.mod[,1])
    aic.1 <- AIC(lin.reg)
    rsq.1 <- summary(lin.reg)$r.squared
    par(mfrow=c(1,1))
    plot(Spp.dat.AFM.mod[,1],Spp.dat.AFM.mod[,2],xlab = "Otolith Weight", ylab = "Ages", main = dat.names[i],pch=21,bg="grey",cex=1.25)
    plot.segmented(out.f,add=TRUE,col="red",lwd=3)
    abline(lin.reg, col="green",lwd=3)
    aic.2 <- AIC(out.f)
    rsq.2 <- summary(out.f)$r.squared
    aic.vec <- c(aic.1,aic.2)
    rsq.vec <- c(rsq.1,rsq.2)
    if (length(bps)==1 )
      {
      print(out.f$psi)
      smartlegend("left","top",paste(c("0 Bp model: AIC=","1 Bp model: AIC=",paste(length(bps))),round(aic.vec,2),c("/ R^2=","/ R^2="),round(rsq.vec,4),sep=" "),lty=1,lwd=2,col=c("green","red"),bty="n")
      }
    if (length(bps)>1 )
    {
      bps.init<-bps
      for(ii in 1:length(bps)){bps.init[ii]<-quantile(Spp.dat.AFM.mod[,1],bps[ii],na.rm=T)}
      tt<-try(try(segmented(x.lm, seg.Z=~OtoWt,psi=bps.init,jt=FALSE),silent=TRUE))
     if(skip.check==T)
     {
       if(class(tt)=="try-error")
       {
         iii<-0
         while(class(tt)=="try-error"|iii<=10)
         {
           tt<-try(segmented(x.lm, seg.Z=~OtoWt,psi=bps.init,jt=FALSE),silent=TRUE)
           iii<-iii+1
         }
       }  
     }
        if(class(tt)!="try-error")
       {
          BP2.models[[i]]<-bp.2<-tt
          plot(bp.2,add=TRUE,col="blue",lwd=3)
          aic.3 <- AIC(bp.2)
          rsq.3 <- summary(bp.2)$r.squared
          print(bp.2$psi)
          print(aic.3)
       } 
       else
         {
           aic.3 <- NA
           rsq.3 <- NA
           print("NA")
         }
      aic.vec <- c(aic.1,aic.2,aic.3)
      rsq.vec <- c(rsq.1,rsq.2,rsq.3)
      if (quad.fit == "F" & add.lengths=="F"){
      smartlegend("left","top",paste(c("0 Bp model: AIC=","1 Bp model: AIC=",paste(length(bps)," Bp model: AIC=",sep="")),round(aic.vec,2),c("/ R^2=","/ R^2=","/ R^2="),round(rsq.vec,4),sep=" "),lty=1,lwd=2,col=c("green","red","blue"),bty="n")
      
      }
      if(quad.fit == "T" & add.lengths=="F"){
        curve(quad$coefficients[[1]]+quad$coefficients[[2]]*x+quad$coefficients[[3]]*x^2,col="purple",lwd=3,add=TRUE)
        aic.4 <- AIC(quad)
        rsq.4 <- summary(quad)$r.squared
        aic.vec <- c(aic.1,aic.2,aic.3,aic.4)
        rsq.vec <- c(rsq.1,rsq.2,rsq.3,rsq.4)
        smartlegend("left","top",paste(c("0 Bp model: AIC=","1 Bp model: AIC=",paste(length(bps)," Bp model: AIC=",sep=""),"Quadratic model: AIC="),round(aic.vec,2),c("/ R^2=","/ R^2=","/ R^2:=","/ R^2="),round(rsq.vec,4),sep=" "),lty=1,lwd=2,col=c("green","red","blue","purple"),bty="n")
      }
      
      if(add.lengths == "T" & quad.fit == "T"){
        pred.zero<-predict(out.zero)
        points(Weights.Spp,pred.zero,col="orange",pch=21,bg="orange",cex=1.25)
        curve(quad$coefficients[[1]]+quad$coefficients[[2]]*x+quad$coefficients[[3]]*x^2,col="purple",lwd=3,add=TRUE)
        aic.4 <- AIC(quad)
        rsq.4 <- summary(quad)$r.squared
        aic.5 <- AIC(out.zero)
        rsq.5 <- summary(out.zero)$r.squared
        if(class(oo)[1]!="try-error"){
          pred.one<-predict(out.one)
          points(Weights.Spp,pred.one,col="yellow",pch=21,bg="yellow",cex=1.25)
          aic.6 <- AIC(out.one)
          rsq.6 <- summary(out.one)$r.squared
        } else{
          aic.6 <- NA
          rsq.6 <- NA
        }
        aic.vec<- c(aic.1,aic.2,aic.3,aic.4,aic.5,aic.6)
        rsq.vec<- c(rsq.1,rsq.2,rsq.3,rsq.4,rsq.5,rsq.6)
        smartlegend("left","top",paste(c("0 Bp model: AIC=","1 Bp model: AIC=",paste(length(bps)," Bp model: AIC=",sep=""),"Quadratic model: AIC=","Zero Breakline model: AIC=","One Breakline model: AIC="),round(aic.vec,2),c("/ R^2=","/ R^2=","/ R^2=","/ R^2=","/ R^2=","/ R^2="),round(rsq.vec,4),sep=" "),lty=1,lwd=2,col=c("green","red","blue","purple","orange","yellow"),bty="n")
      }
    }
      if (comp.plot=="T"){
        meth.1.names<-c("Zero Breakpoint Method","One Breakpoint Method","Two Breakpoint Method","Quadratic Method")
        meth.2.names<-c("Zero Breakline Method", "One Breakline Method")
        if(comp.num[[i]][1]==0){
          meth.1<-predict(lin.reg)
        } else if (comp.num[[i]][1]==1){
          meth.1<-predict(out.f)
        } else if (comp.num[[i]][1]==2){
          meth.1<-predict(bp.2)
        } else {
          meth.1<-predict(quad)
        }
        if(comp.num[[i]][2]==0){
          meth.2<-predict(out.zero)
        } else {
          meth.2<-predict(out.one)
        }
        plot(meth.1,meth.2,col="blue",main=dat.names[i],xlab=meth.1.names[comp.num[[i]][1]+1],ylab=meth.2.names[comp.num[[i]][2]+1])
        abline(a=0,b=1,col="green")
      }
    
    }
    #Create list of modified data, breakpoints, and minimization results. Standardize names.
    #names(Prof.out)<-dat.names
    #if(jitter>0){
      #Dat.Bps.out<-list(Spp.dat.AFM,Spp.bps,Low.breaks,Prof.out,Spp.dat.AFM.mod[,1],Spp.dat.AFM.mod[,2],output.total,pred.ages.out)
      #names(Dat.Bps.out)<-c("Data","Bps","AgeBreaks","Profile","OtoWts","Ages","Jitter Info","Predicted Ages")
    #} else {
     # Dat.Bps.out<-list(Spp.dat.AFM,Spp.bps,Low.breaks,Prof.out,Spp.dat.AFM.mod[,1],Spp.dat.AFM.mod[,2])
      #names(Dat.Bps.out)<-c("Data","Bps","AgeBreaks","Profile","OtoWts","Ages")     
    #}
    #return(Dat.Bps.out)
   names(out.f.all)<-names(pred.ages.out)<-c("All","Females","Males")
   if(length(bps)>1){names(BP2.models)<-c("All","Females","Males")}
   if(jitter>0)
     {
     names(out.f.jitter.all)<-c("All","Females","Males")
     Fxn.out<-c(Dat.Bp.Lm.out,list(out.f.all,out.f.jitter.all,BP2.models,pred.ages.out))
     names(Fxn.out)<-c(names(Dat.Bp.Lm.out),"lms.out","Min.jitter.lms","2Bpts models","Data_for_comps")      
     }

   if(jitter==0)
   {
     Fxn.out<-c(Dat.Bp.Lm.out,list(out.f.all,BP2.models,pred.ages.out))
     names(Fxn.out)<-c(names(Dat.Bp.Lm.out),"lms.out","2Bpts models","Data_for_comps")
   }
  }
  
  if(Bp.find=="F")
  {
    out.f.jitter<-NA
    Lm.out<-list()
    Lm.names<-c("Bp1","Bp2","NoBp")
    for(i in 1:length(Spp.dat.AFM))
    {
      Spp.dat.AFM.mod<-Spp.dat.AFM[[i]][Spp.dat.AFM[[i]][,1]>lowbreaks[i],]
      #Piecewise models
      dat.bp1<-Spp.dat.AFM.mod[Spp.dat.AFM.mod[,1]<=breakpts[i],]
      dat.bp2<-Spp.dat.AFM.mod[Spp.dat.AFM.mod[,1]>breakpts[i],]
      if(intcpt=="Y")
      {
        if(length(dat.bp1)>2) {pw1.lm<-lm(dat.bp1[,2]~dat.bp1[,1])}
        if(length(dat.bp1)<=2) {pw1.lm<-NA}
        if(length(dat.bp2)>2) {pw2.lm<-lm(dat.bp2[,2]~dat.bp2[,1])}
        if(length(dat.bp2)<=2) {pw2.lm<-NA}
        #Linear models
        lm.all<-lm(Spp.dat.AFM.mod[,2]~Spp.dat.AFM.mod[,1])
        Lm.out[[i]]<-list(pw1.lm,pw2.lm,lm.all)
        names(Lm.out[[i]])<-Lm.names
      }
      if(intcpt=="N")
      {
        if(length(dat.bp1)>2) {pw1.lm<-lm(dat.bp1[,2]~dat.bp1[,1]-1)}
        if(length(dat.bp1)<=2) {pw1.lm<-NA}
        if(length(dat.bp2)>2) {pw2.lm<-lm(dat.bp2[,2]~dat.bp2[,1])}
        if(length(dat.bp2)<=2) {pw2.lm<-NA}
        #Linear models
        lm.all<-lm(Spp.dat.AFM.mod[,2]~Spp.dat.AFM.mod[,1]-1)
        Lm.out[[i]]<-list(pw1.lm,pw2.lm,lm.all)
        names(Lm.out[[i]])<-Lm.names
        names(Lm.out)<-dat.names
      }
    }
    
    Fxn.out<-c(Dat.Bp.Lm.out,list(breakpts,Lm.out,pred.ages.out))
    names(Fxn.out)<-c(names(Dat.Bp.Lm.out),"Breakpoints","lms.out","Data_for_comps")
  }
  return(Fxn.out)
}
#Fit PW and LM models on same plot
OtoAge.plot<-function(model.in,data.choice,col.bpts=c("orange","blue"))
{
  age.break.low.dat<-model.in$Data[[data.choice]][model.in$Data[[data.choice]][,1]<=model.in$Low_wt_breaks_used[data.choice],]
  age.break.hi.dat<-model.in$Data[[data.choice]][model.in$Data[[data.choice]][,1]>model.in$Low_wt_breaks_used[data.choice],]
  dat.bp1<-age.break.hi.dat[age.break.hi.dat[,1]<=model.in$Bps[data.choice],]
  dat.bp2<-age.break.hi.dat[age.break.hi.dat[,1]>model.in$Bps[data.choice],]
  plot(age.break.hi.dat,pch=21,col="black",bg="gray",xlim=c(0,ceiling(max(model.in$Data[[data.choice]][,1]))),ylim=c(0,ceiling(max(model.in$Data[[data.choice]][,2]))),xlab="Otolith weight",ylab="Age")
  points(age.break.low.dat,pch=21,col="black",bg="lightblue")
  abline(v=model.in$Bps[data.choice],col="red",lty=2)
  abline(v=model.in$Low_wt_breaks_used[data.choice],col="red",lty=3)
  if(dim(dat.bp1)[1]!=0){lines(dat.bp1[,1],model.in$LMs[[data.choice]]$Bp1$fitted,col=col.bpts[1], lwd=2)}
  if(dim(dat.bp2)[1]!=0){lines(dat.bp2[,1],model.in$LMs[[data.choice]]$Bp2$fitted,col=col.bpts[2],lwd=2)}
  if(dim(age.break.hi.dat)[1]!=0){lines(age.break.hi.dat[,1],model.in$LMs[[data.choice]]$NoBp$fitted,col="black",lwd=2)}
}

LM.comp.plots<-function(model.in,data.choice,add.plot=F,col.bpts=c("orange","blue"))
{
  age.break.low.dat<-model.in$Data[[data.choice]][model.in$Data[[data.choice]][,1]<=model.in$Low_wt_breaks_used[data.choice],]
  age.break.hi.dat<-model.in$Data[[data.choice]][model.in$Data[[data.choice]][,1]>model.in$Low_wt_breaks_used[data.choice],]
  dat.bp1<-age.break.hi.dat[age.break.hi.dat[,1]<=model.in$Bps[data.choice],]
  dat.bp2<-age.break.hi.dat[age.break.hi.dat[,1]>model.in$Bps[data.choice],]
  if(add.plot==F)
  {
    if(dim(dat.bp1)[1]!=0){plot(dat.bp1[,1],model.in$LMs[[data.choice]]$Bp1$fitted,col=col.bpts[1], lwd=2,type="l",xlim=c(0,ceiling(max(model.in$Data[[data.choice]][,1]))),ylim=c(0,ceiling(max(model.in$Data[[data.choice]][,2]))),xlab="Otolith weight",ylab="Age")}
    if(dim(dat.bp2)[1]!=0){lines(dat.bp2[,1],model.in$LMs[[data.choice]]$Bp2$fitted,col=col.bpts[2],lwd=2)}
    #if(dim(age.break.hi.dat)[1]!=0){plot(age.break.hi.dat[,1],model.in$LMs[[data.choice]]$NoBp$fitted,col="black",lwd=2,type="l")}
  }
  if(add.plot==T)
  {
    if(dim(dat.bp1)[1]!=0){lines(dat.bp1[,1],model.in$LMs[[data.choice]]$Bp1$fitted,col=col.bpts[1], lwd=2,lty=2)}
    if(dim(dat.bp2)[1]!=0){lines(dat.bp2[,1],model.in$LMs[[data.choice]]$Bp2$fitted,col=col.bpts[2],lwd=2,lty=2)}
    #if(dim(age.break.hi.dat)[1]!=0){lines(age.break.hi.dat[,1],model.in$LMs[[data.choice]]$NoBp$fitted,col="black",lwd=2)}
  }  
}

#Get numbers
Spp.samp.size<-function(dat.list.in)
{
  Spp.N<-matrix(NA,nrow=length(dat.list.in),ncol=3)
  colnames(Spp.N)<-c("A","F","M")
  rownames(Spp.N)<-names(dat.list.in)
  for(i in 1:length(dat.list.in))
  {Spp.N[i,]<-c(dim(dat.list.in[[i]]$Data$All)[1],dim(dat.list.in[[i]]$Data$Females)[1],dim(dat.list.in[[i]]$Data$Males)[1])}
  return(Spp.N)
}

#Get R^2 values
Spp.r2<-function(dat.list.in)
{
  rtwos<-as.data.frame(matrix(NA,nrow=2*length(dat.list.in),ncol=5))
  colnames(rtwos)<-c("Species","Sex","PW1","PW2","LM")
  genders<-c("Female","Male")
  spp.names<-names(dat.list.in)
  x<-1
  for(i in 1:length(dat.list.in))
  {
    for(g in 1:2)
    {
      for(ii in 1:3)
      {
        rtwos[x,1]<-spp.names[i]
        rtwos[x,2]<-genders[g]
        if(is.na(dat.list.in[[i]]$LMs[[g+1]][[ii]])==FALSE){rtwos[x,2+ii]<-summary(dat.list.in[[i]]$LMs[[g+1]][[ii]])$r.squared}
      }
      x<-1+x
    }
  }
  return(rtwos)
}


#Calcualte predicted ages
#low.breaks: matrix of values with breaks for the low ages (col 1= wt; col 2= age)
#intercept: If > 0, use non-zero intercept
Pred.ages.oto<-function(dat.list.in,low.breaks,intercept=999)
{
  Pred.ages<-list()
  for(i in 1:length(dat.list.in))
  {
    AFM.oto.age.temp<-list()
    for(ii in 1:3)
    {
      #Define low breaks
      low.ages.ind<-dat.list.in[[i]]$Data[[ii]][,1]<=low.breaks[[i]][nrow(low.breaks[[i]]),1]
      if(length(dat.list.in[[i]]$Data[[ii]][low.ages.ind,])>2){low.ages.dat<-dat.list.in[[i]]$Data[[ii]][low.ages.ind,]}
      if(length(dat.list.in[[i]]$Data[[ii]][low.ages.ind,])<=2){low.ages.dat<-matrix(dat.list.in[[i]]$Data[[ii]][low.ages.ind,],nrow=1,ncol=2)}
      low.ages<-matrix(NA,nrow=nrow(low.ages.dat),ncol=2)
      bp.ages.ind<-dat.list.in[[i]]$Data[[ii]][,1]>low.breaks[[i]][nrow(low.breaks[[i]]),1]
      if(length(dat.list.in[[i]]$Data[[ii]][bp.ages.ind,])>2){bp.ages.dat<-dat.list.in[[i]]$Data[[ii]][bp.ages.ind,]}
      if(length(dat.list.in[[i]]$Data[[ii]][bp.ages.ind,])<=2){bp.ages.dat<-matrix(dat.list.in[[i]]$Data[[ii]][bp.ages.ind,],nrow=1,ncol=2)}
      bp.ages<-matrix(NA,nrow=nrow(bp.ages.dat),ncol=2)
      #Predict the low end breaks
      for(iii in 1:nrow(low.breaks[[i]]))
      {
        if(iii==1)
        {
          low.ages[low.ages.dat[,1]>0&low.ages.dat[,1]<=low.breaks[[i]][iii,1],]<-low.breaks[[i]][iii,2]
          
        }
        
        else{low.ages[low.ages.dat[,1]>low.breaks[[i]][iii-1,1]&low.ages.dat[,1]<=low.breaks[[i]][iii,1],]<-low.breaks[[i]][iii,2]}
      }
      #Predict the piecewise breaks
      bp1.temp<-bp.ages.dat[,1]<=dat.list.in[[i]]$Bps[1]
      bp2.temp<-bp.ages.dat[,1]>dat.list.in[[i]]$Bps[1]
      if(intercept==0)
      {
        if(length(dat.list.in[[i]]$LMs[[ii]]$Bp1)>1){bp.ages[bp1.temp,1]<-dat.list.in[[i]]$LMs[[ii]]$Bp1$coef*bp.ages.dat[bp1.temp,1]}
        if(length(dat.list.in[[i]]$LMs[[ii]]$Bp2)>1){bp.ages[bp2.temp,1]<-dat.list.in[[i]]$LMs[[ii]]$Bp2$coef*bp.ages.dat[bp2.temp,1]}
      }
      if(intercept!=0)
      {
        if(length(dat.list.in[[i]]$LMs[[ii]]$Bp1)>1){bp.ages[bp1.temp,1]<-dat.list.in[[i]]$LMs[[ii]]$Bp1$coef[2]*bp.ages.dat[bp1.temp,1]+dat.list.in[[i]]$LMs[[ii]]$Bp1$coef[1]}
        if(length(dat.list.in[[i]]$LMs[[ii]]$Bp2)>1){bp.ages[bp2.temp,1]<-dat.list.in[[i]]$LMs[[ii]]$Bp2$coef[2]*bp.ages.dat[bp2.temp,1]+dat.list.in[[i]]$LMs[[ii]]$Bp2$coef[1] }
      }
      #Predict from linear model
      if(intercept==0){bp.ages[,2]<-dat.list.in[[i]]$LMs[[ii]]$NoBp$coef*bp.ages.dat[,1]}
      if(intercept!=0){bp.ages[,2]<-dat.list.in[[i]]$LMs[[ii]]$NoBp$coef[2]*bp.ages.dat[,1]+dat.list.in[[i]]$LMs[[ii]]$NoBp$coef[1]}
      low.bp.ages<-matrix(NA,nrow=nrow(dat.list.in[[i]]$Data[[ii]]),ncol=2)
      low.bp.ages[low.ages.ind]<-low.ages
      low.bp.ages[bp.ages.ind]<-bp.ages
      names(low.bp.ages)<-c("Bp","noBp")
      AFM.oto.age.temp[[ii]]<-as.data.frame(cbind(dat.list.in[[i]]$Data[[ii]],low.bp.ages,dat.list.in[[i]][[5]][[ii]]))
      colnames(AFM.oto.age.temp[[ii]])<-c("OtoWt","Age","PW_age","Lm_age","Length")
    }
    names(AFM.oto.age.temp)<-c("A","F","M")
    Pred.ages[[i]]<-AFM.oto.age.temp
  }
  names(Pred.ages)<-names(dat.list.in)
  return(Pred.ages)
}

#Piecewise models
PReg.obj <- function(x,x.in,y.in,rngSplit=c(75,180),fixSplt=NULL,verbose=F) {
  #piecewise regression, assuming continuous function at breakpoint.
  #Therefore, same intercept for both equations
  #note: b1 and b2 are slopes and alpha is the y-intercept. splt is the breakpoint.
  splt <- x[1]
  alpha <- x[2]
  b1 <- x[3]
  b2 <- x[4]
  #nobs = number observed
  nobs <- length(x.in)
  #creates an empty vecter with 'nobs' number of elements.
  predY <- rep(NA,nobs)
  #take lower half of data below breakpoint guess
  ind <- x.in <= splt
  
  #plug data into line equation
  predY[ind] <- alpha+b1*(x.in[ind]-splt)
  #difference of squares
  f.left <- sum((predY[ind]-y.in[ind])^2)
  #same as above, but with upper half of data
  ind <- x.in > splt
  predY[ind] <- alpha+b2*(x.in[ind]-splt)
  f.right <- sum((predY[ind]-y.in[ind])^2)
  
  #Sum all the differences of squares (this is to be minimized)
  f <- (f.left+f.right)
  if(verbose){cat(f.left,f.right,f,"\n")}
  
  #if breakpoint is out of domain, return arbitrarily large number as a flag, and so optim will avoid
  if(splt<rngSplit[1] | splt>rngSplit[2]) {
    return(1e10)
  }
  return((nobs/2)*log(f/nobs))    #the likelihood function (assuming normal errors and constant variance)
}



PRegLikeProf.fn <- function(x,x.in,y.in,rngSplit,steppin) {
  PRegLP.obj <- function(x,breakpt,x.in,y.in)
  {
    #piecewise regression, assuming continuous function at breakpoint.
    #Therefore, same intercept for both equations
    #Likelihood profile over splt
    splt <- breakpt
    alpha <- x[1]
    b1 <- x[2]
    b2 <- x[3]
    nobs <- length(x.in)
    predY <- rep(NA,nobs)
    ind <- x.in <= splt
    predY[ind] <- alpha+b1*(x.in[ind]-splt)
    f.left <- sum((predY[ind]-y.in[ind])^2)
    ind <- x.in > splt
    predY[ind] <- alpha+b2*(x.in[ind]-splt)
    f.right <- sum((predY[ind]-y.in[ind])^2)
    f <- (f.left+f.right)
    return((nobs/2)*log(f/nobs))    #the likelihood function (assuming normal errors and constant variance)
  }
  nobs <- length(x.in)
  splits <- seq(rngSplit[1],rngSplit[2],by=steppin)
  out <- matrix(NA,nrow=length(splits),ncol=5,dimnames=list(NULL,c("split","alpha","b1","b2","f")))
  for(i in 1:length(splits))
  {
    opt <- optim(x[-1],PRegLP.obj,breakpt=splits[i],x.in=x.in,y.in=y.in)
    out[i,] <- c(splits[i],opt$par,opt$value)
  }
  as.data.frame(out)
}

#Plot ages comparisons (1:1)
#gen.type:1=All; 2=Females; 3=Males
#age.comp.type:3=PW; 4=LM
#rd.type:1= none; 2=round; 3=ceiling
plot.age2age<-function(spp.dat,spp.dat.nam,gen.type=1,age.comp.type=3,rd.type=1,col.dots="black" ,lab.nam=c("T","T"),add.plot="F")
{
  name.nums<-c(1:length(names(spp.dat)))
  name.num<-name.nums[names(spp.dat)==spp.dat.nam]
  xlab.nam<-ylab.nam<-""
  if(lab.nam[1]=="T"){xlab.nam<-"Age from otolith annuli"}
  if(lab.nam[2]=="T"){ylab.nam<-"Age from otolith weight"}
  if(add.plot=="F")
  {
    if(rd.type==1)
    {
      plot(Spp.age.pred[[name.num]][[gen.type]][,2],Spp.age.pred[[name.num]][[gen.type]][,age.comp.type],xlim=c(0,max(Spp.age.pred[[name.num]][[gen.type]][,c(-1,-5)])),ylim=c(0,max(Spp.age.pred[[name.num]][[gen.type]][,c(-1,-5)])),xlab=xlab.nam, ylab=ylab.nam,pch=21,bg=col.dots)
      abline(lm(Spp.age.pred[[name.num]][[gen.type]][,age.comp.type]~Spp.age.pred[[name.num]][[gen.type]][,2]),lty=2,lwd=1,col=col.dots)
    }
    if(rd.type==2)
    {
      plot(Spp.age.pred[[name.num]][[gen.type]][,2],round(Spp.age.pred[[name.num]][[gen.type]][,age.comp.type]),xlim=c(0,max(Spp.age.pred[[name.num]][[gen.type]][,c(-1,-5)])),ylim=c(0,max(Spp.age.pred[[name.num]][[gen.type]][,c(-1,-5)])),xlab=xlab.nam, ylab=ylab.nam,"Age from otolith weight",pch=21,bg=col.dots)
      abline(lm(round(Spp.age.pred[[name.num]][[gen.type]][,age.comp.type])~Spp.age.pred[[name.num]][[gen.type]][,2]),lty=2,lwd=1,col=col.dots)
    }
    if(rd.type==3)
    {
      plot(Spp.age.pred[[name.num]][[gen.type]][,2],ceiling(Spp.age.pred[[name.num]][[gen.type]][,age.comp.type]),xlim=c(0,max(Spp.age.pred[[name.num]][[gen.type]][,c(-1,-5)])),ylim=c(0,max(Spp.age.pred[[name.num]][[gen.type]][,c(-1,-5)])),xlab=xlab.nam, ylab=ylab.nam,"Age from otolith weight",pch=21,bg=col.dots)
      abline(lm(ceiling(Spp.age.pred[[name.num]][[gen.type]][,age.comp.type])~Spp.age.pred[[name.num]][[gen.type]][,2]),lty=2,lwd=1,col=col.dots)
    }
    abline(a=0,b=1,col="black",lwd=2)
  }
  if(add.plot=="T")
  {
    if(rd.type==1)
    {
      points(Spp.age.pred[[name.num]][[gen.type]][,2],Spp.age.pred[[name.num]][[gen.type]][,age.comp.type],pch=21,bg=col.dots)
      abline(lm(Spp.age.pred[[name.num]][[gen.type]][,age.comp.type]~Spp.age.pred[[name.num]][[gen.type]][,2]),lty=2,lwd=1,col=col.dots)
    }
    if(rd.type==2)
    {
      points(Spp.age.pred[[name.num]][[gen.type]][,2],round(Spp.age.pred[[name.num]][[gen.type]][,age.comp.type]),pch=21,bg=col.dots)
      abline(lm(round(Spp.age.pred[[name.num]][[gen.type]][,age.comp.type])~Spp.age.pred[[name.num]][[gen.type]][,2]),lty=2,lwd=1,col=col.dots)
    }
    if(rd.type==3)
    {
      points(Spp.age.pred[[name.num]][[gen.type]][,2],ceiling(Spp.age.pred[[name.num]][[gen.type]][,age.comp.type]),pch=21,bg=col.dots)
      abline(lm(ceiling(Spp.age.pred[[name.num]][[gen.type]][,age.comp.type])~Spp.age.pred[[name.num]][[gen.type]][,2]),lty=2,lwd=1,col=col.dots)
    }
    abline(a=0,b=1,col="black",lwd=2)
  }
}

### VBGF functions
VBGF<-function(Linf,k,t0,ages)
{
  Lengths_exp<-Linf*(1-exp(-k*(ages-t0)))
  return(Lengths_exp)
}

VBGF.fit<-function(p,obs,return.type=2)
{
  if(return.type==1)
  {
    exp.lts<-VBGF(p[1],p[2],p[3],obs[,1])
    return(exp.lts)
  }
  if(return.type==2)
  {
    exp.lts<-VBGF(p[1],p[2],p[3],obs[,1])
    sigma<-sqrt(sum((obs[,2]-exp.lts)^2)/length(obs[,2]))
    neglogsum<-sum(-log(dnorm(exp.lts,obs[,2],sigma)))
    return(neglogsum)
  }
}

VBGF.fit.plot<-function(AG.age.in,AG.lt.in,max.x,max.y,pch.col="red",title.in="",add.plot="F",print.vbgf="T")
{
  fitvbgf<-nls(lts~Linf*(1-exp(-k*(ages-t0))),data=list(ages=AG.age.in,lts=AG.lt.in),start=list(Linf=max(AG.lt.in,na.rm=TRUE),k=0.1,t0=0),control = list(reltol=0.00000000001))
  AG.fit.par<-summary(fitvbgf)$par
  exp.AG<-VBGF(AG.fit.par[1],AG.fit.par[2],AG.fit.par[3],c(0:max(AG.age.in)))
  if(add.plot=="F")
  {
    plot(AG.age.in,AG.lt.in,xlab="Age (years)", ylab="Length (cm)",xlim=c(0,max.x),ylim=c(0,max.y),pch=21,bg=pch.col,col="black",main=title.in)
    lines(c(0:max(AG.age.in)),exp.AG,lwd=2)
  }
  if(add.plot=="T")
  {
    points(AG.age.in,AG.lt.in,pch=21,bg=pch.col,col="black",main=title.in)
    lines(c(0:max(AG.age.in)),exp.AG,lwd=2,lty=2)
  }
  if(print.vbgf=="T")
  {
    text(0.75*max.x,0.25*max.y,labels=bquote(paste(L[infinity],"= ",.(round(AG.fit.par[1],2)))))
    text(0.75*max.x,0.2*max.y,labels=paste("k= ",round(AG.fit.par[2],2)))
    text(0.75*max.x,0.15*max.y,labels=bquote(paste(t[0],"= ",.(round(AG.fit.par[3],2)))))
  }
  return(AG.fit.par)
}
###

vbgf.pars<-function(spp.dat.in, max.x,max.y)
{
  gender<-c("A", "F", "M")
  age.type<-c("OtA","PwA","LmA")
  num.type<-c("","Rd","Cl")
  temp.out<-list(matrix(NA,nrow=21,ncol=3),matrix(NA,nrow=21,ncol=3))
  r.names<-rep(NA,21)
  colnames(temp.out[[1]])<-colnames(temp.out[[2]])<-c("Linf","k","t0")
  x<-1
  for(i in 1:length(gender))
  {
    for(ii in 1:length(age.type))
    {
      if(ii==1)
      {
        temp.out.parms<-VBGF.fit.plot(spp.dat.in[[i]][[ii+1]],spp.dat.in[[i]][[5]],round(max(spp.dat.in[[i]][[ii+1]])/10)*10,round(max(spp.dat.in[[i]][[5]],na.rm=TRUE)/10)*10,"gray","")
        temp.out[[1]][x,]<-temp.out.parms[,1]
        temp.out[[2]][x,]<-temp.out.parms[,2]
        r.names[x]<-paste(gender[i],"_",age.type[ii],sep="")
        x<-x+1
      }
      if(ii>1)
      {
        temp.out.parms<-VBGF.fit.plot(spp.dat.in[[i]][[ii+1]],spp.dat.in[[i]][[5]],round(max(spp.dat.in[[i]][[ii+1]])/10)*10,round(max(spp.dat.in[[i]][[5]],na.rm=TRUE)/10)*10,"gray","")
        temp.out[[1]][x,]<-temp.out.parms[,1]
        temp.out[[2]][x,]<-temp.out.parms[,2]
        r.names[x]<-paste(gender[i],"_",age.type[ii],"_",num.type[1],sep="")
        x<-x+1
        temp.out.parms<-VBGF.fit.plot(round(spp.dat.in[[i]][[ii+1]]),spp.dat.in[[i]][[5]],round(max(spp.dat.in[[i]][[ii+1]])/10)*10,round(max(spp.dat.in[[i]][[5]],na.rm=TRUE)/10)*10,"gray","")
        temp.out[[1]][x,]<-temp.out.parms[,1]
        temp.out[[2]][x,]<-temp.out.parms[,2]
        r.names[x]<-paste(gender[i],"_",age.type[ii],"_",num.type[2],sep="")
        x<-x+1
        temp.out.parms<-VBGF.fit.plot(ceiling(spp.dat.in[[i]][[ii+1]]),spp.dat.in[[i]][[5]],round(max(spp.dat.in[[i]][[ii+1]])/10)*10,round(max(spp.dat.in[[i]][[5]],na.rm=TRUE)/10)*10,"gray","")
        temp.out[[1]][x,]<-temp.out.parms[,1]
        temp.out[[2]][x,]<-temp.out.parms[,2]
        r.names[x]<-paste(gender[i],"_",age.type[ii],"_",num.type[3],sep="")
        x<-x+1
      }
    }
  }
  rownames(temp.out[[1]])<-rownames(temp.out[[2]])<-r.names
  names(temp.out)<-c("Est","Std")
  return(temp.out)
}

vbgf.par.comp.plot<-function(spp.par.in,rows2use=c(1:21),vbgf.num,labelx="Y",dot.col=c("black","orange","orange","orange","blue","blue","blue"))
{
  par.names<-expression(L[inf],"k",t[0])
  plotCI(spp.par.in$Est[rows2use,vbgf.num],uiw=spp.par.in$Std[rows2use,vbgf.num]*1.96,ylim=c(0, max(spp.par.in$Est[rows2use,vbgf.num]+spp.par.in$Est[rows2use,vbgf.num]*0.25)),pch=21,pt.bg=dot.col,xlab="",ylab=par.names[vbgf.num],axes=F)
  box()
  if(labelx=="Y"){axis(1,at=c(1:dim(spp.par.in$Est)[1]),labels=rownames(spp.par.in$Est)[rows2use],las=2,cex=0.75)}
  if(labelx=="N"){axis(1,at=c(1:dim(spp.par.in$Est)[1]),labels=rep("",dim(spp.par.in$Est)[1]))}
  axis(2)
  if(length(rows2use)/7==3){abline(v=c(7.5,14.5),lty=2)}
  if(length(rows2use)/7==2){abline(v=c(7.5),lty=2)}
}

#Relative errors
Spp.RE<-function(spp.par.in,vbgf.num=1)
{
  F_u<-(spp.par.in$Est[c(9,12),vbgf.num]-spp.par.in$Est[8,vbgf.num])/spp.par.in$Est[8,vbgf.num]
  F_uw<-((spp.par.in$Est[c(9,12),vbgf.num]+spp.par.in$Std[c(9,12),vbgf.num]*1.96)-spp.par.in$Est[8,vbgf.num])/spp.par.in$Est[8,vbgf.num]
  M_u<-(spp.par.in$Est[c(16,19),vbgf.num]-spp.par.in$Est[15,vbgf.num])/spp.par.in$Est[15,vbgf.num]
  M_uw<-((spp.par.in$Est[c(16,19),vbgf.num]+spp.par.in$Std[c(16,19),vbgf.num]*1.96)-spp.par.in$Est[15,vbgf.num])/spp.par.in$Est[15,vbgf.num]
  FM.out<-rbind(F_u,F_uw,M_u,M_uw)
  colnames(FM.out)<-c("PW","LM")
  return(FM.out)
}
#############################################
#############################################
#############################################
Punt_age_prep<-function(Spp.Ages.in)
{
  age.in<-table(round(Spp.Ages.in))
  age.out<-matrix(NA,nrow=length(age.in),ncol=5)
  age.out[,c(1,4)]<-0
  age1<-c(1:dim(age.in)[1])
  age2<-c(1:dim(age.in)[2])
  x<-1
  for(i in 1:length(age1))
  {
    for(ii in 1:length(age2))
    {
      age.out[x,c(2,3,5)]<-c(age1[i],age2[ii],age.in[i,ii])
      x<-x+1
    }
  }
  age.return<-subset(age.out,age.out[,5]>0)
  colnames(age.return)<-c("0","Age1","Age2","0","N")
  return(age.return)
}

#plot SD from Punt program
plot.AErr.comp<-function(dat.in,species.in,col.in="orange")
{
  dat.plot<-subset(dat.in,Species==species.in)
  plot(subset(dat.plot,Reader==1)$SD,subset(dat.plot,Reader==2)$SD,xlim=c(0,round(max(c(subset(dat.plot,Reader==1)$SD,subset(dat.plot,Reader==2)$SD)))),ylim=c(0,round(max(c(subset(dat.plot,Reader==1)$SD,subset(dat.plot,Reader==2)$SD)))),xlab="",ylab="",pch=21,bg=col.in,cex=1.25)
  abline(a=0,b=1,lwd=2,col="red")
}

#New two-breakpoint function for piecewise regression with continuity at breakpoints

PReg.obj.2<-function(x,x.in,y.in,rngSplit=c(75,180),fixSplt=NULL,verbose=F){
  bkpt1<-x[1] #first breakpoint initial parameter
  bkpt2<-x[2] #second breakpoint initial parameter
  y.int.1<-x[3] #y-intercept initial parameters
  m1<-x[4] #first slope initial parameter
  m2<-x[5] #second slope initial parameter
  m3<-x[6] #third slope initial parameter
  
  nobs <- length(x.in) #nobs = number observed
  predY <- rep(NA,nobs) #empty vector to be filled with model y-values
  
  ind1 <- x.in <= bkpt1 #retrieve indexes before first breakpoint
  predY[ind1] <- y.int.1+m1*(x.in[ind1]-bkpt1) #indecies before first breakpoint
  f.left <- sum((predY[ind1]-y.in[ind1])^2) #residual difference of squares
  
  ind2 <- (x.in > bkpt1) & (x.in <= bkpt2)
  predY[ind2] <- y.int.1+m2*(x.in[ind2]-bkpt1) #repeat of above, but with middle data
  f.mid <- sum((predY[ind2]-y.in[ind2])^2)
  
  ind3 <- x.in > bkpt2
  predY[ind3] <- y.int.1+m3*(x.in[ind3]-bkpt2)+m2*(bkpt2-bkpt1)
  f.right <- sum((predY[ind3]-y.in[ind3])^2)
  
  f <- (f.left+f.mid+f.right)
  if(verbose){cat(f.left,f.mid,f.right,f,"\n")}
  
  #if breakpoint is out of domain, return arbitrarily large number as a flag, and so optim will avoid
  if(bkpt1<rngSplit[1] | bkpt1>rngSplit[2] | bkpt2<rngSplit[1] | bkpt2>rngSplit[2]) {
    return(1e10)
  }
  
  return((nobs/2)*log(f/nobs))
}

match.f <- function (file, table, findex = 1, tindex = findex, tcol = NULL, 
                     round. = T, digits=0) 
{
  paste.col <- function(x) {
    if (is.null(dim(x))) 
      return(paste(as.character(x)))
    out <- paste(as.character(x[, 1]))
    for (i in 2:ncol(x)) {
      out <- paste(out, as.character(x[, i]))
    }
    out
  }
  if (is.null(dim(file))) {
    dim(file) <- c(length(file), 1)
  }
  if (round.) {
    for (i in findex) {
      if (is.numeric(file[, i])) 
        file[, i] <- round(file[, i], digits)
    }
    for (i in tindex) {
      if (is.numeric(table[, i])) 
        table[, i] <- round(table[, i], digits)
    }
  }
  if (is.null(tcol)) 
    tcol <- dimnames(table)[[2]][!(dimnames(table)[[2]] %in% 
                                     tindex)]
  cbind(file, table[match(paste.col(file[, findex]), paste.col(table[, 
                                                                     tindex])), tcol, drop = F])
}


# ------------------------------- Latest load() function -------------------------------------------------------------------------------------------

load <- function (file, str. = TRUE, list.len = 15, nrow = 5, ...) 
{
  ls.ext <- function(file, str. = TRUE, list.len, nrow) {
    local({
      base::load(file)
      if (str. == TRUE) {
        Names <- base::ls()
        for (i in Names) {
          OBJ <- eval(parse(text = i))
          cat("\n", i, ":\n\n", sep = "")
          str(OBJ, list.len = list.len)
          cat("\n")
          if(is.matrix(OBJ) | is.data.frame(OBJ)) { 
            print(head(OBJ, nrow)); flush.console(); cat("\n") 
            print(dim(OBJ)); flush.console(); cat("\n")
          }
        }
        rm(i, Names)
        invisible(base::ls())
      }
      else base::ls()
    })
  }
  base::load(file, .GlobalEnv, ...)
  ls.ext(file, str. = str., list.len = list.len, nrow = nrow)
}
