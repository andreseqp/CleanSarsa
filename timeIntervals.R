# --------------------------------- Time intervals -------------------------------------------------#

projDir<-"d:/quinonesa/learning_models_c++/Sarsa/"
simsDir<-"s:/quinonesa/Simulations/Basic_sarsa/"

# libraries ---------------------------------------------------------------------------------------
source('d:/quinonesa/Dropbox/R_files/posPlots.R')
source(paste(projDir,"aesth_par.R",sep=""))
source(paste(projDir,"loadData.R",sep = ""))
library('plotrix')

# Load Data ---------------------------------------------------------------------------------------


# Define data to be loaded 

(listPar<-rep("factRew",2))
(listVal<-c(1,2))
param<-getParam(simsDir,listparam = listPar,values = listVal)

#diffJsons(param[1],param[3])

list.files(simsDir,recursive = TRUE,pattern = "Olle_/")

# Load interval data for FIA from the raw data
FIAtimeInt<-do.call(
  rbind,lapply(
    getFilelist(simsDir,listPar,listVal)$FIA,
    file2timeInter,interV=1001))

# Load FIA data from processed file

# getFilelist(simsDir,listPar,listVal)$FIA

# FIAtimeInt<-do.call(
#   rbind,lapply(getFilelist(projDir,listPar,listVal)$FIA,fread))



# Load interval data for PIA from the raw data
PIAtimeInt<-do.call(
  rbind,lapply(
    getFilelist(simsDir,listPar,listVal)$PIA,
    file2timeInter,interV=1001))

# Load PIA data from processed file

# PIAtimeInt<-do.call(
#   rbind,lapply(getFilelist(genDir,listPar,listVal)$PIA,fread))

# Load DP data from the raw data
DPdataProb<-do.call(rbind,
                    lapply(getFilelist(simsDir,listPar,listVal)$DP,
                           file2lastDP))

# Load DP data from processed file

# DPdataprob<-do.call(
#   rbind,lapply(getFilelist(genDir,listPar,listVal)$DP,fread))

# Plot the dynamics of VR choice -----------------------------------------------------------

extpar<-listPar[1]

FIAIntstats<-FIAtimeInt[,.(meanProb=mean(Prob.RV.V),
                           upIQR=fivenum(Prob.RV.V)[4],
                           lowIQR=fivenum(Prob.RV.V)[2])
                        ,by=.(Interv,Neta,Outbr,Tau,Gamma,get(extpar))]
setnames(FIAIntstats,'get',extpar)
PIAIntstats<-PIAtimeInt[,.(meanProb=mean(Prob.RV.V),
                           upIQR=fivenum(Prob.RV.V)[4],
                           lowIQR=fivenum(Prob.RV.V)[2])
                        ,by=.(Interv,Neta,Outbr,Tau,Gamma,get(extpar))]
setnames(PIAIntstats,'get',extpar)

par(plt=posPlot(numplotx = 2,idplotx = 1)+c(-0.05,-0.05,0,0),yaxt='s',las=1)
with(FIAIntstats[Neta==0.5&factRew==2],{
  plotCI(x=Interv,y=meanProb,
         ui = upIQR,li=lowIQR,
         pch=16,xlab='',ylab='',
         col=colboxes[match(Gamma,unique(Gamma))],
         sfrac=0.002,cex.axis=1.3,ylim=c(0,1))
  lines(x=c(0,max(Interv)),y=c(0.5,0.5),col='grey')
})

with(DPdataProb,  
     {matlines(x = t(matrix(rep(max(FIAtimeInt$Interv)*c(0.75,1),
                                each=length(RV.V)),length(RV.V))),
               y=t(matrix(rep(probRV.V,2),length(RV.V))),
               lwd=2,lty = "dashed",
               col=colboxes[match(Gamma,unique(Gamma))])})

legend('topright',
       legend=unique(FIAIntstats[,Gamma])[order(unique(FIAIntstats[,Gamma]))],
              col=colboxes,pch=15,
              title=eval(extpar),cex=1.5,ncol=3)

par(plt=posPlot(numplotx = 2,idplotx = 2)+c(-0.05,-0.05,0,0),
    new=TRUE,yaxt='s',xpd=TRUE)
with(PIAIntstats,{
  plotCI(x=Interv,y=meanProb,
         ui = upIQR,li=lowIQR,
         pch=16,xlab='',ylab='',
         col=colboxes[match(Gamma,unique(Gamma))],
         sfrac=0.002,cex.axis=1.3,yaxt='n')
  lines(x=c(0,max(Interv)),y=c(0.5,0.5),col='grey')
  axis(side=4,cex.axis=1.3)
})

png(filename = paste(projDir,eval(extpar),".png",sep=""))

dev.off()

# Replicating Olle's simulations -----------------------------------------------

colOlle<-c("red","green")
colOlle2<-c("blue","black")

cexpar<-1.5

png(filename = "d:/quinonesa/Dropbox/Neuchatel/Olle/Sarsa_Olle_par.png",
    width = 800,height = 800)

par(plt=posPlot(numplotx = 1,idplotx = 1),yaxt='s',las=1)
with(FIAIntstats[Neta==0.5&factRew==2],{
  plotCI(x=Interv,y=meanProb,
         ui = upIQR,li=lowIQR,
         pch=16,xlab='Time',ylab='Prob. V over R',cex.lab=2,
         col=colOlle[match(Gamma,unique(Gamma))],
         sfrac=0.002,cex.axis=1.3,ylim=c(0,1),cex=cexpar)
  lines(x=c(0,max(Interv)),y=c(0.5,0.5),col='grey')
})
par(new=TRUE)
with(FIAIntstats[Neta==0&factRew==1],{
  plotCI(x=Interv,y=meanProb,
         ui = upIQR,li=lowIQR,
         pch=16,xlab='',ylab='',
         col=colOlle2[match(Gamma,unique(Gamma))],
         sfrac=0.002,cex.axis=1.3,ylim=c(0,1),cex=cexpar)
  lines(x=c(0,max(Interv)),y=c(0.5,0.5),col='grey')
})

with(DPdataProb[Neta==0.5&factRew==2],  
     {matlines(x = t(matrix(rep(max(FIAtimeInt$Interv)*c(0.75,1),
                                each=length(RV.V)),length(RV.V))),
               y=t(matrix(rep(probRV.V,2),length(RV.V))),
               lwd=2,lty = "dashed",cex=1.2,
               col=colOlle[match(Gamma,unique(Gamma))])})

with(DPdataProb[Neta==0&factRew==1],  
     {matlines(x = t(matrix(rep(max(FIAtimeInt$Interv)*c(0.75,1),
                                each=length(RV.V)),length(RV.V))),
               y=t(matrix(rep(probRV.V,2),length(RV.V))),
               lwd=2,lty = "dashed",cex=1.2,
               col=colOlle2[match(Gamma,unique(Gamma))])})

legend('bottomright',
       legend=c("Punishment and future", "punishment",
                "future","no punishment no future"),
       col=c(colOlle,colOlle2),pch=15,cex=1.5,ncol=1)

dev.off()




  
  