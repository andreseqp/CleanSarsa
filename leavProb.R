# ------------------ Effect of leaving probability ------------------------ #

# Directories --------------------------------------------------------------

projDir<-"d:/quinonesa/learning_models_c++/actCrit/"
simsDir<-"s:/quinonesa/Simulations/actCrit/"


# libraries ----------------------------------------------------------------
source('d:/quinonesa/Dropbox/R_files/posPlots.R')
source(paste(projDir,"aesth_par.R",sep=""))
source(paste(projDir,"loadData.R",sep = ""))
library('plotrix')
# library('lme4')


# Load data ------------------------------------------------------------


# Define data to be loaded 

(listPar<-c("LeavingP"))

(listVal<-c(""))

FIAlastQuart<-do.call(rbind,lapply(
  getFilelist(simsDir,listPar,listVal)$FIA,file2lastProp,0.75,'Vlp'))

FIA.stats<-FIAlastQuart[,.(meanProb=mean(Prob.RV.V),
                              upIQR=fivenum(Prob.RV.V)[4],
                              lowIQR=fivenum(Prob.RV.V)[2])
                           ,by=.(Neta,Gamma,pR,pV,Outbr,Vlp)]


FIAraw<-loadRawData(simsDir,"FIA",listparam = listPar,values = listVal)
param<-getParam(simsDir,listparam = listPar,values = listVal)


FIAagg<-FIAraw[, as.list(unlist(lapply(.SD, function(x) 
  list(mean = mean(x),IQ.h = fivenum(x)[4],IQ.l=fivenum(x)[2])))),
  by=.(Age,Alpha,Gamma,Tau,Neta,Outbr,pR,pV), 
  .SDcols=c('ThetaV','ThetaR','RV','VV','RR','R0','V0','00_')]

setnames(FIAagg,'get',extpar)


# Plots -----------------------------------------------------------------------

FIA.stats[,posit:=ifelse(Gamma==0&Neta==0,0,
                         ifelse(Gamma==0.8&Neta==0,0.01,
                                ifelse(Gamma==0&Neta==1,0.02,0.03)))]

png("d:/quinonesa/Dropbox/Neuchatel/Figs/Actor_critic/Fig3.png",width = 1200,
    height = 1200)

par(plt=posPlot(),las=1)
with(FIA.stats,{
  plotCI(x = Vlp+posit,
         y = meanProb,ui = upIQR
         ,li = lowIQR,
         col=colboxes[ifelse(Gamma==0.8,
                             ifelse(Neta==1,1,2),
                             ifelse(Neta==1,3,4))],
         pch=16,xlab="",ylab="",
         sfrac=0.008,yaxt='s',
         cex.axis=1.3,cex=2,cex.lab=3)
  mtext('Probability of V over R',2,line = 4, cex=3,las=0)
  mtext(expression(l[v]),1,line = 5,las=1,cex=4)
})

legend('topleft',
       legend=c("neg. reward + future", "future",
                "neg. reward","no neg. reward + no future"),
       col=colboxes,pch=15,cex=1.5,ncol=1)

dev.off()
