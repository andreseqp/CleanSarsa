# --------------------------------- Time intervals -------------------------------------------------#
library(here)
projDir<-here()
simsDir<-"e:/NeuchSims/SARSA/AbundLvpTau0.5_/"

# libraries ---------------------------------------------------------------------------------------
source(here("..",'/R_files/posPlots.R'))
source(paste(projDir,"/aesth_par.R",sep=""))
source(paste(projDir,"/loadData.R",sep = ""))
library('plotrix')

# Load Data ---------------------------------------------------------------------------------------


# Define data to be loaded 

(listPar<-c("pR","pV"))
(listVal<-c("0.2","0.2"))
(param<-getParam(simsDir,listparam = listPar,values = listVal))

#diffJsons(param[1],param[3])

tail(list.files(simsDir,recursive = TRUE,pattern = ""),100)

interv<-1001

# Load interval data for FIA from the raw data
FIAtimeInt<-do.call(
  rbind,lapply(
    getFilelist(simsDir,listPar,listVal)$FIA[41:44],
    file2timeInter,interV=interv))



# Load FIA data from processed file

# getFilelist(simsDir,listPar,listVal)$FIA

# FIAtimeInt<-do.call(
  # rbind,lapply(getFilelist(projDir,listPar,listVal)$FIA,fread))



# Load interval data for PIA from the raw data
PIAtimeInt<-do.call(
  rbind,lapply(
    getFilelist(simsDir,listPar,listVal)$PIA,
    file2timeInter,interV=interv))

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

# Calculate stats --------------------------------------------------------------

# extpar<-listPar[1]

FIAIntstats<-FIAtimeInt[,.(meanProb=mean(Prob.RV.V),
                           upIQR=fivenum(Prob.RV.V)[4],
                           lowIQR=fivenum(Prob.RV.V)[2])
                        ,by=.(Interv,Neta,Outbr,Tau,Gamma)]
# setnames(FIAIntstats,'get',extpar)

PIAIntstats<-PIAtimeInt[Alpha==0.005,.(meanProb=mean(Prob.RV.V),
                           upIQR=fivenum(Prob.RV.V)[4],
                           lowIQR=fivenum(Prob.RV.V)[2])
                        ,by=.(Interv,Neta,Outbr,Tau,Gamma)]
# setnames(PIAIntstats,'get',extpar)

# Plot the dynamics of VR choice -----------------------------------------------------------

FIAIntstats[,posit:=ifelse(Gamma==0&Neta==0,0,
                           ifelse(Gamma==0.8&Neta==0,0.1,
                                  ifelse(Gamma==0&Neta==1,0.2,0.3)))]
PIAIntstats[,posit:=ifelse(Gamma==0&Neta==0,0,
                           ifelse(Gamma==0.8&Neta==0,0.1,
                                  ifelse(Gamma==0&Neta==1,0.2,0.3)))]


# png(filename = "d:/quinonesa/Dropbox/Neuchatel/Figs/Sarsa/FigSup1.png",
    # width = 1600,height = 800)

par(plt=posPlot(numplotx = 2,idplotx = 1),yaxt='s',las=1)
with(FIAIntstats,{
  plotCI(x=Interv+posit,y=meanProb,
         ui = upIQR,li=lowIQR,
         pch=16,xlab='',ylab='',
         col=colboxes[ifelse(Gamma==0.8,
                             ifelse(Neta==1,1,2),
                             ifelse(Neta==1,3,4))],
         sfrac=0.002,cex.axis=1.3,ylim=c(0,1))
  mtext(side = 2,text = "frequency of V over R",las=0,cex=2,line=4)
  lines(x=c(0,max(Interv)),y=c(0.5,0.5),col='grey')
  text(x=par('usr')[2]*0.1,y=par('usr')[3]+0.9*(par('usr')[4]-par('usr')[3]),
       labels='FAA',cex=1.5)
})

# with(DPdataProb,  
#      {matlines(x = t(matrix(rep(max(FIAtimeInt$Interv)*c(0.75,1),
#                                 each=length(RV.V)),length(RV.V))),
#                y=t(matrix(rep(probRV.V,2),length(RV.V))),
#                lwd=2,lty = "dashed",
#                col=colboxes[ifelse(Gamma==0.8,
#                                    ifelse(Neta==1,1,2),
#                                    ifelse(Neta==1,3,4))])})

par(plt=posPlot(numplotx = 2,idplotx = 2),
    new=TRUE,yaxt='s',xpd=TRUE)

with(PIAIntstats,{
  plotCI(x=Interv+posit,y=meanProb,
         ui = upIQR,li=lowIQR,
         pch=16,xlab='',ylab='',
         col=colboxes[ifelse(Gamma==0.8,
                             ifelse(Neta==1,1,2),
                             ifelse(Neta==1,3,4))],
         sfrac=0.002,cex.axis=1.3,yaxt='n',ylim=c(0.4,1))
  lines(x=c(0,max(Interv)),y=c(0.5,0.5),col='grey')
  text(x=par('usr')[2]*0.1,y=par('usr')[3]+0.9*(par('usr')[4]-par('usr')[3]),
       labels='PAA',cex=1.5)
})


legend('topright',
       legend=c("neg. reward + future", "future",
                "neg. reward","no neg. reward + no future"),
       col=colboxes,pch=15,cex=1.5,ncol=1)

dev.off()

# Plot the dynamics of VR choice -----------------------------------------------------------



par(plt=posPlot(numplotx = 2,idplotx = 1)+c(-0.05,-0.05,0,0),yaxt='s',las=1)
with(FIAIntstats[(Neta==0)],{
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
yaxRangy<-c('s','n')

ylabRang<-c("Prob. V over R","")

png(filename = "d:/quinonesa/Dropbox/Neuchatel/Olle/Sarsa_Olle_par.png",
    width = 1600,height = 800)
count<-0
for(tau in unique(FIAIntstats$Tau)){

  count<-count+1
  par(plt=posPlot(numplotx = 2,idplotx = count),yaxt=yaxRangy[count],las=1,
      new=c(FALSE,TRUE)[count])
  with(FIAIntstats[(Neta==0.5&factRew==2)&Tau==tau],{
    plotCI(x=Interv,y=meanProb,
           ui = upIQR,li=lowIQR,
           pch=16,xlab='Time',ylab=ylabRang[count],cex.lab=2,
           col=colOlle[match(Gamma,unique(Gamma))],
           sfrac=0.002,cex.axis=1.3,ylim=c(0,1),cex=cexpar)
    lines(x=c(0,max(Interv)),y=c(0.5,0.5),col='grey')
  })
  par(new=TRUE)
  with(FIAIntstats[(Neta==0&factRew==1)&Tau==tau],{
    plotCI(x=Interv,y=meanProb,
           ui = upIQR,li=lowIQR,
           pch=16,xlab='',ylab='',
           col=colOlle2[match(Gamma,unique(Gamma))],
           sfrac=0.002,cex.axis=1.3,ylim=c(0,1),cex=cexpar)
    lines(x=c(0,max(Interv)),y=c(0.5,0.5),col='grey')
  })
  
  with(DPdataProb[(Neta==0.5&factRew==2)&Tau==tau],  
       {matlines(x = t(matrix(rep(max(FIAtimeInt$Interv)*c(0.75,1),
                                  each=length(RV.V)),length(RV.V))),
                 y=t(matrix(rep(probRV.V,2),length(RV.V))),
                 lwd=2,lty = "dashed",cex=1.2,
                 col=colOlle[match(Gamma,unique(Gamma))])})
  
  with(DPdataProb[(Neta==0&factRew==1)&Tau==tau],  
       {matlines(x = t(matrix(rep(max(FIAtimeInt$Interv)*c(0.75,1),
                                  each=length(RV.V)),length(RV.V))),
                 y=t(matrix(rep(probRV.V,2),length(RV.V))),
                 lwd=2,lty = "dashed",cex=1.2,
                 col=colOlle2[match(Gamma,unique(Gamma))])})
  text(x=par('usr')[1]+0.15*(par('usr')[2]-par('usr')[1]),
       y=par('usr')[3]+0.1*(par('usr')[4]-par('usr')[3]),
       labels = bquote(tau==.(tau)),cex = 2)
}
legend('bottomright',
       legend=c("Punishment and future", "punishment",
                "future","no punishment no future"),
       col=c(colOlle,colOlle2),pch=15,cex=1.3,ncol=1)

dev.off()




  
  