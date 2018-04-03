# --------------------------------- Time intervals -------------------------------------------------#

projDir<-"d:/quinonesa/learning_models_c++/functionAprox/"
simsDir<-"s:/quinonesa/Simulations/functionAprox/General"

# libraries ---------------------------------------------------------------------------------------
source('d:/quinonesa/Dropbox/R_files/posPlots.R')
source(paste(projDir,"aesth_par.R",sep=""))
source(paste(projDir,"loadData.R",sep = ""))
library('plotrix')

# Load Data ---------------------------------------------------------------------------------------


# Define data to be loaded 

(listPar<-c(rep("outb",3),"tau","gamma","neta"))
(listVal<-c(0,0.1,0.2,10,0.8,0))
param<-getParam(simsDir,listparam = listPar,values = listVal)

diffJsons(param[1],param[3])

list.files(simsDir,recursive = TRUE,pattern ="outb")

# Load interval data for FIA from the raw data
FIAtimeInt<-do.call(
  rbind,lapply(
    getFilelist(simsDir,listPar,listVal)$FIA[c(4,8,9)],
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

FIAIntstats<-FIAtimeInt[,.(meanProb=mean(Type_choice.mean),
                           upIQR=fivenum(Type_choice.mean)[4],
                           lowIQR=fivenum(Type_choice.mean)[2])
                        ,by=.(Interv,Neta,Outbr,Tau,Gamma,get(extpar))]
setnames(FIAIntstats,'get',extpar)
PIAIntstats<-PIAtimeInt[,.(meanProb=mean(Type_choice.mean),
                           upIQR=fivenum(Type_choice.mean)[4],
                           lowIQR=fivenum(Type_choice.mean)[2])
                        ,by=.(Interv,Neta,Outbr,Tau,Gamma,get(extpar))]
setnames(PIAIntstats,'get',extpar)

par(plt=posPlot(numplotx = 2,idplotx = 1)+c(-0.05,-0.05,0,0),yaxt='s',las=1)
with(FIAIntstats,{
  plotCI(x=Interv,y=meanProb,
         ui = upIQR,li=lowIQR,
         pch=16,xlab='',ylab='',
         col=colboxes[match(get(extpar),unique(get(extpar)))],
         sfrac=0.002,cex.axis=1.3)
  lines(x=c(0,max(Interv)),y=c(0.5,0.5),col='grey')
})
with(DPdataProb,  
     {matlines(x = t(matrix(rep(max(FIAtimeInt$Interv)*c(0.75,1),
                                each=length(RV.V)),length(RV.V))),
               y=t(matrix(rep(probRV.V,2),length(RV.V))),
               lwd=2,lty = "dashed",
               col=colboxes[match(get(extpar),unique(get(extpar)))])})

legend('bottomright',
       legend=unique(FIAIntstats[,get(extpar)])[order(unique(FIAIntstats[,get(extpar)]))],
              col=colboxes,pch=15,
              title=eval(extpar),cex=1.5,ncol=3)

par(plt=posPlot(numplotx = 2,idplotx = 2)+c(-0.05,-0.05,0,0),
    new=TRUE,yaxt='s',xpd=TRUE)
with(PIAIntstats,{
  plotCI(x=Interv,y=meanProb,
         ui = upIQR,li=lowIQR,
         pch=16,xlab='',ylab='',
         col=colboxes[match(get(extpar),unique(get(extpar)))],
         sfrac=0.002,cex.axis=1.3,yaxt='n')
  lines(x=c(0,max(Interv)),y=c(0.5,0.5),col='grey')
  axis(side=4,cex.axis=1.3)
})

png(filename = paste(projDir,eval(extpar),".png",sep=""))

dev.off()

param[1]


  
  