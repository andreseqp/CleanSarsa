######################## Abundance analysis ######################################

# Load libraries and external functions -----------------------------------

projDir<-"d:/quinonesa/learning_models_c++/Sarsa/"
simsDir<-"s:/quinonesa/Simulations/Basic_sarsa/"


source('d:/quinonesa/Dropbox/R_files/posPlots.R')
source(paste(projDir,"aesth_par.R",sep=""))
source(paste(projDir,"loadData.R",sep = ""))
source('D:/quinonesa/Dropbox/R_files/posPlots.R')
source('D:/quinonesa/Dropbox/R_files/vioplot.R')
source('D:/quinonesa/Dropbox/R_files/ternaryAEQP.R')
source(paste(projDir,"data2interp.R",sep=""))
library('plotrix')
library('akima')
library("vcd")

# Delete data reset environment ---------------------------------------------------------

rm(list = ls())

rm(list=ls()[grepl('data.frame', sapply(ls(), function(x) class(get(x))))])


# Useful plotting tricks ------------------------------------------------------------------

axisRangy<-c('s','n','n')
axisRangx<-c('n','n','s')
adjyaxis<-c(rep(-0.04,2),rep(-0.02,2))
lettRang<-matrix(c('A','B','C','D','E','F','G','H','I'),ncol=3,byrow = TRUE)
titY<-matrix(c('','Probability of V over R','',rep('',6)),ncol=3)
titX<-matrix(c(rep('',7),'Learning trial',''),ncol=3,byrow = TRUE)


# Load Data FIA --------------------------------------------------------------------------------------

setwd(simsDir)

(listPar<-c("AbundanceLp"))

(listVal<-c(0))

param<-getParam(simsDir,listparam = listPar,values = listVal)

FIAfirstReach<-do.call(rbind,lapply(
  getFilelist(simsDir,listPar,listVal)$FIA,
                                    loadDataFirstReach,0.6))

FIAfrStats<-FIAfirstReach[,.(firstReach=mean(firstReach),
                             Prob.RV.V=mean(Prob.RV.V)),
                          by=.(Neta,Gamma,pR,pV,Outbr)]

FIAfrStats$notProb<-round(1-FIAfrStats$pR-FIAfrStats$pV,1)


FIAlastQuarData<-do.call(rbind,lapply(
  getFilelist(simsDir,listPar,listVal)$FIA,file2lastProp,0.75))


FIA.stats<-FIAlastQuarData[,.(meanProb=mean(Prob.RV.V),
                                  upIQR=fivenum(Prob.RV.V)[4],
                                  lowIQR=fivenum(Prob.RV.V)[2])
                               ,by=.(Neta,Gamma,pR,pV)]

FIA.stats$notProb<-round(1-FIA.stats$pR-FIA.stats$pV,1)



# Data interpolations ---------------------------------------------------------

FIAinterpData<-AbundData2interp(FIAlastQuarData[Neta==0&Gamma==0.8],
                                Var2int = "Prob.RV.V")

FIAinterpData.Neg<-AbundData2interp(FIAlastQuarData[Neta==1&Gamma==0],
                                    Var2int = "Prob.RV.V")

FIAinterpDataSpeed<-AbundData2interp(FIAfirstReach[Neta==0&Gamma==0.8],
                                     Var2int = "firstReach")

FIAinterpDataSpeed.Neg<-AbundData2interp(FIAfirstReach[Neta==1&Gamma==0],
                                     Var2int = "firstReach")

# Plot real data prob ----------------------------------------------------------

plot.new()
with(FIA.stats[Neta==0&Gamma==0.8],{
  ternaryplot(cbind(pR,pV,notProb),
            col = paletteMeans(100)[findInterval(meanProb,seq(min(meanProb),
                                                              max(meanProb),
                                                              length=100))],
            main="",cex=0.8);
  color.bar.aeqp(paletteMeans(100),min =round(min(meanProb),2),
            max = round(max(meanProb),2),nticks = 5,numplotx = 5,numploty = 3,
            idplotx = 4,idploty = 3)
})

with(FIA.stats[Neta==1&Gamma==0],{
  ternaryplot(cbind(pR,pV,notProb),
              col = paletteMeans(100)[findInterval(meanProb,seq(min(meanProb),
                                                                max(meanProb),
                                                                length=100))],
              main="",cex=0.8);
  color.bar.aeqp(paletteMeans(100),min =round(min(meanProb),2),
                 max = round(max(meanProb),2),nticks = 5,numplotx = 5,numploty = 3,
                 idplotx = 4,idploty = 3)
})

# png(paste(dirfig,"triplex_tau10_neta0_Outbr0.png",sep=""),width=1000,height=1000)

# Plot Interpolated data -------------------------------------------------------

with(FIAinterpData,{
  ternaryplotAEQP(cbind(resProb,visProb,notProb),
              col = paletteMeans(100)[findInterval(Prob.RV.V,
                                                   seq(min(Prob.RV.V),
                                                       max(Prob.RV.V),
                                                       length=100))],
              main="",cex=1,dimnames = c("Resident","Visitor","Absence"),
              border = "white", labels = "outside",labels_rot = c(0,0,0),
              cex.lab = 3.5,cex.grid = 2);
  color.bar.aeqp(paletteMeans(100),min =round(min(Prob.RV.V),2),
                 max = round(max(Prob.RV.V),2),nticks = 5,
                 title = "Probability of \n V over R",cex.tit = 2,numplotx = 5,
                 numploty = 5,idplotx = 5,idploty = 4)
})


with(FIAinterpData.Neg,{
  ternaryplotAEQP(cbind(resProb,visProb,notProb),
                  col = paletteMeans(100)[findInterval(Prob.RV.V,
                                                       seq(min(Prob.RV.V),
                                                           max(Prob.RV.V),
                                                           length=100))],
                  main="",cex=1,dimnames = c("Resident","Visitor","Absence"),
                  border = "white", labels = "outside",labels_rot = c(0,0,0),
                  cex.lab = 3.5,cex.grid = 2);
  color.bar.aeqp(paletteMeans(100),min =round(min(Prob.RV.V),2),
                 max = round(max(Prob.RV.V),2),nticks = 5,
                 title = "Probability of \n V over R",cex.tit = 2,numplotx = 5,
                 numploty = 5,idplotx = 5,idploty = 4)
})

# variation panel -----------------------------------------------------------------------

paletteVar <- colorRampPalette(c('#d8b365','#f5f5f5','#5ab4ac'),alpha=TRUE)

npoints<-100

interpDataVar<-with(stateLastQuarData[Neta==0&Gamma==0.8],
                    {interp(x=pR,y=pV,
                            z=Prob.RV.V,
                            duplicate = "user",
                            dupfun = {function(x) fivenum(x)[4]-fivenum(x)[2]},
                                           nx=npoints,ny=npoints)})
bound <-0.8

state.stats<-stateLastQuarData[,.(meanProb=mean(Prob.RV.V),
                                  upIQR=fivenum(Prob.RV.V)[4],
                                  lowIQR=fivenum(Prob.RV.V)[2])
                               ,by=.(Neta,Gamma,resProb,visProb,Outbr)]



str(interpDataVar)

interpDataVarTrans<-data.table(matrix(0,nrow = npoints*npoints,ncol = 4))

names(interpDataVarTrans)<-c("resProb","visProb","IQR","notProb")
for (i in 1:npoints) {
  for (j in 1:npoints) {
    interpDataVarTrans[(i-1)*npoints+j,resProb:=interpDataVar$x[i]]
    interpDataVarTrans[(i-1)*npoints+j,visProb:=interpDataVar$y[j]]
    interpDataVarTrans[(i-1)*npoints+j,IQR:=interpDataVar$z[i,j]]
  }
  
}

interpDataVarTrans[,4]<-1-interpDataVarTrans[,1]-interpDataVarTrans[,2]

interpDataVarTrans<-interpDataVarTrans[resProb+visProb<0.9]

interpDataVar.Neg<-with(stateLastQuarData[Neta==1&Gamma==0],
                        {interp(x=pR,y=pV,z=Prob.RV.V,
                                duplicate = "user",
                                dupfun = 
                                  {function(x) fivenum(x)[4]-fivenum(x)[2]},
                                                       nx=npoints,ny=npoints)})

interpDataVarTrans.Neg<-data.table(matrix(0,nrow = npoints*npoints,ncol = 4))

names(interpDataVarTrans.Neg)<-c("resProb","visProb","IQR","notProb")
for (i in 1:npoints) {
  for (j in 1:npoints) {
    interpDataVarTrans.Neg[(i-1)*npoints+j,resProb:=interpDataVar.Neg$x[i]]
    interpDataVarTrans.Neg[(i-1)*npoints+j,visProb:=interpDataVar.Neg$y[j]]
    interpDataVarTrans.Neg[(i-1)*npoints+j,IQR:=interpDataVar.Neg$z[i,j]]
  }
  
}

interpDataVarTrans.Neg[,4]<-1-interpDataVarTrans.Neg[,1]-interpDataVarTrans.Neg[,2]

interpDataVarTrans.Neg<-interpDataVarTrans.Neg[resProb+visProb<0.9]

with(state.stats[Neta==1&Gamma==0],{
  ternaryplot(cbind(pR,pV,notProb),
              col = paletteVar(100)[findInterval(upIQR-lowIQR,
                                                 seq(min(upIQR-lowIQR),
                                                     max(upIQR-lowIQR),
                                                     length=100))],
              main="",cex=0.8);
  color.bar.aeqp(paletteVar(100),min =round(min(upIQR-lowIQR),3),
            max = round(max(upIQR-lowIQR),3),
            nticks = 3,numplotx = 15,numploty = 8,idplotx =8,idploty = 8)
})


# Plot  real Speed data -----------------------------------------------------------------------


with(FIAfrStats[Gamma==0.8&Neta==0],{
  ternaryplot(cbind(pR,pV,notProb),
              col = paletteVar(100)[findInterval(firstReach,
                                                 seq(min(firstReach),
                                                     max(firstReach),
                                                     length=100))],
              main="",cex=0.8);
  color.bar.aeqp(paletteVar(100),min =round(min(firstReach/501),3),
                 max = round(max(firstReach/501),3),
                 nticks = 3,numplotx = 15,numploty = 8,idplotx =8,idploty = 8)
})

with(FIAfrStats[Gamma==0&Neta==1],{
  ternaryplot(cbind(pR,pV,notProb),
              col = paletteVar(100)[findInterval(firstReach,
                                                 seq(min(firstReach),
                                                     max(firstReach),
                                                     length=100))],
              main="",cex=0.8);
  color.bar.aeqp(paletteVar(100),min =round(min(firstReach/501),3),
                 max = round(max(firstReach/501),3),
                 nticks = 3,numplotx = 15,numploty = 8,idplotx =8,idploty = 8)
})

# Plot interpolated speed data ---------------------------------------------------------

with(na.omit(FIAinterpDataSpeed),{
  ternaryplot(cbind(resProb,visProb,notProb),
              col = paletteVar(100)[findInterval(firstReach,
                                                 seq(min(firstReach),
                                                     max(firstReach),
                                                     length=100))],
              main="",cex=0.8);
  color.bar.aeqp(paletteVar(100),min =round(min(firstReach/501),3),
                 max = round(max(firstReach/501),3),
                 nticks = 3,numplotx = 15,numploty = 8,idplotx =8,idploty = 8)
})

with(na.omit(FIAinterpDataSpeed.Neg),{
  ternaryplot(cbind(resProb,visProb,notProb),
              col = paletteVar(100)[findInterval(firstReach,
                                                 seq(min(firstReach),
                                                     max(firstReach),
                                                     length=100))],
              main="",cex=0.8);
  color.bar.aeqp(paletteVar(100),min =round(min(firstReach/501),3),
                 max = round(max(firstReach/501),3),
                 nticks = 3,numplotx = 15,numploty = 8,idplotx =8,idploty = 8)
})

# Figure ------------------------------------------------------------------------


png(paste("d:/quinonesa/Dropbox/Neuchatel/Figs/Sarsa/",
          "triplex_pLV0.png",sep=""),
    width=1000,height=1000)

cex.lab.par<-1.2

colorbreaksMeans<-seq(min(c(FIAinterpData$Prob.RV.V,
                                FIAinterpData.Neg$Prob.RV.V)),
                 max(c(FIAinterpData$Prob.RV.V,
                       FIAinterpData.Neg$Prob.RV.V)),length=100)

colorbreaksSpeed<-seq(min(c(FIAinterpDataSpeed$firstReach,
                                FIAinterpDataSpeed.Neg$firstReach),
                          na.rm = TRUE),
                      max(c(FIAinterpDataSpeed$firstReach,
                            FIAinterpDataSpeed.Neg$firstReach),
                          na.rm=TRUE),length=100)

colorbreaksSpeed<-seq(min(FIAinterpDataSpeed$firstReach,
                          na.rm = TRUE),
                      max(FIAinterpDataSpeed$firstReach,
                          na.rm=TRUE),length=100)

plot.new()
with(FIAinterpData,{
  ternaryplotAEQP(rbind(cbind(resProb,visProb,notProb)
                        # ,cbind(unique(FIAagg$pR)[1:3],rep(unique(FIAagg$pV),3)
                        #       ,1-unique(FIAagg$pR)[1:3]-rep(unique(FIAagg$pV),3))
                        ),
                  col = c(paletteMeans(100)[
                    findInterval(Prob.RV.V,colorbreaksMeans)],
                    rep("black",3)),main="",cex=0.4,
                  dimnames = c("Resident","Visitor","Absence"),
                  border = "white",labels = "outside",labels_rot = c(0,0,0),
                  cex.lab = cex.lab.par,cex.grid = 1,
                  numplotx = 2,numploty = 2,idplotx = 1,idploty = 1)
})

# text(x = par('usr')[4]-0.9*(par('usr')[4]-par('usr')[3]),
#      y = par('usr')[2]-0.9*(par('usr')[2]-par('usr')[1]),labels = 'A')

with(na.omit(FIAinterpDataSpeed),{
  ternaryplotAEQP(rbind(cbind(resProb,visProb,notProb)
                        # ,cbind(unique(FIAagg$pR)[1:3],rep(unique(FIAagg$pV),3)
                        #       ,1-unique(FIAagg$pR)[1:3]-rep(unique(FIAagg$pV),3))
                        ),
                  col = c(paletteVar(100)[findInterval(firstReach,
                                                     colorbreaksSpeed)],
                          rep("black",3)),
                  main="",cex=0.4,dimnames = c("Resident","Visitor","Absence"),
                  border = "white",labels = "outside",labels_rot = c(0,0,0),
                  cex.lab = cex.lab.par,cex.grid = 1,
                  numplotx = 2,numploty = 2,idplotx = 1,idploty = 2,new=FALSE);
  # color.bar.aeqp(paletteVar(100),min =round(min(FirstReach),2),
  #           max = round(max(FirstReach),2),nticks = 3,
  #           title = "Probability of \n V over R",cex.tit = 2,numplotx = 10,
  #           numploty = 10,idplotx = 8,idploty = 8)
})


# text(x = par('usr')[4]-0.9*(par('usr')[4]-par('usr')[3]),
#      y = par('usr')[2]-0.9*(par('usr')[2]-par('usr')[1]),labels = 'B')

with(FIAinterpData.Neg,{
  ternaryplotAEQP(cbind(resProb,visProb,notProb),
                  col = paletteMeans(100)[findInterval(Prob.RV.V,
                                                       colorbreaksMeans)],
                  main="",dimnames = c("Resident","Visitor","Absence"),
                  border = "white",labels = "outside",labels_rot = c(0,0,0),
                  cex.lab = cex.lab.par,cex.grid = 1,cex=0.4,newpage = FALSE,
                  numplotx = 2,numploty = 2,idplotx = 2,idploty = 1);
  # color.bar(rgb.palette(100),min =round(min(meanProb),2),max = round(max(meanProb),2),nticks = 5,
  #           title = "Probability of \n V over R",cex.tit = 2)
})
# text(x = par('usr')[4]-0.9*(par('usr')[4]-par('usr')[3]),
#      y = par('usr')[2]-0.9*(par('usr')[2]-par('usr')[1]),labels = 'A')

with(FIAinterpDataSpeed.Neg,{
  ternaryplotAEQP(cbind(resProb,visProb,notProb),
                  col = paletteVar(100)[findInterval(firstReach,colorbreaksSpeed)],
                  main="",dimnames = c("Resident","Visitor","Absence"),
                  border = "white",labels = "outside",labels_rot = c(0,0,0),
                  cex.lab = cex.lab.par,cex.grid = 1,cex=0.5,
                  numplotx = 2,numploty = 2,idplotx = 2,idploty = 2,
                  newpage = FALSE);
  # color.bar(rgb.palette(100),min =round(min(IQR),2),max = round(max(IQR),2),nticks = 3,
  #           title = "Probability of \n V over R",cex.tit = 2)
})
# text(x = par('usr')[4]-0.9*(par('usr')[4]-par('usr')[3]),
#      y = par('usr')[2]-0.9*(par('usr')[2]-par('usr')[1]),labels = 'B')

with(rbind(FIAinterpData,FIAinterpData.Neg),{
  par(new=TRUE)
  color.bar.aeqp(paletteMeans(100),min =round(min(Prob.RV.V),2),
                 max = round(max(Prob.RV.V),2),nticks = 3,
            title = "Probability \n of V over R",cex.tit = 0.8,
            numplotx = 15,numploty = 8,idplotx =7,idploty = 7)})


with(rbind(FIAinterpDataSpeed),{ #,FIAinterpDataSpeed.Neg
  color.bar.aeqp(paletteVar(100),min =round(min(firstReach/1001,na.rm = TRUE),2),
                 max = round(max(firstReach/1001,na.rm = TRUE),2),nticks = 3,
            title = "Time to 0.7",cex.tit = 0.8,
            numplotx = 15,numploty = 8,idplotx =7,idploty = 2)})


dev.off()

# Load Data PIA--------------------------------------------------------------------------------------


setwd(simsDir)

PIALastQuarData<-do.call(rbind,lapply(getFilelist(simsDir,listPar,listVal)$PIA
                                      ,file2lastProp,0.6))

PIA.stats<-PIALastQuarData[,.(meanProb=mean(Prob.RV.V),
                                  upIQR=fivenum(Prob.RV.V)[4],
                                  lowIQR=fivenum(Prob.RV.V)[2])
                               ,by=.(Neta,Gamma,pR,pV,Outbr)]

PIA.stats$notProb<-round(1-PIA.stats$pR-PIA.stats$pV,1)

# Interpolate data PIA -------------------------------------------------------

PIAinterpData<-AbundData2interp(PIALastQuarData[Gamma==0.8&Neta==0],
                                Var2int = "Prob.RV.V")
  
PIAinterpData.Neg<-AbundData2interp(PIALastQuarData[Gamma==0&Neta==1],
                                    Var2int = "Prob.RV.V")

# Plot real data PIA -----------------------------------------------------------


par(yaxt="s")
with(PIA.stats[Neta==0&Gamma==0.8],{
  ternaryplot(cbind(pR,pV,notProb),
              col = paletteMeans(100)[findInterval(meanProb,seq(min(meanProb),
                                                                max(meanProb),
                                                                length=100))],
              main="",cex=0.8);
  color.bar.aeqp(paletteMeans(100),min =round(min(meanProb),2),
                 max = round(max(meanProb),2),
            nticks = 4,numplotx = 15,numploty = 5,idplotx = 12,idploty = 3)
})


# Plot interpolated data PIA --------------------------------------------------

# png(paste(dirfig,"triplex_tau10_neta0_Outbr0.png",sep=""),width=1000,height=1000)

png(paste("d:/quinonesa/Dropbox/Neuchatel/Figs/Actor_critic/",
          "PIA_triplex_Lp_0.png",sep=""),
    width=1000,height=1000)

cex.lab.par<-1.2

colorbreaksMeansPIA<-seq(min(c(PIAinterpData$Prob.RV.V,
                            PIAinterpData.Neg$Prob.RV.V)),
                      max(c(PIAinterpData$Prob.RV.V,
                            PIAinterpData.Neg$Prob.RV.V)),length=100)

plot.new()
with(PIAinterpData,{
  ternaryplotAEQP(cbind(resProb,visProb,notProb),
                  col = paletteMeans(100)[findInterval(Prob.RV.V,
                                                       colorbreaksMeansPIA)],
                  main="",cex=1,dimnames = c("Resident","Visitor","Absence"),
                  border = "white",labels = "outside",labels_rot = c(0,0,0),
                  cex.lab = 3.5,cex.grid = 2,
                  numplotx = 2,numploty = 2,idplotx = 1,idploty = 1);
  # color.bar.aeqp(paletteMeans(100),
  #                min =round(min(meanProb),2),max = round(max(meanProb),2),
  #                nticks = 5,title = "Probability of \n V over R",
  #                cex.tit = 2,numplotx = 15,numploty = 5,
  #                idplotx = 12,idploty = 3)
})

par(new=TRUE)
with(PIAinterpData.Neg,{
  ternaryplotAEQP(cbind(resProb,visProb,notProb),
                  col = paletteMeans(100)[findInterval(Prob.RV.V,
                                                       colorbreaksMeansPIA)],
                  main="",cex=1,dimnames = c("Resident","Visitor","Absence"),
                  border = "white",labels = "outside",labels_rot = c(0,0,0),
                  cex.lab = 3.5,cex.grid = 2,newpage=FALSE,
                  numplotx = 2,numploty = 2,idplotx = 2,idploty = 2);
  # color.bar.aeqp(paletteMeans(100),min =round(min(meanProb),2),max = round(max(meanProb),2),nticks = 5,
  #                title = "Probability of \n V over R",cex.tit = 2,
  #                numplotx = 15,numploty = 5,idplotx = 12,idploty = 3)
})


with(rbind(PIAinterpData,PIAinterpData.Neg),{
  par(new=TRUE)
  color.bar.aeqp(paletteMeans(100),min =round(min(Prob.RV.V),2),
                 max = round(max(Prob.RV.V),2),nticks = 3,
                 title = "Probability \n of V over R",cex.tit = 0.8,
                 numplotx = 15,numploty = 8,idplotx =13,idploty = 7)})

dev.off()

# Figure PIA ------------------------------------------------------------------

# png(paste(dirfig,"triplex_tau10_neta0_Outbr0.png",sep=""),width=1000,height=1000)

cex.lab.par<-1.2

png(paste(dirfig,"SuplFig1_4triplex_PIA.png",sep=''),width = 1000,height = 1000)

colorbreaksMeans<-seq(min(cbind(interpDataTrans$meanProb,interpDataTrans.Neg$meanProb)),
                      max(cbind(interpDataTrans$meanProb,interpDataTrans.Neg$meanProb)),length=100)

colorbreaksIQR<-seq(min(cbind(interpDataVarTrans$IQR,interpDataVarTrans.Neg$IQR)),
                    max(cbind(interpDataVarTrans$IQR,interpDataVarTrans.Neg$IQR)),length=100)

plot.new()
with(interpDataTrans,{
  ternaryplotAEQP(cbind(resProb,visProb,notProb),
                  col = paletteMeans(100)[findInterval(meanProb,colorbreaksMeans)],
                  main="",cex=1,dimnames = c("Resident","Visitor","Absence"),border = "white",
                  labels = "outside",labels_rot = c(0,0,0),cex.lab = cex.lab.par,cex.grid = 1,
                  numplotx = 2,numploty = 2,idplotx = 1,idploty = 1)
})

# text(x = par('usr')[4]-0.9*(par('usr')[4]-par('usr')[3]),
#      y = par('usr')[2]-0.9*(par('usr')[2]-par('usr')[1]),labels = 'A')



with(interpDataVarTrans,{
  ternaryplotAEQP(cbind(resProb,visProb,notProb),
                  col = paletteVar(100)[findInterval(IQR,colorbreaksIQR)],
                  main="",cex=1,dimnames = c("Resident","Visitor","Absence"),border = "white",
                  labels = "outside",labels_rot = c(0,0,0),cex.lab = cex.lab.par,cex.grid = 1,
                  numplotx = 2,numploty = 2,idplotx = 1,idploty = 2,newpage = FALSE);
  # color.bar(rgb.palette(100),min =round(min(IQR),2),max = round(max(IQR),2),nticks = 3,
  #           title = "Probability of \n V over R",cex.tit = 2)
})


# text(x = par('usr')[4]-0.9*(par('usr')[4]-par('usr')[3]),
#      y = par('usr')[2]-0.9*(par('usr')[2]-par('usr')[1]),labels = 'B')

with(interpDataTrans.Neg,{
  ternaryplotAEQP(cbind(resProb,visProb,notProb),
                  col = paletteMeans(100)[findInterval(meanProb,colorbreaksMeans)],
                  main="",cex=1,dimnames = c("Resident","Visitor","Absence"),border = "white",
                  labels = "outside",labels_rot = c(0,0,0),cex.lab = cex.lab.par,cex.grid = 1,
                  numplotx = 2,numploty = 2,idplotx = 2,idploty = 1,newpage = FALSE);
  # color.bar(rgb.palette(100),min =round(min(meanProb),2),max = round(max(meanProb),2),nticks = 5,
  #           title = "Probability of \n V over R",cex.tit = 2)
})
# text(x = par('usr')[4]-0.9*(par('usr')[4]-par('usr')[3]),
#      y = par('usr')[2]-0.9*(par('usr')[2]-par('usr')[1]),labels = 'A')

with(interpDataVarTrans.Neg,{
  ternaryplotAEQP(cbind(resProb,visProb,notProb),
                  col = paletteVar(100)[findInterval(IQR,colorbreaksIQR)],
                  main="",cex=1,dimnames = c("Resident","Visitor","Absence"),border = "white",
                  labels = "outside",labels_rot = c(0,0,0),cex.lab = cex.lab.par,cex.grid = 1,
                  numplotx = 2,numploty = 2,idplotx = 2,idploty = 2,newpage = FALSE);
  # color.bar(rgb.palette(100),min =round(min(IQR),2),max = round(max(IQR),2),nticks = 3,
  #           title = "Probability of \n V over R",cex.tit = 2)
})
# text(x = par('usr')[4]-0.9*(par('usr')[4]-par('usr')[3]),
#      y = par('usr')[2]-0.9*(par('usr')[2]-par('usr')[1]),labels = 'B')

with(rbind(interpDataTrans,interpDataTrans.Neg),{
  par(new=TRUE)
  color.bar(paletteMeans(100),min =round(min(meanProb),2),max = round(max(meanProb),2),nticks = 3,
            title = "Probability \n of V over R",cex.tit = 0.8,
            numplotx = 15,numploty = 8,idplotx =7,idploty = 7)})


with(rbind(interpDataVarTrans,interpDataVarTrans.Neg),{
  color.bar(paletteVar(100),min =round(min(IQR),2),max = round(max(IQR),2),nticks = 3,
            title = "IQR size",cex.tit = 0.8,
            numplotx = 15,numploty = 8,idplotx =7,idploty = 2)})


dev.off()


# Individual runs exploration FIA------------------------------------------------

listParRuns<-c("AbundanceLp","gamma","neta")
listValRuns<-c(0,0.8,0)

sum(as.numeric(gsub("[[:alpha:]]",
     strsplit(getFilelist(simsDir,listParRuns,listValRuns)$FIA[50],"_")[[1]][8:9],
replacement="")))

tmp3<-do.call(rbind,
              lapply(getFilelist(simsDir,listParRuns,listValRuns)$FIA, 
                  function(x){
                    if(sum(as.numeric(gsub("[[:alpha:]]",
                                           strsplit(x,"_")[[1]][7:8],
                                        replacement="")))==0.8){
                      fread(x)
                    }}))
        
tmp3<-do.call(rbind,lapply(
  getFilelist(simsDir,listParRuns,listValRuns)$FIA,
  function(x){
    if(sum(as.numeric(gsub("[[:alpha:]]",
                            strsplit(x,"_")[[1]][8:9],
                            replacement="")))==0.8){
      tmp<-fread(x)
      tmp[,':='(pV=as.numeric(gsub("[[:alpha:]]",
                                 grep("pV",
                                      strsplit(x,"_")[[1]],
                                      value=TRUE),replacement = "")),
              pR=as.numeric(gsub("[[:alpha:]]",
                                 grep("pR",
                                      strsplit(x,"_")[[1]],
                                      value=TRUE),replacement = "")))]
    return(tmp)
    }
    }
  ))

tmp3[,':='(Prob.RV.V=soft_max(RV.V,RV.R,Tau),notProb=1-pR-pV)]

tmp3agg<-tmp3[Age<10000, as.list(unlist(lapply(.SD, function(x) 
  list(mean = mean(x),IQ.h = fivenum(x)[4],IQ.l=fivenum(x)[2])))),
  by=.(Age,Alpha,Gamma,Tau,Neta,Outbr,pR,pV), 
  .SDcols=c('Prob.RV.V','RV.V','RV.R','VV.V','RR.R','R0.R','V0.V','OO.O')]

png(paste("d:/quinonesa/Dropbox/Neuchatel/Figs/Sarsa/",
          "Runs_visLeav_0_pA0.2.png",sep=""),
    width=1000,height=1000)

par(plt=posPlot(numplotx = 1,numploty = 2,idplotx = 1,idploty = 2),yaxt="s",
    xaxt="n")

pVlines<-dcast(tmp3agg,formula = Age~pV,value.var = "Prob.RV.V.mean")
names(pVlines)[2:dim(pVlines)[2]]<-paste0("pv_",
                                          names(pVlines)[2:dim(pVlines)[2]])
with(pVlines,{
  matplot(cbind(pv_0.3,pv_0.4,pv_0.5),type='l',
          ylab = "p(V)",lty=1,las=2)
})

legend("bottomright",legend=c(expression(p[v]==0.3),
                              expression(p[v]==0.4),
                              expression(p[v]==0.5)),
       col=1:5,pch=20,pt.cex = 3,
       ncol=1)
text(x=rep(80,2),y=c(0.65,0.6),labels = c("Sarsa-FIA",expression(gamma==0.8)))

par(plt=posPlot(numplotx = 3,numploty = 2,idplotx = 1,idploty = 1),
    new=TRUE,yaxt="s",xaxt="s")
with(tmp3agg[pV==0.3],{
  matplot(cbind(RV.V.mean,RV.V.mean,VV.V.mean,
                RR.R.mean,V0.V.mean,R0.R.mean),type='l',
          col=colorValues,lty=1,ylab="Values",las=2)
  text(x=15,y=4.5,labels = bquote(p[v]==.(pV)))
})
legend("bottomright",legend=c("RV.V","RV.R","VV.V","RR.R",
                "V0.V","R0.R"),col=colorValues,pch=19,
       ncol=2,cex=1)

par(plt=posPlot(numplotx = 3,numploty = 2,idplotx = 2,idploty = 1),
    new=TRUE,yaxt="n",xaxt="s")
with(tmp3agg[pV==0.4],{
  matplot(cbind(RV.V.mean,RV.V.mean,VV.V.mean,
                RR.R.mean,V0.V.mean,R0.R.mean),type='l',
          col=colorValues,ylab = "",lty=1)
  text(x=15,y=4.5,labels = bquote(p[v]==.(pV)))
})

par(plt=posPlot(numplotx = 3,numploty = 2,idplotx = 3,idploty = 1),
    new=TRUE,yaxt="n",xaxt="s")
with(tmp3agg[pV==0.5],{
  matplot(cbind(RV.V.mean,RV.V.mean,VV.V.mean,
                RR.R.mean,V0.V.mean,R0.R.mean),type='l',
          col=colorValues,ylab = "",lty=1)
  text(x=15,y=4.5,labels = bquote(p[v]==.(pV)))
})

dev.off()

par(plt=posPlot(numplotx = 3,idploty = 1,idplotx = 1,numploty = 2)
    ,xaxt="s",new=TRUE,yaxt='s')
with(tmp3agg[pR==0.],{
  plot(logist(ThetaV.mean,ThetaR.mean)~Age,type="l",ylim = c(0,1),
          ylab="Theta")
})

par(plt=posPlot(numplotx = 3,idploty = 1,idplotx = 2,numploty = 2)
    ,xaxt="s",new=TRUE,yaxt='n')
with(tmp3agg[pR==0.2],{
  plot(logist(ThetaV.mean,ThetaR.mean)~Age,type="l",ylab="",ylim = c(0,1))
})

par(plt=posPlot(numplotx = 3,idploty = 1,idplotx = 3,numploty = 2)
    ,xaxt="s",new=TRUE,yaxt='n')
with(tmp3agg[pR==0.1],{
  plot(logist(ThetaV.mean,ThetaR.mean)~Age,type="l",ylab="",ylim = c(0,1))
})


# Individual runs exploration PIA------------------------------------------------

sum(as.numeric(gsub("[[:alpha:]]",strsplit(getFilelist(simsDir,listPar,listVal)$PIA[38],
                             "_")[[1]][7:8],
                                replacement="")))

tmp3<-do.call(rbind,
              lapply(getFilelist(simsDir,listPar,listVal)$PIA, 
                     function(x){
                       if(round(sum(as.numeric(gsub("[[:alpha:]]",
                                              strsplit(x,"_")[[1]][7:8],
                                              replacement=""))),1)==0.9){
                         fread(x)
                       }}))



tmp3agg<-tmp3[, as.list(unlist(lapply(.SD, function(x) 
  list(mean = mean(x),IQ.h = fivenum(x)[4],IQ.l=fivenum(x)[2])))),
  by=.(Age,Alpha,Gamma,Tau,Neta,Outbr,pR,pV), 
  .SDcols=c('ThetaV','ThetaR','Resident','Visitor','Absence')]

unique(tmp3agg[,pR])

pRrange<-seq(0.3,0.1,length=3)

par(plt=posPlot(numplotx = 3,numploty = 2,idplotx = 1,idploty = 2),yaxt="s",
    xaxt="n")
with(tmp3agg[pR==pRrange[1]],{
  matplot(cbind(Resident.mean,Visitor.mean,Absence.mean),type='l',
          col=colorValues,ylab = "Value",lty=1)
})

legend("bottomright",legend=c("resident","visitor",
                              "absence"),
       col=colorValues,pch=20,pt.cex = 3,title="States values",
       ncol=1)

par(plt=posPlot(numplotx = 3,numploty = 2,idplotx = 2,idploty = 2),
    new=TRUE,yaxt="n",xaxt="n")
with(tmp3agg[pR==pRrange[2]],{
  matplot(cbind(Resident.mean,Visitor.mean,Absence.mean),type='l',
          col=colorValues,ylab = "",lty=1)
})

par(plt=posPlot(numplotx = 3,numploty = 2,idplotx = 3,idploty = 2),
    new=TRUE,yaxt="n",xaxt="n")
with(tmp3agg[pR==pRrange[3]],{
  matplot(cbind(Resident.mean,Visitor.mean,Absence.mean),type='l',
          col=colorValues,ylab = "",lty=1)
})


par(plt=posPlot(numplotx = 3,idploty = 1,idplotx = 1,numploty = 2)
    ,xaxt="s",new=TRUE,yaxt='s')
with(tmp3agg[pR==pRrange[1]],{
  plot(logist(ThetaV.mean,ThetaR.mean)~Age,type="l",ylim = c(0.3,0.7),
       ylab="Theta")
  lines(x = c(0,20000),y = c(0.5,0.5),col="grey")
})

par(plt=posPlot(numplotx = 3,idploty = 1,idplotx = 2,numploty = 2)
    ,xaxt="s",new=TRUE,yaxt='n')
with(tmp3agg[pR==pRrange[2]],{
  plot(logist(ThetaV.mean,ThetaR.mean)~Age,type="l",ylab="",ylim = c(0.3,0.7))
  lines(x = c(0,20000),y = c(0.5,0.5),col="grey")
})

par(plt=posPlot(numplotx = 3,idploty = 1,idplotx = 3,numploty = 2)
    ,xaxt="s",new=TRUE,yaxt='n')
with(tmp3agg[pR==pRrange[3]],{
  plot(logist(ThetaV.mean,ThetaR.mean)~Age,type="l",ylab="",ylim = c(0.3,0.7))
  lines(x = c(0,20000),y = c(0.5,0.5),col="grey")
})
