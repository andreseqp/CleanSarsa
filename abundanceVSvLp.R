######################## Abundance analysis ######################################

# Load libraries and external functions -----------------------------------

library(here)

source(here("aesth_par.R"))
source(here("loadData.R"))
source(here("..","..","R_files","posPlots.R"))
source(here("..","..","R_files","vioplot.R"))
source(here("..","..","R_files",'ternaryAEQP.R'))
source(here("data2interp.R"))
library('plotrix')
library('akima')
library("vcd")

# Delete data reset environment ---------------------------------------------------------

rm(list = ls())

rm(list=ls()[grepl('data.frame', sapply(ls(), function(x) class(get(x))))])


# Load Data FIA --------------------------------------------------------------------------------------


(listPar<-rep("Vlp",10))

(listVal<-seq(0,0.9,0.1))



# param<-getParam(here("Simulations","AbundLvp"),listparam = listPar,values = listVal)

FIAfirstReach<-do.call(rbind,lapply(
  getFilelist(here("Simulations","AbundLvp"),listPar,listVal)$FIA,
                                    loadDataFirstReach,0.7))

FIAfrStats<-FIAfirstReach[,.(firstReach=mean(firstReach),
                             Prob.RV.V=mean(Prob.RV.V)),
                          by=.(Neta,Gamma,pR,pV,Outbr)]

FIAfrStats$notProb<-round(1-FIAfrStats$pR-FIAfrStats$pV,1)



FIAlastQuarData<-do.call(rbind,lapply(
  getFilelist(here("Simulations","AbundLvp_"),
              listPar,listVal)$p1,file2lastProp,0.75,outPar="Vlp",
  genfold="AbundLvp_"))

FIAlastQuarData[,pA:=1-pR-pV]

FIA.stats<-FIAlastQuarData[,.(meanProb=mean(Prob.RV.V),
                                  upIQR=fivenum(Prob.RV.V)[4],
                                  lowIQR=fivenum(Prob.RV.V)[2])
                               ,by=.(Neta,Gamma,pR,pV,Vlp)]

FIA.stats$pA<-round(1-FIA.stats$pR-FIA.stats$pV,1)



# Data interpolations ---------------------------------------------------------

FIAinterpData<-AbundLeavData2interp(FIAlastQuarData[Neta==0&Gamma==0.8],
                                Var2int = "Prob.RV.V",npoints = 1000)

FIAinterpData.Neg<-AbundLeavData2interp(FIAlastQuarData[Neta==1&Gamma==0],
                                    Var2int = "Prob.RV.V",npoints = 1000)


# Plot real data prob ----------------------------------------------------------

png(here("Simulations","AbundLvp_","levelPlotReal.png"),
    width = 1000 , height = 700)

par(yaxt='s')
plot.new()
par(plt=c(posPlot(numplotx = 5,idplotx = 1)[1],
          posPlot(numplotx = 5,idplotx = 2)[2:4]))
with(FIA.stats[Neta==0&Gamma==0.8],{
  plot((1-pA),Vlp,col = paletteMeans(100)[findInterval(meanProb,seq(min(meanProb),
                                                              max(meanProb),
                                                              length=100))],
            main="",cex=0.8,pch=20,xlab="Cleaner abundance",
       ylab="Visitor leaving probability");
  # color.bar.aeqp(paletteMeans(100),min =round(min(meanProb),2),
  #           max = round(max(meanProb),2),nticks = 5,numplotx = 5,numploty = 3,
  #           idplotx = 4,idploty = 3)
})

par(plt=c(posPlot(numplotx = 5,idplotx = 3)[1],
          posPlot(numplotx = 5,idplotx = 4)[2:4]),new=TRUE)
with(FIA.stats[Neta==1&Gamma==0],{
  plot((1-pA),Vlp,col = paletteMeans(100)[findInterval(meanProb,seq(min(meanProb),
                                                                    max(meanProb),
                                                                    length=100))],
       main="",cex=0.8,pch=20,xlab="Cleaner abundance",
       ylab="",yaxt="n");
  color.bar.aeqp(paletteMeans(100),min =round(min(meanProb),2),
                 max = round(max(meanProb),2),nticks = 5,numplotx = 10,numploty = 1,
                 idplotx = 9,idploty = 1,locAxis = 4)
})

dev.off()

# png(paste(dirfig,"triplex_tau10_neta0_Outbr0.png",sep=""),width=1000,height=1000)

# Plot Interpolated data -------------------------------------------------------

png(here("Simulations","AbundLvp_","levelPlotInterp.png"),
    width = 1000 , height = 700)

colorbreaksMeans<-seq(0.45,1,length=100)

plot.new()
par(plt=c(posPlot(numplotx = 5,idplotx = 1)[1],
          posPlot(numplotx = 5,idplotx = 2)[2:4]))
with(FIAinterpData,{
  plot((1-pAbs),VLeavProb,
              col = paletteMeans(100)[findInterval(Prob.RV.V,
                                                   seq(min(Prob.RV.V),
                                                       max(Prob.RV.V),
                                                       length=100))],
              main="",cex=1,cex.lab = 2,pch=20,
       xlab="Cleaner abundance",
       ylab="Visitor leaving probability")
#   color.bar.aeqp(paletteMeans(100),min =round(min(Prob.RV.V),2),
#                  max = round(max(Prob.RV.V),2),nticks = 5,
#                  title = "Probability of \n V over R",cex.tit = 2,numplotx = 5,
#                  numploty = 5,idplotx = 5,idploty = 4)
})
# dev.off()

# png(paste("d:/quinonesa/Dropbox/Neuchatel/Results/actCrit/",listPar[1],"_",
#           listVal[1],"punish.png",sep=""),
#     width = 700 , height = 700)
par(plt=c(posPlot(numplotx = 5,idplotx = 3)[1],
          posPlot(numplotx = 5,idplotx = 4)[2:4]),new=TRUE)
with(FIAinterpData.Neg,{
  plot((1-pAbs),VLeavProb,
              col = paletteMeans(100)[findInterval(Prob.RV.V,
                                                   seq(min(Prob.RV.V),
                                                       max(Prob.RV.V),
                                                       length=100))],
              main="",cex=1,cex.lab = 2,pch=20,xlab="Cleaner abundance",
       ylab="",yaxt="n");
  
  color.bar.aeqp(paletteMeans(100),min =min(colorbreaksMeans),
                 max = max(colorbreaksMeans),nticks = 3,
                 title = "",cex.tit = 1,numplotx = 10,
                 numploty = 1,idplotx = 9,idploty = 1,locAxis = 4)
})

dev.off()

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


png(paste("d:/quinonesa/Dropbox/Neuchatel/Figs/Actor_critic/",
          "triplex_pL",listVal[1],".png",sep=""),
    width=1000,height=1000)


cex.lab.par<-1.8

colorbreaksMeans<-seq(0.40,1,length=100)
  # seq(min(c(FIAinterpData$Prob.RV.V,
  #                               FIAinterpData.Neg$Prob.RV.V)),
  #                max(c(FIAinterpData$Prob.RV.V,
  #                      FIAinterpData.Neg$Prob.RV.V)),length=100)

colorbreaksSpeed<-seq(800,20000,length=100)
  # seq(min(c(FIAinterpDataSpeed$firstReach,
  #                               FIAinterpDataSpeed.Neg$firstReach),
  #                         na.rm = TRUE),
  #                     max(c(FIAinterpDataSpeed$firstReach,
  #                           FIAinterpDataSpeed.Neg$firstReach),
  #                         na.rm=TRUE),length=100)

# colorbreaksSpeed<-seq(min(c(FIAinterpDataFRTrans$FirstReach),
#                           na.rm = TRUE),
#                       max(c(FIAinterpDataFRTrans$FirstReach),
#                           na.rm=TRUE),length=100)

plot.new()
with(FIAinterpData,{
  ternaryplotAEQP(cbind(resProb,visProb,notProb),
                  col = c(paletteMeans(100)[
                    findInterval(Prob.RV.V,colorbreaksMeans)],
                    rep("black",3)),main="",cex=0.4,
                  dimnames = c("Resident","Visitor","Absence")[c(2,3,1)],
                  dimnames_position = "edge",
                  border = "white",labels = "outside",labels_rot = c(0,0,0),
                  cex.lab = cex.lab.par,cex.grid = 1,grid_color = "black",
                  labels_color = "black",
                  numplotx = 2,numploty = 2,idplotx = 1,idploty = 1)
  
  
  ## code to add points to the triplex
  # rbind(cbind(resProb,visProb,notProb),
  #       cbind(unique(FIAagg$pR)[1:3],rep(unique(FIAagg$pV),3)
  #             ,1-unique(FIAagg$pR)[1:3]-rep(unique(FIAagg$pV),3)))
})

# text(x = par('usr')[4]-0.9*(par('usr')[4]-par('usr')[3]),
#      y = par('usr')[2]-0.9*(par('usr')[2]-par('usr')[1]),labels = 'A')

with(na.omit(FIAinterpDataSpeed),{
  ternaryplotAEQP(cbind(resProb,visProb,notProb),
                  col = c(paletteVar(100)[findInterval(firstReach,
                                                     colorbreaksSpeed)],
                          rep("black",3)),
                  main="",cex=0.4,dimnames = c("Resident",
                                               "Visitor","Absence")[c(2,3,1)],
                  border = "white",labels = "outside",labels_rot = c(0,0,0),
                  dimnames_position = "edge",
                  cex.lab = cex.lab.par,cex.grid = 1,grid_color = "black",
                  labels_color = "black",newpage = FALSE,
                  numplotx = 2,numploty = 2,idplotx = 1,idploty = 2);
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
                  main="",dimnames = c("Resident","Visitor","Absence")[c(2,3,1)],
                  dimnames_position = "edge",
                  border = "white",labels = "outside",labels_rot = c(0,0,0),
                  cex.lab = cex.lab.par,cex.grid = 1,grid_color = "black",
                  labels_color = "black",cex=0.4,newpage = FALSE,
                  numplotx = 2,numploty = 2,idplotx = 2,idploty = 1);
  # color.bar(rgb.palette(100),min =round(min(meanProb),2),max = round(max(meanProb),2),nticks = 5,
  #           title = "Probability of \n V over R",cex.tit = 2)
})
# text(x = par('usr')[4]-0.9*(par('usr')[4]-par('usr')[3]),
#      y = par('usr')[2]-0.9*(par('usr')[2]-par('usr')[1]),labels = 'A')

with(FIAinterpDataSpeed.Neg,{
  ternaryplotAEQP(cbind(resProb,visProb,notProb),
                  col = paletteVar(100)[findInterval(firstReach,colorbreaksSpeed)],
                  main="",dimnames = c("Resident","Visitor","Absence")[c(2,3,1)],
                  dimnames_position = "edge",
                  border = "white",labels = "outside",labels_rot = c(0,0,0),
                  cex.lab = cex.lab.par,cex.grid = 1,grid_color = "black",
                  labels_color = "black",cex=0.5,
                  numplotx = 2,numploty = 2,idplotx = 2,idploty = 2,
                  newpage = FALSE);
  # color.bar(rgb.palette(100),min =round(min(IQR),2),max = round(max(IQR),2),nticks = 3,
  #           title = "Probability of \n V over R",cex.tit = 2)
})
# text(x = par('usr')[4]-0.9*(par('usr')[4]-par('usr')[3]),
#      y = par('usr')[2]-0.9*(par('usr')[2]-par('usr')[1]),labels = 'B')

with(rbind(FIAinterpData,FIAinterpData.Neg),{
  par(new=FALSE)
  color.bar.aeqp(paletteMeans(100),min =min(colorbreaksMeans),
                 max = max(colorbreaksMeans),nticks = 3,
            title = "Probability  of Visitors \n chosen over residents",
            cex.tit = 1.2,
            numplotx = 15,numploty = 8,idplotx =7,idploty = 7)})
# using max and min in the scale
# with(rbind(FIAinterpData,FIAinterpData.Neg),{
#   par(new=TRUE)
#   color.bar.aeqp(paletteMeans(100),min =round(min(Prob.RV.V),2),
#                  max = round(max(Prob.RV.V),2),nticks = 3,
#                  title = "Probability \n of V over R",cex.tit = 0.8,
#                  numplotx = 15,numploty = 8,idplotx =7,idploty = 7)})

with(rbind(FIAinterpDataSpeed,FIAinterpDataSpeed.Neg),{
  color.bar.aeqp(paletteVar(100),min = round(min(colorbreaksSpeed)/1001),
                 max = round(max(colorbreaksSpeed)/1001),nticks = 3,
            title = "Time to 0.7",cex.tit = 0.8,
            numplotx = 15,numploty = 8,idplotx =7,idploty = 2)})

# using max and min in the scale
# with(rbind(FIAinterpDataSpeed,FIAinterpDataSpeed.Neg),{
#   color.bar.aeqp(paletteVar(100),min =round(min(firstReach/1001,na.rm = TRUE),2),
#                  max = round(max(firstReach/1001,na.rm = TRUE),2),nticks = 3,
#                  title = "Time to 0.7",cex.tit = 0.8,
#                  numplotx = 15,numploty = 8,idplotx =7,idploty = 2)})



dev.off()



# Load Data PIA--------------------------------------------------------------------------------------


setwd(simsDir)

(listPar<-c("AbundanceLpr"))

(listVal<-c(1))

PIALastQuarData<-do.call(rbind,lapply(getFilelist(simsDir,listPar,listVal)$PIA
                                      ,file2lastProp,0.75))

PIA.stats<-PIALastQuarData[,.(meanProb=mean(Prob.RV.V),
                                  upIQR=fivenum(Prob.RV.V)[4],
                                  lowIQR=fivenum(Prob.RV.V)[2])
                               ,by=.(Neta,Gamma,pR,pV,Outbr)]

PIA.stats$notProb<-round(1-PIA.stats$pR-PIA.stats$pV,1)

PIAfirstReach<-do.call(rbind,lapply(
  getFilelist(simsDir,listPar,listVal)$PIA,
  loadDataFirstReach,0.7))

PIAfrStats<-PIAfirstReach[,.(firstReach=mean(firstReach),
                             Prob.RV.V=mean(Prob.RV.V)),
                          by=.(Neta,Gamma,pR,pV,Outbr)]

PIAfrStats$notProb<-round(1-PIAfrStats$pR-PIAfrStats$pV,1)

# Interpolate data PIA -------------------------------------------------------

PIAinterpData<-AbundData2interp(PIALastQuarData[Gamma==0.8&Neta==0],
                                Var2int = "Prob.RV.V")
  
PIAinterpData.Neg<-AbundData2interp(PIALastQuarData[Gamma==0&Neta==1],
                                    Var2int = "Prob.RV.V")

PIAinterpDataSpeed<-AbundData2interp(PIAfirstReach[Gamma==0.8&Neta==0],
                                Var2int = "firstReach")


PIAinterpDataSpeed.Neg<-AbundData2interp(PIAfirstReach[Gamma==0&Neta==1],
                                    Var2int = "firstReach")

# Plot real data PIA -----------------------------------------------------------

plot.new()
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

par(yaxt="s")
with(PIAfrStats[Neta==0&Gamma==0.8],{
  ternaryplot(cbind(pR,pV,notProb),
              col = paletteVar(100)[findInterval(firstReach,seq(min(firstReach),
                                                                max(firstReach),
                                                                length=100))],
              main="",cex=0.8);
  color.bar.aeqp(paletteVar(100),min =round(min(firstReach),2),
                 max = round(max(firstReach),2),
                 nticks = 4,numplotx = 15,numploty = 5,idplotx = 12,idploty = 3)
})


par(yaxt="s")
with(PIA.stats[Neta==1&Gamma==0],{
  ternaryplot(cbind(pR,pV,notProb),
              col = paletteMeans(100)[findInterval(meanProb,seq(min(meanProb),
                                                                max(meanProb),
                                                                length=100))],
              main="",cex=0.8);
  color.bar.aeqp(paletteMeans(100),min =round(min(meanProb),2),
                 max = round(max(meanProb),2),
                 nticks = 4,numplotx = 15,numploty = 5,idplotx = 12,idploty = 3)
})

par(yaxt="s")
with(PIAfrStats[Neta==1&Gamma==0],{
  ternaryplot(cbind(pR,pV,notProb),
              col = paletteVar(100)[findInterval(firstReach,seq(min(firstReach),
                                                                max(firstReach),
                                                                length=100))],
              main="",cex=0.8);
  color.bar.aeqp(paletteVar(100),min =round(min(firstReach),2),
                 max = round(max(firstReach),2),
                 nticks = 4,numplotx = 15,numploty = 5,idplotx = 12,idploty = 3)
})


# Plot interpolated data PIA --------------------------------------------------

# png(paste(dirfig,"triplex_tau10_neta0_Outbr0.png",sep=""),width=1000,height=1000)

png(paste("d:/quinonesa/Dropbox/Neuchatel/Figs/Actor_critic/",
          "PIA_triplex_InitVal0.png",sep=""),
    width=1000,height=1000)

cex.lab.par<-1.2

colorbreaksMeansPIA<-seq(min(c(PIAinterpData$Prob.RV.V,
                            PIAinterpData.Neg$Prob.RV.V)),
                      max(c(PIAinterpData$Prob.RV.V,
                            PIAinterpData.Neg$Prob.RV.V)),length=100)

# colorbreaksSpeedPIA<-seq(min(c(PIAfrStats$firstReach,
#                                PIAinterpDataSpeed.Neg$firstReach),na.rm = TRUE),
#                          max(c(PIAfrStats$firstReach,
#                                PIAinterpDataSpeed.Neg$firstReach),na.rm = TRUE)
#                          ,length=100)


plot.new()
with(PIAinterpData,{
  ternaryplotAEQP(cbind(resProb,visProb,notProb),
                  col = paletteMeans(100)[findInterval(Prob.RV.V,
                                                       colorbreaksMeansPIA)],
                  main="",cex=1,dimnames = c("Resident","Visitor","Absence"),
                  border = "white",labels = "outside",labels_rot = c(0,0,0),
                  cex.lab = cex.lab.par,cex.grid = 2,
                  numplotx = 2,numploty = 2,idplotx = 1,idploty = 1);
  # color.bar.aeqp(paletteMeans(100),
  #                min =round(min(meanProb),2),max = round(max(meanProb),2),
  #                nticks = 5,title = "Probability of \n V over R",
  #                cex.tit = 2,numplotx = 15,numploty = 5,
  #                idplotx = 12,idploty = 3)
})

# par(new=TRUE)
# with(na.omit(PIAfrStats[Neta==0&Gamma==0.8]),{
#   ternaryplotAEQP(cbind(pR,pV,notProb),
#                   col = c(paletteVar(100)[findInterval(firstReach,
#                                                        colorbreaksSpeedPIA)],
#                           rep("black",3)),
#                   main="",cex=0.4,dimnames = c("Resident","Visitor","Absence"),
#                   border = "white",labels = "outside",labels_rot = c(0,0,0),
#                   cex.lab = cex.lab.par,cex.grid = 1,
#                   numplotx = 2,numploty = 2,idplotx = 1,idploty = 2,new=FALSE);
# })  

par(new=TRUE)
with(PIAinterpData.Neg,{
  ternaryplotAEQP(cbind(resProb,visProb,notProb),
                  col = paletteMeans(100)[findInterval(Prob.RV.V,
                                                       colorbreaksMeansPIA)],
                  main="",cex=1,dimnames = c("Resident","Visitor","Absence"),
                  border = "white",labels = "outside",labels_rot = c(0,0,0),
                  cex.lab = cex.lab.par,cex.grid = 2,newpage=FALSE,
                  numplotx = 2,numploty = 2,idplotx = 2,idploty = 1);
  # color.bar.aeqp(paletteMeans(100),min =round(min(meanProb),2),max = round(max(meanProb),2),nticks = 5,
  #                title = "Probability of \n V over R",cex.tit = 2,
  #                numplotx = 15,numploty = 5,idplotx = 12,idploty = 3)
})

# with(PIAinterpDataSpeed.Neg,{
#   ternaryplotAEQP(cbind(resProb,visProb,notProb),
#                   col = paletteVar(100)[findInterval(firstReach,
#                                                      colorbreaksSpeedPIA)],
#                   main="",dimnames = c("Resident","Visitor","Absence"),
#                   border = "white",labels = "outside",labels_rot = c(0,0,0),
#                   cex.lab = cex.lab.par,cex.grid = 1,cex=0.5,
#                   numplotx = 2,numploty = 2,idplotx = 2,idploty = 2,
#                   newpage = FALSE);
#   # color.bar(rgb.palette(100),min =round(min(IQR),2),max = round(max(IQR),2),nticks = 3,
#   #           title = "Probability of \n V over R",cex.tit = 2)
# })


with(rbind(PIAinterpData,PIAinterpData.Neg),{
  par(new=TRUE)
  color.bar.aeqp(paletteMeans(100),min =round(min(Prob.RV.V),2),
                 max = round(max(Prob.RV.V),2),nticks = 3,
                 title = "Probability \n of V over R",cex.tit = 0.8,
                 numplotx = 15,numploty = 8,idplotx =7,idploty = 7)})

# with(PIAinterpDataSpeed.Neg,{
#   color.bar.aeqp(paletteVar(100),min =round(min(firstReach/1001,na.rm = TRUE),2),
#                  max = round(max(firstReach/1001,na.rm = TRUE),2),nticks = 3,
#                  title = "Time to 0.7",cex.tit = 0.8,
#                  numplotx = 15,numploty = 8,idplotx =7,idploty = 2)})

dev.off()

# Figure PIA ------------------------------------------------------------------

# png(paste(dirfig,"triplex_tau10_neta0_Outbr0.png",sep=""),width=1000,height=1000)

cex.lab.par<-1.2

png(paste(dirfig,"SuplFig1_4triplex_PIA_InitVal_var.png",sep=''),
    width = 1000,height = 1000)

colorbreaksMeans<-seq(min(cbind(PIAinterpData$Prob.RV.V,
                                PIAinterpData.Neg$Prob.RV.V)),
                      max(cbind(PIAinterpData$Prob.RV.V,
                                PIAinterpData.Neg$Prob.RV.V)),length=100)

# colorbreaksIQR<-seq(min(cbind(PIAinterpDataSpeed$,interpDataVarTrans.Neg$IQR)),
#                     max(cbind(interpDataVarTrans$IQR,interpDataVarTrans.Neg$IQR)),length=100)

plot.new()
with(PIAinterpData,{
  ternaryplotAEQP(cbind(resProb,visProb,notProb),
                  col = paletteMeans(100)[findInterval(Prob.RV.V,colorbreaksMeans)],
                  main="",cex=1,dimnames = c("Resident","Visitor","Absence"),border = "white",
                  labels = "outside",labels_rot = c(0,0,0),cex.lab = cex.lab.par,cex.grid = 1,
                  numplotx = 2,numploty = 2,idplotx = 1,idploty = 1)
})

# text(x = par('usr')[4]-0.9*(par('usr')[4]-par('usr')[3]),
#      y = par('usr')[2]-0.9*(par('usr')[2]-par('usr')[1]),labels = 'A')



with(PIAinterpData.Neg,{
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


############# From now on code that old code that does not work anymore ########

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





