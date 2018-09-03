####### Cobined plot for figure 4 #############################################
# Include in one figure tripleces and dot plots with error bars 

# Load libraries and external functions -----------------------------------

projDir<-"d:/quinonesa/learning_models_c++/Sarsa//"
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

# Load Data for interpolations --------------------------------------------------------------------------------------

setwd(simsDir)


(listPar<-c("AbundanceB"))

(listVal<-c(""))

param<-getParam(simsDir,listparam = listPar,values = listVal)



FIAlastQuarData<-do.call(rbind,lapply(
  getFilelist(simsDir,listPar,listVal)$FIA,file2lastProp,0.75))


FIA.stats<-FIAlastQuarData[,.(meanProb=mean(Prob.RV.V),
                              upIQR=fivenum(Prob.RV.V)[4],
                              lowIQR=fivenum(Prob.RV.V)[2])
                           ,by=.(Neta,Gamma,pR,pV)]

FIA.stats$notProb<-round(1-FIA.stats$pR-FIA.stats$pV,1)


# Load data for dot plot -------------------------------------------------------


FIAlast<-do.call(rbind,lapply(getFilelist(simsDir,listPar,listVal)$FIA, 
                              function(x){
                                if(as.numeric(gsub("[[:alpha:]]",
                                                   strsplit(x,"_")[[1]][8],
                                                   replacement=""))==
                                   as.numeric(gsub("[[:alpha:]]",
                                                   strsplit(x,"_")[[1]][9],
                                                   replacement=""))){
                                  file2lastProp(x,0.9)
                                }}))




FIAlast[,pA:=1-pR-pV]

FIA.statsdot<-FIAlast[,.(meanProb=mean(Prob.RV.V),
                         upIQR=fivenum(Prob.RV.V)[4],
                         lowIQR=fivenum(Prob.RV.V)[2])
                      ,by=.(Neta,Gamma,pR,pV,Outbr,pA)]


# Data interpolations ---------------------------------------------------------

FIAinterpData<-AbundData2interp(FIAlastQuarData[Neta==0&Gamma==0.8],
                                Var2int = "Prob.RV.V")

FIAinterpData.Neg<-AbundData2interp(FIAlastQuarData[Neta==1&Gamma==0],
                                    Var2int = "Prob.RV.V")

# Panel a ---------------------------------------------------------------------

FIA.statsdot[,posit:=ifelse(Gamma==0&Neta==0,0,
                            ifelse(Gamma==0.8&Neta==0,0.01,
                                   ifelse(Gamma==0&Neta==1,0.02,0.03)))]

png("d:/quinonesa/Dropbox/Neuchatel/Figs/Sarsa/Fig4_panelA.png",
    width = 700 , height = 1200)

par(plt=posPlot(numplotx = 1,numploty = 1,1,1),las=1)
with(FIA.statsdot,{
  plotCI(x = (1-pA)+posit,
         y = meanProb,ui = upIQR
         ,li = lowIQR,
         col=colboxes[ifelse(Gamma==0.8,
                             ifelse(Neta==1,1,2),
                             ifelse(Neta==1,3,4))],
         pch=16,xlab="",ylab="",
         sfrac=0.008,yaxt='s',
         cex.axis=1.3,cex=2,cex.lab=3)
  mtext('Proportion of visitors chosen over residents',2,line = 4, cex=3,las=0)
  mtext('Overall client abundance',1,line = 5,las=1,cex=4)
  # text(x=par('usr')[1]+0.05*(par('usr')[2]-par('usr')[1]),
  #      y=par('usr')[3]+0.95*(par('usr')[4]-par('usr')[3]),
  #      labels='A',cex=2)
})

legend("topright",
       legend=c("neg. reward + future", "future",
                "neg. reward","no neg. reward + no future"),
       col=colboxes,pch=15,cex=1.5,ncol=1)

dev.off()

# Panel B - future triplex -----------------------------------------------------

png(paste("d:/quinonesa/Dropbox/Neuchatel/Figs/Sarsa/",
          "triplex_panelB.png",sep=""),
    width=500,height=500,units ="px")

cex.lab.par<-1.8

colorbreaksMeans<-seq(0.45,0.75,length=100)

colorbreaksMeans<-seq(min(c(FIAinterpData$Prob.RV.V,
                            FIAinterpData.Neg$Prob.RV.V)),
                      max(c(FIAinterpData$Prob.RV.V,
                            FIAinterpData.Neg$Prob.RV.V)),length=100)

plot.new()
with(FIAinterpData,{
  ternaryplotAEQP(rbind(cbind(resProb,visProb,notProb),
                        cbind(FIA.statsdot[Gamma==0.8&Neta==0]$pR,
                              FIA.statsdot[Gamma==0.8&Neta==0]$pV,
                              FIA.statsdot[Gamma==0.8&Neta==0]$pA)),
                  col = c(paletteMeans(100)[
                    findInterval(Prob.RV.V,colorbreaksMeans)],
                    rep("black",5)),main="",cex=0.4,
                  dimnames = c("Resident","Visitor","Absence")[c(2,3,1)],
                  dimnames_position = "edge",
                  border = "white",labels = "outside",labels_rot = c(0,0,0),
                  cex.lab = cex.lab.par,cex.grid = 1,grid_color = "black",
                  labels_color = "black",
                  numplotx = 1,numploty = 1,idplotx = 1,idploty = 1)
})
dev.off()

# Panel C - triplex negative reward --------------------------------------------

png(paste("d:/quinonesa/Dropbox/Neuchatel/Figs/Sarsa/",
          "triplex_panelC.png",sep=""),
    width=500,height=500,units ="px")

cex.lab.par<-1.8

plot.new()
with(FIAinterpData.Neg,{
  ternaryplotAEQP(rbind(cbind(resProb,visProb,notProb),
                        cbind(FIA.statsdot[Gamma==0&Neta==1]$pR,
                              FIA.statsdot[Gamma==0&Neta==1]$pV,
                              FIA.statsdot[Gamma==0&Neta==1]$pA)),
                  col = c(paletteMeans(100)[
                    findInterval(Prob.RV.V,colorbreaksMeans)],
                    rep("black",5)),main="",cex=0.4,
                  dimnames = c("Resident","Visitor","Absence")[c(2,3,1)],
                  dimnames_position = "edge",
                  border = "white",labels = "outside",labels_rot = c(0,0,0),
                  cex.lab = cex.lab.par,cex.grid = 1,grid_color = "black",
                  labels_color = "black",
                  numplotx = 1,numploty = 1,idplotx = 1,idploty = 1)
})
dev.off()

# ColorScale -------------------------------------------------------------------

png(paste("d:/quinonesa/Dropbox/Neuchatel/Figs/Sarsa/",
          "triplex_panel_colorSca.png",sep=""),
    width=400,height=400,units ="px")

plot.new()
with(rbind(FIAinterpData,FIAinterpData.Neg),{
  par(new=FALSE)
  color.bar.aeqp(paletteMeans(100),min =min(colorbreaksMeans),
                 max = max(colorbreaksMeans),nticks = 3,
                 title = "Probability  of Visitors \n chosen over residents",
                 cex.tit = 1.2,
                 numplotx = 2,numploty = 2,idplotx =2,idploty = 1)})
dev.off()
