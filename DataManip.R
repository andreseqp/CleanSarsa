# Data manipulation -----------------------------------------------------------------

# Directories --------------------------------------------------------------
simsDir<-"s:/quinonesa/Simulations/functionAprox/General/"
projDir<-"d:/quinonesa/learning_models_c++/functionAprox/"


# libraries ----------------------------------------------------------------
source('d:/quinonesa/Dropbox/R_files/posPlots.R')
source(paste(projDir,"loadData.R",sep = ""))


# Check if directories exist if they don't create them

list.dirs(simsDir)

list.dirs(projDir)

parList<-rep("outb_",3)
valList<-c(0,0.1,0.2)


check_create.dir(scriptDir,parList,valList)

listfolders<-list.dirs(simsDir,full.names = TRUE)

listfolders<-grep(parList[1],listfolders,value=TRUE)

lapply(as.list(listfolders), function(x){
  setwd(listfolders)
  fwrite(x =do.call(rbind,lapply(
               getFilelist(x)$FIA,
               file2timeInter,interV=1001)),
         file = paste(scriptDir,strsplit(x,split = '/')[[1]][7],
                      "FIA.timeInt.1001.txt",sep="/")
         ,sep = '\t')
})



lapply(as.list(listfolders), function(x){
  setwd(listfolders)
  fwrite(x =do.call(rbind,lapply(
    getFilelist(x)$PIA,
    file2timeInter,interV=1001)),
    file = paste(scriptDir,strsplit(x,split = '//')[[1]][2],
                 "PIA.timeInt.1001.txt",sep="/")
    ,sep = '\t')
})

lapply(as.list(listfolders), function(x){
  setwd(listfolders)
  fwrite(x =do.call(rbind,lapply(
    getFilelist(x)$DP,
    file2lastDP)),
    file = paste(scriptDir,strsplit(x,split = '//')[[1]][2],
                 "DP.last.txt",sep="/")
    ,sep = '\t')
})





