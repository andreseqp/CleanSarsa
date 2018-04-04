# ------------------------ generate json files  ------------------------------------------ #

library("jsonlite")

projDir<-"D:\\quinonesa\\learning_models_c++\\Sarsa"

simsDir<-"S:/quinonesa/Simulations/Basic_sarsa"

exedir<-paste(projDir,'/./Sarsa.exe',sep='')

fileName<-"parameters.json"


#test<-fromJSON(paste(codedir,"\\test.json",sep=""))

param<-list(totRounds=30000,ResReward=10,VisReward=10,ResProb=0.2,VisProb=0.2,
            ResProbLeav=0,VisProbLeav=1,negativeRew=-10,experiment=FALSE,
            inbr=0,outbr=0,trainingRep=30,forRat=0.0,
            alphaT=0.01,printGen=1,seed=1, gammaRange=c(0,0.8),
            tauRange=c(5,10),netaRange=c(0,0.5),
            folder=simsDir)

param<-list(totRounds=10000,ResReward=1,VisReward=1,ResProb=0.3,VisProb=0.3,
            ResProbLeav=0,VisProbLeav=1,negativeRew=-0.5,experiment=FALSE,
            inbr=0,outbr=0,trainingRep=30,forRat=0.0,
            alphaT=0.01,printGen=1,seed=1, gammaRange=c(0,0.8),
            tauRange=c(0.6667),netaRange=c(0,0.5),
            folder=simsDir)

setwd(simsDir)


check_create.dir<-function(folder,param,values){
  listfolders<-paste(param,values,"_",sep = "")  
  currFolders<-lapply(listfolders,dir.exists)
  if(sum(currFolders>0)){
    warning("At least one of the folders already exists \n Please check",
            immediate. = TRUE)
    print(cbind(listfolders,currFolders))
    ans<-readline("Want to continue?")
    if(substr(ans, 1, 1) == "y"){
      lapply(listfolders,dir.create)
      return(listfolders)
    }
    else{
      return(listfolders)
    }
  }else{
    lapply(listfolders,dir.create)
    return(listfolders)
  }
}

rang<-c(1,2)

listfolders<-check_create.dir(simsDir,rep("factRew",2),rang)

for (i in 1:2) {
  param$folder<-paste(simsDir,'/',listfolders[i],'/',sep='')
  param$ResReward<-param$ResReward*rang[i]
  param$VisReward<-param$VisReward*rang[i]
  param$negativeRew<-param$negativeRew*rang[i]
  outParam<-toJSON(param,auto_unbox = TRUE,pretty = TRUE)
  if(file.exists(paste(param$folder,fileName,sep = ''))){
    currFile<-fromJSON(paste(param$folder,fileName,sep = ''))
    if(sum(unlist(currFile)!=unlist(param))>0){
      warning("You are erasing old files!! n\ Check first!!!",immediate. = TRUE)
      print("OLD value")
      print(unlist(currFile)[unlist(currFile)!=unlist(param)])
      print("NEW value")
      print(unlist(param)[unlist(currFile)!=unlist(param)])
      ans<-readline("Want to continue?")
      if(substr(ans, 1, 1) == "y"){
        write(outParam,paste(simsDir,listfolders[i],fileName,sep="\\"))
      }
    }
  }
  else{
    write(outParam,paste(simsDir,listfolders[i],fileName,sep="\\"))
  }
  # system(paste(exedir,
  #   gsub("\\","/",paste(simsdir,listfolders[i],fileName,sep="\\"),fixed=TRUE)
  #   ,sep = " "))
}
gsub(pattern = "\\",replacement = "/",simsdir,fixed=TRUE)

# system(paste(exedir,
#              gsub("\\","/",paste(simsdir,listfolders[1],fileName,sep="\\"),fixed=TRUE)
#              ,sep = " "))

