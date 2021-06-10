# ------------------------ generate json files  ------------------------------------------ #

library("jsonlite")
library(here)


exedir<-paste(projDir,'/./Sarsa.exe',sep='')

fileName<-"parameters.json"


#test<-fromJSON(paste(codedir,"\\test.json",sep=""))

scenario<-"AbundLvpTauRang"

param<-list(totRounds=20000,ResReward=1,VisReward=1,
            ResProb=c(0.2),
            VisProb=c(0.2),
            ResProbLeav=0,VisProbLeav=1,negativeRew=-0.5,experiment=FALSE,
            inbr=0,outbr=0,trainingRep=30,forRat=0.0,numlearn = 1,
            propfullPrint = 0.7,
            alphaT=0.01,printGen=10000,seed=1, gammaRange=c(0,0.8),
            tauRange=c(0.5,1,2),netaRange=c(0,1),alphaThRange=c(0.01),
            folderL=paste0(here("Simulations"),scenario,"_/"))

clustfolder="/hpcfs/home/a.quinones/CleanSarsa/AbundLvp_/"

clustfolderNeu<-paste0("/home/ubuntu/SARSA/",scenario,"_/")


check_create.dir<-function(folder,param,values){
  setwd(folder)
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

rangLeav<-seq(0,0.3,length.out = 9)
rangAbund<-seq(0.1,0.9,length=9)

check_create.dir(here("Simulations"),param = rep(scenario,1),
                 values = c(""))

listfolders<-check_create.dir(here("Simulations",paste0(scenario,"_")),
                                    param = rep("Vlp",length(rangLeav)),
                              values = rangLeav)

for (i in 1:length(rangLeav)) {
  for(j in 1:length(rangAbund)){
    # if(i==11) param$printGen<-1
    # else param$printGen<-1000
    param$folderL<-paste0(here("Simulations",paste0(scenario,"_"),listfolders[i]),"/")
    param$folder<-paste0(clustfolderNeu,listfolders[i],"/")
    param$ResProb<-c((1-rangAbund[j])/2)
    param$VisProb<-c((1-rangAbund[j])/2)
    param$VisProbLeav<-rangLeav[i]
    outParam<-toJSON(param,auto_unbox = TRUE,pretty = TRUE)
    fileName<-paste("parameters",j,".json",sep="")
  if(file.exists(paste(param$folder,fileName,sep = ''))){
    currFile<-fromJSON(paste(param$folderL,fileName,sep = ''))
    if(sum(unlist(currFile)!=unlist(param))>0){
      # warning("You are erasing old files!! n\ Check first!!!",immediate. = TRUE)
      # print("OLD value")
      # print(unlist(currFile)[unlist(currFile)!=unlist(param)])
      # print("NEW value")
      # print(unlist(param)[unlist(currFile)!=unlist(param)])
      # ans<-readline("Want to continue?")
      # if(substr(ans, 1, 1) == "y"){
        write(outParam,paste(param$folderL,fileName,sep = ""))
        # jobfile(param$folderL,listfolders[i],jobid = j)
      # }
    }
  }
  else{
    write(outParam,paste(param$folderL,fileName,sep = ""))
    # jobfile(param$folderL,listfolders[i],jobid = j)
  }
  # system(paste(exedir,
  #   gsub("\\","/",paste(simsDir,listfolders[i],fileName,sep="\\"),fixed=TRUE)
  #   ,sep = " "))
  }
}
gsub(pattern = "\\",replacement = "/",simsdir,fixed=TRUE)

# system(paste(exedir,
#              gsub("\\","/",paste(simsdir,listfolders[1],fileName,sep="\\"),fixed=TRUE)
#              ,sep = " "))


## Automatically produce job files ---------------------------------------------

jobfile<-function(folder,jobname,timelim="10:00:00",
                  part="short",jobid=""){
  bashafile<-list(line0="#!/bin/bash",
                  jobname="#SBATCH --job-name=",
                  partit="#SBATCH -p ",nodes="#SBATCH -N 1",
                  cpus="#SBATCH --cpus-per-task=1", mem="#SBATCH --mem=2000",
                  time="#SBATCH --time=",
                  mailu="#SBATCH --mail-user=a.quinones@uniandes.edu.co",
                  mailt="#SBATCH --mail-type=END",
                  outp= paste0("#SBATCH -o ","/hpcfs/home/a.quinones/CleanSarsa/AbundLvp_/",
                               jobname,"/TEST_job.o%j"),
                  gethost="host=`/bin/hostname`",
                  getdate="date=/bin/date",
                  module="module load gcc/4.9.4",
                  exec=paste0("/hpcfs/home/a.quinones/CleanSarsa/./cleaner ",
                              "/hpcfs/home/a.quinones/CleanSarsa/AbundLvp_/",
                              jobname,"/parameters",jobid,".json"),
                  printhost="echo \"Run  at: \"$host",
                  printdate="echo  \"Run  on: \"$date")
  
  bashafile$jobname<-paste0(bashafile$jobname,jobname)
  bashafile$time<-paste0(bashafile$time,timelim)
  bashafile$partit<-paste0(bashafile$partit,part)
  if(file.exists(paste0(folder,"jobfile",jobid,".sh"))){
    unlink(paste0(folder,"jobfile",jobid,".sh"))
  }
  conJobfile<-lapply(bashafile, write, 
                     file=file(paste0(folder,"jobfile",jobid,".sh"),"wb"), append=T)
}


