# ------------------------ Load data ------------------------------------------#

library("data.table")
library("jsonlite")
library("rlist")


getFilelist<-# reads de list of files and filters it according to a list of parameters
             # and values of interest
  function(folder, # folder where the files are
           listparam=NULL, # list strings providing the parameters 
                           #of interest
           values=NULL, # list of values matching the list in 
                      # listparam
           fullNam=FALSE # Boolean. Append full path name
           ){
  posAgen<-c("PIA","FIA","DP","p1")
  listRaw<-list.files(folder,recursive = TRUE,full.names = fullNam)
  fullList<-vector("list",4)
  names(fullList)<-posAgen
  if(length(listparam)!=length(values)){
    # parameter and value lists must be the same length
    warning("Parameter list and values don't match",immediate. = TRUE)
  }
  else{
      if(is.null(listparam)){
        # if there is no list, use all the data
        paramList<-listRaw
  }
    else{
      # if there is list, filter the data
      #regExpList<-paste0(listparam,values,"_",sep="")
      for (param in unique(listparam)){
        valsparam<-values[grep(param,listparam)]
        listRaw<-do.call(list.append,lapply(paste(param,valsparam,"_",sep=""),
                                           grep,x=listRaw,value=TRUE))
      }
      
      paramList<-listRaw
    }
      for(agent in posAgen){
        # Create a list with the lists separated by type of agents
        listAgent<-grep(agent,paramList,value = TRUE)
        fullList[[agent]]<-listAgent
      }
      return(fullList)
    }
}




loadRawData<-function(folder,agent,listparam,values){
  setwd(folder)
  fullList<-getFilelist(folder,listparam,values)
  DT<-do.call(rbind,lapply(fullList[[agent]],fread))
  DT$option<-ifelse((DT$Type_choice==1 & DT$Type_discard==0) | 
                   (DT$Type_choice==0 & DT$Type_discard==1),"RV",NA)
  DT$option<-ifelse((DT$Type_choice==0 & DT$Type_discard==0),"RR",DT$option)
  DT$option<-ifelse((DT$Type_choice==1 & DT$Type_discard==1),"VV",DT$option)
  DT$option<-ifelse((DT$Type_choice==0 & DT$Type_discard==2) | 
                   (DT$Type_choice==2 & DT$Type_discard==0),"R0",DT$option)
  DT$option<-ifelse((DT$Type_choice==1 & DT$Type_discard==2) | 
                   (DT$Type_choice==2 & DT$Type_discard==1),"V0",DT$option)
  DT$option<-ifelse((DT$Type_choice==2 & DT$Type_discard==2),"00",DT$option)
  
  return(DT)
}


getParam<-function(folder,agent,listparam=NULL,values=NULL){
  setwd(folder)
  irrelPar<-c("gamma","tau","neta","pR","pV")
  listRaw<-list.files(folder,recursive = TRUE)
  jsonsList<-grep(".json",listRaw,value = TRUE)
  indRelPar<-seq(length(listparam))
  for(param in irrelPar){
    indRelPar<-grep(param,listparam,invert = TRUE)
    listparam<-listparam[indRelPar]
    values<-values[indRelPar]
  }
  if(length(listparam)!=length(values)){
    warning("Parameter list and values don't match",immediate. = TRUE)
  }
  else{
    if(is.null(listparam)){
      finalList<-jsonsList
    }
    else{
      for (param in unique(listparam)){
        valsparam<-values[grep(param,listparam)]
        jsonsList<-do.call(list.append,lapply(paste(param,valsparam,"_",sep=""),
                                              grep,x=jsonsList,value=TRUE))
      }
    }
  }
  jsons<-do.call(list,lapply(jsonsList,fromJSON))
  return(jsons)
}

file2timeInter<-function(filename,interV,maxAge=-2){
  extPar<-strsplit(filename,split ="_/")[[1]][1]
  parVal<-as.numeric(gsub("[[:alpha:]]",extPar,replacement = ''))
  extPar<-gsub("[[:digit:]]",extPar,replacement = '')
  tmp<-fread(filename,nrows = maxAge+1)
  tmp$fullRVoptions<-(tmp$Client1==1& tmp$Client2==0) | 
    (tmp$Client1==0 & tmp$Client2==1)
  tmptimeInter<-
    tmp[fullRVoptions==TRUE,.(Prob.RV.V=mean(Choice)),
      by=.(Interv=floor(Age/interV),Training,Alpha,Gamma,Tau,Neta)]
  if(length(extPar)>0){
    tmptimeInter[,eval(extPar):=parVal]
  }
  return(tmptimeInter)
}

file2lastDP<-function(filename){
  extPar<-strsplit(filename,split ="_/")[[1]][1]
  parVal<-as.numeric(gsub("[[:alpha:]]",extPar,replacement = ''))
  extPar<-gsub("[[:digit:]]",extPar,replacement = '')
  tmp<-fread(filename)
  tmpProbsDP<-tmp[Time==max(Time),
                  .(probRV.V=soft_max(RV.V,RV.R,Tau),RV.V,RV.R),
                  by=.(Alpha,Gamma,Tau,Neta,Outbr)]
  if(length(extPar)>0){
    tmpProbsDP[,eval(extPar):=parVal]
  }
  return(tmpProbsDP)
}

file2lastDP<-function(filename)
{
  extPar<-strsplit(filename,split ="_/")[[1]][1]
  parVal<-as.numeric(gsub("[[:alpha:]]",extPar,replacement = ''))
  extPar<-gsub("[[:digit:]]",extPar,replacement = '')
  tmp<-fread(filename)
  tmpProbsDP<-tmp[Time==max(Time),
                  .(probRV.V=soft_max(RV.V,RV.R,Tau),RV.V,RV.R),
                  by=.(Alpha,Gamma,Tau,Neta,Outbr)]
  if(length(extPar)>0){
    tmpProbsDP[,eval(extPar):=parVal]
  }
  return(tmpProbsDP)
}

file2lastProp<-function(filename,prop,outPar=NULL,genfold=NULL,full.path=FALSE) {
  if(length(outPar)>0){
    extPar<-grep(outPar,strsplit(filename,"_/")[[1]],
                 value=TRUE)
    parVal<-as.numeric(gsub("[[:alpha:]]",extPar,replacement = ''))
    extPar<-gsub("[[:digit:]]",extPar,replacement = '')
  }
  if(is.null(genfold)) 
    if(full.path)   tmp<-fread(filename)
    else tmp<-fread(here("Simulations",filename))
  else tmp<-fread(here("Simulations",genfold,filename))
  tmp[,':='(pV=as.numeric(gsub("[[:alpha:]]",
                               grep("pV",
                                    strsplit(filename,"_")[[1]],
                                    value=TRUE),replacement = "")),
            pR=as.numeric(gsub("[[:alpha:]]",
                               grep("pR",
                                    strsplit(filename,"_")[[1]],
                                    value=TRUE),replacement = ""))
  )]
  tmp$fullRVoptions<-(tmp$Client1==1& tmp$Client2==0) |
    (tmp$Client1==0 & tmp$Client2==1)
  tmp<-tmp[Age>max(Age)*prop]
  lastchunk<-tmp[fullRVoptions==TRUE,
                 list(Prob.RV.V=mean(Choice)),
                 by=.(Training,Gamma,Neta,pR,pV,Outbr)]
  if(length(outPar)>0){
    lastchunk[,eval(outPar):=parVal]
  }
  return(lastchunk)
}

soft_max<-function(x,y,t){
  return(exp(x/t)/(exp(x/t)+exp(y/t)))
}

diffJsons<-function(json1,json2){
  print("JSON.1")
  print(unlist(json1)[unlist(json1)!=unlist(json2)])
  print("JSON.2")
  print(unlist(json2)[unlist(json1)!=unlist(json2)])
}

loadDataFirstReach<-function(filename,bound){
  tmp<-fread(filename)
  tmp$fullRVoptions<-(tmp$Client1==1& tmp$Client2==0) | (tmp$Client1==0 & tmp$Client2==1)
  tmp[,':='(ReachedCut=soft_max(x = RV.V,y = RV.R,t=Tau)>bound,
            pV=as.numeric(gsub("[[:alpha:]]",
                               grep("pV",
                                    strsplit(filename,"_")[[1]],
                                    value=TRUE),replacement = "")),
            pR=as.numeric(gsub("[[:alpha:]]",
                               grep("pR",
                                    strsplit(filename,"_")[[1]],
                                    value=TRUE),replacement = ""))
      )]
  
  tmpFirstR<-tmp[ReachedCut==TRUE,.(firstReach=min(Age),
                                    Prob.RV.V=soft_max(
                                      x = RV.V[Age==min(Age)],
                                      y = RV.R[Age==min(Age)],
                                      t = Tau[Age==min(Age)])),
                 by=.(Training,Gamma,Neta,pR,pV,Outbr)]
  return(tmpFirstR)
}


# Density fucntion for the beta-binomial distribution
dbb <- function(x, N, u, v, log = FALSE) {
  logval <- lbeta(x+u, N-x+v) - lbeta(u,v) + lchoose(N,x)
  if (log) {
    ret <- logval
  } else {
    ret <- exp(logval)
  }
  ret
}

## Random numbers from the beta binomial distribution
rbb <- function(n, N, u, v) {
  p <- rbeta(n, u, v)
  rbinom(n, N, p)
}


## Function to check if new directories already exist---------------------------
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
      here()
      return(listfolders)
    }
    else{
      here()
      return(listfolders)
    }
  }else{
    lapply(listfolders,dir.create)
    here()
    return(listfolders)
  }
}




