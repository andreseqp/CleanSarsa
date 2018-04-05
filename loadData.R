# ------------------------ Load data ------------------------------------------#

library("data.table")
library("jsonlite")
library("rlist")


getFilelist<-# reads de list of files and filters it according to a list of parameters
             # and values of interest
  function(folder, # folder where the files are
           listparam=NULL, # list strings providing the parameters 
                           #of interest
           values=NULL # list of values matching the list in 
                      # listparam
           ){
  posAgen<-c("PIA","FIA","DP")
  listRaw<-list.files(folder,recursive = TRUE)
  fullList<-vector("list",3)
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
  irrelPar<-c("gamma","tau","neta")
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
      by=.(Interv=floor(Age/interV),Training,Alpha,Gamma,Tau,Neta,Outbr)]
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


soft_max<-function(x,y,t){
  return(exp(x/t)/(exp(x/t)+exp(y/t)))
}

diffJsons<-function(json1,json2){
  print("JSON.1")
  print(unlist(json1)[unlist(json1)!=unlist(json2)])
  print("JSON.2")
  print(unlist(json2)[unlist(json1)!=unlist(json2)])
}


