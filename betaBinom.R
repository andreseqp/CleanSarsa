library(data.table)
library(TailRank)
library(ggplot2)
library(dplyr)
source(here("AccFunct.R"))



dbb(c(0,1,2),2,1,1)

rbeta
rchisq


n.dim<-1
probs<-c(0.2,0.5,0.8)


pvar.binom<-data.table(value=c(rbb(n = 1000,N = n.dim,u = 1,v = 1),
                               rbinom(n = 1000,size = n.dim,prob = probs[2]),
                               rbinom(n = 1000,size = n.dim,prob = probs[3])),
                       p.val=as.factor(rep(probs,each=1000)))

ggplot(data = pvar.binom,aes(x=value,fill=p.val))+
  geom_histogram(position = "dodge")+
  theme_classic()


us<-c(0.01)
pvar.betabinom<-data.table(value=as.factor(c(rbb(n = 10000,N = n.dim,
                                                 u = us[1],v = 1),
                               rbb(n = 10000,N = n.dim,u = us[2],v = 1),
                               rbb(n = 10000,N = n.dim,u = us[3],v = 1))),
                       u=as.factor(rep(us,each=10000)))

ggplot(data = pvar.betabinom,aes(x=value,fill=u))+
  geom_bar(position = "dodge")+
  theme_classic()

test<-fread(list.files(here("Simulations","test_"),full.names = TRUE))

test[20000,.(Val.02+Val.11,Val.01+Val.12)]




n.dim<-1
pvar.betabinom.clients<-data.table(value=as.factor(c(rbb(n = 10000,N = n.dim,u = 1,v = 0.01),
                                             rbb(n = 10000,N = n.dim,u = 0.2,v = 1))),
                           client=as.factor(rep(c(0,1),each=10000)))
ggplot(data=pvar.betabinom.clients,aes(x=as.factor(client),fill=as.factor(value)))+
  geom_bar(position = "stack")+
  theme_classic()
  
