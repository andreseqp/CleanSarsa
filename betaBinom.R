library(data.table)
library(TailRank)
library(ggplot2)
n.dim<-4
probs<-c(0.2,0.5,0.8)
pvar.binom<-data.table(value=c(TailRank::rbb(n = 1000,N = n.dim,u = 1,v = 1),
                               rbinom(n = 1000,size = n.dim,prob = probs[2]),
                               rbinom(n = 1000,size = n.dim,prob = probs[3])),
                       p.val=as.factor(rep(probs,each=1000)))

ggplot(data = pvar.binom,aes(x=value,fill=p.val))+
  geom_histogram(position = "dodge")+
  theme_classic()


us<-c(0.1,1,10)
pvar.betabinom<-data.table(value=c(TailRank::rbb(n = 10000,N = n.dim,u = us[1],v = 1),
                               TailRank::rbb(n = 10000,N = n.dim,u = us[2],v = 1),
                               TailRank::rbb(n = 10000,N = n.dim,u = us[3],v = 1)),
                       u=as.factor(rep(us,each=10000)))

ggplot(data = pvar.betabinom,aes(x=value,fill=u))+
  geom_histogram(position = "dodge")+
  theme_classic()
