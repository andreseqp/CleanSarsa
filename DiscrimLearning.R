## Testing the model of discrimination learning 
library(here)
library(data.table)
library(ggplot2)
library(dplyr)
library("cowplot")
source(here("AccFunct.R"))

# load test data
test<-fread(list.files(here("Simulations","test_"),full.names = TRUE)[2])

interv<-1001

simsDir<-here("Simulations","test_")

# Load interval data for FIA from the raw data
PAAtimeInt<-do.call(
  rbind,lapply(
    getFilelist(simsDir,fullNam = TRUE)$PIA[1],
    file2timeInter,interV=interv))

PAAIntstats<-PAAtimeInt[,.(meanProb=mean(Prob.RV.V),
                           upIQR=fivenum(Prob.RV.V)[4],
                           lowIQR=fivenum(Prob.RV.V)[2])
                        ,by=.(Interv,Neta,Tau,Gamma)]

test[,`:=`(Val.vis=Val.02+Val.11,Val.res=Val.01+Val.12)]

Values.long<-melt(test,id.vars = c("Age","Training"),measure.vars = patterns("Val.*"),
                  variable.name = "Features")

ggplot(data = Values.long[Age%%50==0],aes(x=Age,y=value,color=Features))+
  stat_summary(fun.data = function(x) {
    xFive<-fivenum(x)
    return(data.frame(y=xFive[3],ymax=xFive[4],ymin=xFive[2]))})+
  # geom_path()+
  geom_hline(yintercept = c(1,2),color="black")+
  # facet_wrap(~Training)+
  theme_classic()

Alphas.long<-melt(test,id.vars = c("Age","Training"),measure.vars = patterns("alpha.*"),
                  variable.name = "Stimulus")

ggplot(data = Alphas.long[Age%%500==0],aes(x=Age,y=value,color=Stimulus))+
  stat_summary(fun.data = function(x) {
    xFive<-fivenum(x)
    return(data.frame(y=xFive[3],ymax=xFive[4],ymin=xFive[2]))})+
  # geom_step()+
  # facet_wrap(~Training)+
  theme_classic()


plot_grid(
  ggplot(data=test,aes(x=as.factor(Client1),fill=as.factor(Stim1.0)))+
    geom_bar(position = "stack")+
    theme_classic() ,
  ggplot(data=test,aes(x=as.factor(Client1),fill=as.factor(Stim1.1)))+
    geom_bar(position = "stack")+
    theme_classic()
)


plot_grid(
  ggplot(data=test,aes(x=as.factor(Client2),fill=as.factor(Stim2.0)))+
    geom_bar(position = "stack")+
    theme_classic() ,
  ggplot(data=test,aes(x=as.factor(Client2),fill=as.factor(Stim2.1)))+
    geom_bar(position = "stack")+
    theme_classic()
)

ggplot(PAAtimeInt,aes(x=Interv,y=Prob.RV.V))+
  stat_summary(fun.data = function(x) {
      xFive<-fivenum(x)
      return(data.frame(y=xFive[3],ymax=xFive[4],ymin=xFive[2]))
})+ylim(0,1)+
  geom_hline(yintercept = soft_max(1,2,0.5))+
  theme_classic()
