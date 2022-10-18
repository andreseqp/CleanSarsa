## Testing the model of discrimination learning 
library(here)
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library("cowplot")
source(here("AccFunct.R"))

scenario<-"PartialInfo2_2_"

# load test data
list_files<-list.files(here("Simulations",scenario),full.names = TRUE) %>%
  grep(pattern = "PIA",value = TRUE)

rawData<-do.call(rbind,lapply(list_files,FUN = function(file){
  tmpData<-fread(file)
  relPar<-strsplit(file,"_")[[1]] %>% grep(pattern = "AttMech",value = TRUE)
  parVal<-gsub("[^[:digit:]]", "", relPar)  %>%  as.numeric
  tmpData[,AttMech:=parVal]
  return(tmpData)
}))


interv<-501

simsDir<-here("Simulations",scenario)

# Load interval data for FIA from the raw data
PAAtimeInt<-do.call(
  rbind,lapply(
    getFilelist(simsDir,fullNam = TRUE)$PIA,
    file2timeInter,"AttMech",interV=interv))


rawData[,(grep("Val.*",names(rawData),value = TRUE)):=
          lapply(.SD,as.numeric),.SDcols=patterns("Val.*")]

rawData[,`:=`(Val.Clien1=fcase(Stim1.0==0,Val.00,
                               Stim1.0==1,Val.01,Stim1.0==2,Val.02)+
                fcase(Stim1.1==0,Val.10,
                      Stim1.1==1,Val.11,Stim1.1==2,Val.12),
              Val.Clien2=fcase(Stim2.0==0,Val.00,
                               Stim2.0==1,Val.01,Stim2.0==2,Val.02)+
                fcase(Stim2.1==0,Val.10,
                      Stim2.1==1,Val.11,Stim2.1==2,Val.12))]




# rawData[,`:=`(Val.vis=Val.02+Val.11,Val.res=Val.01+Val.12)]

val.vars<-names(rawData) %>% grep(pattern = "Val.*",value = TRUE)


rawData[,AttMech:=factor(AttMech,levels = c(0,1,2),
                      labels = c("no mechanism","Mackintosh","Pearce-Hall"))]

rawData[,`:=`(Client1=factor(Client1,levels = c(0,1,2),labels = c("resident",
                                                       "visitor","absence")),
              Client2=factor(Client2,levels = c(0,1,2),labels = c("resident",
                                                      "visitor","absence")),
              Stim1.0=factor(Stim1.0),Stim2.0=factor(Stim2.0),
              Stim1.1=factor(Stim1.1),Stim2.1=factor(Stim2.1)
)]

Values.long<-melt(rawData,id.vars = c("Age","Training","AttMech"),
                  measure.vars = val.vars[1:6],
                  variable.name = "Features")


png(here("Simulations",scenario,"clientTypesVal.png"))
ggplot(rawData[Age%%500==0],
       aes(y=Val.Clien1,x=Age,color=Client1,
                               shape=AttMech))+
  geom_point()+
  geom_hline(yintercept = c(0,1,2),color="grey")+
  theme_classic()
dev.off()

png(here("Simulations",scenario,"features_val.png"))
ggplot(data = Values.long[Age%%500==0],
       aes(x=Age,y=value,color=Features))+
  stat_summary(fun.data = function(x) {
    xFive<-fivenum(x)
    return(data.frame(y=xFive[3],ymax=xFive[4],ymin=xFive[2]))},
    position = position_dodge(width=200))+
  # geom_path()+
  geom_hline(yintercept = c(1,2),color="black")+
  facet_wrap(~AttMech)+
  theme_classic()
dev.off()

png(here("Simulations",scenario,"Mack_feature.png"))
ggplot(data = Values.long[AttMech=="Mackintosh" & Age%%500==0],
       aes(x=Age,y=value,color=Features))+
  geom_path()+ylim(0,3)+
  geom_hline(yintercept = c(1,2),color="black")+
  facet_wrap(~Training)+labs(title = "Mackintosh")+
  theme_classic()
dev.off()

png(here("Simulations",scenario,"alphaDyn.png"))
Alphas.long<-melt(rawData,id.vars = c("Age","Training","AttMech"),
                  measure.vars = patterns("alpha.*"),
                  variable.name = "Stimulus")

ggplot(data = Alphas.long[Age%%500==0],aes(x=Age,y=value,color=Stimulus))+
  stat_summary(fun.data = function(x) {
    xFive<-fivenum(x)
    return(data.frame(y=xFive[3],ymax=xFive[4],ymin=xFive[2]))},
    geom = "pointrange",position = position_dodge(width=200))+
  # +
  # geom_step()+
  facet_wrap(~AttMech)+
  theme_classic()
dev.off()

png(here("Simulations",scenario,"Mack_alpha_dyn.png"))
ggplot(data = Alphas.long[AttMech=="Mackintosh" & Age%%500==0],
       aes(x=Age,y=value,color=Stimulus))+
  geom_step()+
  facet_wrap(~Training)+
  facet_wrap(~Training)+labs(title = "Mackintosh")+
  theme_classic()
dev.off()

png(here("Simulations",scenario,"PH_alpha_dyn.png"))
ggplot(data = Alphas.long[AttMech=="Pearce-Hall" ],
       aes(x=Age,y=value,color=Stimulus))+
  geom_step()+
  geom_hline(yintercept = 1)+
  facet_wrap(~Training)+
  facet_wrap(~Training)+labs(title ="Pearce-Hall")+
  theme_classic()
  
dev.off()


png(here("Simulations",scenario,"Client_features.png"),width = 966,height = 524)
plot_grid(
  ggplot(data=rawData,aes(x=Client1,fill=Stim1.0))+
    geom_bar(position = "stack")+ylab("frequency")+xlab("Client type")+
    guides(fill=guide_legend(title="Stimulus 1"))+
    theme_classic(),
  ggplot(data=rawData,aes(x=Client1,fill=Stim1.1))+
    geom_bar(position = "stack")+
    geom_bar(position = "stack")+ylab("frequency")+xlab("Client type")+
    guides(fill=guide_legend(title="Stimulus 2"))+
    theme_classic()
)
dev.off()

plot_grid(
  ggplot(data=rawData,aes(x=Client2,fill=Stim2.0))+
    geom_bar(position = "stack")+
    theme_classic() ,
  ggplot(data=rawData,aes(x=Client2,fill=Stim2.1))+
    geom_bar(position = "stack")+
    theme_classic()
)

PAAtimeInt[,AttMech:=factor(AttMech,levels = c(0,1,2),
                            labels = c("no mechanism","Mackintosh","Pearce-Hall"))]
png(here("Simulations",scenario,"RVchoice.png"))

ggplot(PAAtimeInt,aes(x=Interv,y=Prob.RV.V,color=AttMech))+
  stat_summary(fun.data = function(x) {
      xFive<-fivenum(x)
      return(data.frame(y=xFive[3],ymax=xFive[4],ymin=xFive[2]))
  })+
  ylim(0,1)+
  geom_hline(yintercept = soft_max(1,2,0.5))+
  theme_classic()
dev.off()
