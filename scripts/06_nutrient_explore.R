# How does the sponge Ircinia felix influence seagrass bed primary producers? #

# this script will examining nutrient data #


# packages----
if(!require(tidyverse))install.packages("tidyverse");library(tidyverse)
if(!require(car))install.packages("car");library(car)
if(!require(lmerTest))install.packages("lmerTest");library(lmerTest)

theme_set(theme_bw())
# load data----
source("scripts/03_reimport.R")#imports all the data sets
mn<-readxl::read_xlsx("Original_data/missing_nutrients.xlsx")
sg_nuts<-bind_rows(sg_nuts,mn)%>%
  select(-mg)%>%
  group_by(treatment,plot,dist,sampling)%>%
  summarize(PN=sum(PN,na.rm = TRUE),PC=sum(PC,na.rm = TRUE),PP=sum(PP,na.rm = TRUE))
# find missing data----

nm2<-sg_nuts%>%
  mutate(pn=ifelse(PN==0,1,0),
            pp=ifelse(PP==0,1,0),
            pc=ifelse(PC==0,1,0))%>%
  group_by(treatment,
           dist,
           sampling)%>%
  summarize(pn=sum(pn),
            pp=sum(pp),
            pc=sum(pc))%>%
  pivot_longer(pn:pc,names_to = "nutrient",values_to = "nmiss")

# making plots to make this easier

ggplot(nm2,aes(x=nutrient,y=nmiss,fill=treatment))+
  geom_bar(stat="identity",position=position_dodge())+
  facet_grid(sampling~dist)
  

#now no treatment is missing more than 1 sample at each sampling/dist combination

# now I'm going to look at the distributions of the variables

# Because as the experiment goes on fish abundance may be correlated with nutrients 
# I'm going to look at the fish data and how fish abundance correlates with nutrient values

f2<-fish %>%
  group_by(treatment,plot,sampling)%>%
  summarize(abund=sum(abundance))
# join this with nutrient data

sgn<-sg_nuts %>%
  mutate(cn=PC/PN,
         cp=PC/PP,
         np=PN/PP)%>%
  pivot_longer(PN:np,names_to = "nut",values_to = "nvalue")%>%
  left_join(f2)%>%
  mutate(abund=ifelse(is.na(abund),0,abund),)%>%
  filter(!is.na(nvalue))

# before I plot up I'm going to look at the distributions of the data
ggplot(data=sgn,aes(nvalue))+
  geom_histogram()+
  facet_wrap(~nut,scales="free")
# there are some outliers in PC and PN- look at what these are
filter(sgn,nut=="PC" & nvalue>55)
filter(sgn,nut=="PN" & nvalue>2.5)
# carbon and nitrogen are calculated from the same sample and plot 6 (sampling 2 dist 0)
# is an outlier for both- makes me think something besides seagrass got in there.
# I'm going to remove this point before going further. 

sgn_no<-sgn %>%
  filter(plot==6 & dist ==0 & sampling==2 & nut %in% c("PC","PN"))

sgn<-sgn %>%
  anti_join(sgn_no)
# now plot up 

ggplot(data=sgn,aes(x=abund,y=nvalue,color=treatment))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_grid(nut~dist,scales="free")

#there is one plot with a lot more fish than all the others- going to look at 
# things without that plot
ggplot(data=sgn%>% filter(abund<20),aes(x=abund,y=nvalue,color=treatment))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_grid(nut~dist,scales="free")
# doesn't actually look like there's much going on there. 

# going to do the same thing with growth

sgg<-sg_grow%>%
  mutate(gpd=total.growth.mm2/days)%>%
  group_by(treatment,plot,dist,sampling)%>%
  summarize(mgd=mean(gpd),sdgd=sd(gpd))

sgn<-sgn %>%
  left_join(sgg)%>%
  mutate(mnths=case_when(
    sampling==1~0,
    sampling==2~1,
    sampling==3~5,
    sampling==4~12,
    sampling==5~17))

ggplot(data=sgn,aes(x=mgd,y=nvalue,color=treatment))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~nut,scales="free")

# well there doesn't seem to be much there either.

# now I'm going to look at the change in nutrients from initial conditions

sv<-filter(sgn, sampling == 1) %>% 
  rename(snv = nvalue)%>%select(treatment,plot,dist,snv,nut)

dsgn<-left_join(sgn,sv)%>%
  mutate(delta=nvalue-snv)%>%filter(sampling!=1)

ggplot(data=dsgn%>%
         filter(dist==0 &nut %in% c("PC","PN","PP"))%>%
         group_by(treatment,dist,sampling,nut)%>%
         summarize(mv=mean(delta),sev=sd(delta)/sqrt(n())),
       aes(y=mv,color=treatment,x=as.factor(sampling)))+
  geom_point(position=position_dodge(width=.3),size=2)+
  geom_errorbar(aes(ymin=mv-sev,ymax=mv+sev,color=treatment),width=.1,position=position_dodge(width=0.3))+
  facet_grid(rows = vars(nut),scales="free")

# looks like we should look at %C, %N, and %P
pn05<-aov(nvalue~treatment,data=sgn%>%filter(sampling==1&nut=="PN"&dist==0.5))
summary(pn05)
TukeyHSD(pn05)
ipn<-lmer(nvalue~treatment+as.factor(dist)+(1|plot),data=sgn%>%filter(sampling==1&nut=="PN"))
summary(ipn)

ipp<-lmer(nvalue~treatment+(1|plot),data=sgn%>%filter(sampling==1&nut=="PP"))
summary(ipp)

ipc<-lmer(nvalue~treatment+dist+(1|plot),data=sgn%>%filter(sampling==1&nut=="PC"))
summary(ipc)

pn<-aov(delta~treatment,data=dsgn %>%filter(sampling%in%c(4)&dist==0&nut=="PN"))
summary(pn)
TukeyHSD(pn)

pp<-aov(delta~treatment,data=dsgn %>%filter(sampling%in%c(4)&dist==0&nut=="PP"))
summary(pp)

pc<-aov(delta~treatment,data=dsgn %>%filter(sampling%in%c(4)&dist==0&nut=="PC"))
summary(pc)

