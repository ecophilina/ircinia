# How does the sponge Ircinia felix influence seagrass bed primary producers? #
# This script will examining nutrient data #
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

#make seagrass nutriet data frame with nvalues filled in
sgn<-sg_nuts %>%
  filter(PC >0 & PN >0 & PP >0)%>%
  mutate(cn=PC/PN,
         cp=PC/PP,
         np=PN/PP)%>%
  pivot_longer(PN:np,names_to = "nut",values_to = "nvalue")%>%
  filter(!is.na(nvalue))

# I'm going to look at the distributions of the data
ggplot(data=sgn,aes(nvalue))+
  geom_histogram()+
  facet_wrap(~nut,scales="free")
# there are some outliers in PC and PN- look at what these are
filter(sgn,nut=="PC" & nvalue>55)
filter(sgn,nut=="PN" & nvalue>2.5)
# carbon and nitrogen are calculated from the same sample and plot 6 (sampling 2 dist 0), 
# plot 7 sampling 2 dist 0, and plot 8 sampling 2 dist 1
# are outliers for both- makes me think something besides seagrass got in there.
# I'm going to remove these points before going further. 

sgn_no<-bind_rows(sgn %>%
                    filter(plot==7 & sampling==2 & dist==0),
                  sgn %>%
                    filter(plot==8 & sampling==2 & dist==1))

sgn<-sgn %>%
  anti_join(sgn_no)

# now I'm going to look at the change in nutrients from initial conditions

sv<-filter(sgn, sampling == 1) %>% 
  rename(snv = nvalue)%>%select(treatment,plot,dist,snv,nut)

dsgn<-left_join(sgn,sv)%>%
  mutate(delta=nvalue-snv,
    season=case_when(
      sampling==2~"Summer",
      sampling==3~"Winter",
      sampling==4~"Summer",
      sampling==5~"Winter"),
    yr=case_when(
      sampling==2~1,
      sampling==3~1,
      sampling==4~2,
      sampling==5~2),
    mnths=case_when(
      sampling==1~0,
      sampling==2~1,
      sampling==3~5,
      sampling==4~12,
      sampling==5~17))%>%filter(sampling!=1)
# look at sample sizes here
ggplot(dsgn%>%
  group_by(treatment,sampling,dist,nut)%>%
  filter(nut %in% c("PN","PP","PC"))%>%
  summarize(ns=n()))+
  geom_bar(aes(x=treatment,y=ns,fill=nut),position=position_dodge(),stat="identity")+
  facet_grid(sampling~dist)
# check for season effect
Nseason<-lmer(nvalue~season+(1|plot), data = dsgn%>%
    filter(nut=="PN")) 
(Nseason<-Anova(Nseason, type = "III"))

Pseason<-lmer(nvalue~season+(1|plot), data = dsgn%>%
    filter(nut=="PP")) 
(Pseason<-Anova(Pseason, type = "III"))

Cseason<-lmer(nvalue~season+(1|plot), data = dsgn%>%
    filter(nut=="PC")) 
(Cseason<-Anova(Cseason, type = "III"))
# YES for all

dsgn%>% filter(dist==0 & nut %in% c("PC","PN","PP"))%>%
  # group_by(treatment,dist,sampling,nut,season,yr,mnths)%>%
  # summarize(mv=mean(delta),sd=sd(delta)) %>%
  # ggplot(aes(y=mv,color=treatment,x=as.factor(mnths)))+
  # geom_point(position=position_dodge(width=.3),size=2)+
  # geom_errorbar(aes(ymin=mv-sd,ymax=mv+sd,color=treatment),width=.1,position=position_dodge(width=0.3))+
  ggplot(aes(y=delta,color=treatment,x=as.factor(mnths)))+
  geom_point(position=position_dodge(width=.3),size=1.5,alpha=0.7)+
  xlab("Months into experiment") +
  scale_color_viridis_d(
    option="A",
    begin=0, end=0.6,
    name="",
    labels=c("Control","Structure","Sponge"))+
  facet_grid(rows = vars(nut),cols = vars(season),scales="free") +
  ggtitle("At center of plot")

dsgn %>% filter(dist==0.5 &nut %in% c("PC","PN","PP"))%>%
    ggplot(aes(y=delta,color=treatment,x=as.factor(mnths)))+
    geom_point(position=position_dodge(width=.3),size=1.5,alpha=0.7)+
    xlab("Months into experiment") +
  scale_color_viridis_d(
    option="A",
    begin=0, end=0.6,
    name="",
    labels=c("Control","Structure","Sponge"))+
  facet_grid(rows = vars(nut),cols = vars(season),scales="free") +
  ggtitle("At 0.5m from center of plot")

dsgn %>% filter(dist==1 &nut %in% c("PC","PN","PP"))%>%
  ggplot(aes(y=delta,color=treatment,x=as.factor(mnths)))+
  geom_point(position=position_dodge(width=.3),size=1.5,alpha=0.7)+
  xlab("Months into experiment") +
  scale_color_viridis_d(
    option="A",
    begin=0, end=0.6,
    name="",
    labels=c("Control","Structure","Sponge"))+
  facet_grid(rows = vars(nut),cols = vars(season),scales="free") +
  ggtitle("At 1m from center of plot")

dsgn %>% filter(dist==2 &nut %in% c("PC","PN","PP"))%>%
  ggplot(aes(y=delta,color=treatment,x=as.factor(mnths)))+
  geom_point(position=position_dodge(width=.3),size=1.5,alpha=0.7)+
  xlab("Months into experiment") +
  scale_color_viridis_d(
    option="A",
    begin=0, end=0.6,
    name="",
    labels=c("Control","Structure","Sponge"))+
  facet_grid(rows = vars(nut),cols = vars(season),scales="free") +
  ggtitle("At 2m from center of plot")

# because we don't have a before sample during winter, we focus on summer
# look at change after 1 year
#in %N
n.lmerb<-lmer(nvalue~treatment*sampling+(1|plot),
              data=sgn%>%
                filter(nut=="PN" & dist==0 & sampling %in% c(1,4)))
(nbs<-summary(n.lmerb))
# relevel with real as the base
# first make treatment a factor
sgn <- sgn %>% ungroup() 
sgn$treatment<-as.factor(sgn$treatment)

sgn$sampling<-as.factor(sgn$sampling)

n.lmerr<-lmer(nvalue~treatment*sampling+(1|plot),
              data=sgn%>%
                filter(nut=="PN" & dist==0 & sampling %in% c(1,4))%>%
                mutate(treatment=relevel(treatment,ref="real")))
(nrs<-summary(n.lmerr))

# one year later seagrass in plots with sponges have higher %N than either treatment
# even though they started out lower than both and significantly lower than blank.
# although this increase is not statistically significant. 
# seagrass in blank plots did decrease significantly- check fake

n.lmerf<-lmer(nvalue~treatment*sampling+(1|plot),
              data=sgn%>%
                filter(nut=="PN" & dist==0 & sampling %in% c(1,4))%>%
                mutate(treatment=relevel(treatment,ref="fake")))
(nfs<-summary(n.lmerf))
# fake decreased but not significantly.

# now doing percent phosphorus
p.lmerb<-lmer(nvalue~treatment*sampling+(1|plot),
              data=sgn%>%
                filter(nut=="PP" & dist==0 & sampling %in% c(1,4)))
pbs<-summary(p.lmerb)
# no significant change in %P over time in blank plots (but it decreased), but a significant increase in 
# real compared to blank one year later
# relevel with real as the base
p.lmerr<-lmer(nvalue~treatment*sampling+(1|plot),
              data=sgn%>%
                filter(nut=="PP" & dist==0 & sampling %in% c(1,4))%>%
                mutate(treatment=relevel(treatment,ref="real")))
(prs<-summary(p.lmerr))
# one year later seagrass in plots with sponges have higher %P than blank, but not fake
# % P increased in real plots over the year,
# although this increase is not statistically significant. 
# check trend in fake

p.lmerf<-lmer(nvalue~treatment*sampling+(1|plot),
              data=sgn%>%
                filter(nut=="PP" & dist==0 & sampling %in% c(1,4))%>%
                mutate(treatment=relevel(treatment,ref="fake")))
pfs<-summary(p.lmerf)
# %P decreased in fake plots, but again, not significantly

# look at %C now

c.lmerb<-lmer(nvalue~treatment*sampling+(1|plot),
              data=sgn%>%
                filter(nut=="PC" & dist==0 & sampling %in% c(1,4)))
cbs<-summary(c.lmerb)
# real was significantly lower than blank at the beginning
# decrease in %C over the year in blank, but its not significant
# barely not significant increase in %C in sponge plots after 1 year.
# relevel with real as the base
c.lmerr<-lmer(nvalue~treatment*sampling+(1|plot),
              data=sgn%>%
                filter(nut=="PC" & dist==0 & sampling %in% c(1,4))%>%
                mutate(treatment=relevel(treatment,ref="real")))
crls<-summary(c.lmerr)
# increase in %C in sponge plots, but again not significant. There is a significant 
# difference between real and fake after a year

c.lmerf<-lmer(nvalue~treatment*sampling+(1|plot),
              data=sgn%>%
                filter(nut=="PC" & dist==0 & sampling %in% c(1,4))%>%
                mutate(treatment=relevel(treatment,ref="fake")))
cfs<-summary(c.lmerf)
# %C decreased significantly in fake plots after 1 year.

# check distributions
sgn1 <- sgn %>% filter(nut=="PC" & dist==0 & sampling %in% c(1,4))
hist(sgn1$nvalue, breaks=30)
sgn2 <- sgn %>% filter(nut=="PP" & dist==0 & sampling %in% c(1,4))
hist(sgn2$nvalue, breaks=30)
sgn3 <- sgn %>% filter(nut=="PN" & dist==0 & sampling %in% c(1,4))
hist(sgn3$nvalue, breaks=30)


# looking at change after 5 months
n.lmerr5<-lmer(nvalue~treatment*sampling+(1|plot),
              data=sgn%>%
                filter(nut=="PN" & dist==0 & sampling %in% c(1,3))%>%
                mutate(treatment=relevel(treatment,ref="real")))
(nrs5<-summary(n.lmerr5))
# already an increase in %N in real plots

# now doing percent phosphorus
p.lmerr5<-lmer(nvalue~treatment*sampling+(1|plot),
              data=sgn%>%
                filter(nut=="PP" & dist==0 & sampling %in% c(1,5))%>%
                mutate(treatment=relevel(treatment,ref="real")))
(prs5<-summary(p.lmerr5))
# one year later seagrass in plots with sponges have higher %P than blank, but not fake
# % P increased in real plots over the year,
# although this increase is not statistically significant. 
# check trend in fake

# look at %C now

c.lmerr5<-lmer(nvalue~treatment*sampling+(1|plot),
              data=sgn%>%
                filter(nut=="PC" & dist==0 & sampling %in% c(1,3))%>%
                mutate(treatment=relevel(treatment,ref="real")))
(crls5<-summary(c.lmerr5))
# increase in %C in sponge plots, but again not significant. There is a significant 
# difference between real and fake after a year


