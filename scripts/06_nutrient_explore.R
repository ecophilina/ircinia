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
  summarize(fish.abund=sum(abundance))
# join this with nutrient data

sgn<-sg_nuts %>%
  filter(PC >0 & PN >0 & PP >0)%>%
  mutate(cn=PC/PN,
         cp=PC/PP,
         np=PN/PP)%>%
  pivot_longer(PN:np,names_to = "nut",values_to = "nvalue")%>%
  left_join(f2)%>%
  mutate(fish.abund=ifelse(is.na(fish.abund),0,fish.abund),)%>%
  filter(!is.na(nvalue))

# before I plot up I'm going to look at the distributions of the data
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

sgn_no<-sgn %>%
  filter(nut=="PC" & nvalue > 60)

sgn<-sgn %>%
  anti_join(sgn_no)
# now plot up 

ggplot(data=sgn,aes(x=log(fish.abund+1),y=nvalue,color=treatment))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_grid(nut~dist,scales="free")

#there is one plot with a lot more fish than all the others- going to look at 
# things without that plot
ggplot(data=sgn%>% filter(fish.abund<20),aes(x=fish.abund,y=nvalue,color=treatment))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_grid(nut~dist,scales="free")
# doesn't actually look like there's much going on there. 

#exclude immediate vacinity of sponges and still no convincing patterns
ggplot(data=sgn%>% filter(fish.abund<20 & !(treatment=="real" & dist==0)),aes(x=fish.abund,y=nvalue))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~nut,scales="free")

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
  mutate(delta=nvalue-snv,
    dist2=case_when(
      dist==0~"close",
      dist==0.5~"middle dist",
      dist>0.5~"no sponge")
    )%>%filter(sampling!=1)

# check all distances? Focusing on summer for simplicity
ggplot(data=dsgn%>%
    filter(
      sampling %in% c(2,4) & # summer only
      nut %in% c("PC","PN","PP")),
    aes(y=delta,color=treatment,x=as.factor(sampling)))+
  geom_point(position=position_dodge(width=.3),size=2, alpha = 0.5)+
  facet_grid(cols=vars(dist2), rows = vars(nut),scales="free")

# data super noisy, so let's focus only on immediate vacinity of sponges (dist == 0)

ggplot(data=dsgn%>%
         filter(dist==0 &nut %in% c("PC","PN","PP"))%>%
         group_by(treatment,dist,sampling,nut)%>%
         summarize(mv=mean(delta),sd=sd(delta)),
       aes(y=mv,color=treatment,x=as.factor(sampling)))+
  geom_point(position=position_dodge(width=.3),size=2)+
  geom_errorbar(aes(ymin=mv-sd,ymax=mv+sd,color=treatment),width=.1,position=position_dodge(width=0.3))+
  facet_grid(rows = vars(nut),scales="free")

# investigating if there are seasonal differences within blank plots

sgn<-sgn%>%
  mutate(season=case_when(
    sampling==1~"summer",
    sampling==2~"summer",
    sampling==4~"summer",
    sampling==3~"winter",
    sampling==5~"winter"))

ggplot(sgn %>%
         filter(treatment=="blank"))+
  geom_violin(aes(x=season,y=nvalue))+
  facet_wrap(~nut,scales = "free")

pns<-lmer(nvalue~season+(1|plot),
          data=sgn%>%
            filter(treatment=="blank" & nut=="PN"))
summary(pns)
# no difference due to season in %N in blank plots
pps<-lmer(nvalue~season+(1|plot),
          data=sgn%>%
            filter(treatment=="blank" & nut=="PP"))
summary(pps)
# no difference due to season in %P in blank plots
pcs<-lmer(nvalue~season+(1|plot),
          data=sgn%>%
            filter(treatment=="blank" & nut=="PC"))
summary(pcs)
# but there is a significant effect for carbon, but plot explains basically no variance

# look at change after 1 year
#in %N
n.lmerb<-lmer(nvalue~treatment*sampling+(1|plot),
             data=sgn%>%
               filter(nut=="PN" & dist==0 & sampling %in% c(1,4)))
nbs<-summary(n.lmerb)
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
nfs<-summary(n.lmerf)
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
prs<-summary(p.lmerr)
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
