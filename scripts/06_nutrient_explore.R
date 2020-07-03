# How does the sponge Ircinia felix influence seagrass bed primary producers? #

# this script will examining nutrient data #
# run 03_reimport.R fist #

# packages----
if(!require(tidyverse))install.packages("tidyverse");library(tidyverse)

# find missing data----

nut.miss<-sg_nuts%>%
  group_by(treatment,
           plot,
           dist,
           sampling)%>%
  summarize(pn=ifelse(is.na(PN),1,0),
            pp=ifelse(is.na(PP),1,0),
            pc=ifelse(is.na(PC),1,0),
            miss=pn+pp+pc)# this dataset has a 1 anywhere there's missing data 
# there should be 5 measures per treatment per distance per sampling. So I'm going to group by treatment, distance, and sampling
# to see if there are some with 0 bits of missing data across the board. If not, find which ones have the least amount of missing data

nm2<-nut.miss%>%
  group_by(treatment,
           dist,
           sampling)%>%
  summarize(pn=sum(pn),
            pp=sum(pp),
            pc=sum(pc))%>%
  pivot_longer(4:6,names_to = "nutrient",values_to = "nmiss")

# making plots to make this easier

ggplot(nm2,aes(x=nutrient,y=nmiss,fill=treatment))+
  geom_bar(stat="identity",position=position_dodge())+
  facet_grid(sampling~dist)

# inconveniently not much is all the way there. But it looks like we could look at % phosphorus at most distances
# and samplings
# sampling 5 is the best (after sampling 1) for data availability, but this is in the winter. If we go with
# sampling 1 and 4 we could look at the 0 distance... if we go with 0.5 we could leave out sampling 4. I kind of
# vote for looking at patterns at 0 between 1 and 4 and 5 and then looking at patterns associated with
# distance from the sponge at sampling 5 first I'm going to look at mean values over distance and time.

# ----dist----
ggplot(sg_nuts%>%
         pivot_longer(6:8,names_to="nuts",values_to = "percent")%>%
         group_by(treatment,dist,sampling,nuts)%>%
         summarize(mn=mean(percent,na.rm = TRUE)))+
  geom_point(aes(x=as.factor(dist),y=mn,color=treatment))+
  facet_grid(nuts~sampling,scales = "free")
  
# there don't seem to be logical consistent patterns in the nutrient data on first glance. I'm going to look
# at ratios (i.e. C:N, C:P, and N:P)

# --- ratios---
nr<-sg_nuts %>%
  mutate(cn=PC/PN,
         cp=PC/PP,
         np=PN/PP)%>%
  filter(!is.na(cn)|!is.na(cp)|!is.na(np))
# look at carbon:phosphorus
ggplot(nr)+
  geom_boxplot(aes(y=cp,fill=treatment))+
  facet_grid(dist~sampling)
# don't see any interesting patterns here
# look at carbon:nitrogen
ggplot(nr)+
  geom_boxplot(aes(y=cn,fill=treatment))+
  facet_grid(dist~sampling)
# don't see any interesting patterns here
# look at nitrogen:phosphorus
ggplot(nr)+
  geom_boxplot(aes(y=np,fill=treatment))+
  facet_grid(dist~sampling)
# don't see any interesting patterns here

# originally I looked at the change in % nutrients from sampling 1 and 4- there is a pattern there, but after
# being removed from the experiment for so long and looking at the data more thoroughly I don't know if I
# actually buy it. Maybe if I look at the change from initial conditions for each sampling point it 
# will make me feel better about that?

# calculating change from initial conditions

dn<-sg_nuts%>%
  pivot_wider(-mg,names_from = sampling,values_from = c(PN,PC,PP))%>%
  mutate(pn12=(PN_2-PN_1)/PN_1,
         pn13=(PN_3-PN_1)/PN_1,
         pn14=(PN_4-PN_1)/PN_1,
         pn15=(PN_5-PN_1)/PN_1,
         pp12=(PP_2-PP_1)/PP_1,
         pp13=(PP_3-PP_1)/PP_1,
         pp14=(PP_4-PP_1)/PP_1,
         pp15=(PP_5-PP_1)/PP_1,
         pc12=(PC_2-PC_1)/PC_1,
         pc13=(PC_3-PC_1)/PC_1,
         pc14=(PC_4-PC_1)/PC_1,
         pc15=(PC_5-PC_1)/PC_1)%>%
  select(treatment,plot,dist,pn12,pn13,pn14,pn15,pp12,pp13,pp14,pp15,pc12,pc13,pc14,pc15)%>%
  pivot_longer(-1:-3,names_to = "m",values_to = "dn")%>%
  separate(m,into = c("nut","sp"),sep = 2)

ggplot(dn)+
  geom_boxplot(aes(x=as.factor(dist),y=dn,fill=treatment))+
  facet_grid(nut~sp)

# nope that didn't make me feel better. 
# my conclusion is that its hard to explain the nutrient data without the fish data (fish pee a lot and it
# influences seagrass nutrient content). So, if we're keeping this manuscript to the primary producers
# then I think we have to leave the nutrient data out. 
