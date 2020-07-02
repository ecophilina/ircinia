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
# distance from the sponge at sampling 5

# ---- onevfourat0 ----
ggplot(sg_nuts %>%
         filter(sampling==1|sampling==4)%>%
         filter(dist==0 & !is.na(PN)))+
  geom_boxplot(aes(x=as.factor(sampling),y=PN,fill=treatment))
# looks like there's some differences but a lot of variation next plot to look at how individual plots change
# first going to calculate the change in PN (going to do PC and PP at the same time)
dnuts<-sg_nuts%>%
  filter(sampling==1|sampling==4)%>%
  filter(dist==0)%>%
  pivot_wider(-mg,names_from=sampling,values_from=c(PN,PP,PC))%>%
  mutate(dn=(PN_4-PN_1)/PN_1,dc=(PC_4-PC_1)/PC_1,dp=(PP_4-PP_1)/PP_1)%>%
  select(-4:-9)%>%
  pivot_longer(4:6,names_to = "nuts",values_to = "delta")

ggplot(dnuts,aes(y=delta,fill=treatment))+
  geom_boxplot()+
  facet_wrap(~nuts,scales = "free")


# ---- onev5at0 ----
ggplot(sg_nuts %>%
         filter(sampling==1|sampling==5)%>%
         filter(dist==0 & !is.na(PN)))+
  geom_boxplot(aes(x=as.factor(sampling),y=PN,fill=treatment))
# looks like there's some differences but a lot of variation next plot to look at how individual plots change
# first going to calculate the change in PN (going to do PC and PP at the same time)
dnuts5<-sg_nuts%>%
  filter(sampling==1|sampling==5)%>%
  filter(dist==0)%>%
  pivot_wider(-mg,names_from=sampling,values_from=c(PN,PP,PC))%>%
  mutate(dn=(PN_5-PN_1)/PN_1,dc=(PC_5-PC_1)/PC_1,dp=(PP_5-PP_1)/PP_1)%>%
  select(-4:-9)%>%
  pivot_longer(4:6,names_to = "nuts",values_to = "delta")

ggplot(dnuts5,aes(y=delta,fill=treatment))+
  geom_boxplot()+
  facet_wrap(~nuts,scales = "free")
# pattern isn't the same at 5 as it was at 4.

# ----dist----
ggplot(sg_nuts%>%
         filter(sampling==5|sampling==1)%>%
         pivot_longer(6:8,names_to="nuts",values_to = "percent"))+
  geom_boxplot(aes(x=as.factor(dist),y=percent,fill=treatment))+
  facet_grid(nuts~sampling,scales = "free")
  
