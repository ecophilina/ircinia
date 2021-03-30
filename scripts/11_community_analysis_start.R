# organize data for community analysis
library(tidyverse)

# bring in data
source("scripts/03_reimport.R")  

# prep primary producer data
sg<-sg_shoot%>%
  mutate(sg.sd.global=mean(SD))%>%
  group_by(treatment,plot,sampling,sg.sd.global)%>%
  summarize(sg.sd=mean(SD))%>%
  mutate(
    sg.sd=sg.sd-sg.sd.global, # comment this line out if not centering
    sampling=case_when(
      sampling==1~0,
      sampling==2~1,
      sampling==3~5,
      sampling==4~12,
      sampling==5~17)) 

sggrow<-sg_grow%>%
  filter(dist %in% c(0,0.5))%>%
  mutate(grow.global=mean(total.growth.mm2/days))%>%
  group_by(treatment,plot,sampling,grow.global)%>%
  summarize(grow=mean(total.growth.mm2/days))%>%
  mutate(
    grow=grow-grow.global,
    sampling=case_when(
      sampling==1~0,
      sampling==2~1,
      sampling==3~5,
      sampling==4~12,
      sampling==5~17))


# Algae. No need to filter for has algae vs doesn't have algae because no plots had a 0 count for algae :)
alg<-algae%>%
  filter(!taxa %in% c("brown.cyanobacteria","green.cyanobacteria","dictyota"))%>%
  select(-taxa)%>%
  group_by(treatment,plot,sampling)%>%
  summarize(abund=sum(abundance))%>%
  ungroup()%>%
  distinct()%>%
  mutate(abund.global=mean(abund))%>%
  mutate(
    abund=abund-abund.global,
    sampling=case_when(
      sampling==1~0,
      sampling==2~1,
      sampling==3~5,
      sampling==4~12,
      sampling==5~17))

productivity<-left_join(sg,sggrow)%>%
  left_join(alg)

# create community datasets

alg.com.full<-algae%>%
  filter(!taxa %in% c("brown.cyanobacteria","green.cyanobacteria","dictyota"))%>%
  mutate(sampling=case_when(
    sampling==1~0,
    sampling==2~1,
    sampling==3~5,
    sampling==4~12,
    sampling==5~17),
    season=case_when(
      sampling==0~"summer",
      sampling==1~"summer",
      sampling==5~"winter",
      sampling==12~"summer",
      sampling==17~"winter"),
    yr=case_when(
      sampling==0~0,
      sampling==1~1,
      sampling==5~1,
      sampling==12~2,
      sampling==17~2))%>%
  left_join(productivity)%>%
  pivot_wider(names_from=taxa,values_from=abundance,values_fill=0)

inv.com.full<-inverts%>%
  mutate(sampling=case_when(
    sampling==1~0,
    sampling==2~1,
    sampling==3~5,
    sampling==4~12,
    sampling==5~17),
    season=case_when(
      sampling==0~"summer",
      sampling==1~"summer",
      sampling==5~"winter",
      sampling==12~"summer",
      sampling==17~"winter"),
    yr=case_when(
      sampling==0~0,
      sampling==1~1,
      sampling==5~1,
      sampling==12~2,
      sampling==17~2))%>%
  left_join(productivity)%>%
  pivot_wider(names_from=taxa,values_from=abundance,values_fill=0)

fish.com.full<-fish%>%
  mutate(sampling=case_when(
    sampling==1~0,
    sampling==2~1,
    sampling==3~5,
    sampling==4~12,
    sampling==5~17),
    season=case_when(
      sampling==0~"summer",
      sampling==1~"summer",
      sampling==5~"winter",
      sampling==12~"summer",
      sampling==17~"winter"),
    yr=case_when(
      sampling==0~0,
      sampling==1~1,
      sampling==5~1,
      sampling==12~2,
      sampling==17~2))%>%
  left_join(productivity)%>%
  pivot_wider(names_from=taxa,values_from=abundance,values_fill=0)

alg.env<-alg.com.full %>%
  select(treatment,plot,sampling,season,yr,sg.sd.global,
         sg.sd,grow.global,grow,abund,abund.global)

alg.com<-alg.com.full %>%
  select(-treatment,-plot,-sampling,-season,-yr,-sg.sd.global,
         -sg.sd,-grow.global,-grow,-abund,-abund.global)

alg.uni<-alg.env%>%
  mutate(spr=vegan::specnumber(alg.com),
         div=vegan::diversity(alg.com,index = "shannon"),
         j=div/log(spr),
         j=ifelse(is.na(j),0,j))

inv.env<-inv.com.full %>%
  select(treatment,plot,sampling,season,yr,sg.sd.global,
         sg.sd,grow.global,grow,abund,abund.global)

inv.com<-inv.com.full %>%
  select(-treatment,-plot,-sampling,-season,-yr,-sg.sd.global,
         -sg.sd,-grow.global,-grow,-abund,-abund.global)

inv.uni<-inv.env%>%
  mutate(spr=vegan::specnumber(inv.com),
         div=vegan::diversity(inv.com,index = "shannon"),
         j=div/log(spr),
         j=ifelse(is.na(j),0,j))

fish.env<-fish.com.full %>%
  select(treatment,plot,sampling,season,yr,sg.sd.global,
         sg.sd,grow.global,grow,abund,abund.global)

fish.com<-fish.com.full %>%
  select(-treatment,-plot,-sampling,-season,-yr,-sg.sd.global,
         -sg.sd,-grow.global,-grow,-abund,-abund.global)

fish.uni<-fish.env%>%
  mutate(spr=vegan::specnumber(fish.com),
         div=vegan::diversity(fish.com,index = "shannon"),
         j=div/log(spr),
         j=ifelse(is.na(j),0,j))


