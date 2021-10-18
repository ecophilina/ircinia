# How does the sponge Ircinia felix influence seagrass bed primary producers? #
# This script will examining nutrient data #
# packages----
if(!require(tidyverse))install.packages("tidyverse");library(tidyverse)
if(!require(car))install.packages("car");library(car)
if(!require(lmerTest))install.packages("lmerTest");library(lmerTest)
if(!require(htmlTable))install.packages("htmlTable");library(htmlTable)

theme_set(theme_bw())
# load data----
source("MarineBiology_DOI_10.1007/03_reimport.R")#imports all the data sets
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
filter(sgn,nut=="PN" & nvalue>2.5)
ggplot(data=sgn,aes(nvalue))+
  geom_histogram()+
  facet_wrap(~nut,scales="free")

# create a dataset with offset values and season
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

dsgn$treatment<-factor(dsgn$treatment)

# look at sample sizes to see what makes the most sense
(nuttable <- dsgn %>%
    filter(nut %in% c("PP", "PC", "PN")) %>%
    mutate(
      treatment = case_when(
        treatment == "blank" ~ "Control",
        treatment == "fake" ~ "Structure",
        treatment == "real" ~ "Sponge")) %>%
    group_by(treatment, yr, season, dist, nut) %>%
    summarize(ns = n()) %>%
    pivot_wider(
      names_from = c(treatment, nut),
      values_from = ns) %>%
    ungroup() %>%
    select(
      `Distance (m)` = dist,
      Control_PC,
      Control_PN,
      Control_PP,
      Structure_PC,
      Structure_PN,
      Structure_PP,
      Sponge_PC,
      Sponge_PN,
      Sponge_PP) %>%
   # distinct() %>%
    htmlTable(
      rgroup = c("Summer", "Winter", "Summer", "Winter"),
      n.rgroup = c(4, 4, 4, 4),
      tspanner = c("Year 1", "Year 2"),
      n.tspanner = c(8, 8),
      cgroup = c("", "Control", "Structure", "Sponge"),
      n.cgroup = c(1, 3, 3, 3),
      header=c("Distance (m)","%C","%N","%P","%C","%N","%P","%C","% N","% P")))

# now do a model for nitrogen, at 0 distance
# plot this again to make sure I remember what it looks like
ggplot(dsgn %>% 
         filter(dist==0)%>%
         filter(nut=="PN")%>%
         filter(yr==2))+
  geom_jitter(aes(x=season,y=delta,color=treatment),width=.1)

pnb <- lmer(nvalue ~ treatment*season +(1 | plot),
            offset = snv,
            data = dsgn %>% 
              ungroup()%>%
              filter(dist==0)%>%
              filter(nut=="PN")%>%
              filter(yr==2)%>%
              #filter(season=="Summer")%>%
              mutate(treatment = relevel(treatment, ref = "real")))
summary(pnb)

pnb2 <- lmer(nvalue ~ treatment*sampling +(1 | plot),
            data = sgn %>% 
              ungroup()%>%
              filter(dist==0)%>%
              filter(nut=="PN")%>%
              filter(sampling %in% c(1,4))%>%
              #filter(season=="Summer")%>%
              mutate(treatment = relevel(treatment, ref = "real")))
summary(pnb2)

#Phosphorus
ggplot(dsgn %>% 
         filter(dist==0)%>%
         filter(nut=="PP")%>%
         filter(yr==2))+
  geom_jitter(aes(x=season,y=delta,color=treatment),width=.1)

ppb <- lmer(nvalue ~ treatment*season +(1 | plot),
            offset = snv,
            data = dsgn %>% 
              ungroup()%>%
              filter(dist==0)%>%
              filter(nut=="PP")%>%
              filter(yr==2)%>%
              #filter(season=="Summer")%>%
              mutate(treatment = relevel(treatment, ref = "real")))
summary(ppb)

ppb2 <- lmer(nvalue ~ treatment*sampling +(1 | plot),
             data = sgn %>% 
               ungroup()%>%
               filter(dist==0)%>%
               filter(nut=="PP")%>%
               filter(sampling %in% c(1,4))%>%
               #filter(season=="Summer")%>%
               mutate(treatment = relevel(treatment, ref = "real")))
summary(ppb2)

#Carbon
ggplot(dsgn %>% 
         filter(dist==0)%>%
         filter(nut=="PC")%>%
         filter(yr==2))+
  geom_jitter(aes(x=season,y=delta,color=treatment),width=.1)

PCb <- lmer(nvalue ~ treatment*season +(1 | plot),
            offset = snv,
            data = dsgn %>% 
              ungroup()%>%
              filter(dist==0)%>%
              filter(nut=="PC")%>%
              filter(yr==2)%>%
              #filter(season=="Summer")%>%
              mutate(treatment = relevel(treatment, ref = "real")))
summary(PCb)

PCb2 <- lmer(nvalue ~ treatment*sampling +(1 | plot),
             data = sgn %>% 
               ungroup()%>%
               filter(dist==0)%>%
               filter(nut=="PC")%>%
               filter(sampling %in% c(1,4))%>%
               #filter(season=="Summer")%>%
               mutate(treatment = relevel(treatment, ref = "real")))
summary(PCb2)
