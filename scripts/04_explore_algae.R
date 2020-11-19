library(tidyverse)

source("scripts/03_reimport.R")#imports all the data sets
if(!require(lmerTest))install.packages("lmerTest");library(lmerTest)



# ---- algaeplot ----
a2<-algae %>%
  pivot_wider(names_from = taxa,values_from = abundance,values_fill = list(abundance=0))%>%
  mutate(total=rowSums(select(.,-treatment,-plot,-sampling)))%>% # whenever you are using a non-tidyverse function within a 
  #tidyverse pipeline you have to use the "." to tell it to use the dataset you had been working with.
  pivot_longer(4:11,names_to = "taxa",values_to = "abundance")%>%
  mutate(prop=abundance/total,
         season=case_when(
           sampling==1~"summer",
           sampling==2~"summer",
           sampling==4~"summer",
           sampling==3~"winter",
           sampling==5~"winter"),
         plot=paste(treatment,plot),
         taxa=ifelse(taxa=="cladocephalus","udotea",taxa))
# now I'm going to make my plot
(c1<-ggplot(data=a2)+
    geom_bar(aes(x=sampling,y=abundance,fill=taxa),position="stack",stat="identity")+
    facet_wrap(~plot)+
    geom_rect(aes(xmin=2.5,xmax=3.5,ymin=-.01,ymax=1.01),alpha=0.009,size=1,color="black")+
    geom_rect(aes(xmin=4.5,xmax=5.5,ymin=-.01,ymax=1.01),alpha=0.009,size=1,color="black"))

# looking at abundance by species over time

ggplot(data=a2)+
  geom_point(aes(x=sampling,y=abundance,color=treatment))+
  geom_smooth(aes(x=sampling,y=abundance,color=treatment),method = "lm")+
  facet_wrap(~taxa,scales = "free_y")

a3<-a2%>%
  filter(!taxa %in% c("brown.cyanobacteria","dictyota","green.cyanobacteria"))

ggplot(data=a3)+
  geom_point(aes(x=sampling,y=abundance,color=treatment))+
  geom_smooth(aes(x=sampling,y=abundance,color=treatment),method = "lm")+
  facet_wrap(~taxa,scales = "free_y")

# looks like there's likely a good story here.

a4<-a3 %>%
  group_by(treatment, sampling,plot,season)%>%
  summarize(abund=sum(abundance))

ggplot(data=a4)+
  geom_point(aes(x=sampling,y=abund,color=treatment))+
  geom_smooth(aes(x=sampling,y=abund,color=treatment),method = "lm")

# overall algal abundance increases over time

# try a global model with taxa as a random effect.

a3.1<-a3%>%
  filter(sampling==1)%>%
  rename(start.abund=abundance)%>%
  select(treatment, plot,taxa,start.abund)

a5<-a3 %>%
  filter(sampling!=1)%>%
  left_join(a3.1)

a5$treatment<-factor(a5$treatment)


alm1<-glmer(abundance~treatment*sampling+ (1|plot)+(treatment:sampling|taxa),
           data=a5,offset=log(start.abund+1),family=poisson)
alm1.r<-glmer(abundance~treatment*sampling+ (1|plot)+(sampling|taxa),
            data=a5%>%
              mutate(treatment=relevel(treatment, ref = "real")),
            offset=log(start.abund+1),family=poisson)
alm1.f<-glmer(abundance~treatment*sampling+ (1|plot)+(sampling|taxa),
            data=a5%>%
              mutate(treatment=relevel(treatment, ref = "fake")),
            offset=log(start.abund+1),family=poisson)
summary(alm1)
summary(alm1.r)
ranef(alm1)

# overall algal abundance increases over time
#halimeda
ahm1<-glmer(abundance~treatment*sampling+ (1|plot),
            data=a5%>%
              filter(taxa=="halimeda"),
            offset=log(start.abund+1),family=poisson)
ahm1.r<-glmer(abundance~treatment*sampling+ (1|plot),
            data=a5%>%
              filter(taxa=="halimeda")%>%
              mutate(treatment=relevel(treatment, ref = "real")),
            offset=log(start.abund+1),family=poisson)
ahm1.f<-glmer(abundance~treatment*sampling+ (1|plot),
            data=a5%>%
              filter(taxa=="halimeda")%>%
              mutate(treatment=relevel(treatment, ref = "fake")),
            offset=log(start.abund+1),family=poisson)
summary(ahm1)
summary(ahm1.r)
summary(ahm1.f)

# acetabularia
aam1<-glmer(abundance~treatment*sampling+ (1|plot),
            data=a5%>%
              filter(taxa=="acetabularia"),
            offset=log(start.abund+1),family=poisson)
aam1.r<-glmer(abundance~treatment*sampling+ (1|plot),
              data=a5%>%
                filter(taxa=="acetabularia")%>%
                mutate(treatment=relevel(treatment, ref = "real")),
              offset=log(start.abund+1),family=poisson)
aam1.f<-glmer(abundance~treatment*sampling+ (1|plot),
              data=a5%>%
                filter(taxa=="acetabularia")%>%
                mutate(treatment=relevel(treatment, ref = "fake")),
              offset=log(start.abund+1),family=poisson)
summary(aam1)
summary(aam1.r)
summary(aam1.f)

# laurencia
alm1<-glmer(abundance~treatment*sampling+ (1|plot),
            data=a5%>%
              filter(taxa=="laurencia"),
            offset=log(start.abund+1),family=poisson)
alm1.r<-glmer(abundance~treatment*sampling+ (1|plot),
              data=a5%>%
                filter(taxa=="laurencia")%>%
                mutate(treatment=relevel(treatment, ref = "real")),
              offset=log(start.abund+1),family=poisson)
alm1.f<-glmer(abundance~treatment*sampling+ (1|plot),
              data=a5%>%
                filter(taxa=="laurencia")%>%
                mutate(treatment=relevel(treatment, ref = "fake")),
              offset=log(start.abund+1),family=poisson)
summary(alm1)
summary(alm1.r)
summary(alm1.f)

