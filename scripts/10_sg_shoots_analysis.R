# How does the sponge Ircinia felix influence seagrass bed primary producers? #

# this script will examine seagrass shoot density #
# run 03_reimport.R fist #

# load data and necessary packages----
source("scripts/03_reimport.R")#imports all the data sets
if(!require(lmerTest))install.packages("lmerTest");library(lmerTest)
#make things into factors
sgsd0<-sg_shoot%>%
  mutate(dist2=factor(dist,
                      level=c("0-1","1-2","2-3"),labels = c("near","mid","far")),
         treatment=factor(treatment))
#calcualte the mean of the pre-experiment sampling to add as an offset in analysis
sgsd1<-sgsd0%>%
  filter(sampling==1)%>%
  group_by(treatment,plot,dist2)%>%
  mutate(mtsd=mean(T.SD),mshsd=mean(SH.SD),msd=mean(SD))%>%
  ungroup()
# look at differences in the beginning
sd1lm<-lmer(SD~treatment*dist2+(1|plot),data=sgsd1)
sd1aov<-anova(sd1lm)

tsd1lm<-lmer(T.SD~treatment*dist2+(1|plot),data=sgsd1)
anova(tsd1lm)
shsd1lm<-lmer(SH.SD~treatment*dist2+(1|plot),data=sgsd1)
anova(shsd1lm)

sgsd1<-sgsd1%>%
  select(-sampling,-T.SD,-SH.SD,-SD)%>%
  distinct()


#join back to main dataset
sgsd<-sgsd0%>%
  filter(sampling!=1)%>%
  left_join(sgsd1)%>%
  mutate(yr=case_when(
    sampling==2~1,
    sampling==3~1,
    sampling==4~2,
    sampling==5~2),
    season=case_when(
      sampling==2~"summer",
      sampling==3~"winter",
      sampling==4~"summer",
      sampling==5~"winter"))

#looking at overall shoot densities
ggplot(data=sgsd,aes(y=SD-msd,x=sampling))+
  geom_jitter(aes(color=treatment))+
  geom_smooth()+
  geom_hline(yintercept = 0)+
  facet_grid(treatment~dist2)

sdlmb<-lmer(SD~treatment * yr+treatment*season+(1|plot/dist2),
           offset=msd,
           data=sgsd)#%>%filter(dist2!="near"))

sdb.sum<-summary(sdlmb)
sd.aov<-anova(sdlmb)

sdlmf<-lmer(SD~treatment * yr+treatment*season+(1|plot/dist2),
            offset=msd,
            data=sgsd%>%
              #filter(dist2!="near")%>%
              mutate(treatment=relevel(treatment, ref = "fake")))

sdf.sum<-summary(sdlmf)


sdlmr<-lmer(SD~treatment * yr+treatment*season+(1|plot/dist2),
            offset=msd,
            data=sgsd%>%
              #filter(dist2!="near")%>%
              mutate(treatment=relevel(treatment, ref = "real")))

sdr.sum<-summary(sdlmr)


#looking at Thalassia shoot densities
ggplot(data=sgsd,aes(y=T.SD-mtsd,x=sampling))+
  geom_jitter(aes(color=treatment))+
  geom_smooth()+
  geom_hline(yintercept = 0)+
  facet_grid(treatment~dist2)

tsdlmb<-lmer(T.SD~treatment * yr+treatment*season+(1|plot/dist2),
            offset=mtsd,
            data=sgsd)#%>%
            #filter(dist2!="near"))

tsdb.sum<-summary(tsdlmb)
tsd.aov<-anova(tsdlmb)

tsdlmf<-lmer(T.SD~treatment * yr+treatment*season+(1|plot/dist2),
            offset=mtsd,
            data=sgsd%>%
              #filter(dist2!="near")%>%
              mutate(treatment=relevel(treatment, ref = "fake")))

tsdf.sum<-summary(tsdlmf)


tsdlmr<-lmer(T.SD~treatment * yr+treatment*season+(1|plot/dist2),
            offset=mtsd,
            data=sgsd%>%
 #             filter(dist2!="near")%>%
              mutate(treatment=relevel(treatment, ref = "real")))

tsdr.sum<-summary(tsdlmr)

#looking at syringodium/halodule shoot densities
ggplot(data=sgsd,aes(y=SH.SD-mshsd,x=sampling))+
  geom_jitter(aes(color=as.factor(plot)))+
  geom_smooth()+
  geom_hline(yintercept = 0)+
  facet_grid(treatment~dist2)

shsdlmb<-lmer(SH.SD~treatment * yr+treatment*season+(1|plot/dist2),
             offset=mshsd,
             data=sgsd)#%>%
               #filter(dist2!="near"))

shsdb.sum<-summary(shsdlmb)
shsd.aov<-anova(shsdlmb)

shsdlmf<-lmer(SH.SD~treatment * yr+treatment*season+(1|plot/dist2),
             offset=mshsd,
             data=sgsd%>%
#               filter(dist2!="near")%>%
               mutate(treatment=relevel(treatment, ref = "fake")))

shsdf.sum<-summary(shsdlmf)


shsdlmr<-lmer(SH.SD~treatment * yr+treatment*season+(1|plot/dist2),
             offset=mshsd,
             data=sgsd%>%
               filter(dist2!="near")%>%
               mutate(treatment=relevel(treatment, ref = "real")))

shsdr.sum<-summary(shsdlmr)


