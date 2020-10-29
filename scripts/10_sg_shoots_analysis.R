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
(sd1aov<-anova(sd1lm))
# distance is significant at the start and this looks to be be because of the far distance. To remove pre-existing bias we're 
# going to see if removing the far distance makes us have no significant differences before the start of the experiment.
# this will also work out to be consistent with the growth analysis.
# this effect was likely do to some plots being near the edge of the seagrass bed.

# look at differences in the beginning without far distance
sd1lm2<-lmer(SD~treatment*dist2+(1|plot),data=sgsd1%>%filter(dist2!="far"))
(sd1aov2<-anova(sd1lm2))
# for overall shoot density this gets rid of the significant differences. Make sure this holds true for each spp group individually

tsd1lm<-lmer(T.SD~treatment*dist2+(1|plot),data=sgsd1%>%filter(dist2!="far"))
anova(tsd1lm)

shsd1lm<-lmer(SH.SD~treatment*dist2+(1|plot),data=sgsd1%>%filter(dist2!="far"))
anova(shsd1lm)


#join back to main dataset to create the dataset for analysis.

sgsd1<-sgsd1%>%
  select(-sampling,-T.SD,-SH.SD,-SD)%>%
  distinct()

sgsd<-sgsd0%>%
  filter(sampling!=1)%>%
  left_join(sgsd1)%>%
  filter(dist2!="far")%>%
  mutate(yr=case_when(
    sampling==2~1,
    sampling==3~1,
    sampling==4~2,
    sampling==5~2),
    season=case_when(
      sampling==2~"summer",
      sampling==3~"winter",
      sampling==4~"summer",
      sampling==5~"winter"),
    samp2=case_when(
      sampling==2~1,
      sampling==3~5,
      sampling==4~12,
      sampling==5~17))

#looking at overall shoot densities
ggplot(data=sgsd,aes(y=SD-msd,x=samp2))+
  geom_jitter(aes(color=treatment))+
  geom_smooth()+
  geom_hline(yintercept = 0)+
  facet_grid(treatment~dist2)

# now start analysis

# overall shoot density
# run the model with each treatment as the reference level to help with interpretation
# for all groupings of seagrass shoot density the model we decided on was treatment interacting with months into the experiment and then treatment
# interacting with distance and having plot as a random factor. We included months into the experiment as the time variable
# because there is no a priori reason to think shoot density would vary with season.

sdlmb<-lmer(SD~treatment*samp2 + treatment*dist2+(1|plot),
           offset=msd,
           data=sgsd)
sdlmf<-lmer(SD~treatment*samp2 + dist2*treatment+(1|plot),
            offset=msd,
            data=sgsd%>%
            mutate(treatment=relevel(treatment, ref = "fake")))
sdlmr<-lmer(SD~treatment*samp2 + dist2*treatment+(1|plot),
            offset=msd,
            data=sgsd%>%
              mutate(treatment=relevel(treatment, ref = "real")))
# now look at results, first overall anova then each linear model
(sd.aov<-anova(sdlmb))
(sdb.sum<-summary(sdlmb))
(sdf.sum<-summary(sdlmf))
(sdr.sum<-summary(sdlmr))

# over time, there is a decrease in overall shoot density in control and structure control plots (only significant in struct. control)
# there is a non-significant increase in overall shoot density in sponge plots. While the increase is not significant the 
# pattern is significantly different than both the control and structure control.

# Look at overall shoot density at the end of the experiment.
sgsd.end<-sgsd%>%
  filter(sampling==5)

sd.endb<-lmer(SD~treatment+(1|plot),
              offset=msd,
              data=sgsd.end%>%
                filter(dist2=="near"))

(sd.end.aov<-anova(sd.endb))
summary(sd.endb)
# not significant


#looking at Thalassia shoot densities
ggplot(data=sgsd,aes(y=T.SD-mtsd,x=samp2))+
  geom_jitter(aes(color=treatment))+
  geom_smooth()+
  geom_hline(yintercept = 0)+
  facet_grid(treatment~dist2)

# run the model with each treatment as the reference level to help with interpretation

tsdlmb<-lmer(T.SD~treatment * samp2+dist2*treatment+(1|plot),
            offset=mtsd,
            data=sgsd)
tsdlmf<-lmer(T.SD~treatment *  samp2+dist2*treatment+(1|plot),
            offset=mtsd,
            data=sgsd%>%
            mutate(treatment=relevel(treatment, ref = "fake")))
tsdlmr<-lmer(T.SD~treatment *  samp2+dist2*treatment+(1|plot),
            offset=mtsd,
            data=sgsd%>%
            mutate(treatment=relevel(treatment, ref = "real")))
# now look at results, first overall anova then each linear model
(tsd.aov<-anova(tsdlmb))
(tsdb.sum<-summary(tsdlmb))
(tsdf.sum<-summary(tsdlmf))
(tsdr.sum<-summary(tsdlmr))

# conclusion here is that all treatments lost Thalassia shoot density over time. However, this is not significant for sponge
# treatments. Conversely, there is an increase in Thalassia shoot density at mid-distances in all treatments, but its 
# only significant in the real sponge treatment. 

# # could simplify if we need overall estimate of Thalassia density change
# 
# tsdlmb<-lmer(T.SD~treatment + samp2 + dist2 + (1|plot),
#   offset=mtsd,
#   data=sgsd)
# (tsd.sum<-summary(tsdlmb))

#looking at syringodium/halodule shoot densities
ggplot(data=sgsd,aes(y=SH.SD-mshsd,x=samp2))+
  geom_jitter(aes(color=as.factor(plot)))+
  geom_smooth()+
  geom_hline(yintercept = 0)+
  facet_grid(treatment~dist2)
# run the model with each treatment as the reference level to help with interpretation
shsdlmb<-lmer(SH.SD~treatment *  samp2+dist2*treatment+(1|plot),
             offset=mshsd,
             data=sgsd)

shsdlmf<-lmer(SH.SD~treatment  *  samp2+dist2*treatment+(1|plot),
             offset=mshsd,
             data=sgsd%>%
             mutate(treatment=relevel(treatment, ref = "fake")))

shsdlmr<-lmer(SH.SD~treatment * samp2+dist2*treatment+(1|plot),
             offset=mshsd,
             data=sgsd%>%
             mutate(treatment=relevel(treatment, ref = "real")))

# now look at results, first overall anova then each linear model
(shsd.aov<-anova(shsdlmb))

(shsdb.sum<-summary(shsdlmb))
(shsdf.sum<-summary(shsdlmf))
(shsdr.sum<-summary(shsdlmr))

# Syringodium and Halodule (SH) increase slightly, but non-significantly over time in the control. There is a significant decrease in the mid
# distance in the control- this appears to be driven by a single plot.
# SH decrease slightly, but not significantly over time in structure control plots. There is no effect of distance
# in structure control plots
# There is a significant increase in SH in sponge plots over time and this is significantly different than
# both the control and structure control. There is no effect of distance in sponge plots.




