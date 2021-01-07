# How does the sponge Ircinia felix influence seagrass bed primary producers? #

# this script will examine seagrass shoot density #
# run 03_reimport.R fist #

# load data and necessary packages----
source("scripts/03_reimport.R")#imports all the data sets

if(!require(DHARMa))install.packages("DHARMa");library(DHARMa)
if(!require(lmerTest))install.packages("lmerTest");library(lmerTest)
if(!require(glmmTMB))install.packages("glmmTMB");library(glmmTMB)

# custom function for checking for overdispersion
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}


# make things into factors and create a continuous distance variable
sgsd0<-sg_shoot%>%
  mutate(dist2=factor(dist,
                      level=c("0-1","1-2","2-3"),labels = c("near","mid","far")),
         dist1=case_when(dist2=="near"~0,
           dist2=="mid"~1,
           dist2=="far"~2),
         treatment=factor(treatment))

#calcualte the mean of the pre-experiment sampling to add as an offset in analysis
sgsd1<-sgsd0%>%
  filter(sampling==1)%>%
  group_by(treatment,plot,dist2)%>%
  mutate(mtsd=mean(T.SD),mshsd=mean(SH.SD),msd=mean(SD))%>%
  ungroup()
# look at differences in the beginning

# check distributions
# raw counts
hist(sgsd1$SD, breaks = 20)
hist(sgsd1$T.SD, breaks = 20)
hist(sgsd1$SH.SD, breaks = 20)

# change in counts
hist(sgsd1$SD-sgsd1$msd, breaks = 20)
hist(sgsd1$T.SD-sgsd1$mtsd, breaks = 20)
hist(sgsd1$SH.SD-sgsd1$mshsd, breaks = 20)

# some look normalish so try a range of model types


### test for initial differences between treatments and distance as factor
# residuals ok for first but terrible for last...
sd1lm<-lmer(SD~treatment*dist2+(1|plot), data=sgsd1) 
tsd1lm<-lmer(T.SD~treatment*dist2+(1|plot), data=sgsd1) 
shsd1lm<-lmer(SH.SD~treatment*dist2+(1|plot), data=sgsd1) 
plot(sd1lm)
plot(tsd1lm)
plot(shsd1lm)

shsd1lm<-glmer(SH.SD~treatment*dist2+(1|plot), data=sgsd1, family = poisson) # overdispersed
shsd1lm<-glmer.nb(SH.SD~treatment*dist2+(1|plot), data=sgsd1, family = quasipoisson) # not converging
overdisp_fun(shsd1lm)

# this model both solves overdispersion and converges for all 
sd1glm<-glmmTMB(SD~treatment*dist2+(1|plot), data=sgsd1, family = nbinom1) 
tsd1glm<-glmmTMB(T.SD~treatment*dist2+(1|plot), data=sgsd1, family = nbinom1)
shsd1glm<-glmmTMB(SH.SD~treatment*dist2+(1|plot), data=sgsd1, family = nbinom1)

# # residual checks look great
tsd1glm_simres <- simulateResiduals(tsd1glm)
plot(tsd1glm_simres)
testOverdispersion(tsd1glm_simres)

sd1glm_simres <- simulateResiduals(shsd1glm)
plot(sd1glm_simres)
testOverdispersion(sd1glm_simres)


(sd1aov<-glmmTMB:::Anova.glmmTMB(sd1glm, type = "III"))
(tsd1aov<-glmmTMB:::Anova.glmmTMB(tsd1glm, type = "III"))
(shsd1aov<-glmmTMB:::Anova.glmmTMB(shsd1glm, type = "III"))

# (sd1sum<-summary(sd1glm))
# (tsd1sum<-summary(tsd1glm))
# (shsd1sum<-summary(shsd1glm))

# no significant differences between treatments or distances
# weak, non-sig lower T density at far distance in sponge treatments
# likely due to some plots being near the edge of the seagrass bed
# similar non-sig lower SH density in sponge treatments

# same pattern for continuous distance, so will use that in models
tsd1lm<-glmmTMB(T.SD~treatment*dist1+(1|plot), data=sgsd1, family = nbinom1) 
(tsd1aov<-glmmTMB:::Anova.glmmTMB(tsd1lm, type = "III"))
(tsd1aov<-summary(tsd1lm))

shsd1lm<-glmmTMB(SH.SD~treatment*dist1+(1|plot),data=sgsd1, family = nbinom1)
(shsd1aov<-glmmTMB:::Anova.glmmTMB(shsd1lm, type = "III"))
(shsd1aov<-summary(shsd1lm))


# #looking at overall shoot densities
# ggplot(data=sgsd,aes(y=SD-msd,x=samp2))+
#   geom_jitter(aes(color=season))+
#   geom_smooth(method = "lm", colour = "black")+
#   geom_hline(yintercept = 0)+
#   facet_grid(treatment~dist2)
# 
# #looking at overall shoot densities
# ggplot(data=sgsd%>% filter(dist2=="near"),
#   aes(y=SD-msd,x=samp2))+
#   geom_jitter(aes(color=season))+
#   geom_smooth(method = "lm", colour = "darkgrey")+
#   geom_hline(yintercept = 0)+
#   facet_grid(rows=vars(treatment))
 

#looking at Thalassia shoot density change
ggplot(data=sgsd,aes(y=T.SD-mtsd,x=samp2))+
  # geom_jitter(aes(color=treatment))+
  geom_jitter(aes(color=season))+
  geom_smooth(method = "lm", colour = "black")+
  geom_hline(yintercept = 0)+
  facet_grid(treatment~dist2)
# facet_grid(row=vars(treatment))
# appears to decline everywhere, but maybe less so in real near

#looking at syringodium/halodule shoot density change
ggplot(data=sgsd,aes(y=SH.SD-mshsd,x=samp2))+
  # geom_jitter(aes(color=treatment))+
  geom_jitter(aes(color=season))+
  geom_smooth(method = "lm", colour = "black")+
  geom_hline(yintercept = 0)+
  facet_grid(treatment~dist2)
# increaseing for real only, especially in near

# actual syringodium/halodule shoot densities
ggplot(data=sgsd,aes(y=SH.SD,x=samp2))+
  # geom_jitter(aes(color=treatment))+
  geom_boxplot(aes(group = as.factor(samp2), color=season))+
  geom_jitter(aes(color=season))+
  geom_smooth(method = "lm", colour = "black")+
  geom_hline(yintercept = 0)+
  facet_grid(treatment~dist2)



# # overall shoot density
# # run the model with each treatment as the reference level to help with interpretation
# # for all groupings of seagrass shoot density the model we decided on was treatment interacting with months into the experiment and then treatment
# # interacting with distance and having plot as a random factor. We included months into the experiment as the time variable
# # because there is no a priori reason to think shoot density would vary with season.
# 
# sdlmb1<-glmmTMB(SD ~ treatment*samp2 + treatment*dist1 + samp2*dist1 +
#     (1|plot),
#   offset=log(msd),
#   data=sgsd %>%
#     mutate(treatment=relevel(treatment, ref = "real")), 
#   REML = F,
#   family = nbinom1) 
# 
# sdlmb_simres <- simulateResiduals(sdlmb1) #, refit = TRUE
# testDispersion(sdlmb_simres) #, "greater"
# plot(sdlmb_simres)
# (sdlmbaov<-summary(sdlmb1))
# 
# # 
# # # because of offset resulting in negative and possitive changes in counts, gaussian actually works better for SD and T, but not for SH
# # sdlmb2<-glmmTMB(SD ~ treatment*samp2 + treatment*dist1 +
# #     (dist1|plot),
# #   offset=msd,
# #   data=sgsd %>%
# #     mutate(treatment=relevel(treatment, ref = "real")),
# #   REML = F,
# #   family = gaussian)
# # 
# # sdlmb_simres2 <- simulateResiduals(sdlmb2) #, refit = TRUE
# # testDispersion(sdlmb_simres2) #, "greater"
# # plot(sdlmb_simres2)
# # (sdlmbaov<-summary(sdlmb2))
# 
# # over time, there is a decrease in overall shoot density in control and structure control plots (only significant in struct. control)
# # there is a non-significant increase in overall shoot density in sponge plots. While the increase is not significant the 
# # pattern is significantly different than both the control and structure control.


# Why not just look at shoot density at the end of the experiment?

range(sgsd$mtsd) # no zeros means log(mtsd) could work
range(sgsd$mshsd) # but zeros means log(mshsd +1) needed

sgsd.end<-sgsd%>%
  filter(sampling==5)

tsd.endb<-glmmTMB(T.SD~treatment*dist1+(1|plot),
  offset=log(mtsd+1),
  data=sgsd.end%>% mutate(treatment=relevel(treatment, ref = "real")), 
  REML = F,
  family = nbinom1)

(tsd.end.aov<-glmmTMB:::Anova.glmmTMB(tsd.endb, type = "III"))
summary(tsd.endb)
# no treatment effects for T

shsd.endb<-glmmTMB(SH.SD~treatment*dist1+(1|plot),
  offset=(mshsd+1),
  data=sgsd.end%>% mutate(treatment=relevel(treatment, ref = "real")), 
  REML = F,
  family = nbinom1)

(shsd.end.aov<-glmmTMB:::Anova.glmmTMB(shsd.endb, type = "III"))
summary(shsd.endb)
# treatment interacts with dist for SH



# join back to main dataset to create the dataset for analysis.
sgsd1<-sgsd1%>%
  select(-sampling,-T.SD,-SH.SD,-SD)%>%
  distinct()

sgsd<-sgsd0%>%
  filter(sampling!=1)%>%
  left_join(sgsd1)%>%
  # filter(dist2=="near")%>%
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
      sampling==5~17),
    delta = (SD-msd),
    ratio = delta/msd
  )

# confirm no season effect
tsdseason<-glmmTMB(T.SD~season+yr+(1|plot), data=sgsd, family = nbinom1) 
(tsdseason<-glmmTMB:::Anova.glmmTMB(tsdseason, type = "III"))
# (tsdseason2<-summary(tsdseason))

shsdseason<-glmmTMB(SH.SD~season+yr+(1|plot),data=sgsd, family = nbinom1)
(shsdseason<-glmmTMB:::Anova.glmmTMB(shsdseason, type = "III"))
# (shsdseason2<-summary(shsdseason))

############
# FINAL MODELS

tsdlmr<-glmmTMB(T.SD ~ 
    treatment*samp2*dist1 +
    # treatment*samp2 + treatment*dist1 + samp2*dist1 + # still no sig treatment effects
    (1|plot),
  offset=log(mtsd+1),
  data=sgsd %>% mutate(treatment=relevel(treatment, ref = "real")), 
  REML = F, 
  family = nbinom1)
  # family = gaussian) # improves residuals but still no sig. treatment effects

(tsd.aov<-glmmTMB:::Anova.glmmTMB(tsdlmr, type = "III"))
(tsdr.sum<-summary(tsdlmr))

# check residuals: not amazing, but not terrible
tsdlmr_simres <- simulateResiduals(tsdlmr) 
testDispersion(tsdlmr_simres) 
plot(tsdlmr_simres)


# other intercepts for getting coefs
tsdlmb<-glmmTMB(T.SD ~ treatment*samp2*dist1 +
    (1|plot),
  offset=log(mtsd+1),
  data=sgsd, 
  REML = F, family = nbinom1) 
(tsdb.sum<-summary(tsdlmb))

tsdlmf<-glmmTMB(T.SD ~ treatment*samp2*dist1 +
    (1|plot),
  offset=log(mtsd+1),
  data=sgsd %>% mutate(treatment=relevel(treatment, ref = "fake")), 
  REML = F, family = nbinom1) 
(tsdf.sum<-summary(tsdlmf))


# conclusion here is still not treatment effect but there is some evidence of an overall loss in Thalassia shoot density over time
# could simplify if we need overall estimate of Thalassia density change
tsdlm<-glmmTMB(T.SD ~ samp2  +
    # treatment*dist1 + 
    # (1|plot,
     (1|plot/dist2),
  offset=log(mtsd+1),
  data=sgsd, 
  REML = F, 
  family = nbinom1)
# family = gaussian) # improves residuals but still no sig. treatment effects

(tsd.sum<-summary(tsdlm))

# tsdlm_simres <- simulateResiduals(tsdlm)
# testDispersion(tsdlm_simres)
# plot(tsdlm_simres)



# Knowing that Syringodium and Halodule (SH) differ by treatment by end of experiment. 
# How does this develop through time and with distance from sponge?

shsdlmr<-glmmTMB(SH.SD ~ 
    treatment * samp2 * dist1 + # squared term never sig
    # treatment * poly(samp2,2) * dist1 + # squared term never sig
    # treatment*samp2 + treatment*dist1 + #samp2*dist1 +
    (1|plot),
  offset=log(mshsd+1),
  data=sgsd %>%
    mutate(treatment=relevel(treatment, ref = "real")), 
  REML = F,
  family = nbinom1) 

(shsd.aov<-glmmTMB:::Anova.glmmTMB(shsdlmr, type = "III"))
(shsdr.sum<-summary(shsdlmr))

shsdlmr_simres <- simulateResiduals(shsdlmr) 
testDispersion(shsdlmr_simres) 
plot(shsdlmr_simres)

# other intercepts for getting coefs
shsdlmb<-glmmTMB(SH.SD ~ 
    treatment*samp2*dist1 +
    (1|plot),
  offset=log(mshsd+1),
  data=sgsd, REML = F, family = nbinom1) 
(shsdb.sum<-summary(shsdlmb))

shsdlmf<-glmmTMB(SH.SD ~ 
    treatment*samp2*dist1 +
    # treatment*samp2 + treatment*dist1 + 
    # (dist1|plot),
    (1|plot),
  offset=log(mshsd+1),
  data=sgsd %>%
    mutate(treatment=relevel(treatment, ref = "fake")), 
  REML = F, family = nbinom1) 
(shsdf.sum<-summary(shsdlmf))



shsdr.sum
shsdb.sum
shsdf.sum

# Syringodium and Halodule (SH) increase significantly in sponge plots over time and this is significantly different than both the control and structure control. 
# There is now an effect of distance in sponge plots.
# Further plots start out with more SH, but show less of an increase through time.
# Control and structural control shows no change through time and both only differ sig from real.






