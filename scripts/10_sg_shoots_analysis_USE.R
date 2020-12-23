# Cleaned up seagrass shoot density analysis

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

#Thalassia shoot density analysis
tsdlmr<-glmmTMB(T.SD ~ 
                  treatment*samp2*dist1 +
                (1|plot),
                offset=log(mtsd+1),
                data=sgsd %>% mutate(treatment=relevel(treatment, ref = "real")), 
                REML = F, 
                family = nbinom1)

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


# Knowing that Syringodium and Halodule (SH) differ by treatment by end of experiment. 
# How does this develop through time and with distance from sponge?

shsdlmr<-glmmTMB(SH.SD ~ 
                   treatment * samp2 * dist1 + 
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



