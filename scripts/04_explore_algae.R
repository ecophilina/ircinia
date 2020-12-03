library(tidyverse)

source("scripts/03_reimport.R")#imports all the data sets
if(!require(lmerTest))install.packages("lmerTest");library(lmerTest)
if(!require(DHARMa))install.packages("DHARMa");library(DHARMa)
if(!require(glmmTMB))install.packages("glmmTMB");library(glmmTMB)


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


# restarting all of this trying out the appropriate model

# Philina's custom function for checking for overdispersion
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

# check distributions of things
hist(a3$total, breaks = 20) # total algal abundance
hist((a3%>%
       filter(taxa=="halimeda"))$abundance, breaks = 20)
hist((a3%>%
        filter(taxa=="laurencia"))$abundance, breaks = 20)
hist((a3%>%
        filter(taxa=="acetabularia"))$abundance, breaks = 20)
hist((a3%>%
        filter(taxa=="udotea"))$abundance, breaks = 20)
hist((a3%>%
        filter(taxa=="penicillus"))$abundance, breaks = 20)

# try out more appropriate model

# trying the one global model first
alm1<-glmmTMB(abundance~treatment*sampling+ (1|plot)+(treatment:sampling|taxa),
            data=a5,
            offset=log(start.abund+1),
            REML=F,
            family=poisson)
#this model doesn't converge no matter the family you choose using glmmTMB

# going with species-specific models. 
#First I'm going to try a model with only total abundance of algae

a6<-a3 %>%
  select(treatment, plot,sampling,total,season)%>%
  distinct()

a6.1<-a6%>%
  filter(sampling==1)%>%
  rename(start.total=total)%>%
  select(-sampling,-season)

a6<-a6 %>%
  filter(sampling!=1)%>%
  left_join(a6.1)%>%
  mutate(samp2=case_when(
    sampling==2~1,
    sampling==3~5,
    sampling==4~12,
    sampling==5~17))

a6$treatment<-factor(a6$treatment)

atot<-glmmTMB(total~treatment*samp2+(1|plot),
              data=a6%>%
                mutate(treatment=relevel(treatment, ref = "real")),
              offset=log(start.total),
              REML = FALSE,
              family=nbinom1)
# look at residuals
atot_simres <- simulateResiduals(atot)
testDispersion(atot_simres)
plot(atot_simres)
# residuals look good

summary(atot)

# after changing sampling to represent the actual # of months into experiment
# algae increase significantly in the real treatment over time, but there's
# no significant difference between treatments.

a5<-a5%>%
  mutate(samp2=case_when(
    sampling==2~1,
    sampling==3~5,
    sampling==4~12,
    sampling==5~17))

# now look at halimeda alone
hal<-glmmTMB(abundance~treatment*samp2+(1|plot),
              data=a5%>%
               filter(taxa=="halimeda")%>%
                mutate(treatment=relevel(treatment, ref = "real")),
              offset=log(start.abund),
              REML = FALSE,
              family=nbinom1)
# look at residuals
hal_simres <- simulateResiduals(hal)
testDispersion(hal_simres)
plot(hal_simres)
# residuals look good

summary(hal)
# nothing significant
# note gaussian residuals look good too, and sampling is almost sig. but not there

# look at acetabularia- might not be enough there for models to converge

ace<-glmmTMB(abundance~treatment*samp2+(1|plot),
             data=a5%>%
               filter(taxa=="acetabularia")%>%
               mutate(treatment=relevel(treatment, ref = "real")),
             offset=log(start.abund+1),
             REML = FALSE,
             family=nbinom1)
# look at residuals
ace_simres <- simulateResiduals(ace)
testDispersion(ace_simres)
plot(ace_simres)
# residuals look good

summary(ace)
# model did converge- significant increase over time in real, but not different
# than other two treatments. But note that change over time is not significant in other two

# now look at Laurencia
lau<-glmmTMB(abundance~treatment*samp2+(1|plot),
             data=a5%>%
               filter(taxa=="laurencia")%>%
               mutate(treatment=relevel(treatment, ref = "real")),
             offset=log(start.abund+1),
             REML = FALSE,
             family=nbinom1())
# look at residuals
lau_simres <- simulateResiduals(lau)
testDispersion(lau_simres)
plot(lau_simres)
# residuals look good

summary(lau)
# no significant change over time and no difference between treatments for real.
# note that fake and blank do significantly decrease over time

# now penicillus

pen<-glmmTMB(abundance~treatment*samp2+(1|plot),
             data=a5%>%
               filter(taxa=="penicillus")%>%
               mutate(treatment=relevel(treatment, ref = "real")),
             offset=log(start.abund+1),
             REML = FALSE,
             family=poisson)
#nbinom1 didn't converge
# look at residuals
pen_simres <- simulateResiduals(pen)
testDispersion(pen_simres)
plot(pen_simres)
# residuals look good

summary(pen)
# penicillus increases over time, no significant differences with other treatments
# but note: model doesn't converge with blank or fake as the reference level and there is
# ridiculous error associated with the fake sampling this is a result of 
# penicilus mostly not occurring in the fake and blank... perhaps a better 
# approach for this taxa is just to see if it increases in real plots over time?

# finally udotea

udot<-glmmTMB(abundance~treatment*samp2+(1|plot),
             data=a5%>%
               filter(taxa=="udotea")%>%
               mutate(treatment=relevel(treatment, ref = "real")),
             offset=log(start.abund+1),
             REML = FALSE,
             family=nbinom1)
# look at residuals
udot_simres <- simulateResiduals(udot)
testDispersion(udot_simres)
plot(udot_simres)
# residuals look good

summary(udot)
# nothing significant going on


# my thoughts- since the by species info isnt that much more interesting than 
# the total algae story we just look at total algae and live with not having 
# a totally significant result for one set of primary producers.

# reminder what the total algae results look like

summary(atot)

# potential figure

a6<-a6%>%
  mutate(delta.tot=total-start.total)
a7<-a6%>%
  group_by(treatment,samp2)%>%
  summarize(m.tot=mean(delta.tot),sd.tot=sd(delta.tot))
ylab<-expression(paste(Delta," macroalgae abundance"))

(afig <- ggplot() +
    geom_hline(aes(yintercept = 0)) +
    geom_point(
      data = a6,
      aes(x = samp2, y = delta.tot, color = treatment),
      position = position_dodge(1.8)) +
    geom_point(data = a7,
      aes(x = samp2, y = m.tot, color = treatment),
      size = 8,
      position = position_dodge(+1.8),
      alpha = .5) +
    geom_errorbar(
      data = a7,
      aes(x = samp2,ymin = m.tot - sd.tot,ymax = m.tot + sd.tot,color = treatment),
      width = 1,
      position = position_dodge(1.8)) +
    theme_bw() +
    xlab("Months into experiment") +
    ylab(ylab) +
    scale_color_brewer(type = "qual",
      palette = "Set2",
      name = "",
      labels = c("Control", "Structure control", "Sponge")) +
    theme_bw() +
    theme(panel.grid = element_blank(),
      panel.spacing = unit(0, "lines"),
      strip.background = element_blank(),
      legend.position = "top",
      legend.text = element_text(size = 14),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      strip.text = element_text(size = 14)))
ggsave("figures/algal_abundance.jpg")
