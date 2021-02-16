# this script will explore the invert data

library(tidyverse)
if(!require(vegan))install.packages("vegan"); library(vegan)

# bring in data
source("scripts/03_reimport.R")  

# prep primary producer data
sg<-sg_shoot%>%
  group_by(treatment,plot,sampling)%>%
  summarize(sg.sd=mean(SD))%>%
  mutate(sampling=case_when(
    sampling==1~0,
    sampling==2~1,
    sampling==3~5,
    sampling==4~12,
    sampling==5~17))

alg<-algae%>%
  group_by(treatment,plot,sampling)%>%
  summarize(alg=sum(abundance))%>%
  mutate(sampling=case_when(
    sampling==1~0,
    sampling==2~1,
    sampling==3~5,
    sampling==4~12,
    sampling==5~17))

sggrow<-sg_grow%>%
  filter(dist %in% c(0,0.5))%>%
  group_by(treatment,plot,sampling)%>%
  summarize(grow=mean(total.growth.mm2/days))%>%
  mutate(sampling=case_when(
    sampling==1~0,
    sampling==2~1,
    sampling==3~5,
    sampling==4~12,
    sampling==5~17))


# inverts

i2<-inverts%>%
  pivot_wider(names_from=taxa,values_from=abundance,values_fill=0)

i2<-i2%>%
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
      sampling==17~2))


i0<-i2%>%filter(yr == 0)
i2<-i2%>%filter(yr != 0)

# Organize data
i.env<-i2 %>% select(treatment, plot, yr, sampling, season)%>%
  left_join(alg)%>%
  left_join(sg)%>%
  left_join((sggrow))
i.com<-i2 %>% select(-treatment, -plot, -yr, -sampling, -season)

i.com$dummy<-1

# pre-experiment data
i.env0<-i0 %>% select(treatment, plot, yr, sampling, season)%>%
  left_join(alg)%>%
  left_join(sg)%>%
  left_join((sggrow))
i.com0<-i0 %>% select(-treatment, -plot, -yr, -sampling, -season)

i.com0$dummy<-1

# use hellinger: square root of method to standardize species data
i.com.pa<-decostand(i.com,"pa")

i.com.hel<-decostand(i.com.pa,"hellinger")

i.pca<-rda(i.com.hel)

i.scores<-data.frame(scores(i.pca,1:3)$sites)%>%
  bind_cols(i.env)

ggplot(data=i.scores)+
  geom_jitter(aes(x=PC2,y=PC3,
                 color=treatment),size=2,alpha=.65, 
    width = 0.03, height = 0.03)+
  scale_color_viridis_d(option="B",end=.8)+
  facet_wrap(~sampling)

#start examining statistical relationship

# first look at initial pre-experiment
i.com.pa.0<-decostand(i.com0,"pa")
i.com.hel0<-decostand(i.com.pa.0,"hellinger")

i.com.pa<-decostand(i.com,"pa")
i.com.hel<-decostand(i.com.pa,"hellinger")

# look at just treatment at 0 and later
i0.rda.null<-rda(i.com.hel0~1)
i0.rda<-rda(i.com.hel0~treatment,data=i.env0)
RsquareAdj(i0.rda)
(invert.initial<-anova(i0.rda.null,i0.rda))

# treatment does not explain a significant amount of variance between plots initially

# now start building more and more complex models for experiment data
#add in interactions
trt<-as.factor(i.env$treatment)
seas<-i.env$season
samp<-i.env$sampling
yr<-as.factor(i.env$yr)
plts<-as.factor(i.env$plot)

# start with the simplest model that makes sense - there are two of these - season and year
# with an interaction with treatment. Look at these with and without plots

# first question is an interaction between samp* treatment better than intercept only
# first build intercept only
i.rda.null<-rda(i.com.hel~1)

# now make samp*treat matrix
tr.samp.mat<-data.frame(model.matrix(~ samp*trt + plts, 
                                       contrasts=list(trt="contr.helmert", seas="contr.helmert")))[,-1]
i.rda.samp.treat.plot<-rda(i.com.hel~.,data=tr.samp.mat)

# now check to see if its better than intercept
anova(i.rda.null,i.rda.samp.treat.plot)

# it is better than null
# does including plot help
tr.samp.mat.np<-data.frame(model.matrix(~ samp*trt, 
                                     contrasts=list(trt="contr.helmert", seas="contr.helmert")))[,-1]
i.rda.samp.treat<-rda(i.com.hel~.,data=tr.samp.mat.np)

# check to see if plot should be included here
anova(i.rda.samp.treat.plot,i.rda.samp.treat)
# no difference as far as anova - look at rsquare
RsquareAdj(i.rda.samp.treat.plot)$adj.r.squared
RsquareAdj(i.rda.samp.treat)$adj.r.squared

# r squared better without plot
# current best model is just samp*treat

# now look at next simplest "base" model - this is treatment*year*season

tr.s.yr.mat<-data.frame(model.matrix(~ yr*seas*trt + plts, 
  contrasts=list(trt="contr.helmert", seas="contr.helmert",yr="contr.helmert")))[,-1]
tr.s.yr.mat.nop<-data.frame(model.matrix(~ yr*seas*trt, 
  contrasts=list(trt="contr.helmert", seas="contr.helmert",yr="contr.helmert")))[,-1]

# look at whether or not model without plot is better than current best model
i.rda.yr.s.treat<-rda(i.com.hel~.,tr.s.yr.mat.nop)

# models aren't nested so have to look at adjusted R2
RsquareAdj(i.rda.yr.s.treat)$adj.r.squared
RsquareAdj(i.rda.samp.treat)$adj.r.squared

# now yr*seas*treat is best model - see if plots make a difference here
i.rda.yr.s.treat.plot<-rda(i.com.hel~.,tr.s.yr.mat)
anova(i.rda.yr.s.treat,i.rda.yr.s.treat.plot)

# plot doesn't improve things here - how about for r2
RsquareAdj(i.rda.yr.s.treat)$adj.r.squared
RsquareAdj(i.rda.yr.s.treat.plot)$adj.r.squared

#tiny bit better but could arguably still leave plot out.

# now how about a series of two-way interactions vs the three

tr.s.yr.mat2<-data.frame(model.matrix(~ yr*seas + yr*trt + seas*trt +  plts, 
  contrasts=list(trt="contr.helmert", seas="contr.helmert", yr="contr.helmert")))[,-1]
tr.s.yr.mat2.nop<-data.frame(model.matrix(~ yr*seas + yr*trt + seas*trt, # without plot
  contrasts=list(trt="contr.helmert", seas="contr.helmert", yr="contr.helmert")))[,-1]

# look to see if 3-way interaction improves matters
i.rda.yr.s.treat2<-rda(i.com.hel~.,tr.s.yr.mat2.nop)
anova(i.rda.yr.s.treat,i.rda.yr.s.treat2)

# model with the 3-way interaction is better - confirm with r2
RsquareAdj(i.rda.yr.s.treat)$adj.r.squared
RsquareAdj(i.rda.yr.s.treat2)$adj.r.squared

# is season important
tr.yr.mat<-data.frame(model.matrix(~ yr*trt, # without plot
  contrasts=list(trt="contr.helmert",yr="contr.helmert")))[,-1]

i.rda.yr.treat<-rda(i.com.hel~.,tr.yr.mat)

anova(i.rda.yr.s.treat,i.rda.yr.treat)

# sampling is important. 
# now what about year
tr.seas.mat<-data.frame(model.matrix(~ seas*trt, # without plot
  contrasts=list(trt="contr.helmert",seas="contr.helmert")))[,-1]

i.rda.seas.treat<-rda(i.com.hel~.,tr.seas.mat)
anova(i.rda.yr.s.treat,i.rda.seas.treat)

#yep. Now finally treatmemnt
yr.seas.mat<-data.frame(model.matrix(~ yr*seas, # without plot
   contrasts=list(yr="contr.helmert",seas="contr.helmert")))[,-1]

i.rda.seas.yr<-rda(i.com.hel~.,yr.seas.mat)

anova(i.rda.yr.s.treat,i.rda.seas.yr)

#yep treatment is important

# now look into whether or not productivity measures explain community patterns better

i.env.prod<-i.env[,6:8]
i.rda.prod<-rda(i.com.hel~.,i.env.prod)

# does this model do better than the null
anova(i.rda.null,i.rda.prod)

# yes it does
# now does it do better than the treatment model
RsquareAdj(i.rda.yr.s.treat)$adj.r.squared
RsquareAdj(i.rda.prod)$adj.r.squared

# not on its own, no. What if we include different measures of productivity in our treatment model
# from now on I'm referring to the year*season*treatment as i.rda.best
i.rda.best<-rda(i.com.hel~.,tr.s.yr.mat.nop)
alg.e<-i.env$alg
sg<-i.env$sg.sd
sggrow<-i.env$grow

b.alg.mat<-data.frame(model.matrix(~ yr*seas*trt+alg.e, 
                                   contrasts=list(trt="contr.helmert", seas="contr.helmert",yr="contr.helmert")))[,-1]
b.sg.mat<-data.frame(model.matrix(~ yr*seas*trt+sg, 
                                   contrasts=list(trt="contr.helmert", seas="contr.helmert",yr="contr.helmert")))[,-1]
b.sgg.mat<-data.frame(model.matrix(~ yr*seas*trt+sggrow, 
                                   contrasts=list(trt="contr.helmert", seas="contr.helmert",yr="contr.helmert")))[,-1]
b.alg.sg.mat<-data.frame(model.matrix(~ yr*seas*trt+alg.e+sg, 
                                   contrasts=list(trt="contr.helmert", seas="contr.helmert",yr="contr.helmert")))[,-1]
b.alg.sgg.mat<-data.frame(model.matrix(~ yr*seas*trt+alg.e+sggrow, 
                                   contrasts=list(trt="contr.helmert", seas="contr.helmert",yr="contr.helmert")))[,-1]
b.sg.sgg.mat<-data.frame(model.matrix(~ yr*seas*trt+sg+sggrow, 
                                   contrasts=list(trt="contr.helmert", seas="contr.helmert",yr="contr.helmert")))[,-1]
b.alg.sg.sgg.mat<-data.frame(model.matrix(~ yr*seas*trt+alg.e+sg+sggrow, 
                                   contrasts=list(trt="contr.helmert", seas="contr.helmert",yr="contr.helmert")))[,-1]

# start with all of them added
i.rda.bprod<-rda(i.com.hel~.,b.alg.sg.sgg.mat)

anova(i.rda.best,i.rda.bprod)

# no difference between models - look at r2
RsquareAdj(i.rda.best)$adj.r.squared
RsquareAdj(i.rda.bprod)$adj.r.squared

# ever so slight increase in adjusted R2
# look at just sg growth

i.rda.bgrow<-rda(i.com.hel~.,b.sgg.mat)

anova(i.rda.best,i.rda.bgrow)

# no difference between models - look at r2
RsquareAdj(i.rda.best)$adj.r.squared
RsquareAdj(i.rda.bgrow)$adj.r.squared

# ever so slight increase in adjusted R2
# look at sggrow and algae
i.rda.balggrow<-rda(i.com.hel~.,b.alg.sgg.mat)

anova(i.rda.best,i.rda.balggrow)

# no difference between models - look at r2
RsquareAdj(i.rda.best)$adj.r.squared
RsquareAdj(i.rda.balggrow)$adj.r.squared

# ever so slight increase in adjusted R2
# look at sggrow and sg
i.rda.bsggrow<-rda(i.com.hel~.,b.sg.sgg.mat)

anova(i.rda.best,i.rda.bsggrow)

# no difference between models - look at r2
RsquareAdj(i.rda.best)$adj.r.squared
RsquareAdj(i.rda.bsggrow)$adj.r.squared
RsquareAdj(i.rda.bprod)$adj.r.squared

# including sg density gets us the best r2 so far - so looking at just sg density

i.rda.bsg<-rda(i.com.hel~.,b.sg.mat)

anova(i.rda.best,i.rda.bsg)

# no difference between models - look at r2
RsquareAdj(i.rda.best)$adj.r.squared
RsquareAdj(i.rda.bsg)$adj.r.squared

# best model at the moment: is seas*samp*year + sg + sggrow, but this model is 
#not significantly better than seas*samp*year

#best model at the moment
plot(i.rda.bsggrow, scaling = 3, display = c("sp", "cn"))

# note the order matters for adonis2 with by="terms" and the by="margin" doesn't seem to work
i.mod <- adonis2(i.com.hel~., data = b.sg.sgg.mat, method="euclidean", by="terms") 
i.mod

#this is where I stopped. I think the moral of the story is that the interaction between 
# treatment, season, and year has an affect on the community
# now we have to figure out what that is.

tr.s.yr3<-data.frame(model.matrix(~ yr*seas*trt+plts, 
  contrasts=list(trt="contr.helmert", seas="contr.helmert", yr="contr.helmert")))[,-1]

i.mod.best <- adonis2(i.com.hel~., data=tr.s.yr3, method="euclidean", sqrt.dist =T) 
i.mod.best

#bring in productivity data
i.env1<-left_join(i.env,sg)%>%
  left_join(alg)

#check to make sure data is still in the same order
summary(i.env[,1:4]==i.env1[,1:4])

# I'm worried that the helmert contrasts don't actually test difference between fake and real
# so attempt this part without blank data so that effect of cage applies to all plots
i.env2a <- i.env1 %>% filter(treatment!="blank")

# bind back the community data for just the retained plots
i3 <- left_join(i.env2a, i2)

# separate env and com again to sure order matches in this new dataset 
i.env2 <- i3 %>% select(treatment, plot, yr, sampling, season, sg.sd, alg)
i.com2 <- i3 %>% select(-treatment, -plot, -yr, -sampling, -season, -sg.sd, -alg)


i.com2$dummy<-1
com.dist2<-vegdist(i.com2,"bray")
i.mds2<-metaMDS(com.dist2,trymax = 100)
# plot(i.mds2)

# rda
i.com.hel2<-decostand(i.com2,"hellinger")


sg.sd<-i.env2$sg.sd
alg.ab<-i.env2$alg

trt<-as.factor(as.character(i.env2$treatment))
seas<-i.env2$season
samp<-i.env2$sampling
yr<-as.factor(i.env2$yr)
plts<-as.factor(i.env2$plot)

# rerun treatment model to make sure data matches exactly
tr.s.yr.mat3<-data.frame(model.matrix(~ yr*seas*trt+plts, 
  contrasts=list(trt="contr.helmert", seas="contr.helmert", yr="contr.helmert")))[,-1]
i.rda.yr3<-rda(i.com.hel2~.,data=tr.s.yr.mat3)

# just add linear effect of primary producers

tr.s.yr.mat3.prod<-data.frame(model.matrix(~ trt*seas*yr+sg.sd+alg.ab+plts, 
  contrasts=list(trt="contr.helmert", seas="contr.helmert", yr="contr.helmert")))[,-1]
# tr.s.yr.mat3.sg<-data.frame(model.matrix(~ trt*seas*yr+plts+sg.sd, 
#   contrasts=list(trt="contr.helmert", seas="contr.helmert", yr="contr.helmert")))[,-1]
# tr.s.yr.mat3.alg<-data.frame(model.matrix(~ trt*seas*yr+plts+alg.ab, 
#   contrasts=list(trt="contr.helmert", seas="contr.helmert", yr="contr.helmert")))[,-1]

i.rda.yr3.prod<-rda(i.com.hel2~.,data=tr.s.yr.mat3.prod)
anova(i.rda.yr3.prod,i.rda.yr3)

anova(i.rda.yr3.prod)
adonis2(i.com.hel2~., data = tr.s.yr.mat3.prod, method = "euclidean",by =NULL) 

i.mod.prod <- adonis2(i.com.hel2~., data = tr.s.yr.mat3.prod, method = "euclidean",by =
    "terms") 
i.mod.prod

# if above was sig we would check each type of producer separately
# i.rda.yr3.sg<-rda(i.com.hel~.,data=tr.s.yr.mat3.sg)
# i.rda.yr3.alg<-rda(i.com.hel~.,data=tr.s.yr.mat3.alg)

# add seasonal interaction with total of each type of primary producer
# or interaction between types
tr.s.yr.mat3.prod2<-data.frame(model.matrix(~ yr*seas*trt+sg.sd*alg.ab*seas+plts, 
  contrasts=list(trt="contr.helmert", seas="contr.helmert", yr="contr.helmert")))[,-1]
i.rda.prod.full<-rda(i.com.hel2~.,data=tr.s.yr.mat3.prod2)
anova(i.rda.prod.full,i.rda.yr3)

# chance order of variables to test effect of treatment with producers in model already
prod.yr.seas.trt<-data.frame(model.matrix(~ sg.sd*alg.ab*seas + yr*seas*trt + plts, 
  contrasts=list(trt="contr.helmert", seas="contr.helmert", yr="contr.helmert")))[,-1]

i.mod.full <- adonis2(i.com.hel2~., data=prod.yr.seas.trt, method="euclidean", by="terms") 
i.mod.full

# producers don't eliminate the effect of treatment

# # unnecessary deep dive into other possible model structures for producers
# # replace treatment with total of each type of primary producer
# tr.s.yr.prod3<-data.frame(model.matrix(~ seas*yr*sg.sd+seas*yr*alg.ab+plts, 
#   contrasts=list(seas="contr.helmert", yr="contr.helmert")))[,-1]
# 
# tr.s.yr.sg3<-data.frame(model.matrix(~ seas*yr*sg.sd+plts, 
#   contrasts=list(seas="contr.helmert", yr="contr.helmert")))[,-1]
# tr.s.yr.alg3<-data.frame(model.matrix(~ seas*yr*alg.ab+plts, 
#   contrasts=list(seas="contr.helmert", yr="contr.helmert")))[,-1]
# 
# i.rda.prod3<-rda(i.com.hel2~.,data=tr.s.yr.prod3)
# i.rda.sg3<-rda(i.com.hel2~.,data=tr.s.yr.sg3)
# i.rda.alg3<-rda(i.com.hel2~.,data=tr.s.yr.alg3)
# 
# anova(i.rda.prod3,i.rda.sg3) 
# anova(i.rda.prod3,i.rda.alg3) # algae is almost identical to full
# 
# 
# tr.s.yr.prod2<-data.frame(model.matrix(~ seas*yr+sg.sd*seas+alg.ab*seas+plts, 
#   contrasts=list(seas="contr.helmert", yr="contr.helmert")))[,-1]
# i.rda.prod2<-rda(i.com.hel2~.,data=tr.s.yr.prod2)
# anova(i.rda.prod3,i.rda.prod2)
# 
# tr.s.yr.alg2<-data.frame(model.matrix(~ seas*yr+alg.ab*seas+plts, 
#   contrasts=list(seas="contr.helmert", yr="contr.helmert")))[,-1]
# i.rda.alg2<-rda(i.com.hel2~.,data=tr.s.yr.alg2)
# anova(i.rda.alg2,i.rda.prod2)
# 
# ## any chance types interact with eachother? YES
# tr.s.yr.prod2b<-data.frame(model.matrix(~ yr*seas+plts+sg.sd*alg.ab, 
#   contrasts=list(seas="contr.helmert", yr="contr.helmert")))[,-1]
# tr.s.yr.prod1<-data.frame(model.matrix(~ yr*seas+plts+sg.sd+alg.ab, 
#   contrasts=list(seas="contr.helmert", yr="contr.helmert")))[,-1]
# 
# i.rda.prod2b<-rda(i.com.hel2~.,data=tr.s.yr.prod2b)
# i.rda.prod1<-rda(i.com.hel2~.,data=tr.s.yr.prod1)
# anova(i.rda.prod1,i.rda.prod2b)
# 
# ## doesn't improve with season 3-way
# # tr.s.yr.prod3b<-data.frame(model.matrix(~ yr*seas+plts+sg.sd*alg.ab*seas, 
# #   contrasts=list(seas="contr.helmert", yr="contr.helmert")))[,-1]
# # i.rda.prod3b<-rda(i.com.hel2~.,data=tr.s.yr.prod3b)
# # anova(i.rda.prod3b,i.rda.prod2b)
# 
# 
# # not nested so can't use anova, but R squared less than for treatment model for each alone
# RsquareAdj(i.rda.yr3) #treatment model
# RsquareAdj(i.rda.sg3) #seagrass model
# RsquareAdj(i.rda.alg3) #algae model
# RsquareAdj(i.rda.prod2b) #both producers interacting model comes closest but still pretty far off
# 


#fish 

hist(fish$abundance, breaks = 100)
nozeros <- fish%>%filter(abundance!=0) 
hist(nozeros$abundance, breaks = 20)
# there appears to be only one sample > 5 and it's 22! 
# not sure best approach to addressing this outlier, but for now I'm replacing it with 6

# with zeros
f2<-fish %>% 
  # mutate(abundance = ifelse(abundance>5, 6, abundance)) %>%
  pivot_wider(names_from=taxa,values_from=abundance,values_fill=0)%>%
 mutate(dummy=1)

# #without zeros
# f2<-fish%>%
#   filter(abundance!=0)%>%
#   pivot_wider(names_from=taxa,values_from=abundance,values_fill=0)


f2<-f2%>%
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
      sampling==17~2))
# f.env$plot<-as.factor(f.env$plot)


f2<-left_join(f2,sg)%>%
  left_join(alg)

# NMDS
f.env<-f2 %>% select(treatment, plot, yr, sampling, season, sg.sd, alg)
f.com<-f2 %>% select(-treatment, -plot, -yr, -sampling, -season, -sg.sd, -alg)


com.pa<-decostand(f.com,"pa")
com.dist<-vegdist(com.pa,"bray")

f.mds<-metaMDS(com.dist,trymax = 100)
plot(f.mds)
f.mds.scores<-data.frame(scores(f.mds))%>%
  bind_cols(f.env)

ggplot(data=f.mds.scores)+
  geom_jitter(aes(x=NMDS1,y=NMDS2,
                 color=treatment),size=2,alpha=.5, 
    width = 0.03, height = 0.03)+
  scale_color_viridis_d(option="B",end=.8)+
  facet_wrap(~sampling)

# rda
f.com.hel<-decostand(f.com,"hellinger")

f.pca<-rda(f.com.hel)
f.scores<-data.frame(scores(f.pca,1:3)$sites)%>%
  bind_cols(f.env)

ggplot(data=f.scores)+
  geom_jitter(aes(x=PC1,y=PC2,
    color=treatment),size=2,alpha=.5, 
    width = 0.03, height = 0.03)+
  scale_color_viridis_d(option="B",end=.8)+
  facet_wrap(~sampling)


f.env0<-f2 %>% select(sampling, treatment, plot) %>% mutate(plot=as.factor(plot))
# insufficient degrees of freedom to retain plot as factor
# f.env0<-f2 %>% select(sampling, treatment)


f.rda0<-rda(f.com.hel~.,f.env0)
# summary(f.rda)
step.fa0 <- ordistep(f.rda0,scope = formula(f.rda0),direction = "backward")
# just sampling and treatment retained
anova(f.rda0)
plot(f.rda0, scaling = 3, display = c("sp","wa", "cn"))


f.env1<-f2 %>% select(yr, season, treatment)
f.rda1<-rda(f.com.hel~.,f.env1)
# summary(f.rda)
step.fa1 <- ordistep(f.rda1,scope = formula(f.rda1),direction = "backward")
# all retained
anova(f.rda1)
plot(f.rda1, scaling = 3, display = c("sp", "cn"))

RsquareAdj(f.rda0)
RsquareAdj(f.rda1)


## redo with presence-absence data
f.com.pa<-decostand(f.com,"pa")

f.pca2<-rda(f.com.pa)
f.scores2<-data.frame(scores(f.pca2,1:3)$sites)%>%
  bind_cols(f.env)

ggplot(data=f.scores2)+
  geom_jitter(aes(x=PC1,y=PC2,
    color=treatment),size=2,alpha=.5, 
    width = 0.03, height = 0.03)+
  scale_color_viridis_d(option="B",end=.8)+
  facet_grid(season~yr)


f.env0<-f2 %>% select(plot, sampling, treatment)%>% mutate(plot=as.factor(plot))
f.rda2<-rda(f.com.pa~.,f.env0)
# summary(f.rda2)
step.fp <- ordistep(f.rda2,scope = formula(f.rda2),direction = "backward")


f.env1<-f2 %>% select(plot, yr, season, treatment)%>% mutate(plot=as.factor(plot))
f.rda2<-rda(f.com.pa~.,f.env1)
# summary(f.rda2)
step.fp <- ordistep(f.rda2,scope = formula(f.rda2),direction = "backward")

f.env0<-f2 %>% select(sampling, treatment)
f.rda2a<-rda(f.com.pa~.,f.env0)
f.env1<-f2 %>% select(yr, treatment)
f.rda2b<-rda(f.com.pa~.,f.env1)

RsquareAdj(f.rda2a)
RsquareAdj(f.rda2b)

plot(f.rda2a, scaling = 3, display = c("sp", "cn"))


# add in interactions

f.env<-f2 %>% select(treatment, plot, yr, sampling, season, sg.sd, alg)

trt<-as.factor(f.env$treatment)
samp<-f.env$sampling
plts<-as.factor(f.env$plot)
sg.sd<-f.env$sg.sd
alg.ab<-f.env$alg

f.data.trt<-data.frame(model.matrix(~ samp*trt, 
  contrasts=list(trt="contr.helmert")))[,-1]

f.rda.trt <- rda(f.com.hel~., f.data.trt) 
RsquareAdj(f.rda.trt)
plot(f.rda.trt, scaling=3, display=c("sp", "cn"))

f.mod.trt <- adonis2(f.com.hel~., data = f.data.trt, method="euclidean", by="terms") 
f.mod.trt

# add in producers 
f.data.prod<-data.frame(model.matrix(~ sg.sd + alg.ab + samp*trt, 
  contrasts=list(trt="contr.helmert")))[,-1]

f.rda.prod <- rda(f.com.hel~., f.data.prod) 
RsquareAdj(f.rda.prod)
plot(f.rda.prod, scaling=3, display=c("sp", "cn"))

f.mod.prod <- adonis2(f.com.hel~., data = f.data.prod, method="euclidean", by="terms") 
f.mod.prod

#  producers only

f.rda.prod0 <- rda(f.com.hel~ sg.sd + alg.ab + samp, f.env) 
RsquareAdj(f.rda.prod0)
plot(f.rda.prod0, scaling=3, display=c("sp", "cn"))

f.mod.prod0 <- adonis2(f.com.hel~sg.sd + alg.ab + samp, data = f.env, method="euclidean", by="terms") 
f.mod.prod0


RsquareAdj(f.rda.prod0) # producers only
RsquareAdj(f.rda.trt)# treatment only
RsquareAdj(f.rda.prod) # both 


# appears to be an effect of sponge beyond that of producers
# result nearly identical with or without outlier 

# presence-absence
f.mod1 <- adonis2(f.com.pa~., data = f.data, method="euclidean", by="terms") 
f.mod1

f.mod1.rda <- rda(f.com.pa~., f.data)
plot(f.mod1.rda, scaling=3, display=c("sp", "cn"))


# experiment with other standardizations
f.com.m<-decostand(f.com,"rank")
f.com.m<-decostand(f.com,"frequency")
f.com.m<-decostand(f.com,"max")
f.com.m <- f.com.pa
f.com.m <- f.com.hel
mod <- rda(f.com.m~., f.data)

plot(mod, type="n", scaling="species")
points(mod, pch=1, col="grey", cex=1, scaling="species")
text(mod, dis="cn", scaling="species")
text(mod, "species", col="blue", cex=0.8, scaling="site")
# not sure what makes most sense



##### fish for sponge vs. structural control only ####
f.env2a <- f.env %>% filter(treatment!="blank")

# bind back the community data for just the retained plots
f3 <- left_join(f.env2a, f2)

# separate env and com again to sure order matches in this new dataset 
f.env2 <- f3 %>% select(treatment, plot, yr, sampling, season, sg.sd, alg)
f.com2 <- f3 %>% select(-treatment, -plot, -yr, -sampling, -season, -sg.sd, -alg)

trt<-as.factor(f.env2$treatment)
samp<-f.env2$sampling
plts<-as.factor(f.env2$plot)
sg.sd<-f.env2$sg.sd
alg.ab<-f.env2$alg

f.data2<-data.frame(model.matrix(~ sg.sd + alg.ab + samp*trt, 
  contrasts=list(trt="contr.helmert")))[,-1]

f.com.hel2<-decostand(f.com2,"hellinger")

f.mod2 <- adonis2(f.com.hel2~., data = f.data2, method="euclidean", by="terms") 
f.mod2
# there is still an effect of sponge beyond that of producers and structure control

f.mod2.rda <- rda(f.com.hel2~., f.data2)
plot(f.mod2.rda, scaling=3, display=c("sp", "cn"))


# and for presense-absence
f.data2<-data.frame(model.matrix(~ sg.sd + alg.ab + samp*trt, 
  contrasts=list(trt="contr.helmert")))[,-1]

f.com.pa2<-decostand(f.com2,"pa")

f.mod3 <- adonis2(f.com.pa2~., data = f.data2, method="euclidean", by="terms") 
f.mod3
# there is still an effect of sponge beyond that of producers and structure control

f.mod3.rda <- rda(f.com.pa2~., f.data2)
plot(f.mod3.rda, scaling=3, display=c("sp", "cn"))

