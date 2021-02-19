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
i.com0<-i.com0[,colSums(i.com0)!=0]
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
trt0<-as.factor(i.env0$treatment)
trt0.mat<-data.frame(model.matrix(~ trt0, 
                                  contrasts=list(trt0="contr.helmert")))[,-1]
i0.rda<-rda(i.com.hel0~.,data=trt0.mat)
RsquareAdj(i0.rda)
(invert.initial<-anova(i0.rda.null,i0.rda,permutations=how(nperm=999)))

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
               contrasts=list(trt="contr.helmert")))[,-1]
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
alg.e<-i.env$alg
sg<-i.env$sg.sd
sggrow<-i.env$grow


i.rda.prod<-rda(i.com.hel~.,i.env.prod)

# does this model do better than the null
anova(i.rda.null,i.rda.prod)

# yes it does
# now does it do better than the treatment model
RsquareAdj(i.rda.yr.s.treat)$adj.r.squared
RsquareAdj(i.rda.prod)$adj.r.squared

#What about a model that uses sg grow instead of treatment
sg.s.yr.mat.nop<-data.frame(model.matrix(~ yr*seas*sggrow, 
                                         contrasts=list(seas="contr.helmert",yr="contr.helmert")))[,-1]
i.rda.sgsyr<-rda(i.com.hel~.,sg.s.yr.mat.nop)

# does it do better than the null
anova(i.rda.null,i.rda.sgsyr)

# yes it does - how about compared to model with treatment
RsquareAdj(i.rda.yr.s.treat)$adj.r.squared
RsquareAdj(i.rda.sgsyr)$adj.r.squared

# on its own productivity doesn't do a better job than just treatment

# not on its own, no. What if we include different measures of productivity in our treatment model
# from now on I'm referring to the year*season*treatment as i.rda.best
i.rda.best<-rda(i.com.hel~.,tr.s.yr.mat.nop)

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
plot(i.rda.best, scaling = 3, display = c("sp", "cn"))

# note the order matters for adonis2 with by="terms" and the by="margin" doesn't seem to work
i.mod <- adonis2(i.com.hel~., data = tr.s.yr.mat.nop, method="euclidean", by="terms") 
i.mod

#this is where I stopped. I think the moral of the story is that the interaction between 
# treatment, season, and year has an affect on the community
# now we have to figure out what that is.

# first going to make sure where there's a difference between treatments

# look at 1 month into experiment
i.com1<-i.com[i.env$sampling==1,]
i.com1<-i.com1[,colSums(i.com1)!=0]
i.env1<-i.env[i.env$sampling==1,]

i.com1.pa<-decostand(i.com1,"pa")
i.com1.hel<-decostand(i.com1.pa,"hellinger")


i1.rda.null<-rda(i.com1.hel~1)
trt1<-as.factor(i.env1$treatment)
trt1.mat<-data.frame(model.matrix(~ trt1, 
                                  contrasts=list(trt1="contr.helmert")))[,-1]
i1.rda<-rda(i.com1.hel~.,data=trt1.mat)
RsquareAdj(i1.rda)
(invert.s1<-anova(i1.rda.null,i1.rda))

# look at 5 months into experiment
i.com5<-i.com[i.env$sampling==5,]
i.com5<-i.com5[,colSums(i.com5)!=0]
i.env5<-i.env[i.env$sampling==5,]

i.com5.pa<-decostand(i.com5,"pa")
i.com5.hel<-decostand(i.com5.pa,"hellinger")


i5.rda.null<-rda(i.com5.hel~1)
trt5<-as.factor(i.env5$treatment)
trt5.mat<-data.frame(model.matrix(~ trt5, 
                                  contrasts=list(trt5="contr.helmert")))[,-1]
i5.rda<-rda(i.com5.hel~.,data=trt5.mat)
RsquareAdj(i5.rda)
(invert.s5<-anova(i5.rda.null,i5.rda))

# treatment is significant by 5 months into the experiment. 
# Does this hold in the next two samplings?

# look at 12 months into experiment
i.com12<-i.com[i.env$sampling==12,]
i.com12<-i.com12[,colSums(i.com12)!=0]
i.env12<-i.env[i.env$sampling==12,]

i.com12.pa<-decostand(i.com12,"pa")
i.com12.hel<-decostand(i.com12.pa,"hellinger")


i12.rda.null<-rda(i.com12.hel~1)
trt12<-as.factor(i.env12$treatment)
trt12.mat<-data.frame(model.matrix(~ trt12, 
                                  contrasts=list(trt12="contr.helmert")))[,-1]
i12.rda<-rda(i.com12.hel~.,data=trt12.mat)
RsquareAdj(i12.rda)
(invert.s12<-anova(i12.rda.null,i12.rda))
anova(i12.rda,by="axis",permutations = how(n=999))

# look at 17 months into experiment
i.com17<-i.com[i.env$sampling==17,]
i.com5<-i.com5[,colSums(i.com5)!=0]
i.env17<-i.env[i.env$sampling==17,]

i.com17.pa<-decostand(i.com17,"pa")
i.com17.hel<-decostand(i.com17.pa,"hellinger")


i17.rda.null<-rda(i.com17.hel~1)
trt17<-as.factor(i.env17$treatment)
trt17.mat<-data.frame(model.matrix(~ trt17, 
                                   contrasts=list(trt17="contr.helmert")))[,-1]
i17.rda<-rda(i.com17.hel~.,data=trt17.mat)
RsquareAdj(i17.rda)
(invert.s17<-anova(i17.rda.null,i17.rda))

# treatment stays an important factor from 5 months on. 
# at 5 months which treatments are different from each other?

i.com5.bf<-i.com5[i.env5$treatment!="real",]
i.env5.bf<-i.env5[i.env5$treatment!="real",]

i.com5bf.pa<-decostand(i.com5.bf,"pa")
i.com5bf.hel<-decostand(i.com5bf.pa,"hellinger")

i5bf.null<-rda(i.com5bf.hel~1)

trt5bf<-as.factor(i.env5.bf$treatment)
i5bf.rda<-rda(i.com5bf.hel~treatment,i.env5.bf)
RsquareAdj(i5bf.rda)

(invert.s5bf<-anova(i5bf.null,i5bf.rda))

# blank and fake aren't different

i.com5.br<-i.com5[i.env5$treatment!="fake",]
i.env5.br<-i.env5[i.env5$treatment!="fake",]

i.com5br.pa<-decostand(i.com5.br,"pa")
i.com5br.hel<-decostand(i.com5br.pa,"hellinger")

i5br.null<-rda(i.com5br.hel~1)

trt5br<-as.factor(i.env5.br$treatment)
i5br.rda<-rda(i.com5br.hel~treatment,i.env5.br)
RsquareAdj(i5br.rda)

(invert.s5br<-anova(i5br.null,i5br.rda))

# close but blank and real are not actually different

i.com5.fr<-i.com5[i.env5$treatment!="blank",]
i.env5.fr<-i.env5[i.env5$treatment!="blank",]

i.com5fr.pa<-decostand(i.com5.fr,"pa")
i.com5fr.hel<-decostand(i.com5fr.pa,"hellinger")

i5fr.null<-rda(i.com5fr.hel~1)

trt5fr<-as.factor(i.env5.fr$treatment)
i5fr.rda<-rda(i.com5fr.hel~treatment,i.env5.fr)
RsquareAdj(i5fr.rda)

(invert.s5fr<-anova(i5fr.null,i5fr.rda))

# real and fake are different


# at 12 months which treatments are different from each other?

i.com12.bf<-i.com12[i.env12$treatment!="real",]
i.env12.bf<-i.env12[i.env12$treatment!="real",]

i.com12bf.pa<-decostand(i.com12.bf,"pa")
i.com12bf.hel<-decostand(i.com12bf.pa,"hellinger")

i12bf.null<-rda(i.com12bf.hel~1)

trt12bf<-as.factor(i.env12.bf$treatment)
i12bf.rda<-rda(i.com12bf.hel~treatment,i.env12.bf)
RsquareAdj(i12bf.rda)

(invert.s12bf<-anova(i12bf.null,i12bf.rda))

# blank and fake aren't different

i.com12.br<-i.com12[i.env12$treatment!="fake",]
i.env12.br<-i.env12[i.env12$treatment!="fake",]

i.com12br.pa<-decostand(i.com12.br,"pa")
i.com12br.hel<-decostand(i.com12br.pa,"hellinger")

i12br.null<-rda(i.com12br.hel~1)

trt12br<-as.factor(i.env12.br$treatment)
i12br.rda<-rda(i.com12br.hel~treatment,i.env12.br)
RsquareAdj(i12br.rda)

(invert.s12br<-anova(i12br.null,i12br.rda))

# close but blank and real are not actually different

i.com12.fr<-i.com12[i.env12$treatment!="blank",]
i.env12.fr<-i.env12[i.env12$treatment!="blank",]

i.com12fr.pa<-decostand(i.com12.fr,"pa")
i.com12fr.hel<-decostand(i.com12fr.pa,"hellinger")

i12fr.null<-rda(i.com12fr.hel~1)

trt12fr<-as.factor(i.env12.fr$treatment)
i12fr.rda<-rda(i.com12fr.hel~treatment,i.env12.fr)
RsquareAdj(i12fr.rda)

(invert.s12fr<-anova(i12fr.null,i12fr.rda))

# real and fake are different

# at 17 months which treatments are different from each other?

i.com17.bf<-i.com17[i.env17$treatment!="real",]
i.env17.bf<-i.env17[i.env17$treatment!="real",]

i.com17bf.pa<-decostand(i.com17.bf,"pa")
i.com17bf.hel<-decostand(i.com17bf.pa,"hellinger")

i17bf.null<-rda(i.com17bf.hel~1)

trt17bf<-as.factor(i.env17.bf$treatment)
i17bf.rda<-rda(i.com17bf.hel~treatment,i.env17.bf)
RsquareAdj(i17bf.rda)

(invert.s17bf<-anova(i17bf.null,i17bf.rda))

# blank and fake ARE different

i.com17.br<-i.com17[i.env17$treatment!="fake",]
i.env17.br<-i.env17[i.env17$treatment!="fake",]

i.com17br.pa<-decostand(i.com17.br,"pa")
i.com17br.hel<-decostand(i.com17br.pa,"hellinger")

i17br.null<-rda(i.com17br.hel~1)

trt17br<-as.factor(i.env17.br$treatment)
i17br.rda<-rda(i.com17br.hel~treatment,i.env17.br)
RsquareAdj(i17br.rda)

(invert.s17br<-anova(i17br.null,i17br.rda))

# blank and real ARE different

i.com17.fr<-i.com17[i.env17$treatment!="blank",]
i.env17.fr<-i.env17[i.env17$treatment!="blank",]

i.com17fr.pa<-decostand(i.com17.fr,"pa")
i.com17fr.hel<-decostand(i.com17fr.pa,"hellinger")

i17fr.null<-rda(i.com17fr.hel~1)

trt17fr<-as.factor(i.env17.fr$treatment)
i17fr.rda<-rda(i.com17fr.hel~treatment,i.env17.fr)
RsquareAdj(i17fr.rda)

(invert.s17fr<-anova(i17fr.null,i17fr.rda))

# real and fake are different

# so sponge and structure control diverge from each other early on but neither is different
# from the control until 17 months until all three plot types have a different community

# make figures for this

anova(i0.rda, permutations=how(nperm=999))
# Tests of all canonical axes
anova(i0.rda, by="axis", permutations=how(nperm=999))
init.scores<-scores(i0.rda.null,scaling=1,1:2)

i.env0p<-bind_cols(i.env0,data.frame(init.scores$sites))

hull0 <- i.env0p %>%
  group_by(treatment)%>%
  slice(chull(PC1, PC2))
  # slice(chull(RDA1, RDA2))

spe.good <- goodness(i0.rda.null,display = "species")

spr<-init.scores$species[which(spe.good[,1]>=0.6 | spe.good[,2]>=0.6),]

(i0<-ggplot()+
    geom_hline(aes(yintercept=0),linetype="dashed",color="grey")+
    geom_vline(aes(xintercept=0),linetype="dashed",color="grey")+
    geom_polygon(data = hull0, aes(x=RDA1,y=RDA2,color=treatment,fill=treatment),alpha = 0.1)+
    geom_point(data=i.env0p,aes(x=RDA1,y=RDA2,color=treatment),size=4)+
    theme_bw()+
    theme(panel.grid = element_blank())+
    coord_fixed()+
    ylim(-0.65,0.65)+
    xlim(-0.65,0.65)+  
    scale_color_viridis_d(option="A",begin=0,end=0.6,"")+
    scale_fill_viridis_d(option="A",begin=0,end=0.6,""))





# make figures for month 5
m5.scores<-scores(i5.rda.null,scaling=3,1:2)

i.env5p<-bind_cols(i.env5,data.frame(m5.scores$sites))

hull5 <- i.env5p %>%
  group_by(treatment)%>%
  slice(chull(PC1, PC2))

(i5<-ggplot()+
    geom_hline(aes(yintercept=0),linetype="dashed",color="grey")+
    geom_vline(aes(xintercept=0),linetype="dashed",color="grey")+
    geom_polygon(data = hull5, aes(x=PC1,y=PC2,color=treatment,fill=treatment),alpha = 0.1)+
    geom_point(data=i.env5p,aes(x=PC1,y=PC2,color=treatment),size=4)+
    theme_bw()+
    theme(panel.grid = element_blank())+
    coord_fixed()+
    ylim(-0.65,0.65)+
    xlim(-0.65,0.65)+  
    scale_color_viridis_d(option="A",begin=0,end=0.6,"")+
    scale_fill_viridis_d(option="A",begin=0,end=0.6,""))

# make figures for month 12
m12.scores<-scores(i12.rda.null,scaling=3,1:2)

i.env12p<-bind_cols(i.env12,data.frame(m12.scores$sites))

hull12 <- i.env12p %>%
  group_by(treatment)%>%
  slice(chull(PC1, PC2))

(i12<-ggplot()+
    geom_hline(aes(yintercept=0),linetype="dashed",color="grey")+
    geom_vline(aes(xintercept=0),linetype="dashed",color="grey")+
    geom_polygon(data = hull12, aes(x=PC1,y=PC2,color=treatment,fill=treatment),alpha = 0.1)+
    geom_point(data=i.env12p,aes(x=PC1,y=PC2,color=treatment),size=4)+
    theme_bw()+
    theme(panel.grid = element_blank())+
    coord_fixed()+
    ylim(-0.6,0.6)+
    xlim(-0.6,0.6)+  
    scale_color_viridis_d(option="A",begin=0,end=0.65,"")+
    scale_fill_viridis_d(option="A",begin=0,end=0.65,""))


# make figures for month 17
m17.scores<-scores(i17.rda.null,scaling=3,1:2)

i.env17p<-bind_cols(i.env17,data.frame(m17.scores$sites))

hull17 <- i.env17p %>%
  group_by(treatment)%>%
  slice(chull(PC1, PC2))

(i17<-ggplot()+
    geom_hline(aes(yintercept=0),linetype="dashed",color="grey")+
    geom_vline(aes(xintercept=0),linetype="dashed",color="grey")+
    geom_polygon(data = hull17, aes(x=PC1,y=PC2,color=treatment,fill=treatment),alpha = 0.1)+
    geom_point(data=i.env17p,aes(x=PC1,y=PC2,color=treatment),size=4)+
    theme_bw()+
    theme(panel.grid = element_blank())+
    coord_fixed()+
    ylim(-0.65,0.65)+
    xlim(-0.65,0.65)+  
    scale_color_viridis_d(option="A",begin=0,end=0.6,"")+
    scale_fill_viridis_d(option="A",begin=0,end=0.6,""))


# now look at univariate results

i.env$spr<-specnumber(i.com)
i.env$div<-diversity(i.com)
i.env$plot<-as.factor(i.env$plot)
i.env$treatment<-as.factor(i.env$treatment)
i.env$season<-as.factor(i.env$season)

library(lmerTest)

i.env0uni<-i.env0%>%
  mutate(strt.spr=specnumber(i.com0),
         strt.div=diversity(i.com0),
         plot=factor(plot))%>%
  select(treatment,plot,strt.spr,strt.div)

i.env.uni<-i.env%>%
  filter(sampling!=0)%>%
  left_join(i.env0uni)

i.env.uni$treatment<-as.factor(i.env.uni$treatment)

spr.lmer<-lmer(spr~treatment*as.factor(sampling) + sg.sd+grow+(1|plot)+
                 offset(strt.spr),
                  data = i.env.uni%>%
                 mutate(treatment=relevel(treatment, ref = "real")))
summary(spr.lmer)
anova(spr.lmer)

# looks like species richness slowly increases in sponge plots. Eventually
# it increases in structure control but doesn't start to diverge from
# control until the last sampling


i.sum<-i.env.uni%>%
  group_by(treatment,sampling)%>%
  summarize(spr.m=mean(spr),div.m=mean(div),spr.sd=sd(spr),div.sd=sd(div))

i.sum$treatment<-factor(i.sum$treatment,labels=c("Control","Structure Control","Sponge"))

ggplot()+
  geom_line(data=i.sum,aes(x=sampling,y=spr.m,group=treatment,color=treatment),position=position_dodge(0.5))+
  geom_errorbar(data=i.sum,aes(x=sampling,ymin=spr.m-spr.sd,ymax=spr.m+spr.sd,color=treatment),width=.5,position=position_dodge(0.5))+
  geom_point(data=i.sum,aes(x=sampling,y=spr.m,color=treatment),size=5,position=position_dodge(0.5))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        legend.text = element_text(size=10))+
  scale_color_viridis_d(option="A", begin=0, end=0.6,"")+
  ylab("Species Richness")+
  xlab("Months into Experiment")

ggsave("figures/invert_sp_richness.jpg")

# now look at diversity

div.lmer<-lmer(div~treatment*as.factor(sampling) + sg.sd+grow+(1|plot)+
                 offset(strt.div),
               data = i.env.uni%>%
                 mutate(treatment=relevel(treatment, ref = "real")))
summary(div.lmer)
anova(div.lmer)

ggplot()+
  geom_line(data=i.sum,aes(x=sampling,y=div.m,group=treatment,color=treatment),position=position_dodge(0.5))+
  geom_errorbar(data=i.sum,aes(x=sampling,ymin=div.m-div.sd,ymax=div.m+div.sd,color=treatment),width=.5,position=position_dodge(0.5))+
  geom_point(data=i.sum,aes(x=sampling,y=div.m,color=treatment),size=5,position=position_dodge(0.5))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        legend.text = element_text(size=10))+
  scale_color_viridis_d(option="A", begin=0, end=0.6,"")+
  ylab("Diversity")+
  xlab("Months into Experiment")

ggsave("figures/invert_sp_diversity.jpg")

# the diversity story isn't as clear cut but generally the same as richness




