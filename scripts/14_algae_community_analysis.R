# this script will explore the algae community data

library(tidyverse)
if(!require(vegan))install.packages("vegan"); library(vegan)

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
alg<-algae

#From these two tables I can see acetabularia, dictyota, cladocephalus, and penicillus may be important in showing differences between treatment. Also, there may be a trend in algal abundance decreasing through time with blank and fake treatment while real mostly stayed constant(possible seasonal effect) 
table(alg$treatment,alg$taxa)
table(alg$treatment,alg$sampling)


algae<-algae%>%
  pivot_wider(names_from=taxa,values_from=abundance,values_fill=0)

algae<-algae%>%
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



a0<-algae%>%filter(yr == 0)
a2<-algae%>%filter(yr != 0)

# Organize data


a.env<-a2 %>% select(treatment, plot, yr, sampling, season)%>%
  left_join(sg)%>%
  left_join((sggrow))

a.com<-a2 %>% select(-treatment, -plot, -yr, -sampling, -season)

# pre-experiment data
a.env0<-a0 %>% select(treatment, plot, yr, sampling, season)%>%
  left_join(sg)%>%
  left_join((sggrow))
a.com0<-a0 %>% select(-treatment, -plot, -yr, -sampling, -season)

a.com0<-a.com0[,colSums(a.com0)!=0]

# use hellinger: square root of method to standardize species data
a.com.pa<-decostand(a.com,"pa")

a.com.hel<-decostand(a.com.pa,"hellinger")

a.pca<-rda(a.com.hel)
a.pca

a.scores<-data.frame(scores(a.pca,1:3)$sites)%>%
  bind_cols(a.env)

ggplot(data=a.scores)+
  geom_jitter(aes(x=PC1,y=PC2,
    color=treatment),size=2,alpha=.65, 
    width = 0.03, height = 0.03)+
  scale_color_viridis_d(option="B",end=.8)+
  facet_wrap(~sampling)

ggplot(data=a.scores)+
  geom_jitter(aes(x=PC2,y=PC3,
    color=treatment),size=2,alpha=.65, 
    width = 0.03, height = 0.03)+
  scale_color_viridis_d(option="B",end=.8)+
  facet_wrap(~sampling)

#So it looks like over time treatments are becomings really differnt.
#start examining statistical relationship

# first look at initial pre-experiment
a.com.pa.0<-decostand(a.com0,"pa")
a.com.hel0<-decostand(a.com.pa.0,"hellinger")

a.com.pa<-decostand(a.com,"pa")
a.com.hel<-decostand(a.com.pa,"hellinger")

# look at just treatment at 0 and later
a0.rda.null<-rda(a.com.hel0~1)
trt0<-as.factor(a.env0$treatment)
trt0.mat<-data.frame(model.matrix(~ trt0, 
                                  contrasts=list(trt0="contr.helmert")))[,-1]
a0.rda<-rda(a.com.hel0~.,data=trt0.mat)
RsquareAdj(a0.rda)
(algae.initial<-anova(a0.rda.null,a0.rda,permutations=how(nperm=999)))

# treatment does not explain a significant amount of variance between plots initially

# confirming this with an anosim
(a0.anosim<-anosim(a.com.pa.0,a.env0$treatment,permutations = 999,distance="bray"))
summary(a0.anosim)
# the anosim results confirm this. 

# now start building more and more complex models for experiment data
#add in interactions
trt<-as.factor(a.env$treatment)
seas<-a.env$season
samp<-a.env$sampling
yr<-as.factor(a.env$yr)
plts<-as.factor(a.env$plot)

# start with the simplest model that makes sense - there are two of these - season and year
# with an interaction with treatment. Look at these with and without plots

# first question is an interaction between samp* treatment better than intercept only
# first build intercept only
a.rda.null<-rda(a.com.hel~1)

# now make samp*treat matrix
tr.samp.mat.np<-data.frame(model.matrix(~ samp*trt, 
                                        contrasts=list(trt="contr.helmert")))[,-1]
a.rda.samp.treat<-rda(a.com.hel~.,data=tr.samp.mat.np)

# now check to see if its better than intercept

anova(a.rda.null,a.rda.samp.treat)

# it is better than null
# does including plot help
tr.samp.mat<-data.frame(model.matrix(~ samp*trt + plts, 
                                     contrasts=list(trt="contr.helmert")))[,-1]
a.rda.samp.treat.plot<-rda(a.com.hel~.,data=tr.samp.mat)

# check to see if plot should be included here
anova(a.rda.samp.treat.plot,a.rda.samp.treat)

# no difference as far as anova - look at rsquare
RsquareAdj(a.rda.samp.treat.plot)$adj.r.squared
RsquareAdj(a.rda.samp.treat)$adj.r.squared

# r squared better without plot
# current best model is just samp*treat

# now look at next simplest "base" model - this is treatment*year*season

tr.s.yr.mat<-data.frame(model.matrix(~ yr*seas*trt + plts, 
                                     contrasts=list(trt="contr.helmert", seas="contr.helmert",yr="contr.helmert")))[,-1]
tr.s.yr.mat.nop<-data.frame(model.matrix(~ yr*seas*trt, 
                                         contrasts=list(trt="contr.helmert", seas="contr.helmert",yr="contr.helmert")))[,-1]

# look at whether or not model without plot is better than current best model
a.rda.yr.s.treat<-rda(a.com.hel~.,tr.s.yr.mat.nop)

# models aren't nested so have to look at adjusted R2
RsquareAdj(a.rda.yr.s.treat)$adj.r.squared
RsquareAdj(a.rda.samp.treat)$adj.r.squared

# So it seems like year*season*treatment is the best model - see if plots make a difference here
a.rda.yr.s.treat.plot<-rda(a.com.hel~.,tr.s.yr.mat)
anova(a.rda.yr.s.treat,a.rda.yr.s.treat.plot)

# plot doesn't improve things here - how about for r2
RsquareAdj(a.rda.yr.s.treat)$adj.r.squared
RsquareAdj(a.rda.yr.s.treat.plot)$adj.r.squared # this is the best model at the moment
RsquareAdj(a.rda.samp.treat.plot)$adj.r.squared 


# is the interaction between sampling and treatment important?
tr.samp.mat1<-data.frame(model.matrix(~ samp+trt+plts, 
                                      contrasts=list(trt="contr.helmert")))[,-1]
a.rda.samp.treat.noint<-rda(a.com.hel~.,data=tr.samp.mat1)

anova(a.rda.samp.treat.plot,a.rda.samp.treat.noint)

# there is a significant difference between the models but let's look at r2 just to be safe.
RsquareAdj(a.rda.samp.treat.noint)$adj.r.squared
RsquareAdj(a.rda.samp.treat.plot)$adj.r.squared
RsquareAdj(a.rda.yr.s.treat.plot)$adj.r.squared #This is still the best model

# now look into whether or not productivity measures explain community patterns better

a.env.prod<-a.env[,6:7]
sg<-a.env$sg.sd
sggrow<-a.env$grow


a.rda.prod<-rda(a.com.hel~.,a.env.prod)

# does this model do better than the null
anova(a.rda.null,a.rda.prod)

# yes it does
# now does it do better than the best model so far
RsquareAdj(a.rda.yr.s.treat.plot)$adj.r.squared
RsquareAdj(a.rda.prod)$adj.r.squared

#no

#What about a model that uses sg grow instead of treatment
sg.s.yr.mat.nop<-data.frame(model.matrix(~ yr*seas*sggrow + plts,                                       
          contrasts=list(yr="contr.helmert", seas="contr.helmert")))[,-1]
a.rda.sgsyr<-rda(a.com.hel~.,sg.s.yr.mat.nop)

# does it do better than the null
anova(a.rda.null,a.rda.sgsyr)

# yes it does - how about compared to model with treatment
RsquareAdj(a.rda.yr.s.treat.plot)$adj.r.squared
RsquareAdj(a.rda.sgsyr)$adj.r.squared

# not on its own, no. What if we include different measures of productivity in our treatment model
# from now on I'm referring to the year*season*treatment+plot as a.rda.best
a.rda.best<-rda(a.com.hel~.,tr.s.yr.mat)

b.sg.mat<-data.frame(model.matrix(~ yr*seas*trt + plts +sg, 
                                  contrasts=list(trt="contr.helmert", seas="contr.helmert", yr="contr.helmert" )))[,-1]
b.sgg.mat<-data.frame(model.matrix(~ yr*seas*trt + plts +sggrow, 
                                   contrasts=list(trt="contr.helmert", seas="contr.helmert", yr="contr.helmert")))[,-1]
b.sg.sgg.mat<-data.frame(model.matrix(~ yr*seas*trt+ plts +sg +sggrow, 
                                      contrasts=list(trt="contr.helmert", seas="contr.helmert", yr="contr.helmert")))[,-1]

# start with all of them added
a.rda.bprod<-rda(a.com.hel~.,b.sg.sgg.mat)

anova(a.rda.best,a.rda.bprod)

# anova didn't show significant difference, but I will still look at r2
RsquareAdj(a.rda.best)$adj.r.squared
RsquareAdj(a.rda.bprod)$adj.r.squared

# best model still highest r2... only by 0.01 though

# look at just sg growth

a.rda.bgrow<-rda(a.com.hel~.,b.sgg.mat)

anova(a.rda.best,a.rda.bgrow)

# no sig difference between models - look at r2
RsquareAdj(a.rda.best)$adj.r.squared
RsquareAdj(a.rda.bgrow)$adj.r.squared

# model with just seagrass growth is now slightly higher by 0.002 r2

# look at just sg
a.rda.bsg<-rda(a.com.hel~.,b.sg.mat)

anova(a.rda.best,a.rda.bsg)

# no difference between models - look at r2
RsquareAdj(a.rda.best)$adj.r.squared
RsquareAdj(a.rda.bsg)$adj.r.squared
RsquareAdj(a.rda.bprod)$adj.r.squared
RsquareAdj(a.rda.bgrow)$adj.r.squared

# Conclusion: adding any of the productivity variables didn't change the amount of variance enough to be significantly different than the best model. So we will keep using the best model yr*season*treatment +plot.

# best model at the moment: year*season*treatment + plot

#best model at the moment
plot(a.rda.best, scaling = 3, display = c("sp", "cn"))

# note the order matters for adonis2 with by="terms" and the by="margin" doesn't seem to work
a.mod <- adonis2(a.com.hel~., data = tr.s.yr.mat, method="euclidean", by="terms") 
a.mod


# first going to make sure where there's a difference between treatments using anosim

# look at 1 month into experiment
a.com1<-a.com.pa[a.env$sampling==1,]
a.com1<-a.com1[,colSums(a.com1)!=0]
a.env1<-a.env[a.env$sampling==1,]

(a1.anosim<-anosim(a.com1,a.env1$treatment))
# no difference

# look at 5 months into experiment
a.com5<-a.com.pa[a.env$sampling==5,]
a.com5<-a.com5[,colSums(a.com5)!=0]
a.env5<-a.env[a.env$sampling==5,]

(a5.anosim<-anosim(a.com5,a.env5$treatment))

# Still no significant difference

# look at 12 months into experiment
a.com12<-a.com.pa[a.env$sampling==12,]
a.com12<-a.com12[,colSums(a.com12)!=0]
a.env12<-a.env[a.env$sampling==12,]

(a12.anosim<-anosim(a.com12,a.env12$treatment))
# summary(a12.anosim)
# significantly different. R is larger (closer to 1) and sigificance is less than 0.05

# look at 17 months into experiment
a.com17<-a.com.pa[a.env$sampling==17,]
a.com17<-a.com17[,colSums(a.com17)!=0]
a.env17<-a.env[a.env$sampling==17,]

(a17.anosim<-anosim(a.com17,a.env17$treatment))
summary(a17.anosim)
# Significantly different but R isn't quite as large.

#I think this corresponds with the results of the RDA that year season and treatment are important

# treatment stays an important factor from 12 months on. 
# at 12 months which treatments are different from each other?
#This first anosim compares blank and fake treatment at 12 months

# at 12 months which treatments are different from each other?

a.com12.bf<-a.com12[a.env12$treatment!="real",]
a.env12.bf<-a.env12[a.env12$treatment!="real",]

(a12bf.anosim<-anosim(a.com12.bf,a.env12.bf$treatment))

# blank and fake are significantly different

# a.com12.br<-a.com12[a.env12$treatment!="fake",]
# a.env12.br<-a.env12[a.env12$treatment!="fake",]
# 
# (a12br.anosim<-anosim(a.com12.br,a.env12.br$treatment))
# 
# # blank and real are VERY different

a.com12.fr<-a.com12[a.env12$treatment!="blank",]
a.env12.fr<-a.env12[a.env12$treatment!="blank",]

(a12fr.anosim<-anosim(a.com12.fr,a.env12.fr$treatment))

# real and fake are VERY different

# at 17 months which treatments are different from each other?

a.com17.bf<-a.com17[a.env17$treatment!="real",]
a.env17.bf<-a.env17[a.env17$treatment!="real",]

(a17bf.anosim<-anosim(a.com17.bf,a.env17.bf$treatment))

# blank and fake are not different

# a.com17.br<-a.com17[a.env17$treatment!="fake",]
# a.env17.br<-a.env17[a.env17$treatment!="fake",]
# 
# (a17br.anosim<-anosim(a.com17.br,a.env17.br$treatment))
# 
# # blank and real ARE different

a.com17.fr<-a.com17[a.env17$treatment!="blank",]

a.env17.fr<-a.env17[a.env17$treatment!="blank",]

(a17fr.anosim<-anosim(a.com17.fr,a.env17.fr$treatment))

# real and fake ARE different

# so algal communities in sponge plots diverge from control about a year later 
# control and structure control diverge at 12 months but then don't at 17 months?


# which species are driving this difference?
(a12.simper<-simper(a.com12,a.env12$treatment,permutations = 999))
summary(a12.simper)

(a17.simper<-simper(a.com17,a.env17$treatment,permutations = 999))
summary(a17.simper)

#looks like cladocephalus, laurencia, acetabularia, and penicillus are influential
# and consistently contribute to differences between the groups here

acom17diff<-bind_cols(a.com.pa[,colnames(a.com.pa) %in% c(
    "acetabularia",
    "cladocephalus",
    "laurencia",
    "penicillus")],
                      a.env)%>%
  pivot_longer(1:4,names_to="taxa",values_to="presence")%>%
  group_by(treatment,taxa,yr)%>%
  summarize(n.plots = sum(presence))%>%
  mutate(treatment=factor(treatment,levels=c("blank","fake","real"),labels = c("Control","Structure Control","Sponge")))

ggplot(data=acom17diff)+
  geom_col(aes(x=taxa,y=n.plots,fill=treatment),width=.5,
           position=position_dodge(.5))+
  scale_fill_viridis_d(option="A",begin=0,end=0.6,"")+
  # geom_point(aes(x=taxa, y=n.plots,color=treatment),
  #            position=position_dodge(.5),size=5)+
  scale_color_viridis_d(option="A",begin=0,end=0.6,"")+
  facet_wrap(~yr)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        legend.position = "top")+
  ylab("Number of plots where taxa is present")

#TODO: add time = 0 to this plot?


# make figures for this

init.scores<-scores(rda(a.com.hel0),scaling=1,1:2)

a.env0p<-bind_cols(a.env0,data.frame(init.scores$sites))%>%
  mutate(treatment=factor(treatment, labels=c("Control","Structure Control","Sponge")))

hull0 <- a.env0p %>%
  group_by(treatment)%>%
  slice(chull(PC1, PC2))
#  slice(chull(RDA1, RDA2))
circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

# the denominator of the diameter is the number of principle components where the 
# eigenvalue is >0
circ <- circleFun(center=c(0,0),diameter=sqrt(2/7),npoints = 500)


spr0<-data.frame(init.scores$species)
sprp0<-data.frame(spr0[abs(spr0[,1])>=max(circ[,1])|abs(spr0[,2])>=max(circ[,2]),])
sprp0$taxa<-rownames(sprp0)
sprp0<-sprp0%>%
  filter(taxa %in% c("acetabularia",
                     "cladocephalus",
                     "laurencia",
                     "penicillus"))
taxa<-rownames(sprp0)

# look at how much variation each axis explains to add to axis labels
summary(a0.rda.null)$cont$importance

(a0<-ggplot()+
    ylim(-1.2,1.2)+
    xlim(-1.2,1.2)+
    geom_hline(aes(yintercept=0),linetype="dashed",color="grey")+
    geom_vline(aes(xintercept=0),linetype="dashed",color="grey")+
    geom_path(data = circ,aes(x,y), lty = 2, color = "grey", alpha = 0.7)+
    # geom_segment(data=sprp0,aes(x=0,xend=PC1,y=0,yend=PC2),
    #              arrow = arrow(length = unit(0.025, "npc"), type = "open"),
    #              lwd = .5)+
    # geom_text(data = sprp0,
    #           aes(x = PC1*1.1, y =  PC2*1.1,
    #               label = taxa),
    #           check_overlap = T, size = 3) +
    stat_ellipse(geom="polygon", 
                 aes(x=PC1,y=PC2,fill = treatment),
                 data=a.env0p,
                 alpha = 0.2, 
                 show.legend = FALSE,
                 level = 0.95)+
    #    geom_polygon(data = hull0, aes(x=PC1,y=PC2,color=treatment,fill=treatment),alpha = 0.1)+
    geom_point(data=a.env0p,aes(x=PC1,y=PC2,color=treatment),size=2)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          legend.position = c(.2,0.93),
          legend.background = element_blank())+
    coord_fixed()+
    xlab("PC1 40.09")+
    ylab("PC2 32.57%")+
    scale_color_viridis_d(option="A",begin=0,end=0.6,"")+
    scale_fill_viridis_d(option="A",begin=0,end=0.6,""))

(a0hull<-ggplot()+
    ylim(-1,1)+
    xlim(-1,1)+
    geom_hline(aes(yintercept=0),linetype="dashed",color="grey")+
    geom_vline(aes(xintercept=0),linetype="dashed",color="grey")+
    geom_path(data = circ,aes(x,y), lty = 2, color = "grey", alpha = 0.7)+
    # geom_segment(data=sprp0,aes(x=0,xend=PC1,y=0,yend=PC2),
    #              arrow = arrow(length = unit(0.025, "npc"), type = "open"),
    #              lwd = .5)+
    # geom_text(data = sprp0,
    #           aes(x = PC1*1.1, y =  PC2*1.1,
    #               label = taxa),
    #           check_overlap = T, size = 3) +
    # stat_ellipse(geom="polygon", 
    #              aes(x=PC1,y=PC2,fill = treatment),
    #              data=f.env0p,
    #              alpha = 0.2, 
  #              show.legend = FALSE,
  #              level = 0.95)+
  geom_polygon(data = hull0, aes(x=PC1,y=PC2,color=treatment,fill=treatment),alpha = 0.1)+
    geom_point(data=a.env0p,aes(x=PC1,y=PC2,color=treatment),size=2)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          legend.position = c(.2,0.93),
          legend.background = element_blank())+
    coord_fixed()+
    xlab("PC1 40.09%")+
    ylab("PC2 32.57%")+
    scale_color_viridis_d(option="A",begin=0,end=0.6,"")+
    scale_fill_viridis_d(option="A",begin=0,end=0.6,""))

# these figures are kind of informative. You can see a big divergence of the control 

# make figures for month 1
a.com1.hel<-decostand(a.com1,"hellinger")
a1.rda.null<-rda(a.com1.hel~1)
m1.scores<-scores(a1.rda.null,scaling=1,1:2)

a.env1p<-bind_cols(a.env1,data.frame(m1.scores$sites))

# look at how much variation each axis explains to add to axis labels
summary(a1.rda.null)$cont$importance

(summary(a1.rda.null)$cont$importance)

# how many eigenvalues to use in circ 
N_eigens <- ncol(summary(a1.rda.null)$cont$importance)

circ <- circleFun(center=c(0,0),diameter=sqrt(2/N_eigens),npoints = 500)

hull1 <- a.env1p %>%
  group_by(treatment)%>%
  slice(chull(PC1, PC2))

spr1<-data.frame(m1.scores$species)
sprp1<-data.frame(spr1[abs(spr1[,1])>=max(circ[,1])|abs(spr1[,2])>=max(circ[,2]),])
sprp1$taxa<-rownames(sprp1)
sprp1<-sprp1%>%
  filter(taxa %in% c("acetabularia",
    "cladocephalus",
    "laurencia",
    "penicillus"))
taxa<-rownames(sprp1)



(a1<-ggplot()+
    #    ylim(-1.2,1.2)+
    #    xlim(-1.2,1.2)+
    geom_hline(aes(yintercept=0),linetype="dashed",color="grey")+
    geom_vline(aes(xintercept=0),linetype="dashed",color="grey")+
    geom_path(data = circ,aes(x,y), lty = 2, color = "grey", alpha = 0.7)+
    geom_segment(data=sprp1,aes(x=0,xend=PC1,y=0,yend=PC2),
      arrow = arrow(length = unit(0.025, "npc"), type = "open"),
      lwd = .5)+
    geom_text(data = sprp1,
      aes(x = PC1*1.1, y =  PC2*1.1,
        label = taxa),
      check_overlap = T, size = 3) +
    stat_ellipse(geom="polygon", 
      aes(x=PC1,y=PC2,fill = treatment),
      data=a.env1p,
      alpha = 0.2, 
      show.legend = FALSE,
      level = 0.95)+
    #    geom_polygon(data = hull0, aes(x=PC1,y=PC2,color=treatment,fill=treatment),alpha = 0.1)+
    geom_point(data=a.env1p,aes(x=PC1,y=PC2,color=treatment),size=2)+
    theme_bw()+
    theme(panel.grid = element_blank(),
      legend.position = "none")+
    coord_fixed()+
    xlab("PC1 47.00%")+
    ylab("PC2 27.23%")+
    scale_color_viridis_d(option="A",begin=0,end=0.6,"")+
    scale_fill_viridis_d(option="A",begin=0,end=0.6,""))

(a1hull<-ggplot()+
    ylim(-1,1)+
    xlim(-1,1)+
    geom_hline(aes(yintercept=0),linetype="dashed",color="grey")+
    geom_vline(aes(xintercept=0),linetype="dashed",color="grey")+
    geom_path(data = circ,aes(x,y), lty = 2, color = "grey", alpha = 0.7)+
    geom_segment(data=sprp1,aes(x=0,xend=PC1,y=0,yend=PC2),
      arrow = arrow(length = unit(0.025, "npc"), type = "open"),
      lwd = .5)+
    geom_text(data = sprp1,
      aes(x = PC1*1.1, y =  PC2*1.1,
        label = taxa),
      check_overlap = T, size = 3) +
    # stat_ellipse(geom="polygon", 
    #              aes(x=PC1,y=PC2,fill = treatment),
    #              data=f.env5p,
    #              alpha = 0.2, 
    #              show.legend = FALSE,
    #              level = 0.95)+
    geom_polygon(data = hull1, aes(x=PC1,y=PC2,color=treatment,fill=treatment),alpha = 0.1)+
    geom_point(data=a.env1p,aes(x=PC1,y=PC2,color=treatment),size=2)+
    theme_bw()+
    theme(panel.grid = element_blank(),
      legend.position = "none")+
    coord_fixed()+
    xlab("PC1 47.00%")+
    ylab("PC2 27.23%")+
    scale_color_viridis_d(option="A",begin=0,end=0.6,"")+
    scale_fill_viridis_d(option="A",begin=0,end=0.6,""))



# make figures for month 5
a.com5.hel<-decostand(a.com5,"hellinger")
a5.rda.null<-rda(a.com5.hel~1)
m5.scores<-scores(a5.rda.null,scaling=1,1:2)

a.env5p<-bind_cols(a.env5,data.frame(m5.scores$sites))

# look at how much variation each axis explains to add to axis labels
summary(a5.rda.null)$cont$importance

# how many eigenvalues to use in circ 
N_eigens <- ncol(summary(a5.rda.null)$cont$importance)

circ <- circleFun(center=c(0,0),diameter=sqrt(2/N_eigens),npoints = 500)

hull5 <- a.env5p %>%
  group_by(treatment)%>%
  slice(chull(PC1, PC2))

spr5<-data.frame(m5.scores$species)
sprp5<-data.frame(spr5[abs(spr5[,1])>=max(circ[,1])|abs(spr5[,2])>=max(circ[,2]),])
sprp5$taxa<-rownames(sprp5)
sprp5<-sprp5%>%
  filter(taxa %in% c("acetabularia",
                     "cladocephalus",
                     "laurencia",
                     "penicillus"))
taxa<-rownames(sprp5)



(a5<-ggplot()+
#    ylim(-1.2,1.2)+
#    xlim(-1.2,1.2)+
    geom_hline(aes(yintercept=0),linetype="dashed",color="grey")+
    geom_vline(aes(xintercept=0),linetype="dashed",color="grey")+
    geom_path(data = circ,aes(x,y), lty = 2, color = "grey", alpha = 0.7)+
    geom_segment(data=sprp5,aes(x=0,xend=PC1,y=0,yend=PC2),
                 arrow = arrow(length = unit(0.025, "npc"), type = "open"),
                 lwd = .5)+
    geom_text(data = sprp5,
              aes(x = PC1*1.1, y =  PC2*1.1,
                  label = taxa),
              check_overlap = T, size = 3) +
    stat_ellipse(geom="polygon", 
                 aes(x=PC1,y=PC2,fill = treatment),
                 data=a.env5p,
                 alpha = 0.2, 
                 show.legend = FALSE,
                 level = 0.95)+
    #    geom_polygon(data = hull0, aes(x=PC1,y=PC2,color=treatment,fill=treatment),alpha = 0.1)+
    geom_jitter(data=a.env5p,aes(x=PC1,y=PC2,color=treatment),
      width = 0.1, height = 0.1,
      size=2)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          legend.position = "none")+
    coord_fixed()+
    xlab("PC1 47.00%")+
    ylab("PC2 27.23%")+
    scale_color_viridis_d(option="A",begin=0,end=0.6,"")+
    scale_fill_viridis_d(option="A",begin=0,end=0.6,""))

(a5hull<-ggplot()+
    ylim(-1,1)+
    xlim(-1,1)+
    geom_hline(aes(yintercept=0),linetype="dashed",color="grey")+
    geom_vline(aes(xintercept=0),linetype="dashed",color="grey")+
    geom_path(data = circ,aes(x,y), lty = 2, color = "grey", alpha = 0.7)+
    geom_segment(data=sprp5,aes(x=0,xend=PC1,y=0,yend=PC2),
                 arrow = arrow(length = unit(0.025, "npc"), type = "open"),
                 lwd = .5)+
    geom_text(data = sprp5,
              aes(x = PC1*1.1, y =  PC2*1.1,
                  label = taxa),
              check_overlap = T, size = 3) +
    # stat_ellipse(geom="polygon", 
    #              aes(x=PC1,y=PC2,fill = treatment),
    #              data=f.env5p,
    #              alpha = 0.2, 
    #              show.legend = FALSE,
    #              level = 0.95)+
    geom_polygon(data = hull5, aes(x=PC1,y=PC2,color=treatment,fill=treatment),alpha = 0.1)+
    geom_jitter(data=a.env5p,aes(x=PC1,y=PC2,color=treatment),
      # width = 0.1, 
      height = 0.1,
      size=2)+theme_bw()+
    theme(panel.grid = element_blank(),
          legend.position = "none")+
    coord_fixed()+
    xlab("PC1 47.00%")+
    ylab("PC2 27.23%")+
    scale_color_viridis_d(option="A",begin=0,end=0.6,"")+
    scale_fill_viridis_d(option="A",begin=0,end=0.6,""))

# still not hugely informative

# make figures for month 12
a.com12.hel<-decostand(a.com12,"hellinger")
a12.rda.null<-rda(a.com12.hel~1)

m12.scores<-scores(a12.rda.null,scaling=1,1:2)

a.env12p<-bind_cols(a.env12,data.frame(m12.scores$sites))

# look at how much variation each axis explains to add to axis labels
summary(a12.rda.null)$cont$importance

circ <- circleFun(center=c(0,0),diameter=sqrt(2/5),npoints = 500)

hull12 <- a.env12p %>%
  group_by(treatment)%>%
  slice(chull(PC1, PC2))

sprp12<-data.frame(m12.scores$species)
sprp12$taxa<-rownames(sprp12)
sprp12<-sprp12%>%
  filter(taxa %in% c("acetabularia",
                     "cladocephalus",
                     "laurencia",
                     "penicillus"))
taxa<-rownames(sprp12)


(a12<-ggplot()+
    ylim(-1.2,1.2)+
    xlim(-1.2,1.2)+
    geom_hline(aes(yintercept=0),linetype="dashed",color="grey")+
    geom_vline(aes(xintercept=0),linetype="dashed",color="grey")+
    geom_path(data = circ,aes(x,y), lty = 2, color = "grey", alpha = 0.7)+
    geom_segment(data=sprp12,aes(x=0,xend=PC1,y=0,yend=PC2),
                 arrow = arrow(length = unit(0.025, "npc"), type = "open"),
                 lwd = .5)+
    geom_text(data = sprp12,
              aes(x = PC1*1.1, y =  PC2*1.1,
                  label = taxa),
              check_overlap = T, size = 3) +
    stat_ellipse(geom="polygon", 
                 aes(x=PC1,y=PC2,fill = treatment),
                 data=a.env12p,
                 alpha = 0.2, 
                 show.legend = FALSE,
                 level = 0.95)+
    #    geom_polygon(data = hull0, aes(x=PC1,y=PC2,color=treatment,fill=treatment),alpha = 0.1)+
    geom_jitter(data=a.env12p,aes(x=PC1,y=PC2,color=treatment),
      width = 0.05, height = 0.05,
      size=2)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          legend.position = "none")+
    coord_fixed()+
    xlab("PC1 48.99%")+
    ylab("PC2 21.87%")+
    scale_color_viridis_d(option="A",begin=0,end=0.6,"")+
    scale_fill_viridis_d(option="A",begin=0,end=0.6,""))

(a12hull<-ggplot()+
    ylim(-1,1)+
    xlim(-1,1)+
    geom_hline(aes(yintercept=0),linetype="dashed",color="grey")+
    geom_vline(aes(xintercept=0),linetype="dashed",color="grey")+
    geom_path(data = circ,aes(x,y), lty = 2, color = "grey", alpha = 0.7)+
    geom_segment(data=sprp12,aes(x=0,xend=PC1,y=0,yend=PC2),
                 arrow = arrow(length = unit(0.025, "npc"), type = "open"),
                 lwd = .5)+
    geom_text(data = sprp12,
              aes(x = PC1*1.1, y =  PC2*1.1,
                  label = taxa),
              check_overlap = T, size = 3) +
    # stat_ellipse(geom="polygon", 
    #              aes(x=PC1,y=PC2,fill = treatment),
    #              data=f.env12p,
    #              alpha = 0.2, 
    #              show.legend = FALSE,
    #              level = 0.95)+
    geom_polygon(data = hull12, aes(x=PC1,y=PC2,color=treatment,fill=treatment),alpha = 0.1)+
    geom_jitter(data=a.env12p,aes(x=PC1,y=PC2,color=treatment),
      width = 0.05, height = 0.05,
      size=2)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          legend.position = "none")+
    coord_fixed()+
    xlab("PC1 48.99%")+
    ylab("PC2 21.87%")+
    scale_color_viridis_d(option="A",begin=0,end=0.6,"")+
    scale_fill_viridis_d(option="A",begin=0,end=0.6,""))


# make figures for month 17

a.com17.hel<-decostand(a.com17,"hellinger")
a17.rda.null<-rda(a.com17.hel~1)
m17.scores<-scores(a17.rda.null,scaling=3,1:2)

a.env17p<-bind_cols(a.env17,data.frame(m17.scores$sites))

# look at how much variation each axis explains to add to axis labels
summary(a17.rda.null)$cont$importance

circ <- circleFun(center=c(0,0),diameter=sqrt(2/6),npoints = 500)

hull17 <- a.env17p %>%
  group_by(treatment)%>%
  slice(chull(PC1, PC2))

sprp17<-data.frame(m17.scores$species)
sprp17$taxa<-rownames(sprp17)
sprp17<-sprp17%>%
  filter(taxa %in% c("acetabularia",
                     "cladocephalus",
                     "laurencia",
                     "penicillus"))
taxa<-rownames(sprp17)



(a17<-ggplot()+
    ylim(-1.2,1.2)+
    xlim(-1.2,1.2)+
    geom_hline(aes(yintercept=0),linetype="dashed",color="grey")+
    geom_vline(aes(xintercept=0),linetype="dashed",color="grey")+
    geom_path(data = circ,aes(x,y), lty = 2, color = "grey", alpha = 0.7)+
    geom_segment(data=sprp17,aes(x=0,xend=PC1,y=0,yend=PC2),
                 arrow = arrow(length = unit(0.025, "npc"), type = "open"),
                 lwd = .5)+
    geom_text(data = sprp17,
              aes(x = PC1*1.1, y =  PC2*1.1,
                  label = taxa),
              check_overlap = T, size = 3) +
    stat_ellipse(geom="polygon", 
                 aes(x=PC1,y=PC2,fill = treatment),
                 data=a.env17p,
                 alpha = 0.2, 
                 show.legend = FALSE,
                 level = 0.95)+
    #    geom_polygon(data = hull17, aes(x=PC1,y=PC2,color=treatment,fill=treatment),alpha = 0.1)+
    geom_point(data=a.env17p,aes(x=PC1,y=PC2,color=treatment),size=2)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          legend.position = "none")+
    coord_fixed()+
    xlab("PC1 56.99%")+
    ylab("PC2 20.89%")+
    scale_color_viridis_d(option="A",begin=0,end=0.6,"")+
    scale_fill_viridis_d(option="A",begin=0,end=0.6,""))

(a17hull<-ggplot()+
    ylim(-1,1)+
    xlim(-1,1)+
    geom_hline(aes(yintercept=0),linetype="dashed",color="grey")+
    geom_vline(aes(xintercept=0),linetype="dashed",color="grey")+
    geom_path(data = circ,aes(x,y), lty = 2, color = "grey", alpha = 0.7)+
    geom_segment(data=sprp17,aes(x=0,xend=PC1,y=0,yend=PC2),
                 arrow = arrow(length = unit(0.025, "npc"), type = "open"),
                 lwd = .5)+
    geom_text(data = sprp17,
              aes(x = PC1*1.1, y =  PC2*1.1,
                  label = taxa),
              check_overlap = T, size = 3) +
    # stat_ellipse(geom="polygon", 
    #              aes(x=PC1,y=PC2,fill = treatment),
    #              data=f.env17p,
    #              alpha = 0.2, 
    #              show.legend = FALSE,
    #              level = 0.95)+
    geom_polygon(data = hull17, aes(x=PC1,y=PC2,color=treatment,fill=treatment),alpha = 0.1)+
    geom_point(data=a.env17p,aes(x=PC1,y=PC2,color=treatment),size=2)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          legend.position = "none")+
    coord_fixed()+
    xlab("PC1 56.99%")+
    ylab("PC2 20.89%")+
    scale_color_viridis_d(option="A",begin=0,end=0.6,"")+
    scale_fill_viridis_d(option="A",begin=0,end=0.6,""))

# can go back and look at which species are consistently important

#devtools::install_github("thomasp85/patchwork")
library(patchwork)
a0+a5+a12+a17+plot_layout(widths = 1,heights = 1)

a0hull+a5hull+a12hull+a17hull+plot_layout(widths = 1,heights = 1)

# now look at univariate results

a.env<-algae %>% select(treatment, plot, yr, sampling, season)%>%
  left_join(sg)%>%
  left_join((sggrow))

a.com<-algae %>% select(-treatment, -plot, -yr, -sampling, -season)


a.env$spr<-specnumber(a.com)
a.env$div<-diversity(a.com)
a.env$plot<-as.factor(a.env$plot)
a.env$treatment<-as.factor(a.env$treatment)
a.env$season<-as.factor(a.env$season)


a.env0uni<-a.env0%>%
  mutate(strt.spr=specnumber(a.com0),
         strt.div=diversity(a.com0),
         plot=factor(plot))%>%
  select(treatment,plot,strt.spr,strt.div)

a.env.uni<-a.env%>%
  # filter(sampling!=0)%>%
  left_join(a.env0uni) %>%
  mutate(spr.change = spr - strt.spr)

a.env.uni$treatment<-as.factor(a.env.uni$treatment)

library(glmmTMB)
library(lmerTest)

hist(a.env.uni$spr)

hist(a.env.uni$spr.change)
ggplot(a.env.uni, aes(spr.change)) + geom_histogram() + facet_wrap(~treatment)

spr.lmer<-lmer(spr.change~treatment*yr*season +#sg.sd + grow + 
    (1|plot),
  data = a.env.uni%>%
    mutate(treatment=relevel(treatment, ref = "real")))
summary(spr.lmer)


spr.glmm<-glmmTMB(spr~treatment*yr + season + #sg.sd + grow + 
              (1|plot)+
              offset(log(strt.spr+1)),
              # family = poisson,
              family = nbinom2(link = "log"),
               data = a.env.uni%>%
                 mutate(treatment=relevel(treatment, ref = "real")))

summary(spr.glmm)

if(!require(DHARMa))install.packages("DHARMa");library(DHARMa)

# look at residuals
spr.glmm_simres <- simulateResiduals(spr.glmm)
testDispersion(spr.glmm_simres)
plot(spr.glmm_simres)



# looks like species richness is NOT significantly different between treatments


a.sum<-a.env.uni%>%
  group_by(treatment,sampling)%>%
  summarize(spr.m=mean(spr),div.m=mean(div),spr.sd=sd(spr),div.sd=sd(div))

a.sum$treatment<-factor(a.sum$treatment,labels=c("Control","Structure Control","Sponge"))

ggplot()+
  geom_line(data=a.sum,aes(x=sampling,y=spr.m,group=treatment,color=treatment),position=position_dodge(0.5))+
  geom_errorbar(data=a.sum,aes(x=sampling,ymin=spr.m-spr.sd,ymax=spr.m+spr.sd,color=treatment),width=.5,position=position_dodge(0.5))+
  geom_point(data=a.sum,aes(x=sampling,y=spr.m,color=treatment),size=5,position=position_dodge(0.5))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        legend.text = element_text(size=10))+
  scale_color_viridis_d(option="A", begin=0, end=0.6,"")+
  ylab("Taxa Richness")+
  xlab("Months into Experiment")

ggsave("figures/algae_sp_richness.jpg")

# now look at diversity

div.lmer<-lmer(div~treatment*as.factor(sampling) + sg.sd+grow+(1|plot)+
                 offset(strt.div),
               data = a.env.uni%>%
                 mutate(treatment=relevel(treatment, ref = "real")))
summary(div.lmer)
anova(div.lmer)
#diversity is significantly different among treatments

ggplot()+
  geom_line(data=a.sum,aes(x=sampling,y=div.m,group=treatment,color=treatment),position=position_dodge(0.5))+
  geom_errorbar(data=a.sum,aes(x=sampling,ymin=div.m-div.sd,ymax=div.m+div.sd,color=treatment),width=.5,position=position_dodge(0.5))+
  geom_point(data=a.sum,aes(x=sampling,y=div.m,color=treatment),size=5,position=position_dodge(0.5))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        legend.text = element_text(size=10))+
  scale_color_viridis_d(option="A", begin=0, end=0.6,"")+
  ylab("Diversity")+
  xlab("Months into Experiment")

ggsave("figures/algal_sp_diversity.jpg")



# we probably made this before but it's usefull here for thinking about the above patterns
alg2<-readxl::read_xlsx(here("Original_data","ForFinella_Transplant_data.xlsx"),sheet="Algae")

a.taxa.abun<-alg2 %>% # this tells R that I'd like to work with this data set
  pivot_longer(c(-Treatment,-plot,-sampling),names_to = "taxa",values_to = "abundance")%>% 
  rename(treatment=Treatment)%>%
  mutate(sampling=case_when(
    sampling==1~0,
    sampling==2~1,
    sampling==3~5,
    sampling==4~12,
    sampling==5~17),
    taxa=ifelse(taxa=="cladocephalus","udotea",taxa))%>%
  group_by(treatment,sampling,taxa)%>%
  summarize(abun.m=mean(abundance),abun.sd=sd(abundance))

ggplot()+
  geom_line(data=a.taxa.abun,
    aes(x=sampling,y=abun.m,group=treatment,color=treatment),position=position_dodge(0.5))+
  geom_errorbar(data=a.taxa.abun,aes(x=sampling,ymin=abun.m-abun.sd,ymax=abun.m+abun.sd,color=treatment),width=.5,position=position_dodge(0.5))+
  geom_point(data=a.taxa.abun,aes(x=sampling,y=abun.m,color=treatment),size=5,position=position_dodge(0.5))+
  facet_wrap(~taxa, ncol=1, scales = "free_y") +
  theme_bw()+
  theme(panel.grid = element_blank(),
    axis.text = element_text(size=12),
    axis.title = element_text(size=14),
    legend.text = element_text(size=10))+
  scale_color_viridis_d(option="A", begin=0, end=0.6,"")+
  ylab("Adundance")+
  xlab("Months into Experiment")

ggsave("figures/algal_sp_abundance.jpg")

# seems like there was a shift away from cyanobacteria towards halimeda across the board 
# (but doing best in sponge plots)
# meanwhile several other algae species also took off only in sponge plots


