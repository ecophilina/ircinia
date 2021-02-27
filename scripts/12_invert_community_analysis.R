# this script will explore the invert data

library(tidyverse)
if(!require(vegan))install.packages("vegan"); library(vegan)
# if(!require(ggpointdensity))devtools::install_github("LKremer/ggpointdensity"); library(ggpointdensity)
if(!require(ggrepel))install.packages("ggrepel"); library(ggrepel)

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
  geom_jitter(aes(x=PC1,y=PC2,
    color=treatment),size=2,alpha=.5, 
    width = 0.02, height = 0.02)+
  scale_color_viridis_d(option="B",end=.8)+
  facet_wrap(~sampling) + 
  ggsidekick::theme_sleek()

# PC2 and PC3 
ggplot(data=i.scores)+
  geom_jitter(aes(x=PC2,y=PC3,
                 color=treatment),size=2,alpha=.65, 
    width = 0.02, height = 0.02)+
  scale_color_viridis_d(option="B",end=.8)+
  facet_wrap(~sampling)+ 
  ggsidekick::theme_sleek()

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

# confirming this with an anosim
(i0.anosim<-anosim(i.com.pa.0,i.env0$treatment,permutations = 999,distance="bray"))
summary(i0.anosim)
# the anosim results confirm this. 

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

#no

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


# first going to make sure where there's a difference between treatments using anosim

# look at 1 month into experiment
i.com1<-i.com.pa[i.env$sampling==1,]
i.com1<-i.com1[,colSums(i.com1)!=0]
i.env1<-i.env[i.env$sampling==1,]

(i1.anosim<-anosim(i.com1,i.env1$treatment))
# no difference

# look at 5 months into experiment
i.com5<-i.com.pa[i.env$sampling==5,]
i.com5<-i.com5[,colSums(i.com5)!=0]
i.env5<-i.env[i.env$sampling==5,]

(i5.anosim<-anosim(i.com5,i.env5$treatment))
# Now we're getting more robust differences between treatments

# look at 12 months into experiment
i.com12<-i.com.pa[i.env$sampling==12,]
i.com12<-i.com12[,colSums(i.com12)!=0]
i.env12<-i.env[i.env$sampling==12,]

(i12.anosim<-anosim(i.com12,i.env12$treatment))

# not as different as at month 5 but still significantly different.

# look at 17 months into experiment
i.com17<-i.com.pa[i.env$sampling==17,]
i.com17<-i.com17[,colSums(i.com17)!=0]
i.env17<-i.env[i.env$sampling==17,]

(i17.anosim<-anosim(i.com17,i.env17$treatment))
# and again back to being solidly different. 

#this corresponds with the results of the RDA that treatment and season are important

# treatment stays an important factor from 5 months on. 
# at 5 months which treatments are different from each other?

i.com5.bf<-i.com5[i.env5$treatment!="real",]
i.env5.bf<-i.env5[i.env5$treatment!="real",]

(i5bf.anosim<-anosim(i.com5.bf,i.env5.bf$treatment))

# blank and fake aren't different

i.com5.br<-i.com5[i.env5$treatment!="fake",]
i.env5.br<-i.env5[i.env5$treatment!="fake",]

(i5br.anosim<-anosim(i.com5.br,i.env5.br$treatment))

# close but blank and real are not actually different

i.com5.fr<-i.com5[i.env5$treatment!="blank",]
i.env5.fr<-i.env5[i.env5$treatment!="blank",]

(i5fr.anosim<-anosim(i.com5.fr,i.env5.fr$treatment))

# real and fake are different


# at 12 months which treatments are different from each other?

i.com12.bf<-i.com12[i.env12$treatment!="real",]
i.env12.bf<-i.env12[i.env12$treatment!="real",]

(i12bf.anosim<-anosim(i.com12.bf,i.env12.bf$treatment))

# blank and fake aren't different

i.com12.br<-i.com12[i.env12$treatment!="fake",]
i.env12.br<-i.env12[i.env12$treatment!="fake",]

(i12br.anosim<-anosim(i.com12.br,i.env12.br$treatment))

# blank and real are not different

i.com12.fr<-i.com12[i.env12$treatment!="blank",]
i.env12.fr<-i.env12[i.env12$treatment!="blank",]

(i12fr.anosim<-anosim(i.com12.fr,i.env12.fr$treatment))

# real and fake are different

# at 17 months which treatments are different from each other?

i.com17.bf<-i.com17[i.env17$treatment!="real",]
i.env17.bf<-i.env17[i.env17$treatment!="real",]

(i17bf.anosim<-anosim(i.com17.bf,i.env17.bf$treatment))

# blank and fake are different

i.com17.br<-i.com17[i.env17$treatment!="fake",]
i.env17.br<-i.env17[i.env17$treatment!="fake",]

(i17br.anosim<-anosim(i.com17.br,i.env17.br$treatment))

# blank and real ARE different

i.com17.fr<-i.com17[i.env17$treatment!="blank",]
i.env17.fr<-i.env17[i.env17$treatment!="blank",]

(i17fr.anosim<-anosim(i.com17.fr,i.env17.fr$treatment))

# real and fake are different

# so sponge and structure control diverge from each other early on but neither is different
# from the control until 17 months until all treatments have a different community


# which species are driving this difference?
(i17.simper<-simper(i.com17,i.env17$treatment,permutations = 999))
summary(i17.simper)

# looks like ceriths, oysters, sea cucumbers, anemones, and blue crabs are both significant
# and consistently contribute to differences between the groups here

icom17diff<-bind_cols(i.com17[,colnames(i.com17) %in% c("blue crab",
                                              "Anemone",
                                              "cerith",
                                              "sea cucumber",
                                              "oyster")],
                      i.env17)%>%
  pivot_longer(1:5,names_to="taxa",values_to="presence")%>%
  group_by(treatment,taxa,sampling)%>%
  summarize(n.plots = sum(presence))%>%
  mutate(treatment=factor(treatment,levels=c("blank","fake","real"),labels = c("Control","Structure Control","Sponge")))

ggplot(data=icom17diff)+
  geom_col(aes(x=taxa,y=n.plots,fill=treatment),width=.5,
             position=position_dodge(.5))+
  scale_fill_viridis_d(option="A",begin=0,end=0.6,"")+
  # geom_point(aes(x=taxa, y=n.plots,color=treatment),
  #            position=position_dodge(.5),size=5)+
  scale_color_viridis_d(option="A",begin=0,end=0.6,"")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        legend.position = "top")+
  ylab("Number of plots where taxa is present")
  



# make figures for this

anova(i0.rda.null, permutations=how(nperm=999))
# Tests of all canonical axes
anova(i0.rda.null, by="axis", permutations=how(nperm=999))
init.scores<-scores(rda(i.com.hel0),scaling=1,1:2)

i.env0p<-bind_cols(i.env0,data.frame(init.scores$sites))%>%
  mutate(treatment=factor(treatment, labels=c("Control","Structure Control","Sponge")))

hull0 <- i.env0p %>%
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
circ <- circleFun(center=c(0,0),diameter=sqrt(2/13),npoints = 500)


spr0<-data.frame(init.scores$species)
sprp0<-data.frame(spr0[abs(spr0[,1])>=max(circ[,1])|abs(spr0[,2])>=max(circ[,2]),])
sprp0$taxa<-rownames(sprp0)
sprp0<-sprp0%>%
  filter(taxa %in% c("Anemone","blue crab","cerith","little white snail","mantis shrimp"))
taxa<-rownames(sprp0)

# look at how much variation each axis explains to add to axis labels
summary(i0.rda.null)$cont$importance

(i0<-ggplot()+
    ylim(-1.2,1.2)+
    xlim(-1.2,1.2)+
    geom_hline(aes(yintercept=0),linetype="dashed",color="grey")+
    geom_vline(aes(xintercept=0),linetype="dashed",color="grey")+
    geom_path(data = circ,aes(x,y), lty = 2, color = "grey", alpha = 0.7)+
    geom_segment(data=sprp0,aes(x=0,xend=PC1,y=0,yend=PC2),
                 arrow = arrow(length = unit(0.025, "npc"), type = "open"),
                 lwd = .5)+
    geom_text(data = sprp0,
              aes(x = PC1*1.1, y =  PC2*1.1,
                  label = taxa),
              check_overlap = T, size = 3) +
    
    stat_ellipse(geom="polygon", 
                 aes(x=PC1,y=PC2,fill = treatment),
                 data=i.env0p,
                 alpha = 0.1, 
                 show.legend = FALSE,
                 level = 0.95)+
    stat_ellipse(geom="polygon", aes(x=PC1,y=PC2,fill = treatment),
      data=i.env0p, alpha = 0.2, show.legend = FALSE,
      level = 0.75)+
    stat_ellipse(geom="polygon", aes(x=PC1,y=PC2,fill = treatment),
      data=i.env0p,alpha = 0.1, show.legend = FALSE,
      level = 0.8)+
    stat_ellipse(geom="polygon", aes(x=PC1,y=PC2,fill = treatment), 
      data=i.env0p,alpha = 0.1, show.legend = FALSE,
      level = 0.85)+
    stat_ellipse(geom="polygon", aes(x=PC1,y=PC2,fill = treatment),
      data=i.env0p, alpha = 0.1, show.legend = FALSE,
      level = 0.9)+
#    geom_polygon(data = hull0, aes(x=PC1,y=PC2,color=treatment,fill=treatment),alpha = 0.1)+
    geom_point(data=i.env0p,aes(x=PC1,y=PC2,color=treatment),size=2)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          legend.position = c(.14,0.95),
          legend.background = element_blank())+
     coord_fixed()+
    xlab("PC1 23.74%")+
    ylab("PC2 16.54%")+
    scale_color_viridis_d(option="A",begin=0,end=0.6,"")+
    scale_fill_viridis_d(option="A",begin=0,end=0.6,""))

(i0hull<-ggplot()+
    ylim(-1,1)+
    xlim(-1,1)+
    geom_hline(aes(yintercept=0),linetype="dashed",color="grey")+
    geom_vline(aes(xintercept=0),linetype="dashed",color="grey")+
    geom_path(data = circ,aes(x,y), lty = 2, color = "grey", alpha = 0.7)+
    geom_segment(data=sprp0,aes(x=0,xend=PC1,y=0,yend=PC2),
                 arrow = arrow(length = unit(0.025, "npc"), type = "open"),
                 lwd = .5)+
    geom_text(data = sprp0,
              aes(x = PC1*1.1, y =  PC2*1.1,
                  label = taxa),
              check_overlap = T, size = 3) +
    # stat_ellipse(geom="polygon", 
    #              aes(x=PC1,y=PC2,fill = treatment),
    #              data=i.env0p,
    #              alpha = 0.2, 
    #              show.legend = FALSE,
    #              level = 0.95)+
    geom_polygon(data = hull0, 
         aes(x=PC1,y=PC2,color=treatment,fill=treatment),alpha = 0.1)+
    geom_point(data=i.env0p,aes(x=PC1,y=PC2,color=treatment),size=2)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          legend.position = c(.14,0.95),
          legend.background = element_blank())+
    coord_fixed()+
    xlab("PC1 23.74%")+
    ylab("PC2 16.54%")+
    scale_color_viridis_d(option="A",begin=0,end=0.6,"")+
    scale_fill_viridis_d(option="A",begin=0,end=0.6,""))





# make figures for month 5
i.com5.hel<-decostand(i.com5,"hellinger")
i5.rda.null<-rda(i.com5.hel~1)
m5.scores<-scores(i5.rda.null,scaling=1,1:2)

i.env5p<-bind_cols(i.env5,data.frame(m5.scores$sites))

# look at how much variation each axis explains to add to axis labels
summary(i5.rda.null)$cont$importance

circ <- circleFun(center=c(0,0),diameter=sqrt(2/12),npoints = 500)

hull5 <- i.env5p %>%
  group_by(treatment)%>%
  slice(chull(PC1, PC2))

spr5<-data.frame(m5.scores$species)
sprp5<-data.frame(spr5[abs(spr5[,1])>=max(circ[,1])|abs(spr5[,2])>=max(circ[,2]),])
sprp5$taxa<-rownames(sprp5)
sprp5<-sprp5%>%
  filter(taxa %in% c("Anemone","blue crab","cerith","little white snail","mantis shrimp"))
taxa<-rownames(sprp5)



(i5<-ggplot()+
    ylim(-1,1)+
    xlim(-1,1)+
    geom_hline(aes(yintercept=0),linetype="dashed",color="grey")+
    geom_vline(aes(xintercept=0),linetype="dashed",color="grey")+
    geom_path(data = circ,aes(x,y), lty = 2, color = "grey", alpha = 0.7)+
    geom_segment(data=sprp5,aes(x=0,xend=PC1,y=0,yend=PC2),
                 arrow = arrow(length = unit(0.025, "npc"), type = "open"),
                 lwd = .5)+
    # geom_text(check_overlap = T, 
  geom_text_repel(min.segment.length = 5, #point.padding = 0.2, direction = "both",
      nudge_x = 0.17, nudge_y = 0.07,
      aes(x = PC1, y =  PC2, label = taxa),
      size = 3, data=sprp5) +
    stat_ellipse(geom="polygon", aes(x=PC1,y=PC2,fill = treatment),
      data=i.env5p, alpha = 0.2, show.legend = FALSE,
      level = 0.75)+
    stat_ellipse(geom="polygon", aes(x=PC1,y=PC2,fill = treatment),
      data=i.env5p,alpha = 0.1, show.legend = FALSE,
      level = 0.8)+
    stat_ellipse(geom="polygon", aes(x=PC1,y=PC2,fill = treatment), 
      data=i.env5p,alpha = 0.1, show.legend = FALSE,
      level = 0.85)+
    stat_ellipse(geom="polygon", aes(x=PC1,y=PC2,fill = treatment),
      data=i.env5p, alpha = 0.1, show.legend = FALSE,
      level = 0.9)+
    # stat_ellipse(geom="polygon", aes(x=PC1,y=PC2,fill = treatment),
    #   data=i.env5p, alpha = 0.1, show.legend = FALSE,
    #   level = 0.95)+
    #    geom_polygon(data = hull0, aes(x=PC1,y=PC2,color=treatment,fill=treatment),alpha = 0.1)+
    geom_jitter(data=i.env5p,aes(x=PC1,y=PC2,color=treatment),
      size=2, alpha = 0.75, height = 0.03, width = 0.03)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          legend.position = "none")+
    coord_fixed()+
    xlab("PC1 23.45%")+
    ylab("PC2 19.19%")+
    scale_color_viridis_d(option="A",begin=0,end=0.6,"")+
    scale_fill_viridis_d(option="A",begin=0,end=0.6,""))

(i5hull<-ggplot()+
    ylim(-1,1)+
    xlim(-1,1)+
    geom_hline(aes(yintercept=0),linetype="dashed",color="grey")+
    geom_vline(aes(xintercept=0),linetype="dashed",color="grey")+
    geom_path(data = circ,aes(x,y), lty = 2, color = "grey", alpha = 0.7)+
    # geom_polygon(data = hull12, aes(x=PC1,y=PC2,color=treatment,fill=treatment),
    #   alpha = 0.3, lwd=0.1)+
    stat_density2d(geom="polygon", method="ndensity",
      aes(x=PC1, y=PC2, fill=treatment), alpha = 0.2, bins = 2,
      data=filter(i.env5p, treatment == "blank"), show.legend = FALSE) +
    stat_density2d(geom="polygon", method="ndensity",
      aes(x=PC1, y=PC2, fill = treatment), alpha = 0.2, bins = 2,
      data=filter(i.env5p, treatment=="fake"), show.legend = FALSE) +
    stat_density2d(geom="polygon", method="ndensity",
      aes(x=PC1, y=PC2, fill = treatment), alpha = 0.2, bins = 2,
      data=filter(i.env5p, treatment=="real"), show.legend = FALSE) +
    geom_jitter(data=i.env5p,aes(x=PC1, y=PC2, color=treatment), 
      width = 0.01, height = 0.01, alpha =0.75, size=2)+
    geom_segment(data=sprp5,aes(x=0, xend=PC1, y=0, yend=PC2),
      arrow = arrow(length = unit(0.025, "npc"), type = "open"),
      lwd = .5)+
    geom_text_repel(min.segment.length = 5, #point.padding = 0.2, direction = "both",
      nudge_x = 0.17, nudge_y = 0.07,
      # geom_text(check_overlap = T, 
      #   position=position_jitter(width = -0.3, height = 0.1
      #     # width=ifelse(sprp12$taxa=='little white snail', -0.2, 0),
      #     # height=ifelse(sprp12$taxa=='little white snail', -0.2, 0)
      #     ), 
      aes(x = PC1, y =  PC2, label = taxa), size = 3, data = sprp5) +
    geom_segment(data=sprp5,aes(x=0,xend=PC1,y=0,yend=PC2),
                 arrow = arrow(length = unit(0.025, "npc"), type = "open"),
                 lwd = .5)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          legend.position = "none")+
    coord_fixed()+
    xlab("PC1 23.45%")+
    ylab("PC2 19.19%")+
    scale_color_viridis_d(option="A",begin=0,end=0.6,"")+
    scale_fill_viridis_d(option="A",begin=0,end=0.6,""))

# make figures for month 12
i.com12.hel<-decostand(i.com12,"hellinger")
i12.rda.null<-rda(i.com12.hel~1)

m12.scores<-scores(i12.rda.null,scaling=1,1:2)

i.env12p<-bind_cols(i.env12,data.frame(m12.scores$sites))

# look at how much variation each axis explains to add to axis labels
summary(i12.rda.null)$cont$importance

circ <- circleFun(center=c(0,0),diameter=sqrt(2/13),npoints = 500)

hull12 <- i.env12p %>%
  group_by(treatment)%>%
  slice(chull(PC1, PC2))

sprp12<-data.frame(m12.scores$species)
sprp12$taxa<-rownames(sprp12)
sprp12<-sprp12%>%
  filter(taxa %in% c("Anemone","blue crab","cerith","little white snail","mantis shrimp"))
taxa<-rownames(sprp12)


(i12<-ggplot()+
    ylim(-1,1)+
    xlim(-1,1)+
    geom_hline(aes(yintercept=0),linetype="dashed",color="grey")+
    geom_vline(aes(xintercept=0),linetype="dashed",color="grey")+
    geom_path(data = circ,aes(x,y), lty = 2, color = "grey", alpha = 0.7)+
    geom_segment(data=sprp12,aes(x=0,xend=PC1,y=0,yend=PC2),
                 arrow = arrow(length = unit(0.025, "npc"), type = "open"),
                 lwd = .5)+
    # geom_text(check_overlap = T, 
      geom_text_repel(min.segment.length = 5, #point.padding = 0.2, direction = "both",
        nudge_x = -0.17, nudge_y = -0.07,
              aes(x = PC1, y =  PC2, label = taxa),
              size = 3, data=sprp12) +
    stat_ellipse(geom="polygon", aes(x=PC1,y=PC2,fill = treatment),
      data=i.env12p, alpha = 0.2, show.legend = FALSE,
      level = 0.75)+
    stat_ellipse(geom="polygon", aes(x=PC1,y=PC2,fill = treatment),
      data=i.env12p,alpha = 0.1, show.legend = FALSE,
      level = 0.8)+
    stat_ellipse(geom="polygon", aes(x=PC1,y=PC2,fill = treatment), 
      data=i.env12p,alpha = 0.1, show.legend = FALSE,
      level = 0.85)+
    stat_ellipse(geom="polygon", aes(x=PC1,y=PC2,fill = treatment),
      data=i.env12p, alpha = 0.1, show.legend = FALSE,
      level = 0.9)+
    # stat_ellipse(geom="polygon", aes(x=PC1,y=PC2,fill = treatment),
      # data=i.env12p, alpha = 0.2, show.legend = FALSE,
      # level = 0.95)+
    #    geom_polygon(data = hull0, aes(x=PC1,y=PC2,color=treatment,fill=treatment),alpha = 0.1)+
    geom_point(data=i.env12p,aes(x=PC1,y=PC2,color=treatment),size=2)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          legend.position = "none")+
    coord_fixed()+
    xlab("PC1 24.11%")+
    ylab("PC2 18.83%")+
    scale_color_viridis_d(option="A",begin=0,end=0.6,"")+
    scale_fill_viridis_d(option="A",begin=0,end=0.6,""))

(i12hull<-ggplot()+
    ylim(-1,1)+
    xlim(-1,1)+
    geom_hline(aes(yintercept=0),linetype="dashed",color="grey")+
    geom_vline(aes(xintercept=0),linetype="dashed",color="grey")+
    geom_path(data = circ,aes(x,y), lty = 2, color = "grey", alpha = 0.7)+
    # geom_polygon(data = hull12, aes(x=PC1,y=PC2,color=treatment,fill=treatment),
    #   alpha = 0.3, lwd=0.1)+
    stat_density2d(geom="polygon", method="ndensity",
      aes(x=PC1, y=PC2, fill=treatment), alpha = 0.2, bins = 2,
      data=filter(i.env12p, treatment == "blank"), show.legend = FALSE) +
    stat_density2d(geom="polygon", method="ndensity",
      aes(x=PC1, y=PC2, fill = treatment), alpha = 0.2, bins = 2,
      data=filter(i.env12p, treatment=="fake"), show.legend = FALSE) +
    stat_density2d(geom="polygon", method="ndensity",
      aes(x=PC1, y=PC2, fill = treatment), alpha = 0.2, bins = 2,
      data=filter(i.env12p, treatment=="real"), show.legend = FALSE) +
    geom_jitter(data=i.env12p,aes(x=PC1, y=PC2, color=treatment), 
      width = 0.01, height = 0.01, alpha =0.75, size=2)+
    geom_segment(data=sprp12,aes(x=0, xend=PC1, y=0, yend=PC2),
                 arrow = arrow(length = unit(0.025, "npc"), type = "open"),
                 lwd = .5)+
    geom_text_repel(min.segment.length = 5, #point.padding = 0.2, direction = "both",
      nudge_x = -0.17, nudge_y = -0.07,
    # geom_text(check_overlap = T, 
    #   position=position_jitter(width = -0.3, height = 0.1
    #     # width=ifelse(sprp12$taxa=='little white snail', -0.2, 0),
    #     # height=ifelse(sprp12$taxa=='little white snail', -0.2, 0)
    #     ), 
      aes(x = PC1, y =  PC2, label = taxa), size = 3, data = sprp12) +
    theme_bw()+
    theme(panel.grid = element_blank(),
          legend.position = "none")+
    coord_fixed()+
    xlab("PC1 24.11%")+
    ylab("PC2 18.83%")+
    scale_color_viridis_d(option="A",begin=0,end=0.6,"")+
    scale_fill_viridis_d(option="A",begin=0,end=0.6,""))


# make figures for month 17

i.com17.hel<-decostand(i.com17,"hellinger")
i17.rda.null<-rda(i.com17.hel~1)
m17.scores<-scores(i17.rda.null,scaling=3,1:2)

i.env17p<-bind_cols(i.env17,data.frame(m17.scores$sites))

# look at how much variation each axis explains to add to axis labels
summary(i17.rda.null)$cont$importance

circ <- circleFun(center=c(0,0),diameter=sqrt(2/13),npoints = 500)

hull17 <- i.env17p %>%
  group_by(treatment)%>%
  slice(chull(PC1, PC2))

sprp17<-data.frame(m17.scores$species)
sprp17$taxa<-rownames(sprp17)
sprp17<-sprp17%>%
  filter(taxa %in% c("Anemone","blue crab","cerith","little white snail","mantis shrimp"))
taxa<-rownames(sprp17)

(i17<-ggplot()+
    ylim(-1,1)+
    xlim(-1,1.15)+
    geom_hline(aes(yintercept=0),linetype="dashed",color="grey")+
    geom_vline(aes(xintercept=0),linetype="dashed",color="grey")+
    geom_path(data = circ,aes(x,y), lty = 2, color = "grey", alpha = 0.7)+
  stat_ellipse(geom="polygon", aes(x=PC1,y=PC2,fill = treatment),
    data=i.env17p, alpha = 0.2, show.legend = FALSE,
    level = 0.75)+
  stat_ellipse(geom="polygon", aes(x=PC1,y=PC2,fill = treatment),
    data=i.env17p,alpha = 0.1, show.legend = FALSE,
    level = 0.8)+
  stat_ellipse(geom="polygon", aes(x=PC1,y=PC2,fill = treatment),
    data=i.env17p,alpha = 0.1, show.legend = FALSE,
    level = 0.85)+
  stat_ellipse(geom="polygon", aes(x=PC1,y=PC2,fill = treatment),
    data=i.env17p, alpha = 0.1, show.legend = FALSE,
    level = 0.9)+
  # stat_ellipse(geom="polygon", aes(x=PC1,y=PC2,fill = treatment),
  #   data=i.env17p, alpha = 0.2, show.legend = FALSE,
  #   level = 0.95)+
  # geom_polygon(data = hull17, aes(x=PC1,y=PC2,color=treatment,fill=treatment),
  #   alpha = 0.1)+
  geom_jitter(data=i.env17p, aes(x=PC1, y=PC2, color=treatment), 
    width = 0.01, height = 0.01, alpha =0.75, size=2)+
  geom_segment(data=sprp17,aes(x=0,xend=PC1,y=0,yend=PC2),
    arrow = arrow(length = unit(0.025, "npc"), type = "open"),
    lwd = .5)+
  geom_text_repel(min.segment.length = 5, #point.padding = 0.2, direction = "both",
    nudge_x = -0.15, nudge_y = 0.09,
    aes(x = PC1, y =  PC2, label = taxa),
  # geom_text(check_overlap = T, aes(x = PC1, y =  PC2, label = taxa),
  #   position=position_jitter(width = -0.3, height = 0.1
  #     # width=ifelse(sprp12$taxa=='little white snail', -0.2, 0),
  #     # height=ifelse(sprp12$taxa=='little white snail', -0.2, 0)
  #     ), 
      size = 3, data = sprp17) +
    theme_bw()+
    theme(panel.grid = element_blank(),
      legend.position = "none")+
    coord_fixed()+
    xlab("PC1 26.87%")+
    ylab("PC2 20.53%")+
    scale_color_viridis_d(option="A",begin=0,end=0.6,"")+
    scale_fill_viridis_d(option="A",begin=0,end=0.6,""))


(i17<-ggplot()+
    ylim(-1,1)+
    xlim(-1,1.15)+
    geom_hline(aes(yintercept=0),linetype="dashed",color="grey")+
    geom_vline(aes(xintercept=0),linetype="dashed",color="grey")+
    geom_path(data = circ,aes(x,y), lty = 2, color = "grey", alpha = 0.7)+
    stat_density2d(geom="polygon", method="ndensity",
      aes(x=PC1, y=PC2, fill=treatment), alpha = 0.2, bins = 2,
      data=filter(i.env17p, treatment=="blank"), show.legend = FALSE) +
    stat_density2d(geom="polygon", method="ndensity",
      aes(x=PC1, y=PC2, fill=treatment), alpha = 0.2, bins = 2,
      data=filter(i.env17p, treatment=="fake"), show.legend = FALSE) +
    stat_density2d(geom="polygon", method="ndensity",
      aes(x=PC1, y=PC2, fill=treatment), alpha = 0.2, bins = 2,
      data=filter(i.env17p, treatment=="real"), show.legend = FALSE) +
    geom_jitter(data=i.env17p, aes(x=PC1, y=PC2, color=treatment), 
      width = 0.01, height = 0.01, alpha =0.75, size=2)+
    #    geom_polygon(data = hull17, aes(x=PC1,y=PC2,color=treatment,fill=treatment),alpha = 0.1)+
    geom_segment(data=sprp17,aes(x=0,xend=PC1,y=0,yend=PC2),
      arrow = arrow(length = unit(0.025, "npc"), type = "open"),
      lwd = .5)+
    geom_text_repel(min.segment.length = 5, #point.padding = 0.2, direction = "both",
      nudge_x = -0.15, nudge_y = 0.09,
      aes(x = PC1, y =  PC2, label = taxa),
      # geom_text(check_overlap = T, aes(x = PC1, y =  PC2, label = taxa),
      #   position=position_jitter(width = -0.3, height = 0.1
      #     # width=ifelse(sprp12$taxa=='little white snail', -0.2, 0),
      #     # height=ifelse(sprp12$taxa=='little white snail', -0.2, 0)
      #     ), 
      size = 3, data = sprp17) +
    theme_bw()+
    theme(panel.grid = element_blank(),
          legend.position = "none")+
    coord_fixed()+
    xlab("PC1 26.87%")+
    ylab("PC2 20.53%")+
    scale_color_viridis_d(option="A",begin=0,end=0.6,"")+
    scale_fill_viridis_d(option="A",begin=0,end=0.6,""))


# can go back and look at which species are consistently important
taxa<-bind_rows(sprp5,sprp12,sprp17)%>%
  group_by(taxa)%>%
  summarize(n=n())%>%
  filter(n==3)

devtools::install_github("thomasp85/patchwork")
library(patchwork)
i0+i5+i12+i17+plot_layout(widths = 1,heights = 1)

i0hull+i5hull+i12hull+i17hull+plot_layout(widths = 1,heights = 1)

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




