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


# fish
f.hasfish<-fish %>%
  filter(abundance!=0)

table(f.hasfish$treatment,f.hasfish$taxa)
table(f.hasfish$treatment,f.hasfish$sampling)

f2<-fish%>%
  pivot_wider(names_from=taxa,values_from=abundance,values_fill=0)

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



f0<-f2%>%filter(yr == 0)
f2<-f2%>%filter(yr != 0)

# Organize data


f.env<-f2 %>% select(treatment, plot, yr, sampling, season)%>%
  left_join(alg)%>%
  left_join(sg)%>%
  left_join((sggrow))

f.com<-f2 %>% select(-treatment, -plot, -yr, -sampling, -season)

f.com$dummy<-1

# pre-experiment data
f.env0<-f0 %>% select(treatment, plot, yr, sampling, season)%>%
  left_join(alg)%>%
  left_join(sg)%>%
  left_join((sggrow))
f.com0<-f0 %>% select(-treatment, -plot, -yr, -sampling, -season)

f.com0$dummy<-1
f.com0<-f.com0[,colSums(f.com0)!=0]

# use hellinger: square root of method to standardize species data
f.com.pa<-decostand(f.com,"pa")

f.com.hel<-decostand(f.com.pa,"hellinger")

f.pca<-rda(f.com.hel)

f.scores<-data.frame(scores(f.pca,1:3)$sites)%>%
  bind_cols(f.env)

ggplot(data=f.scores)+
  geom_jitter(aes(x=PC1,y=PC2,
                  color=treatment),size=2,alpha=.65, 
              width = 0.03, height = 0.03)+
  scale_color_viridis_d(option="B",end=.8)+
  facet_wrap(~sampling)

ggplot(data=f.scores)+
  geom_jitter(aes(x=PC2,y=PC3,
    color=treatment),size=2,alpha=.65, 
    width = 0.03, height = 0.03)+
  scale_color_viridis_d(option="B",end=.8)+
  facet_wrap(~sampling)

#start examining statistical relationship

# first look at initial pre-experiment
f.com.pa.0<-decostand(f.com0,"pa")
f.com.hel0<-decostand(f.com.pa.0,"hellinger")

f.com.pa<-decostand(f.com,"pa")
f.com.hel<-decostand(f.com.pa,"hellinger")

# look at just treatment at 0 and later
f0.rda.null<-rda(f.com.hel0~1)
trt0<-as.factor(f.env0$treatment)
trt0.mat<-data.frame(model.matrix(~ trt0, 
                                  contrasts=list(trt0="contr.helmert")))[,-1]
f0.rda<-rda(f.com.hel0~.,data=trt0.mat)
RsquareAdj(f0.rda)
(fish.initial<-anova(f0.rda.null,f0.rda,permutations=how(nperm=999)))

# treatment does not explain a significant amount of variance between plots initially

# confirming this with an anosim
(f0.anosim<-anosim(f.com.pa.0,f.env0$treatment,permutations = 999,distance="bray"))
summary(f0.anosim)
# the anosim results confirm this. 

# now start building more and more complex models for experiment data
#add in interactions
trt<-as.factor(f.env$treatment)
seas<-f.env$season
samp<-f.env$sampling
yr<-as.factor(f.env$yr)
plts<-as.factor(f.env$plot)

# start with the simplest model that makes sense - there are two of these - season and year
# with an interaction with treatment. Look at these with and without plots

# first question is an interaction between samp* treatment better than intercept only
# first build intercept only
f.rda.null<-rda(f.com.hel~1)

# now make samp*treat matrix
tr.samp.mat<-data.frame(model.matrix(~ samp*trt + plts, 
                                     contrasts=list(trt="contr.helmert")))[,-1]
f.rda.samp.treat.plot<-rda(f.com.hel~.,data=tr.samp.mat)

# now check to see if its better than intercept
anova(f.rda.null,f.rda.samp.treat.plot)

# it is better than null
# does including plot help
tr.samp.mat.np<-data.frame(model.matrix(~ samp*trt, 
                                        contrasts=list(trt="contr.helmert")))[,-1]
f.rda.samp.treat<-rda(f.com.hel~.,data=tr.samp.mat.np)

# check to see if plot should be included here
anova(f.rda.samp.treat.plot,f.rda.samp.treat)

# no difference as far as anova - look at rsquare
RsquareAdj(f.rda.samp.treat.plot)$adj.r.squared
RsquareAdj(f.rda.samp.treat)$adj.r.squared

# r squared better without plot
# current best model is just samp*treat

# now look at next simplest "base" model - this is treatment*year*season

tr.s.yr.mat<-data.frame(model.matrix(~ yr*seas*trt + plts, 
                                     contrasts=list(trt="contr.helmert", seas="contr.helmert",yr="contr.helmert")))[,-1]
tr.s.yr.mat.nop<-data.frame(model.matrix(~ yr*seas*trt, 
                                         contrasts=list(trt="contr.helmert", seas="contr.helmert",yr="contr.helmert")))[,-1]

# look at whether or not model without plot is better than current best model
f.rda.yr.s.treat<-rda(f.com.hel~.,tr.s.yr.mat.nop)

# models aren't nested so have to look at adjusted R2
RsquareAdj(f.rda.yr.s.treat)$adj.r.squared
RsquareAdj(f.rda.samp.treat)$adj.r.squared

# unlike for inverts samp*treat+plot is still the best model - see if plots make a difference here
f.rda.yr.s.treat.plot<-rda(f.com.hel~.,tr.s.yr.mat)
anova(f.rda.yr.s.treat,f.rda.yr.s.treat.plot)

# plot doesn't improve things here - how about for r2
RsquareAdj(f.rda.yr.s.treat)$adj.r.squared
RsquareAdj(f.rda.yr.s.treat.plot)$adj.r.squared
RsquareAdj(f.rda.samp.treat.plot)$adj.r.squared # this is the best model at the moment


# is the interaction between sampling and treatment important?
tr.samp.mat1<-data.frame(model.matrix(~ samp+trt+plts, 
                                        contrasts=list(trt="contr.helmert")))[,-1]
f.rda.samp.treat.noint<-rda(f.com.hel~.,data=tr.samp.mat1)

anova(f.rda.samp.treat.plot,f.rda.samp.treat.noint)

# the interaction doesn't actually significantly improve the model. Look at r2
RsquareAdj(f.rda.samp.treat.noint)$adj.r.squared
RsquareAdj(f.rda.samp.treat.plot)$adj.r.squared# this is still the best model


# now look into whether or not productivity measures explain community patterns better

f.env.prod<-f.env[,6:8]
alg.e<-f.env$alg
sg<-f.env$sg.sd
sggrow<-f.env$grow


f.rda.prod<-rda(f.com.hel~.,f.env.prod)

# does this model do better than the null
anova(f.rda.null,f.rda.prod)

# yes it does
# now does it do better than the treatment model
RsquareAdj(f.rda.samp.treat.plot)$adj.r.squared
RsquareAdj(f.rda.prod)$adj.r.squared

#no

#What about a model that uses sg grow instead of treatment
sg.s.yr.mat.nop<-data.frame(model.matrix(~ samp*sggrow + plts))[,-1]
f.rda.sgsyr<-rda(f.com.hel~.,sg.s.yr.mat.nop)

# does it do better than the null
anova(f.rda.null,f.rda.sgsyr)

# yes it does - how about compared to model with treatment
RsquareAdj(f.rda.samp.treat.plot)$adj.r.squared
RsquareAdj(f.rda.sgsyr)$adj.r.squared

# on its own productivity doesn't do a better job than just treatment - but its close

# not on its own, no. What if we include different measures of productivity in our treatment model
# from now on I'm referring to the year*season*treatment as f.rda.best
f.rda.best<-f.rda.samp.treat.plot

b.alg.mat<-data.frame(model.matrix(~ samp*trt+plts+alg.e, 
                                   contrasts=list(trt="contr.helmert")))[,-1]
b.sg.mat<-data.frame(model.matrix(~ samp*trt+plts+sg, 
                                  contrasts=list(trt="contr.helmert")))[,-1]
b.sgg.mat<-data.frame(model.matrix(~ samp*trt+plts+sggrow, 
                                   contrasts=list(trt="contr.helmert")))[,-1]
b.alg.sg.mat<-data.frame(model.matrix(~ samp*trt+plts+alg.e+sg, 
                                      contrasts=list(trt="contr.helmert")))[,-1]
b.alg.sgg.mat<-data.frame(model.matrix(~ samp*trt+plts+alg.e+sggrow, 
                                       contrasts=list(trt="contr.helmert")))[,-1]
b.sg.sgg.mat<-data.frame(model.matrix(~ samp*trt+plts+sg+sggrow, 
                                      contrasts=list(trt="contr.helmert")))[,-1]
b.alg.sg.sgg.mat<-data.frame(model.matrix(~ samp*trt+plts+alg.e+sg+sggrow, 
                                          contrasts=list(trt="contr.helmert")))[,-1]

# start with all of them added
f.rda.bprod<-rda(f.com.hel~.,b.alg.sg.sgg.mat)

anova(f.rda.best,f.rda.bprod)

# no difference between models - look at r2
RsquareAdj(f.rda.best)$adj.r.squared
RsquareAdj(f.rda.bprod)$adj.r.squared

# best model still highest r2

# look at just sg growth

f.rda.bgrow<-rda(f.com.hel~.,b.sgg.mat)

anova(f.rda.best,f.rda.bgrow)

# no difference between models - look at r2
RsquareAdj(f.rda.best)$adj.r.squared
RsquareAdj(f.rda.bgrow)$adj.r.squared

# best model still highest r2

# look at sggrow and algae
f.rda.balggrow<-rda(f.com.hel~.,b.alg.sgg.mat)

anova(f.rda.best,f.rda.balggrow)

# no difference between models - look at r2
RsquareAdj(f.rda.best)$adj.r.squared
RsquareAdj(f.rda.balggrow)$adj.r.squared

# best model still highest r2

# look at sggrow and sg
f.rda.bsggrow<-rda(f.com.hel~.,b.sg.sgg.mat)

anova(f.rda.best,f.rda.bsggrow)

# no difference between models - look at r2
RsquareAdj(f.rda.best)$adj.r.squared
RsquareAdj(f.rda.bsggrow)$adj.r.squared
RsquareAdj(f.rda.bprod)$adj.r.squared

# best model still highest r2

f.rda.bsg<-rda(f.com.hel~.,b.sg.mat)

anova(f.rda.best,f.rda.bsg)

# no difference between models - look at r2
RsquareAdj(f.rda.best)$adj.r.squared
RsquareAdj(f.rda.bsg)$adj.r.squared

# best model at the moment: is samp*treatment + plot

#best model at the moment
plot(f.rda.best, scaling = 3, display = c("sp", "cn"))

# note the order matters for adonis2 with by="terms" and the by="margin" doesn't seem to work
f.mod <- adonis2(f.com.hel~., data = tr.samp.mat, method="euclidean", by="terms") 
f.mod


# first going to make sure where there's a difference between treatments using anosim

# look at 1 month into experiment
f.com1<-f.com.pa[f.env$sampling==1,]
f.com1<-f.com1[,colSums(f.com1)!=0]
f.env1<-f.env[f.env$sampling==1,]

(f1.anosim<-anosim(f.com1,f.env1$treatment))
# no difference

# look at 5 months into experiment
f.com5<-f.com.pa[f.env$sampling==5,]
f.com5<-f.com5[,colSums(f.com5)!=0]
f.env5<-f.env[f.env$sampling==5,]

(f5.anosim<-anosim(f.com5,f.env5$treatment))

# Now we're getting more robust differences between treatments

# look at 12 months into experiment
f.com12<-f.com.pa[f.env$sampling==12,]
f.com12<-f.com12[,colSums(f.com12)!=0]
f.env12<-f.env[f.env$sampling==12,]

(f12.anosim<-anosim(f.com12,f.env12$treatment))

# significantly different. R is larger

# look at 17 months into experiment
f.com17<-f.com.pa[f.env$sampling==17,]
f.com17<-f.com17[,colSums(f.com17)!=0]
f.env17<-f.env[f.env$sampling==17,]

(f17.anosim<-anosim(f.com17,f.env17$treatment))

# different but R isn't quite as large

#this corresponds with the results of the RDA that sampling and treatment are important

# treatment stays an important factor from 5 months on. 
# at 5 months which treatments are different from each other?

f.com5.bf<-f.com5[f.env5$treatment!="real",]
f.env5.bf<-f.env5[f.env5$treatment!="real",]

(f5bf.anosim<-anosim(f.com5.bf,f.env5.bf$treatment))

# blank and fake aren't different

f.com5.br<-f.com5[f.env5$treatment!="fake",]
f.env5.br<-f.env5[f.env5$treatment!="fake",]

(f5br.anosim<-anosim(f.com5.br,f.env5.br$treatment))

# blank and real are VERY different

f.com5.fr<-f.com5[f.env5$treatment!="blank",]
f.env5.fr<-f.env5[f.env5$treatment!="blank",]

(f5fr.anosim<-anosim(f.com5.fr,f.env5.fr$treatment))

# real and fake are different


# at 12 months which treatments are different from each other?

f.com12.bf<-f.com12[f.env12$treatment!="real",]
f.env12.bf<-f.env12[f.env12$treatment!="real",]

(f12bf.anosim<-anosim(f.com12.bf,f.env12.bf$treatment))

# blank and fake aren't different

f.com12.br<-f.com12[f.env12$treatment!="fake",]
f.env12.br<-f.env12[f.env12$treatment!="fake",]

(f12br.anosim<-anosim(f.com12.br,f.env12.br$treatment))

# blank and real are VERY different

f.com12.fr<-f.com12[f.env12$treatment!="blank",]
f.env12.fr<-f.env12[f.env12$treatment!="blank",]

(f12fr.anosim<-anosim(f.com12.fr,f.env12.fr$treatment))

# real and fake are VERY different

# at 17 months which treatments are different from each other?

f.com17.bf<-f.com17[f.env17$treatment!="real",]
f.env17.bf<-f.env17[f.env17$treatment!="real",]

(f17bf.anosim<-anosim(f.com17.bf,f.env17.bf$treatment))

# blank and fake are not different

f.com17.br<-f.com17[f.env17$treatment!="fake",]
f.env17.br<-f.env17[f.env17$treatment!="fake",]

(f17br.anosim<-anosim(f.com17.br,f.env17.br$treatment))

# blank and real ARE different

f.com17.fr<-f.com17[f.env17$treatment!="blank",]
f.env17.fr<-f.env17[f.env17$treatment!="blank",]

(f17fr.anosim<-anosim(f.com17.fr,f.env17.fr$treatment))

# real and fake are different

# so fish communities in sponge plots diverge from control and structure control early.
# control and structure control never diverge


# which species are driving this difference?
(f17.simper<-simper(f.com17,f.env17$treatment,permutations = 999))
summary(f17.simper)

# looks like ceriths, oysters, sea cucumbers, anemones, and blue crabs are both significant
# and consistently contribute to differences between the groups here

fcom17diff<-bind_cols(f.com.pa[,colnames(f.com.pa) %in% c("damselfish","grunt")],
                      f.env)%>%
  pivot_longer(1:2,names_to="taxa",values_to="presence")%>%
  group_by(treatment,taxa,sampling)%>%
  summarize(n.plots = sum(presence))%>%
  mutate(treatment=factor(treatment,levels=c("blank","fake","real"),labels = c("Control","Structure Control","Sponge")))

ggplot(data=fcom17diff)+
  geom_col(aes(x=taxa,y=n.plots,fill=treatment),width=.5,
           position=position_dodge(.5))+
  scale_fill_viridis_d(option="A",begin=0,end=0.6,"")+
  # geom_point(aes(x=taxa, y=n.plots,color=treatment),
  #            position=position_dodge(.5),size=5)+
  scale_color_viridis_d(option="A",begin=0,end=0.6,"")+
  facet_wrap(~sampling)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        legend.position = "top")+
  ylab("Number of plots where taxa is present")




# make figures for this

init.scores<-scores(rda(f.com.hel0),scaling=1,1:2)

f.env0p<-bind_cols(f.env0,data.frame(init.scores$sites))%>%
  mutate(treatment=factor(treatment, labels=c("Control","Structure Control","Sponge")))

hull0 <- f.env0p %>%
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
  filter(taxa %in% c("damselfish","grunt"))
taxa<-rownames(sprp0)

# look at how much variation each axis explains to add to axis labels
summary(f0.rda.null)$cont$importance

(f0<-ggplot()+
    # ylim(-1.2,1.2)+
    # xlim(-1.2,1.2)+
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
                 data=f.env0p,
                 alpha = 0.2, 
                 show.legend = FALSE,
                 level = 0.95)+
    #    geom_polygon(data = hull0, aes(x=PC1,y=PC2,color=treatment,fill=treatment),alpha = 0.1)+
    geom_jitter(data=f.env0p,aes(x=PC1,y=PC2,color=treatment),size=2)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          legend.position = c(.14,0.95),
          legend.background = element_blank())+
    coord_fixed()+
    xlab("PC1 70.22%")+
    ylab("PC2 29.73%")+
    scale_color_viridis_d(option="A",begin=0,end=0.6,"")+
    scale_fill_viridis_d(option="A",begin=0,end=0.6,""))

(f0hull<-ggplot()+
    # ylim(-1,1)+
    # xlim(-1,1)+
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
    geom_point(data=f.env0p,aes(x=PC1,y=PC2,color=treatment),size=2)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          legend.position = c(.14,0.95),
          legend.background = element_blank())+
    coord_fixed()+
    xlab("PC1 23.74%")+
    ylab("PC2 16.54%")+
    scale_color_viridis_d(option="A",begin=0,end=0.6,"")+
    scale_fill_viridis_d(option="A",begin=0,end=0.6,""))

# these figures are not informative 



# make figures for month 5
f.com5.hel<-decostand(f.com5,"hellinger")
f5.rda.null<-rda(f.com5.hel~1)
m5.scores<-scores(f5.rda.null,scaling=1,1:2)

f.env5p<-bind_cols(f.env5,data.frame(m5.scores$sites))

# look at how much variation each axis explains to add to axis labels
summary(f5.rda.null)$cont$importance

circ <- circleFun(center=c(0,0),diameter=sqrt(2/12),npoints = 500)

hull5 <- f.env5p %>%
  group_by(treatment)%>%
  slice(chull(PC1, PC2))

spr5<-data.frame(m5.scores$species)
sprp5<-data.frame(spr5[abs(spr5[,1])>=max(circ[,1])|abs(spr5[,2])>=max(circ[,2]),])
sprp5$taxa<-rownames(sprp5)
sprp5<-sprp5%>%
  filter(taxa %in% c("damselfish","grunt"))
taxa<-rownames(sprp5)



(f5<-ggplot()+
    # ylim(-1.2,1.2)+
    # xlim(-1.2,1.2)+
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
                 data=f.env5p,
                 alpha = 0.2, 
                 show.legend = FALSE,
                 level = 0.95)+
    #    geom_polygon(data = hull0, aes(x=PC1,y=PC2,color=treatment,fill=treatment),alpha = 0.1)+
    geom_point(data=f.env5p,aes(x=PC1,y=PC2,color=treatment),size=2)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          legend.position = "none")+
    coord_fixed()+
    xlab("PC1 47.03%")+
    ylab("PC2 28.73%")+
    scale_color_viridis_d(option="A",begin=0,end=0.6,"")+
    scale_fill_viridis_d(option="A",begin=0,end=0.6,""))

(f5hull<-ggplot()+
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
    geom_point(data=f.env5p,aes(x=PC1,y=PC2,color=treatment),size=2)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          legend.position = "none")+
    coord_fixed()+
    xlab("PC1 23.45%")+
    ylab("PC2 19.19%")+
    scale_color_viridis_d(option="A",begin=0,end=0.6,"")+
    scale_fill_viridis_d(option="A",begin=0,end=0.6,""))

# still not hugely informative

# make figures for month 12
f.com12.hel<-decostand(f.com12,"hellinger")
f12.rda.null<-rda(f.com12.hel~1)

m12.scores<-scores(f12.rda.null,scaling=1,1:2)

f.env12p<-bind_cols(f.env12,data.frame(m12.scores$sites))

# look at how much variation each axis explains to add to axis labels
summary(f12.rda.null)$cont$importance

circ <- circleFun(center=c(0,0),diameter=sqrt(2/13),npoints = 500)

hull12 <- f.env12p %>%
  group_by(treatment)%>%
  slice(chull(PC1, PC2))

sprp12<-data.frame(m12.scores$species)
sprp12$taxa<-rownames(sprp12)
sprp12<-sprp12%>%
  filter(taxa %in% c("damselfish","grunt"))
taxa<-rownames(sprp12)


(f12<-ggplot()+
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
                 data=f.env12p,
                 alpha = 0.2, 
                 show.legend = FALSE,
                 level = 0.95)+
    #    geom_polygon(data = hull0, aes(x=PC1,y=PC2,color=treatment,fill=treatment),alpha = 0.1)+
    geom_point(data=f.env12p,aes(x=PC1,y=PC2,color=treatment),size=2)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          legend.position = "none")+
    coord_fixed()+
    xlab("PC1 54.10%")+
    ylab("PC2 16.34%")+
    scale_color_viridis_d(option="A",begin=0,end=0.6,"")+
    scale_fill_viridis_d(option="A",begin=0,end=0.6,""))

(f12hull<-ggplot()+
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
    geom_point(data=f.env12p,aes(x=PC1,y=PC2,color=treatment),size=2)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          legend.position = "none")+
    coord_fixed()+
    xlab("PC1 24.11%")+
    ylab("PC2 18.83%")+
    scale_color_viridis_d(option="A",begin=0,end=0.6,"")+
    scale_fill_viridis_d(option="A",begin=0,end=0.6,""))


# make figures for month 17

f.com17.hel<-decostand(f.com17,"hellinger")
f17.rda.null<-rda(f.com17.hel~1)
m17.scores<-scores(f17.rda.null,scaling=3,1:2)

f.env17p<-bind_cols(f.env17,data.frame(m17.scores$sites))

# look at how much variation each axis explains to add to axis labels
summary(f17.rda.null)$cont$importance

circ <- circleFun(center=c(0,0),diameter=sqrt(2/15),npoints = 500)

hull17 <- f.env17p %>%
  group_by(treatment)%>%
  slice(chull(PC1, PC2))

sprp17<-data.frame(m17.scores$species)
sprp17$taxa<-rownames(sprp17)
sprp17<-sprp17%>%
  filter(taxa %in% c("damselfish","grunt"))
taxa<-rownames(sprp17)



(f17<-ggplot()+
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
                 data=f.env17p,
                 alpha = 0.2, 
                 show.legend = FALSE,
                 level = 0.95)+
    #    geom_polygon(data = hull17, aes(x=PC1,y=PC2,color=treatment,fill=treatment),alpha = 0.1)+
    geom_point(data=f.env17p,aes(x=PC1,y=PC2,color=treatment),size=2)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          legend.position = "none")+
    coord_fixed()+
    xlab("PC1 26.87%")+
    ylab("PC2 20.53%")+
    scale_color_viridis_d(option="A",begin=0,end=0.6,"")+
    scale_fill_viridis_d(option="A",begin=0,end=0.6,""))

(f17hull<-ggplot()+
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
    geom_point(data=f.env17p,aes(x=PC1,y=PC2,color=treatment),size=2)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          legend.position = "none")+
    coord_fixed()+
    xlab("PC1 26.87%")+
    ylab("PC2 20.53%")+
    scale_color_viridis_d(option="A",begin=0,end=0.6,"")+
    scale_fill_viridis_d(option="A",begin=0,end=0.6,""))

# can go back and look at which species are consistently important

#devtools::install_github("thomasp85/patchwork")
library(patchwork)
f0+f5+f12+f17+plot_layout(widths = 1,heights = 1)

f0hull+f5hull+f12hull+f17hull+plot_layout(widths = 1,heights = 1)

# now look at univariate results

f.env$spr<-specnumber(f.com)
f.env$div<-diversity(f.com)
f.env$plot<-as.factor(f.env$plot)
f.env$treatment<-as.factor(f.env$treatment)
f.env$season<-as.factor(f.env$season)

library(lmerTest)

f.env0uni<-f.env0%>%
  mutate(strt.spr=specnumber(f.com0),
         strt.div=diversity(f.com0),
         plot=factor(plot))%>%
  select(treatment,plot,strt.spr,strt.div)

f.env.uni<-f.env%>%
  filter(sampling!=0)%>%
  left_join(f.env0uni)

f.env.uni$treatment<-as.factor(f.env.uni$treatment)

spr.lmer<-lmer(spr~treatment*sampling + sg.sd+grow+(1|plot)+
                 offset(strt.spr),
               data = f.env.uni%>%
                 mutate(treatment=relevel(treatment, ref = "real")))
summary(spr.lmer)
anova(spr.lmer)

# looks like species richness is just always higher than controls in sponge plots


f.sum<-f.env.uni%>%
  group_by(treatment,sampling)%>%
  summarize(spr.m=mean(spr),div.m=mean(div),spr.sd=sd(spr),div.sd=sd(div))

f.sum$treatment<-factor(f.sum$treatment,labels=c("Control","Structure Control","Sponge"))

ggplot()+
  geom_line(data=f.sum,aes(x=sampling,y=spr.m,group=treatment,color=treatment),position=position_dodge(0.5))+
  geom_errorbar(data=f.sum,aes(x=sampling,ymin=spr.m-spr.sd,ymax=spr.m+spr.sd,color=treatment),width=.5,position=position_dodge(0.5))+
  geom_point(data=f.sum,aes(x=sampling,y=spr.m,color=treatment),size=5,position=position_dodge(0.5))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        legend.text = element_text(size=10))+
  scale_color_viridis_d(option="A", begin=0, end=0.6,"")+
  ylab("Species Richness")+
  xlab("Months into Experiment")

ggsave("figures/fish_sp_richness.jpg")

# now look at diversity

div.lmer<-lmer(div~treatment*as.factor(sampling) + sg.sd+grow+(1|plot)+
                 offset(strt.div),
               data = f.env.uni%>%
                 mutate(treatment=relevel(treatment, ref = "real")))
summary(div.lmer)
anova(div.lmer)

ggplot()+
  geom_line(data=f.sum,aes(x=sampling,y=div.m,group=treatment,color=treatment),position=position_dodge(0.5))+
  geom_errorbar(data=f.sum,aes(x=sampling,ymin=div.m-div.sd,ymax=div.m+div.sd,color=treatment),width=.5,position=position_dodge(0.5))+
  geom_point(data=f.sum,aes(x=sampling,y=div.m,color=treatment),size=5,position=position_dodge(0.5))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        legend.text = element_text(size=10))+
  scale_color_viridis_d(option="A", begin=0, end=0.6,"")+
  ylab("Diversity")+
  xlab("Months into Experiment")

ggsave("figures/fish_sp_diversity.jpg")

# the diversity story isn't as clear cut but generally the same as richness




