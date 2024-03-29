# How does the sponge Ircinia felix influence seagrass bed primary producers? #
# This script will examining nutrient data #
# packages----
if(!require(tidyverse))install.packages("tidyverse");library(tidyverse)
if(!require(car))install.packages("car");library(car)
if(!require(lmerTest))install.packages("lmerTest");library(lmerTest)
if(!require(DHARMa))install.packages("DHARMa");library(DHARMa)
if(!require(glmmTMB))install.packages("glmmTMB");library(glmmTMB)
source("MarineBiology_DOI_10.1007/00_results_functions.R")
theme_set(theme_bw())
# load data----
source("MarineBiology_DOI_10.1007/03_reimport.R")#imports all the data sets
mn<-readxl::read_xlsx("Original_data/missing_nutrients.xlsx")
sg_nuts<-bind_rows(sg_nuts,mn)%>%
  select(-mg)%>%
  group_by(treatment,plot,dist,sampling)%>%
  summarize(PN=sum(PN,na.rm = TRUE),PC=sum(PC,na.rm = TRUE),PP=sum(PP,na.rm = TRUE))

#make seagrass nutriet data frame with nvalues filled in
sgn<-sg_nuts %>%
  filter(PC >0 & PN >0 & PP >0)%>%
  mutate(cn=PC/PN,
         cp=PC/PP,
         np=PN/PP)%>%
  pivot_longer(PN:np,names_to = "nut",values_to = "nvalue")%>%
  filter(!is.na(nvalue))

# I'm going to look at the distributions of the data
ggplot(data=sgn,aes(nvalue))+
  geom_histogram()+
  facet_wrap(~nut,scales="free")
# there are some outliers in PC and PN- look at what these are
filter(sgn,nut=="PC" & nvalue>55)
filter(sgn,nut=="PN" & nvalue>2.5)
# carbon and nitrogen are calculated from the same sample and plot 6 (sampling 2 dist 0), 
# plot 7 sampling 2 dist 0, and plot 8 sampling 2 dist 1
# are outliers for both- makes me think something besides seagrass got in there.
# I'm going to remove these points before going further. 

sgn_no<-bind_rows(sgn %>%
                    filter(plot==6 & sampling==2 & dist==0),
                  sgn %>%
                    filter(plot==7 & sampling==2 & dist==0),
                  sgn %>%
                    filter(plot==8 & sampling==2 & dist==1))

sgn<-sgn %>%
  anti_join(sgn_no)

# now I'm going to look at the change in nutrients from initial conditions

sv<-filter(sgn, sampling == 1) %>% 
  rename(snv = nvalue)%>%select(treatment,plot,dist,snv,nut)

dsgn<-left_join(sgn,sv)%>%
  mutate(delta=nvalue-snv,
    season=case_when(
      sampling==1~"summer",
      sampling==2~"summer",
      sampling==3~"winter",
      sampling==4~"summer",
      sampling==5~"winter"),
    yr=as.factor(case_when(
      sampling==1~0,
      sampling==2~1,
      sampling==3~1,
      sampling==4~2,
      sampling==5~2)),
    mnths=case_when(
      sampling==1~0,
      sampling==2~1,
      sampling==3~5,
      sampling==4~12,
      sampling==5~17))#%>%filter(sampling!=1)

# look at sample sizes here
ggplot(dsgn%>%
  group_by(treatment,sampling,dist,nut)%>%
  filter(nut %in% c("PN","PP","PC"))%>%
  summarize(ns=n()))+
  geom_bar(aes(x=treatment,y=ns,fill=nut),position=position_dodge(),stat="identity")+
  facet_grid(sampling~dist)


# check for season effect
Nseason<-lmer(nvalue~season+(1|plot), data = dsgn%>%
    filter(nut=="PN")) 
(Nsum<-summary(Nseason))

efsize(Nsum,2)

(Nseason<-Anova(Nseason, type = "III"))

Pseason<-lmer(nvalue~season+(1|plot), data = dsgn%>%
    filter(nut=="PP")) 

(Psum<-summary(Pseason))

efsize(Psum,2)

(Pseason<-Anova(Pseason, type = "III"))

Cseason<-lmer(nvalue~season+(1|plot), data = dsgn%>%
    filter(nut=="PC")) 

(Csum<-summary(Cseason))

efsize(Csum,2)

(Cseason<-Anova(Cseason, type = "III"))
# YES for all

# first make treatment a factor to allow releveling 
sgn$treatment<-as.factor(sgn$treatment)
sv$treatment<-as.factor(sv$treatment)
dsgn$treatment<-as.factor(dsgn$treatment)
# check starting conditions at the distances for which there is balanced data
n.lmerb<-lmer(snv~treatment+(1|plot),
  data=sv%>%
    filter(nut=="PN" & dist<=0.5)%>%
    mutate(treatment=relevel(treatment,ref="real")))
# (nstart<-Anova(n.lmerb, type = "III"))
(nstart<-summary(n.lmerb))

p.lmerb<-lmer(snv~treatment+(1|plot),
  data=sv%>%
    filter(nut=="PP" & dist<=0.5)%>%
    mutate(treatment=relevel(treatment,ref="real")))
# (pstart<-Anova(p.lmerb, type = "III"))
(pstart<-summary(p.lmerb))

c.lmerb<-lmer(snv~treatment+(1|plot),
  data=sv%>%
    filter(nut=="PC" & dist<=0.5)%>%
    mutate(treatment=relevel(treatment,ref="real")))
# (cstart<-Anova(c.lmerb, type = "III"))
(cstart<-summary(c.lmerb))

# fake never differs from real at start, only % carbon in blank plots differs sig from real

# make a function to look at residuals
glmm.resids<-function(model){
  t1 <- simulateResiduals(model)
  print(testDispersion(t1))
  plot(t1)
}


dsgn%>% filter(dist<=0 & nut %in% c("PC","PN","PP"))%>%
  # group_by(treatment,dist,sampling,nut,season,yr,mnths)%>%
  # summarize(mv=mean(delta),sd=sd(delta)) %>%
  # ggplot(aes(y=mv,color=treatment,x=as.factor(mnths)))+
  # geom_point(position=position_dodge(width=.3),size=2)+
  # geom_errorbar(aes(ymin=mv-sd,ymax=mv+sd,color=treatment),width=.1,position=position_dodge(width=0.3))+
  ggplot(aes(
    # y=delta,
    y=nvalue,
    color=treatment,x=as.factor(mnths)))+
  geom_point(position=position_dodge(width=.3),size=1.5,alpha=0.7)+
  xlab("Months into experiment") +
  scale_color_viridis_d(
    option="A",
    begin=0, end=0.6,
    name="",
    labels=c("Control","Structure","Sponge"))+
  facet_grid(rows = vars(nut),cols = vars(season),scales="free") +
  ggtitle("At center of plot")

dsgn %>% filter(dist<=0.5 &nut %in% c("PC","PN","PP"))%>%
    ggplot(aes(y=nvalue,color=treatment,x=as.factor(mnths)))+
    geom_point(position=position_dodge(width=.3),size=1.5,alpha=0.7)+
    xlab("Months into experiment") +
  scale_color_viridis_d(
    option="A",
    begin=0, end=0.6,
    name="",
    labels=c("Control","Structure","Sponge"))+
  facet_grid(rows = vars(nut),cols = vars(season),scales="free") +
  ggtitle("Up to 0.5m from center of plot")

# ggsave("figures/allnutrient0.5m.png")

dsgn %>% filter(dist<=1 &nut %in% c("PC","PN","PP"))%>%
  ggplot(aes(y=nvalue,color=treatment,x=as.factor(mnths)))+
  geom_point(position=position_dodge(width=.3),size=1.5,alpha=0.7)+
  xlab("Months into experiment") +
  scale_color_viridis_d(
    option="A",
    begin=0, end=0.6,
    name="",
    labels=c("Control","Structure","Sponge"))+
  facet_grid(rows = vars(nut),cols = vars(season),scales="free") +
  ggtitle("Up to 1m from center of plot")

dsgn %>% filter(dist<=2 & nut %in% c("PC","PN","PP"))%>%
  ggplot(aes(y=nvalue,color=treatment,x=as.factor(mnths)))+
  geom_point(position=position_dodge(width=.3),size=1.5,alpha=0.7)+
  xlab("Months into experiment") +
  scale_color_viridis_d(
    option="A",
    begin=0, end=0.6,
    name="",
    labels=c("Control","Structure","Sponge"))+
  facet_grid(rows = vars(nut),cols = vars(season),scales="free") +
  ggtitle("Entire plots")

# ggsave("figures/allnutrients.png")


# because we don't have a before sample during winter, we focus on summer
# due to unbalanced sample number in some distances and sample times
# look at change after 1 year and include only nearest two distances 

#in %N
# relevel with real as the base
# one year later seagrass in plots with sponges have higher %N than either treatment

# # summer offset version
# n.lmerb.o<-glmmTMB(nvalue~
#     treatment*dist+
#     # treatment*as.factor(dist)+ #use as.factor if including more than 2 levels due to likely non-linearity
#     (1|plot),
#   offset=snv,
#   data=dsgn%>%
#     filter(nut=="PN" &
#         dist<=0.5 & # both these distances have complete samples
#         sampling %in% c(4))) # second summer
# (nbs<-summary(n.lmerb.o))
# glmm.resids(n.lmerb.o) # ugly residuals


# summer BACI version
n.lmerb<-lmer(nvalue~
    treatment*yr*dist+ 
    (1|plot),
  data=dsgn%>%
    filter(nut=="PN" & 
        dist<=0.5 & # both these distances have complete samples
        sampling %in% c(1,4))%>%# second summer
    mutate(treatment=relevel(treatment,ref="blank")))
(nbs<-summary(n.lmerb))
glmm.resids(n.lmerb) # good!

# # summer glmmTMB BACI version similar but less concervative so stick with lmer
# n.lmerb<-glmmTMB(nvalue~
#     treatment*yr*dist+ 
#     (1|plot),
#   data=dsgn%>%
#     filter(nut=="PN" & 
#         dist<=0.5 & # both these distances have complete samples
#         sampling %in% c(1,4)) %>%
#   mutate(treatment=relevel(treatment,ref="blank")))
# (nbs<-summary(n.lmerb))
# glmm.resids(n.lmerb)


n.lmerr<-lmer(nvalue~
    treatment*yr*dist+
    (1|plot),
  data=dsgn%>%
    filter(nut=="PN" & 
        dist<=0.5 & # both these distances have complete samples
        sampling %in% c(1,4))%>% # second summer
    mutate(
      treatment=relevel(treatment,ref="real")))
(nrs<-summary(n.lmerr))

n.lmerf<-lmer(nvalue~    
    treatment*yr*dist+ 
    (1|plot),
  data=dsgn%>%
    filter(nut=="PN" & 
        dist<=0.5 & # both these distances have complete samples
        sampling %in% c(1, 4))%>% # second summer
    mutate(treatment=relevel(treatment,ref="fake")))
(nfs<-summary(n.lmerf))
# fake starts higher and decreases non-significantly in nearest distance, and increases at 0.5 m.

# now doing percent phosphorus

# summer BACI version
p.lmerb<-lmer(nvalue~
    treatment*yr*dist+ 
    # treatment*as.factor(dist)+ #use as.factor if including more than 2 levels due to likely non-linearity
    (1|plot),
  # offset=snv,
  data=dsgn%>%
    filter(nut=="PP" & 
        dist<=0.5 & # both these distances have complete samples
        sampling %in% c(1,4))) # second summer
(pbs<-summary(p.lmerb))
glmm.resids(p.lmerb)

p.lmerr<-lmer(nvalue~
    treatment*yr*dist+ 
    (1|plot),
  # offset=snv,
  data=dsgn%>%
    filter(nut=="PP" & 
        dist<=0.5 & # both these distances have complete samples
        sampling %in% c(1, 4))%>% # second summer
    mutate(
      treatment=relevel(treatment,ref="real")))
(prs<-summary(p.lmerr))

p.lmerf<-lmer(nvalue ~    
    treatment*yr*dist+ 
    (1|plot),
  # offset=snv,
  data=dsgn%>%
    filter(nut=="PP" & 
        dist<=0.5 & # both these distances have complete samples
        sampling %in% c(1, 4))%>% # second summer
    mutate(treatment=relevel(treatment,ref="fake")))
(pfs<-summary(p.lmerf))



# look at %C now

# summer BACI version
c.lmerb<-lmer(nvalue~
    treatment*yr*dist+
    # treatment*as.factor(dist)+ #use as.factor if including more than 2 levels due to likely non-linearity
    (1|plot),
  # offset=snv,
  data=dsgn%>%
    filter(nut=="PC" & 
        dist<=0.5 & # both these distances have complete samples
        sampling %in% c(1,4))) # second summer
(cbs<-summary(c.lmerb))
glmm.resids(c.lmerb) #nquantile deviations detected

c.lmerr<-lmer(nvalue~
    treatment*yr*dist+ 
    (1|plot),
  # offset=snv,
  data=dsgn%>%
    filter(nut=="PC" & 
        dist<=0.5 & # both these distances have complete samples
        sampling %in% c(1,4))%>% # second summer
    mutate(
      treatment=relevel(treatment,ref="real")))
(crs<-summary(c.lmerr))

c.lmerf<-lmer(nvalue~    
    treatment*yr*dist+ 
    (1|plot),
  # offset=snv,
  data=dsgn%>%
    filter(nut=="PC" & 
        dist<=0.5 & # both these distances have complete samples
        sampling %in% c(1,4))%>% # second summer
    mutate(treatment=relevel(treatment,ref="fake")))
(cfs<-summary(c.lmerf))
glmm.resids(c.lmerf)

# sponge plots were significantly lower in %C than blank at the beginning 
# sponge plots remain stable (no sig change) in both time and at both distances, but differs from blank an fake by yr 2
# decrease in %C over the year in blank, but its not significant
# There is a significant difference between real and fake after a year
# %C decreased significantly in fake plots after 1 year.





# check distributions
sgn1 <- sgn %>% filter(nut=="PC" & dist==0 & sampling %in% c(1,4))
hist(sgn1$nvalue, breaks=30)
sgn2 <- sgn %>% filter(nut=="PP" & dist==0 & sampling %in% c(1,4))
hist(sgn2$nvalue, breaks=30)
sgn3 <- sgn %>% filter(nut=="PN" & dist==0 & sampling %in% c(1,4))
hist(sgn3$nvalue, breaks=30)



# Now looking at Nutrients after 5 months - 

n.lmerr5<-lmer(nvalue~treatment*season*dist+(1|plot),
              data=dsgn%>%
                filter(nut=="PN" & dist %in% c(0,0.5) & sampling %in% c(1,3))%>%
                mutate(treatment=relevel(treatment,ref="real")))
(nrs5<-summary(n.lmerr5))

glmm.resids(n.lmerr5)

n.lmerb5<-lmer(nvalue~treatment*season*dist+(1|plot),
               data=dsgn%>%
                 filter(nut=="PN" & dist %in% c(0,0.5) & sampling %in% c(1,3))%>%
                 mutate(treatment=relevel(treatment,ref="blank")))
(nbs5<-summary(n.lmerb5))

n.lmerf5<-lmer(nvalue~treatment*season*dist+(1|plot),
               data=dsgn%>%
                 filter(nut=="PN" & dist %in% c(0,0.5) & sampling %in% c(1,3))%>%
                 mutate(treatment=relevel(treatment,ref="fake")))
(nfs5<-summary(n.lmerf5))


# now doing percent phosphorus
p.lmerr5<-lmer(nvalue~treatment*season*dist+(1|plot),
              data=dsgn%>%
                filter(nut=="PP" & dist<=0.5 & sampling %in% c(1,3))%>%
                mutate(treatment=relevel(treatment,ref="real")))
(prs5<-summary(p.lmerr5))

glmm.resids(p.lmerr5)

p.lmerb5<-lmer(nvalue~treatment*season*dist+(1|plot),
               data=dsgn%>%
                 filter(nut=="PP" & dist<=0.5 & sampling %in% c(1,3))%>%
                 mutate(treatment=relevel(treatment,ref="blank")))
(pbs5<-summary(p.lmerb5))

p.lmerf5<-lmer(nvalue~treatment*season*dist+(1|plot),
               data=dsgn%>%
                 filter(nut=="PP" & dist<=0.5 & sampling %in% c(1,3))%>%
                 mutate(treatment=relevel(treatment,ref="fake")))
(pfs5<-summary(p.lmerf5))

# look at %C now

c.lmerr5<-lmer(nvalue~treatment*season*dist+(1|plot),
              data=dsgn%>%
                filter(nut=="PC" & dist<=0.5 & sampling %in% c(1,3))%>%
                mutate(treatment=relevel(treatment,ref="real")))

glmm.resids(c.lmerr5)

(crls5<-summary(c.lmerr5))

c.lmerb5<-lmer(nvalue~treatment*season*dist+(1|plot),
               data=dsgn%>%
                 filter(nut=="PC" & dist<=0.5 & sampling %in% c(1,3))%>%
                 mutate(treatment=relevel(treatment,ref="blank")))

(cbls5<-summary(c.lmerb5))

c.lmerf5<-lmer(nvalue~treatment*season*dist+(1|plot),
               data=dsgn%>%
                 filter(nut=="PC" & dist<=0.5 & sampling %in% c(1,3))%>%
                 mutate(treatment=relevel(treatment,ref="fake")))

(cfls5<-summary(c.lmerf5))



