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

# get data organized

a3.1<-a3%>%
  filter(sampling==1)%>%
  rename(start.abund=abundance)%>%
  select(treatment, plot,taxa,start.abund)

a5<-a3 %>%
  filter(sampling!=1)%>%
  left_join(a3.1)

a5$treatment<-factor(a5$treatment)


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

a5<-a5%>%
  mutate(
    delta.abund = abundance-start.abund,
    delta.abund2 = exp(log(abundance+1)-log(start.abund+1)),
    year=case_when(
    sampling==2~1,
    sampling==3~1,
    sampling==4~2,
    sampling==5~2),
    start = log(start.abund+1))


# if we need sampling in months
a5<-a5%>%
  mutate(samp2=case_when(
    sampling==2~1,
    sampling==3~5,
    sampling==4~12,
    sampling==5~17))


# trying the one global model first
alm1<-glmmTMB(abundance~treatment*year+treatment*season + 
    offset(start) + 
    (1|plot)+(sampling|taxa),
    data = a5 %>% mutate(treatment=relevel(treatment, ref = "real")),
    REML=F,
    family=nbinom2)

# look at residuals
alm1_simres <- simulateResiduals(alm1)
testDispersion(alm1_simres)
plot(alm1_simres)

# residuals look meh
summary(alm1)

# BUT we now have a significant increase over time in real treatment that is 
# significantly different than the other treatments.

library(ggeffects)
# library(sjstats)
# ggpredict(alm1, "year")

p1 <- ggpredict(alm1, terms = c("year", "season","treatment" )) %>% 
  rename(year = x, season = group, treatment = facet)%>% 
  mutate(
    sampling = case_when(
      year==1&season=="summer"~2,
      year==1&season=="winter"~3,
      year==2&season=="summer"~4,
      year==2&season=="winter"~5
    )
  )

p1$treatment<-factor(p1$treatment,levels=c("blank","fake","real"),labels=c("Control","Structure Control","Sponge"))

p2 <- ggpredict(alm1, terms = c("sampling", "taxa", "treatment"), type = "random")
plot(p2)

# making a more appropriate figure
ylab<-expression(paste(Delta," Macroalgae Abundance"))

a6<-a5%>%
  mutate(delta.abund=abundance-start.abund)

a6$treatment<-factor(a6$treatment,levels=c("blank","fake","real"),labels=c("Control","Structure Control","Sponge"))
my.labels <- paste0(c("S","W","S","W"), "\n", c("Year 1","Year 1","Year 2", "Year 2"))

a7<-a6%>%
  filter(season!="winter")
p3<-ggpredict(alm1, terms = c("year", "treatment" ))%>% 
  rename(year = x,treatment = group)
p3$treatment<-factor(p3$treatment,levels=c("blank","fake","real"),labels=c("Control","Structure Control","Sponge"))



ggplot(data=a6)+
  geom_hline(aes(yintercept=0), linetype = "dashed", colour = "darkgrey")+
  scale_x_continuous(name="",breaks=c(2,3,4,5),label = my.labels)+
  geom_ribbon(data = p3, aes(sampling,
    ymin = conf.low,
    ymax = conf.high, 
    group = season),
    alpha = 0.3
    ) + 
  geom_line(data = p1, aes(sampling, predicted, group = season)) +
  geom_point(aes(x=sampling,y=delta.abund,color=taxa), alpha =0.5,
    position=position_dodge(0.5))+
  facet_wrap(~treatment)+
  theme_bw()+
  theme(panel.grid = element_blank(),
    legend.title.align = 0.5,
    strip.background = element_blank(),
    strip.text = element_text(size=14))+
  ylab("Macroalgal Abundance")+
  scale_color_brewer(palette = "Set2",name="Taxa",labels=c("Acetabularia",
    "Halimeda",
    "Laurencia",
    "Penicillus",
    "Udotea"))

# only summer

ggplot(data=a7)+
  geom_hline(aes(yintercept=0), linetype = "dashed", colour = "darkgrey")+
#  scale_x_continuous(name="",breaks=c(2,3,4,5),label = my.labels)+
  geom_ribbon(data = p3, aes(year,
                             ymin = conf.low,
                             ymax = conf.high, 
                             group = treatment),
              alpha = 0.3
  ) + 
  geom_line(data = p3, aes(year, predicted, group = treatment)) +
  geom_point(aes(x=year,y=delta.abund,color=taxa), alpha =0.5,
             position=position_dodge(0.5))+
  facet_wrap(~treatment)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.title.align = 0.5,
        strip.background = element_blank(),
        strip.text = element_text(size=14))+
  ylab("Macroalgal Abundance")+
  scale_color_brewer(palette = "Set2",name="Taxa",labels=c("Acetabularia",
                                                           "Halimeda",
                                                           "Laurencia",
                                                           "Penicillus",
                                                           "Udotea"))

# try summer with means
a8<-a7 %>%
  group_by(year,treatment,taxa)%>%
  summarize(m.abund=mean(delta.abund),se.abund=sd(delta.abund)/sqrt(5))


ggplot(data=a8)+
  geom_hline(aes(yintercept=0), linetype = "dashed", colour = "darkgrey")+
  scale_x_continuous(name="Year of Experiment",breaks=c(1,2),label = c(1,2))+
  geom_ribbon(data = p3, aes(year,
                             ymin = conf.low,
                             ymax = conf.high, 
                             group = treatment),
              alpha = 0.3
  ) + 
  geom_line(data = p3, aes(year, predicted, group = treatment)) +
  geom_errorbar(aes(x=year,ymin=m.abund-se.abund,ymax=m.abund+se.abund,color=taxa),
                position=position_dodge(0.5),width=.3)+
  geom_point(aes(x=year,y=m.abund,color=taxa), alpha =0.95,
             position=position_dodge(0.5),size=4)+
  facet_wrap(~treatment)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.title.align = 0.5,
        strip.background = element_blank(),
        strip.text = element_text(size=14))+
  ylab(ylab)+
  # scale_color_brewer(palette = "Set2",name="Taxa",labels=c("Acetabularia",
  #                                                          "Halimeda",
  #                                                          "Laurencia",
  #                                                          "Penicillus",
  #                                                          "Udotea"))
  scale_color_viridis_d(name="Taxa",
                        end=.9,
                        labels=c("Acetabularia",
                                             "Halimeda",
                                             "Laurencia",
                                             "Penicillus",
                                             "Udotea"))


ggsave("figures/algae_summer_means_viridis.jpg")

# make figure for supplemental
a6$season<-factor(a6$season,levels=c("summer","winter"),labels=c("Summer","Winter"))
p1$season<-factor(p1$season,levels=c("summer","winter"),labels=c("Summer","Winter"))

ggplot()+
  geom_hline(yintercept=0, linetype = "dashed", colour = "darkgrey")+
  geom_point(aes(x=year,y=delta.abund,color=taxa),
             data=a6,
             position = position_dodge(0.5),
             size=3,
             alpha=.5)+
  geom_line(aes(x=year,y=predicted,group=season),data=p1)+
  geom_ribbon(aes(year,
                  ymin = conf.low,
                  ymax = conf.high, 
                  group = season),
              alpha = 0.3,data=p1)+
  facet_grid(season~treatment,scales = "free")+
  ylab(ylab)+
  scale_x_continuous(breaks=c(1,2),name = "Year of Experiment")+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size=14),
    legend.title.align = 0.5)+
  scale_color_viridis_d(name="Taxa",
                        end=.9,
                        labels=c("Acetabularia",
                                 "Halimeda",
                                 "Laurencia",
                                 "Penicillus",
                                 "Udotea"))


ggsave("figures/algae_for_supplemental.jpg")  

# old global figure
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
    sampling==5~17),
    year=case_when(
      sampling==2~1,
      sampling==3~1,
      sampling==4~2,
      sampling==5~2))

a6$treatment<-factor(a6$treatment)
a6$season<-factor(a6$season)

a6<-a6%>%
  mutate(delta.tot=total-start.total)
a7<-a6%>%
  group_by(treatment,samp2)%>%
  summarize(m.tot=mean(delta.tot),sd.tot=sd(delta.tot))


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
