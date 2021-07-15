library(tidyverse)

source("scripts/03_reimport.R")#imports all the data sets
if(!require(lmerTest))install.packages("lmerTest");library(lmerTest)
if(!require(DHARMa))install.packages("DHARMa");library(DHARMa)
if(!require(glmmTMB))install.packages("glmmTMB");library(glmmTMB)
if(!require(ggeffects))install.packages("ggeffects");library(ggeffects)


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

# looks like there's something here.

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
  mutate(year=case_when(
    sampling==2~1,
    sampling==3~1,
    sampling==4~2,
    sampling==5~2),
    logstart= log(start.abund+1))

# check for season effect
Aseason<-glmmTMB(abundance~season +
    (1|plot) + (1|taxa), data=a5, REML=F,
  family=nbinom2) 
(Aseason<-glmmTMB:::Anova.glmmTMB(Aseason, type = "III"))
# YES

# trying the one global model first
alm1.r<-glmmTMB(abundance~treatment*year+treatment*season + 
    offset(logstart) +
    (1|plot) +
    (year:season|taxa),
    # (sampling|taxa),
    data=a5%>%mutate(treatment=relevel(treatment, ref = "real")),
    REML=F,
    family=nbinom2)

alm1.f<-glmmTMB(abundance~treatment*year+treatment*season +
    offset(logstart) +
    (1|plot) +
    (year:season|taxa),
    # (sampling|taxa),
    data=a5%>%mutate(treatment=relevel(treatment, ref = "fake")),
    REML=F,
    family=nbinom2)

alm1.b<-glmmTMB(abundance~treatment*year+treatment*season +       
    offset(logstart) +
    (1|plot) +
    (year:season|taxa),
    # (sampling|taxa),
    data=a5%>%mutate(treatment=relevel(treatment, ref = "blank")),
    REML=F,
    family=nbinom2)
                  

# look at residuals
alm1_simres <- simulateResiduals(alm1.r)
testDispersion(alm1_simres)
plot(alm1_simres)

# residuals look meh
(alm1r.sum<-summary(alm1.r))
(alm1f.sum<-summary(alm1.f))
(alm1b.sum<-summary(alm1.b))


p1 <- ggpredict(alm1.r, terms = c("year", "season", "treatment" )) %>% 
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

p2 <- ggpredict(alm1.r, terms = c("year", "treatment", "taxa", "season" ), type = "random")
plot(p2)


# BUT we now have a significant increase over time in real treatment that is 
# significantly different than the other treatments.

a5<-a5%>%
  mutate(samp2=case_when(
    sampling==2~1,
    sampling==3~5,
    sampling==4~12,
    sampling==5~17))

ylab<-expression(paste(Delta," macroalgae abundance"))

a5$treatment<-factor(a5$treatment,levels=c("blank","fake","real"),labels=c("Control","Structure Control","Sponge"))
my.labels <- paste0(c("S","W","S","W"), "\n", c("Year 1","Year 1","Year 2", "Year 2"))

a6<-a5%>%
  mutate(delta.abund=abundance-start.abund)


p3<-ggpredict(alm1.r, terms = c("year", "treatment", "season"))%>% 
  rename(year = x,treatment = group, season = facet)%>% 
  mutate(
    sampling = case_when(
      year==1&season=="summer"~2,
      year==1&season=="winter"~3,
      year==2&season=="summer"~4,
      year==2&season=="winter"~5
    )
  )
p3$treatment<-factor(p3$treatment,levels=c("blank","fake","real"),labels=c("Control","Structure Control","Sponge"))

# summer only
a7<-a6%>%
  filter(season!="winter")

# with means
a8<-a7 %>%
  group_by(year,treatment,taxa)%>%
  summarize(m.abund=mean(delta.abund),se.abund=sd(delta.abund)/sqrt(5),sd.abund=sd(delta.abund))


# ggplot(data=a8)+
#   geom_hline(aes(yintercept=0), linetype = "dashed", colour = "darkgrey")+
#   scale_x_continuous(name="Year of Experiment",breaks=c(1,2),label = c(1,2))+
#   # global effects from model
#   geom_ribbon(data = filter(p3, season=="summer"), aes(year,
#     ymin = conf.low, ymax = conf.high, 
#     group = treatment, fill = treatment),
#     alpha = 0.2
#   ) + 
#   geom_line(data = filter(p3, season=="summer"),  
#     aes(year, predicted, 
#       colour = treatment,
#       group = treatment), 
#     lwd = 2) +
#   scale_fill_viridis_d(name=" ", option="A", begin=0, end=0.6) + 
#   scale_colour_viridis_d(name=" ", option="A", begin=0, end=0.6) + 
#   # species means and variability
#   ggnewscale::new_scale_color() +
#   geom_errorbar(aes(x=year,
#     # ymin=m.abund-sd.abund, ymax=m.abund+sd.abund,
#     # ymin=m.abund-se.abund, ymax=m.abund+se.abund,
#     ymin=m.abund-se.abund*1.96, ymax=m.abund+se.abund*1.96, # 95% CI
#     color=taxa,
#     group = taxa
#     ), position=position_dodge(0.5), lwd=0.5, alpha = 0.3)+
#   geom_point(aes(x=year,y=m.abund,
#     color=taxa,
#     # shape=taxa,
#     group = taxa
#     ), position=position_dodge(0.5), size=2, alpha = 0.8) +
#   scale_color_viridis_d(name="Taxa",
#     end=.9,
#     labels=c("Acetabularia",
#       "Halimeda",
#       "Laurencia",
#       "Penicillus",
#       "Udotea")) +
#   ylab(ylab)+
#   coord_cartesian(ylim = c(-2, 23)) +
#   facet_wrap(~treatment)+
#   theme_bw()+
#   theme(panel.grid = element_blank(),
#     legend.title.align = 0.5,
#     strip.background = element_blank(),
#     strip.text = element_text(size=14))
# 
# ggsave("figures/algae_summer_means_viridis_95CI.jpg")
## ggsave("figures/algae_summer_means_viridis_95CI_bw.jpg")
