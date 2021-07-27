# All figures for productivity

# figure 1 - algae abundance
source("scripts/04_algae_analysis.R")
a8$taxa<-factor(a8$taxa, labels=c("Acetabularia",
                                    "Halimeda",
                                    "Laurencia",
                                    "Penicillus",
                                    "Udotea"))
(algaeplot<- ggplot(data=a8)+
  geom_hline(aes(yintercept=0), linetype = "dashed", colour = "darkgrey")+
  scale_x_continuous(name="Year of Experiment",breaks=c(1,2),label = c(1,2))+
  # global effects from model
  geom_ribbon(data = filter(p3, season=="summer"), aes(year,
    ymin = conf.low, ymax = conf.high,
    group = treatment, fill = treatment),
    alpha = 0.2
  ) +
  geom_line(data = filter(p3, season=="summer"),
    aes(year, predicted,
      colour = treatment,
      group = treatment),
    lwd = 2) +
  scale_fill_viridis_d(name=" ", option="A",
    begin=0, end=0.6, guide =F
    ) +
  scale_colour_viridis_d(name=" ", option="A",
    begin=0, end=0.6, guide =F
    ) +
  # species means and variability
  # ggnewscale::new_scale_color() +
  geom_errorbar(aes(x=year,
    # ymin=m.abund-sd.abund, ymax=m.abund+sd.abund,
    # ymin=m.abund-se.abund, ymax=m.abund+se.abund,
    ymin=m.abund-se.abund*1.96, ymax=m.abund+se.abund*1.96, # 95% CI
    # color=taxa,
    group = taxa
    ), position=position_dodge(0.5), lwd=0.5, alpha = 0.3)+
  geom_point(aes(x=year,y=m.abund,
    # color=taxa,
    shape=taxa,
    group = taxa
  ), position=position_dodge(0.5), size=2, alpha = 0.8) +
  #scale_shape_discrete(name="Taxa") +
    scale_shape_manual(values=c(15,5,17,16,6))+
  # scale_color_viridis_d(name="Taxa",
  # end=.9,
  # labels=c("Acetabularia",
  #   "Halimeda",
  #   "Laurencia",
  #   "Penicillus",
  #   "Udotea")) +
  ylab(ylab)+
  # coord_cartesian(ylim = c(-2, 23)) +
  facet_wrap(~treatment)+
  theme_bw()+
  theme(panel.grid = element_blank(),
    legend.title.align = 0.5,
    strip.background = element_blank(),
    strip.text = element_text(size=14)))
(algaeplot<-egg::tag_facet(algaeplot, hjust = 0, y = Inf,
  tag_pool = c('(a)  Control','(b)  Structure', "(c)  Sponge"),
  open = " ", close = " ") )
# # ggsave("figures/algae_summer_means_viridis_SE.jpg", plot = algaeplot, width = 7,height=4)
# # ggsave("figures/algae_summer_means_viridis_95CI.jpg", plot = algaeplot, width = 7,height=4)
#ggsave("figures/algae_summer_95CIbw.jpg", plot = algaeplot, width = 7,height=3.5)
ggsave("figures/algae_summer_95CIbw_shapes.jpg", plot = algaeplot, width = 7,height=3.5)

# figure 2 - seagrass shoot density
source("scripts/07_shoots_analysis.R")
sdplot<-sgsd%>%
  mutate(thalassia=T.SD-mtsd,sh=SH.SD-mshsd)%>%
  select(treatment,samp2,thalassia,sh)%>%
  pivot_longer(thalassia:sh,names_to="seagrass",values_to="shoots")

sdplot2<-sdplot%>%
  mutate(samp=case_when(
    treatment=="blank"~samp2-1,
    treatment=="fake"~samp2,
    treatment=="real"~samp2+1),
    seagrass=factor(seagrass,levels=c("thalassia","sh"))
    )

sdsum<-sdplot2%>%
  group_by(treatment,samp,seagrass)%>%
  summarize(msd=mean(shoots),
    sdsd=sd(shoots),
    sdse=sd(shoots)/sqrt(12*5)) # could be just 5?

sdlab<-expression(paste(Delta," seagrass shoots per m"^2))

sdp<-ggplot()+
  geom_hline(yintercept = 0, linetype = "dashed", colour = "darkgrey") +
  ## dot plot versions
  # # geom_jitter(data=sdplot2,aes(x=samp,shoots,color=treatment),width=.5,alpha=.2)+
  geom_point(data=sdsum,aes(x=samp,y=msd,color=treatment),size=4)+
  geom_errorbar(data=sdsum,aes(x=samp,
    # ymin=msd-sdsd, ymax=msd+sdsd,
    ymin=msd-sdse*1.96, ymax=msd+sdse*1.96, #95% CI
    color=treatment),width=.4)+
  scale_color_viridis_d(
                   option="A",
                   begin=0, end=0.6,
                   name="",
                   labels=c("Control","Structure","Sponge"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        #    legend.position = "none", 
        legend.position = "top",
        strip.background = element_blank(),
        strip.text = element_blank())+
  ylab(sdlab)+
  xlab("Months into experiment")+
  facet_wrap(~seagrass,scales="free_y",nrow=2)

sdp<-egg::tag_facet(sdp, hjust = 0, y = Inf,
  tag_pool = c(
    "(a)", "(b)"
    ## tried adding italic names but it didn't work
    # expression("(a)"~italic(Ttestudinum)), 
    # expression("(b)"~italic(S filliforme)~"and"~italic(H wrightii))
    ), open = " ", close = " ", fontface = 2)
sdp
ggsave(filename="figures/shootdensity-95CI.jpg", plot=sdp, width = 3.5, height=4.5)

# figure 3 - seagrass growth
source("scripts/05_growth_analysis.R")
sgf<-sgg2%>%
  # mutate(dy=case_when(
  #   yr==1 & dist_factor=="closer"~"Yr 1 near",
  #   yr==1 & dist_factor!="closer"~"Yr 1 far",
  #   yr==2 & dist_factor=="closer"~"Yr 2 near",
  #   yr==2 & dist_factor!="closer"~"Yr 2 far"))%>%
  group_by(treatment,dist_factor,season,yr)%>%
  # group_by(treatment,dy,season)%>%
  summarize(mdsg=mean(delta_gpd_st,na.rm=TRUE),
    sddsg=sd(delta_gpd_st,na.rm=TRUE),
    sedsg=sddsg/sqrt(5*2*5)
    )

sgf$season<-factor(sgf$season,levels=c("s","w"),labels=c("Summer","Winter"))

sgf$dist_factor<-factor(sgf$dist_factor,levels=c("closer","farther"),labels=c("Near","Far"))

sglab<-expression(paste(Delta," Seagrass growth (mm "^2,"d "^"-1",")"))

sgplot <- ggplot(sgf, aes(x=as.factor(yr),y=mdsg,group=treatment,color=treatment)) +
  geom_hline(aes(yintercept=0), linetype = "dashed", colour = "darkgrey")+
  geom_point(size=5,position=position_dodge(0.5))+
  geom_path(position=position_dodge(0.5))+
  geom_errorbar(aes(x=as.factor(yr),
    # ymin=mdsg-sddsg,ymax=mdsg+sddsg),
    ymin=mdsg-sedsg*1.96,ymax=mdsg+sedsg*1.96),
    width=.2,position=position_dodge(0.5))+
  facet_grid(season~dist_factor,scales = "free")+
  # scale_color_brewer(type="qual",
  #                    palette="Set2",
  #                    name="",
  #                    labels=c("Control","Structure","Sponge"))+
 scale_color_viridis_d(option="A",
                  begin=0, end=0.6,
                  name="",
                  labels=c("Control","Structure","Sponge"))+
  scale_y_continuous(expand=expansion(add = c(1, 2))) +
  theme_bw()+
  xlab("Year")+
  ylab(sglab)
(sgplot <- egg::tag_facet(sgplot,
  #x = Inf, y = Inf, # hjust = -1.5,
  open = "(", close = ")", fontface = 2)+
  theme(panel.grid=element_blank(),
        panel.spacing = unit(1, "lines"),
        strip.background = element_blank(),
        #legend.position = "none",
        legend.position = "top",
        legend.text = element_text(size=12),
        axis.title=element_text(size=14),
        axis.text = element_text(size=12),
        strip.text = element_text(size=14)))
# ggsave("figures/sggrow.jpg",width=6,height=6.25)
ggsave("figures/sggrow-95CI.jpg", plot=sgplot, width=5.5,height=5.5)


# figure 4 - nutrients
source("scripts/06_nutrient_analysis.R")
pn<-ggplot(sgn %>%
    group_by(treatment, sampling, nut, dist) %>%
    summarize(mn = mean(nvalue), sdn = sd(nvalue),
      se = sdn/sqrt(5)) %>%
    filter(sampling %in% c(1, 4) &
             dist <= 0.5 &
             nut %in% c("PN")) %>% 
    mutate(
      dist = if_else(dist == 0, "0 m", "0.5 m")
    )
  ) +
   geom_line(aes(group=treatment,
                 x=as.factor(sampling),
                 y=mn,
                 color=treatment),
             position = position_dodge(width = .5))+
  geom_point(aes( x = as.factor(sampling),
      y = mn,
      color = treatment),
    size = 4,
    position = position_dodge(width = .5)) +
  geom_errorbar(aes(x = as.factor(sampling),
      # ymin = mn - sdn, ymax = mn + sdn,
      ymin = mn - se*1.96, ymax = mn + se*1.96,
      color = treatment),
    width = .1,
    position = position_dodge(width = .5))+
  scale_y_continuous(expand=expansion(add = c(0.05, 0.1)),
    breaks=c(1.60,1.80,2.00),labels=c("1.60","1.80","2.00")) +
  ylab("% Nitrogen")+
  scale_x_discrete(labels = c(0,12))+
  xlab("Months into the experiment")+
  # scale_color_brewer(type="qual",
  #                     palette="Set2",
  #                     name="",
  #                     labels=c("Control","Structure","Sponge"))+
scale_color_viridis_d(option="A",
                 begin=0, end=0.6,
                 name="",
                 labels=c("Control","Structure","Sponge"))+
  facet_grid(cols = vars(dist)) + 
  theme_bw()+
  theme(panel.grid = element_blank(),
         legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(size=14),
        axis.title.x = element_blank())

pp<-ggplot(sgn %>%
    group_by(treatment, sampling, nut, dist) %>%
    summarize(mn = mean(nvalue), sdn = sd(nvalue),
      se = sdn/sqrt(5)) %>%
    filter(sampling %in% c(1, 4) &
             dist <= 0.5 &
             nut %in% c("PP"))) +
  geom_line(aes(group=treatment,
                 x=as.factor(sampling),
                 y=mn,
                 color=treatment),
             position = position_dodge(width = .5))+
  geom_point(aes( x = as.factor(sampling),
      y = mn,
      color = treatment),
    size = 4,
    position = position_dodge(width = .5)) +
  geom_errorbar(aes(x = as.factor(sampling),
      # ymin = mn - sdn, ymax = mn + sdn,
      ymin = mn - se*1.96, ymax = mn + se*1.96,
      color = treatment),
    width = .1,
    position = position_dodge(width = .5))+  
  scale_y_continuous(expand=expansion(add = c(0.001, 0.005)), 
    breaks=c(0.05,0.06,0.07)) +
  # ylab("")+
  ylab("% Phosphorus")+
  # scale_color_brewer(type="qual",palette="Set2")+
  scale_color_viridis_d(option="A", begin=0, end=0.6)+
  scale_x_discrete(labels = c(0,12))+
  xlab("Months into the experiment")+
  facet_grid(cols = vars(dist)) +
  theme(panel.grid = element_blank(),
    # axis.title.x = element_blank(),
         legend.position = "none")

pc<-ggplot(sgn %>%
    group_by(treatment, sampling, nut, dist) %>%
    summarize(mn = mean(nvalue), sdn = sd(nvalue),
      se = sdn/sqrt(5)) %>%
    filter(sampling %in% c(1, 4) &
             dist <= 0.5 &
             nut %in% c("PC"))) +
  geom_line(aes(group=treatment,
                 x=as.factor(sampling),
                 y=mn,
                 color=treatment),
             position = position_dodge(width = .5))+
  geom_point(aes( x = as.factor(sampling),
      y = mn,
      color = treatment),
    size = 4,
    position = position_dodge(width = .5)) +
  geom_errorbar(aes(x = as.factor(sampling),
      # ymin = mn - sdn, ymax = mn + sdn,
      ymin = mn - se*1.96, ymax = mn + se*1.96,
      color = treatment),
    width = .1,
    position = position_dodge(width = .5))+
  scale_y_continuous(expand=expansion(add = c(0.5, 1.75)), 
    breaks = c(30.0, 33.0, 36.0),labels = c("30.0", "33.0", "36.0")) +
  # ylab("")+
  ylab("% Carbon")+
  # scale_color_brewer(type="qual",palette="Set2")+
  scale_color_viridis_d(option="A", begin=0, end=0.6)+
  scale_x_discrete(labels = c(0,12))+
  xlab("Months into the experiment")+
  facet_grid(cols = vars(dist)) +   
  theme(panel.grid = element_blank(),
    # axis.title.x = element_blank(),
         legend.position = "none")

library(cowplot)
lgnd<-get_legend(pn +
    guides(color = guide_legend(nrow = 1, align_plots = "right")) +
    theme(legend.position = "bottom"))
# hjust = 0, y = Inf,
  # tag_pool = c('(a)  Control','(b)  Structure', "(c)  Sponge"),
  # open = " ", close = " "
pn<-egg::tag_facet(pn
  , tag_pool = c("a", "b")
  # , tag_pool = c('(a)  Nitrogen'), open = " ", close = " ", hjust = 0, y = Inf
  )+ 
  # ylab("") +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size=14),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),axis.title.x = element_blank()) 
pc<-egg::tag_facet(pc 
  , tag_pool = c("c","d"))
  # , tag_pool = c('(c)  Carbon'), open = " ", close = " ", hjust = 0, y = Inf)+
  # theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),
  # axis.title.X = element_blank())+ 
  # ylab("")
pp<-egg::tag_facet(pp
  , tag_pool = c("e", "f")
  # , tag_pool = c('(b)  Phosphorus'), open = " ", close = " ", hjust = 0, y = Inf
  )+ 
  # ylab("Nutrient concentrations (%)") + 
  theme(axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()) 
#without carbon
#n1<-plot_grid(lgnd,pn,pp,labels = c("","A","B"),nrow=3,ncol = 1,rel_heights = c(.1,1,1))
#with carbon
(n2<-plot_grid(
  lgnd,pn,pp,pc,nrow=4,ncol = 1,rel_heights = c(.1,1,1,1.1), align = "right"
  #pn,pc,pp,nrow=3,ncol = 1,rel_heights = c(1,1,1.1)
  ))
# ggsave(filename = "figures/nuts.pdf",plot=n2,width=3,height = 6)
# ggsave(filename = "figures/nuts.jpg",plot=n2,width=3,height = 6)
# ggsave(filename = "figures/nuts-95CI.jpg",plot=n2,width=3.5,height=6.5)
ggsave(filename = "figures/nuts-95CI-w-0.5.jpg",plot=n2,width=6,height=6.5)


# Figures for supplement


# Algae 
alg2<-readxl::read_xlsx(here("Original_data","ForFinella_Transplant_data.xlsx"),sheet="Algae")

firstup <- function(x){
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

a.taxa.abun<-alg2 %>% 
  pivot_longer(c(-Treatment,-plot,-sampling),names_to = "taxa",values_to = "abundance")%>% 
  rename(treatment=Treatment)%>%
  mutate(sampling=case_when(
    sampling==1~0,
    sampling==2~1,
    sampling==3~5,
    sampling==4~12,
    sampling==5~17),
    taxa=ifelse(taxa=="cladocephalus","udotea",taxa))%>%
  filter(!taxa %in% c("brown.cyanobacteria","green.cyanobacteria","dictyota"))%>%
  mutate(taxa = firstup(taxa)) %>%
  group_by(treatment,sampling,taxa)%>%
  summarize(abun.m=mean(abundance),abun.sd=sd(abundance))

(a1 <- ggplot()+
  geom_line(data=a.taxa.abun,
    aes(x=sampling,y=abun.m,group=treatment,color=treatment),position=position_dodge(0.5))+
  geom_errorbar(data=a.taxa.abun,aes(x=sampling,ymin=abun.m-abun.sd,ymax=abun.m+abun.sd,color=treatment),width=.5,position=position_dodge(0.5))+
  geom_point(data=a.taxa.abun,aes(x=sampling,y=abun.m,color=treatment),size=3,position=position_dodge(0.5))+
  facet_wrap(~taxa, ncol=1, scales = "free_y") +
  theme_bw()+
  theme(panel.grid = element_blank(),
    axis.text = element_text(size=12),
    axis.title = element_text(size=14),
    legend.text = element_text(size=10))+
  scale_color_viridis_d(option="A", begin=0, end=0.6,"")+
  ylab("Adundance")+
  xlab("Months into Experiment"))

ggsave("figures/algal_sp_abundance.jpg",plot=a1,width=5,height=7.5)

# several algae species also took off only in sponge plots

