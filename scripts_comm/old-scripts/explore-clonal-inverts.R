# plot clonal inverts 
source("scripts_comm/02_community_data_org.R")

col.inv.spp<-inverts%>%
  filter(taxa %in% col.inverts)%>%
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
      sampling==17~2))%>%
  left_join(productivity)

col.inv.uni<-col.inv.env%>%
  mutate(spr=vegan::specnumber(col.inv.com),
    div=vegan::diversity(col.inv.com,index = "shannon"),
    j=div/log(spr),
    j=ifelse(is.na(j),0,j),
    i.abund=rowSums(col.inv.com),
    wg_tunicate=col.inv.com$`white and green tunicate`
  )


# col.inv.com.full
ggplot(data = col.inv.spp #%>% filter(season == "summer" ) %>% 
  # mutate(treatment = factor(treatment, 
  #                           levels = c("blank", "fake", "real"), 
  #                           labels = c("Control", "Structure Control", "Sponge")))
) + 
  geom_jitter(aes(sg.sd.c, abundance, colour = treatment, alpha = as.factor(sampling)), 
    width = 0.2, height = 0.2, size = 2.5) + 
  # geom_line(aes(a.abund, spr, group = as.factor(plot))) +
  # geom_smooth(method = "lm", aes(sg.sd.c, spr), colour = "black") + 
  scale_alpha_discrete(range = c(0.3,1), name = "Months") + 
  scale_color_viridis_d(option="A", begin=0, end=0.6,"")+
  # xlab("Algae abundance") +
  ylab("Clonal Invertebrate species richness") +
  facet_wrap(~taxa, scales = "free") +
  # coord_cartesian(expand = F) +
  ggsidekick::theme_sleek()
# ggsave("figures/clonal-invert-species-by-sg.png", width = 5, height = 3)

# col.inv.com.full
ggplot(data = col.inv.spp #%>% filter(season == "summer" ) %>% 
  # mutate(treatment = factor(treatment, 
  #                           levels = c("blank", "fake", "real"), 
  #                           labels = c("Control", "Structure Control", "Sponge")))
) + 
  geom_jitter(aes(pp.struct, abundance, colour = treatment, alpha = as.factor(sampling)), 
    width = 0.2, height = 0.2, size = 2.5) + 
  # geom_line(aes(a.abund, spr, group = as.factor(plot))) +
  # geom_smooth(method = "lm", aes(sg.sd.c, spr), colour = "black") + 
  scale_alpha_discrete(range = c(0.3,1), name = "Months") + 
  scale_color_viridis_d(option="A", begin=0, end=0.6,"")+
  # xlab("Algae abundance") +
  ylab("Clonal Invertebrate species richness") +
  facet_wrap(~taxa, scales = "free") +
  # coord_cartesian(expand = F) +
  ggsidekick::theme_sleek()
# ggsave("figures/clonal-invert-spr-by-algae.png", width = 5, height = 3)


ggplot(data = col.inv.uni %>% filter(season == "summer" & sampling > 0) %>% 
  mutate(invasive = if_else(wg_tunicate > 0, 1, 0),
    spr.not.wg = if_else(invasive == 1, spr-1, as.double(spr))
    )
) + 
  geom_jitter(aes(as.factor(invasive), spr.not.wg, colour = treatment, alpha = as.factor(sampling)), 
    width = 0.2, height = 0.2, size = 2.5) + 
 # geom_boxplot(aes(as.factor(invasive), spr.not.wg, colour = treatment)) +
  # geom_smooth(method = "lm", aes(sg.sd.c, spr), colour = "black") + 
  scale_alpha_discrete(range = c(0.3,1), name = "Months") + 
  scale_color_viridis_d(option="A", begin=0, end=0.6,"")+
  # xlab("Algae abundance") +
  ylab("Clonal Invertebrate species richness") +
  # coord_cartesian(expand = F) +
  ggsidekick::theme_sleek()
# ggsave("figures/clonal-invert-spr-by-algae.png", width = 5, height = 3)

ggplot(data = col.inv.uni %>% filter(season == "summer") %>% 
    mutate(invasive = if_else(wg_tunicate > 0, 1, 0),
      spr.not.wg = if_else(invasive == 1, spr-1, as.double(spr))
    )
) + 
  geom_jitter(aes(as.factor(sampling), spr, colour = treatment, alpha = as.factor(sampling)), 
    width = 0.1, height = 0, size = 5) + 
  geom_line(aes(as.factor(sampling), spr, colour = treatment, group = as.factor(plot)), alpha = 0.2, size = 2) +
  scale_alpha_discrete(range = c(0.3,1), name = "Months") + 
  scale_color_viridis_d(option="A", begin=0, end=0.6,"")+
  # xlab("Algae abundance") +
  ylab("Clonal Invertebrate species richness") +
  # coord_cartesian(expand = F) +
  ggsidekick::theme_sleek()
# ggsave("figures/clonal-invert-spr-by-algae.png", width = 5, height = 3)

