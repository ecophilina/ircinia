# How does the sponge Ircinia felix influence seagrass bed primary producers? #

# this script will examine seagrass shoot density #
# run 03_reimport.R fist #

# packages----
if(!require(tidyverse))install.packages("tidyverse");library(tidyverse)
if(!require(ggpubr))install.packages("ggpubr");library(ggpubr)

theme_set(theme_bw())


# notes ----
# what does the shoot data look like
glimpse(sg_shoot)
str(sg_shoot)
head(sg_shoot)

# sg_shoot <- as_tibble(sg_shoot)

# there are four replicates within each plot/dist/sampling combination
# T.SD = T. testudinum 
# SH.SD = S. filliforme and H. wrightii 
# SD = total density

# transform longer
ddat <- sg_shoot %>% pivot_longer(
  T.SD:SD, 
  names_to = "type", values_to = "density"
  ) %>% mutate(
    id = paste0(plot, "_", dist),
    season = case_when(
      sampling==1~"summer",
      sampling==2~"summer",
      sampling==4~"summer",
      sampling==3~"winter",
      sampling==5~"winter")
  ) 

glimpse(ddat)

ddat %>% filter(type== "T.SD") %>%
  ggplot(aes(as.factor(dist), density))+
  geom_point(alpha=0.5) +
  facet_grid(treatment ~ sampling)

ddat %>% filter(type== "SH.SD") %>%
  ggplot(aes(as.factor(dist), density))+
  geom_point(alpha=0.5) +
  facet_grid(treatment ~ sampling)

ddat %>% filter(type== "SD") %>%
  ggplot(aes(as.factor(dist), density))+
  geom_point(alpha=0.5) +
  facet_grid(treatment ~ sampling)

ddat %>% filter(type== "T.SD") %>%
ggplot(aes(as.factor(sampling), density, group = id))+
  geom_point() +
  geom_line() +
  facet_grid(treatment ~ dist)

ddat %>% filter(type== "SH.SD") %>%
ggplot(aes(as.factor(sampling), density, group = id))+
  geom_point() +
  geom_line() +
  facet_grid(treatment ~ dist)


ddat %>% filter(type== "SD") %>%
  ggplot(aes(as.factor(sampling), density))+
  geom_violin() +
  facet_grid(rows=vars(treatment))+
  ggtitle("Total")

ddat %>% filter(type== "T.SD") %>%
  ggplot(aes(as.factor(sampling), density))+
  geom_violin() +
  facet_grid(rows=vars(treatment))+
  ggtitle("T.SD")

ddat %>% filter(type== "SH.SD") %>%
  ggplot(aes(as.factor(sampling), density))+
  geom_violin() +
  facet_grid(rows=vars(treatment))+
  ggtitle("SH.SD")


# replicates are confusing things for now... 
# since they aren't fixed through time will need to by averaged to assess change anyway


ddat2 <- ddat %>% group_by(id, type, sampling, treatment) %>% 
  mutate(raw_density = density, density = mean(raw_density)) %>% 
  select(-raw_density) %>% ungroup() %>% unique() 


ddat2 %>% filter(type== "T.SD") %>%
  ggplot(aes(as.factor(sampling), density, group = id))+
  geom_point() +
  geom_line() +
  facet_grid(treatment ~ dist)+
  ggtitle("T.SD")

ddat2 %>% filter(type== "SH.SD") %>%
  ggplot(aes(as.factor(sampling), density, group = id))+
  geom_point() +
  geom_line() +
  facet_grid(treatment ~ dist) +
  ggtitle("SH.SD")

ddat2 %>% filter(type== "SD") %>%
  ggplot(aes(as.factor(sampling), density, group = id))+
  geom_point() +
  geom_line() +
  facet_grid(treatment ~ dist) +
  ggtitle("Total")


### add change variables
glimpse(ddat2)

start_val <- filter(ddat2, sampling == 1) %>% 
  rename(start_dens = density) %>% select(id, type, start_dens)

ddat3 <- left_join(ddat2, start_val) %>% 
  mutate(delta = density - start_dens) %>% filter(sampling != 1) 

ddat3 %>% filter(type== "SD") %>%
  ggplot(aes(as.factor(sampling), delta, group = id))+
  geom_hline(yintercept = 0, colour= "red") +
  geom_point() +
  geom_line() +
  facet_grid(treatment ~ dist) +
  ggtitle("Total")

ddat3 %>% filter(type== "SH.SD") %>%
  ggplot(aes(as.factor(sampling), delta, group = id))+
  geom_hline(yintercept = 0, colour= "red") +
  geom_point() +
  geom_line() +
  facet_grid(treatment ~ dist) +
  ggtitle("SH.SD")

ddat3 %>% filter(type== "T.SD") %>%
  ggplot(aes(as.factor(sampling), delta, group = id))+
  geom_hline(yintercept = 0, colour= "red") +
  geom_point() +
  geom_line() +
  facet_grid(treatment ~ dist) +
  ggtitle("T.SD")


