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

# sampling 2 is 1 month after sampling 1
# sampling 3 is 5 months after 2
# sampling 4 is 12 months after 2 (7 after 3)
# sampling 5 is 17 months after 2 (5 after 4)

sg_shoot %>% ggplot(aes(T.SD, SH.SD )) + geom_point() + 
  facet_grid(treatment~sampling)

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


m <- aov(density~treatment, data = filter(ddat, dist=="0-1" & sampling == 1))
summary(m)
TukeyHSD(m)


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

start_val <- filter(ddat2, sampling == 2) %>% 
  rename(start_dens = density) %>% select(id, type, start_dens)

ddat3 <- left_join(ddat2, start_val) %>% 
  mutate(delta = density - start_dens) %>% filter(!(sampling %in% c(1,2))) 

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

ddat3 %>% filter(type== "T.SD") %>%
  ggplot(aes((sampling), delta))+
  geom_hline(yintercept = 0, colour= "red") +
  geom_point() +
  geom_smooth(method = lm) +
  facet_grid(treatment ~ dist) +
  ggtitle("T.SD")

m <- lm(delta~treatment, data = filter(ddat3, dist=="0-1" & sampling == 2))
summary(m)
m <- lm(delta~treatment, data = filter(ddat3, dist=="1-2" & sampling == 2))
summary(m)
m <- lm(delta~treatment, data = filter(ddat3, dist=="0-1" & sampling == 3))
summary(m)
m <- lm(delta~treatment, data = filter(ddat3, dist=="1-2" & sampling == 3))
summary(m)



