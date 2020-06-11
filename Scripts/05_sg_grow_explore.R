## this is a project to learn best practices of data management in R ##
## this script is only to load packages and data ##

## @knitr loadpackages
# load packages ####

if(!require(here))install.packages('here');library(here)
if(!require(tidyverse))install.packages('tidyverse');library(tidyverse)
if(!require(readxl))install.packages('readxl');library(readxl)
## this is a project to learn best practices of data management in R ##

## this script will import the reorganized data and do initial data explorations ##

# @knitr reimport

## import reworked data ####
files <- list.files(here("working_data"), pattern = "2020-06-04")
# this bit of code finds all the files in the working data directory that match the date I
# exported the reorganized data you can replace the date with "*.csv" if you want to find all the csv
# files in a folder- or replace the date with any text pattern really.

fd <- data.frame(files = files, dname = files) %>%
  # this line creates a data frame that I'll use to import the files in a couple lines
  separate(dname, into = c("dname", "extra"), sep = "2020")
# this bit separates out of the name of the dataset from the date and the file extension- note
# we still have one column with the full file name

for (i in 1:nrow(fd))assign(fd[i, 2], read.csv(here("working_data", fd[i, 1]))[, -1])
# this bit of code says for every row in the fd data frame assign the file in column 1
# to the name in column 2. You should now have all the reorganized data sets imported.

#growth per day added as a new column
sg_grow$gd <- sg_grow$total.growth.mm2/sg_grow$days

#adding seasons to the plots
sg_grow <- sg_grow %>% 
  mutate(season=case_when(
    sampling==1~"summer",
    sampling==2~"summer",
    sampling==3~"winter",
    sampling==4~"summer",
    sampling==5~"winter" ))

# @knitr dataexplore

bp<-ggplot(data=sg_grow)

#boxplot of growth per day by treatment sorted into sampling categories
bp+
  geom_boxplot(aes(x=treatment,y=gd))+ 
facet_wrap(~sampling)


#boxplot of gpd by treated sorted into seasons and distance
bp+
  geom_boxplot(aes(x=treatment,y=gd))+ 
  facet_grid(dist~season)

#boxplot of gpd by treatment categorized by distance and sampling
bp+
  geom_boxplot(aes(x=treatment,y=gd))+ 
  facet_grid(dist~sampling)

#change in growth over time(sampling) by plot
#changing the data set

bp<-sg_grow %>% filter(dist==0.5) %>% 
  group_by(plot,sampling,treatment) %>%
  summarise(mgd=mean(gd), sdgd=sd(gd)) %>%
  ggplot()

bp+
  geom_point(aes(x=sampling, y=mgd, color=treatment))+
  geom_line(aes(x=sampling, y=mgd, color=treatment, group=plot))


#comparing by distance

bp<-sg_grow %>% 
  group_by(plot,sampling,treatment,dist) %>%
  summarise(mgd=mean(gd), sdgd=sd(gd)) %>%
  ggplot()

bp+
  geom_point(aes(x=sampling, y=mgd, color=treatment))+
  geom_line(aes(x=sampling, y=mgd, color=treatment, group=plot))+
  facet_wrap(~dist)


#Philina's way----
install.packages('ggpubr')

library(ggpubr)
ggpaired(sg_grow, x = "season", y = "gd",
         id = "plot", facet.by = c("treatment") )

2:06
ggpaired(sg_grow, x = "season", y = "gd",
         id = "plot", facet.by = c("treatment")) +
  stat_compare_means(
    aes(label = paste0("p = ", ..p.format..)),
    method = "t.test",
    # paired = TRUE,
    ref.group = NULL)

#another way
ggpaired(sg_grow, x = "season", y = "gd",
         id = "plot", facet.by = c("treatment") )
ggpaired(sg_grow %<, x = "season", y = "gd",
         id = "plot", facet.by = c("treatment")) +
  stat_compare_means(
    aes(label = paste0("p = ", ..p.format..)),
    method = "t.test",
    # paired = TRUE,
    ref.group = NULL)