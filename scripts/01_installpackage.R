## this is a project to learn best practices of data management in R ##
## this script is only to load packages and data ##

## @knitr loadpackages
# load packages ####
# these should be all packages needed to run all scripts in the Ircinia Project. 
#**note** If you are using a new computer, you may need to go to https://cran.r-project.org/bin/windows/Rtools/ website and download RTools, if there are errors when running the scripts.
#If Rstudio still says you do not have RTools use the 2 lines of code below
#writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
#Sys.which("make")

if(!require(here))install.packages('here');library(here)
if(!require(tidyverse))install.packages('tidyverse');library(tidyverse)
if(!require(readxl))install.packages('readxl');library(readxl)
if(!require(vegan))install.packages('vegan');library(vegan)

if(!require(usethis))install.packages('usethis')
#usethis::use_git_config(user.name = "username", user.email = "email")
if(!require(bookdown))install.packages('bookdown')


if(!require(lmerTest))install.packages("lmerTest");library(lmerTest)
if(!require(DHARMa))install.packages("DHARMa");library(DHARMa)
if(!require(glmmTMB))install.packages("glmmTMB");library(glmmTMB)
if(!require(ggeffects))install.packages("ggeffects");library(ggeffects)
if(!require(car))install.packages("car");library(car)
if(!require(ggpubr))install.packages("ggpubr");library(ggpubr)
if(!require(egg))install.packages("egg")
if(!require(cowplot))install.packages("cowplot");library(cowplot)