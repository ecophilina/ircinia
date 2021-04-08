## this script is only to load packages used in community analysis ##
## @knitr loadpackages
# load packages ####
#**note** If you are using a new computer, you may need to go to https://cran.r-project.org/bin/windows/Rtools/ website and download RTools, if there are errors when running the scripts.
#If Rstudio still says you do not have RTools use the 2 lines of code below
#writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
#Sys.which("make")

if(!require(here))install.packages('here')
if(!require(tidyverse))install.packages('tidyverse')
if(!require(readxl))install.packages('readxl')
if(!require(vegan))install.packages('vegan')

if(!require(usethis))install.packages('usethis')
#usethis::use_git_config(user.name = "username", user.email = "email")
if(!require(bookdown))install.packages('bookdown')

if(!require(glmmTMB))install.packages("glmmTMB")
# if(!require(lmerTest))install.packages("lmerTest")
if(!require(DHARMa))install.packages("DHARMa")
if(!require(vegan))install.packages('performance')


if(!require(ggeffects))install.packages("ggeffects")
if(!require(car))install.packages("car")
if(!require(ggpubr))install.packages("ggpubr")
if(!require(egg))install.packages("egg")
if(!require(cowplot))install.packages("cowplot")

# install.packages("devtools")
if(!require(ggnewscale))devtools::install_github("eliocamp/ggnewscale")

install.packages("see", dependencies = TRUE)