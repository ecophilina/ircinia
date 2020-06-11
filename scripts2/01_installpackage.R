## this is a project to learn best practices of data management in R ##
## this script is only to load packages and data ##

## @knitr loadpackages
# load packages ####

if(!require(here))install.packages('here');library(here)
if(!require(tidyverse))install.packages('tidyverse');library(tidyverse)
if(!require(readxl))install.packages('readxl');library(readxl)
if(!require(vegan))install.packages('vegan');library(vegan)

if(!require(usethis))install.packages('usethis')
#usethis::use_git_config(user.name = "username", user.email = "email")
