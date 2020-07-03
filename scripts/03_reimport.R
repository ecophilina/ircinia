## this is a project to learn best practices of data management in R ##

## this script will import the reorganized data and do initial data explorations ##
# packages----
if(!require(here))install.packages("here");library(here)
if(!require(tidyverse))install.packages("tidyverse");library(tidyverse)
# @knitr reimport

## import reworked data ####
files <- list.files(here("working_data"), pattern = "*.csv")
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


