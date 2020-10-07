##### EEMB 144L Intro to R with CAL FIRE data #####
#Tim Lamour
#10/2/2020

# This is a comment, below this will be a command

Library("tidyverse")
library(tidyverse)

install.packages("readxl")              
library(readxl)

##### Load Data #####

#unlike .csv files .xlsx files can have multiple data sheets

excel_sheets("Input_Data/week1/2013_2019_CALFIRE_Redbook.xlsx")

calfire.metadata <- read_excel("Input_Data/week1/2013_2019_CALFIRE_Redbook.xlsx", sheet = metadata
calfire.data <- read_excel("Input_Data/week1/2013_2019_CALFIRE_Redbook.xlsx", sheet = 2)                               
