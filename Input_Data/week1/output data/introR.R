##### Intro to R with the CAL FIRE MAJOR FIRES dataset (2013-2019) #####
# Tim Lamour 
# 10/08/2020

# install.packages("tidyverse")
library(tidyverse)
?tidyverse
# install.packages("dplyr")
library(dplyr)
?dplyr
# install.packages("readxl")
library(readxl)
# install.packages("praise")
library(praise)

#### Load Dataset ####

excel_sheets("Input_Data/week1/2013_2019_CALFIRE_Redbook.xlsx")

metadata <- read_excel("Input_Data/week1/2013_2019_CALFIRE_Redbook.xlsx", sheet = 1)
view(metadata)
data <- read_excel("Input_Data/week1/2013_2019_CALFIRE_Redbook.xlsx", sheet = "Data")
view(data)

# to clear your environments you can type in the following: rm(list = ls())

#### Initial Data Exploration ####

names(data)
?names
dim(data)
class(data)
head(data)
tail(data)
str(data)
glimpse(data)
typeof(data$Total_Acres_Burned)

# single columns can be reffered to using '$' 

max(data$Total_Acres_Burned)
max(data$Structures_Destroyed)
?max
max(data$Structures_Destroyed, na.rm = T)

summary(data)

#### Basic data wrangling (dplyr functions) ####

df1 <- select(data, County_Unit:Controlled_Date, Total_Acres_Burned, Cause:Structures_Damaged)
unique(df1$County_Unit)

df2 <- filter(df1, County_Unit %in% c("SANTA BARBARA", "VENTURA", "LOS ANGELES", "SAN DIEGO", "ORANGE", "VENTURA/SANTA BARBARA") & Total_Acres_Burned >= 500)

# | = "or" 
# == = "equals"/"matches" , %in% c()
# & = "and" 

df3 <- arrange(df2, desc(Total_Acres_Burned))

df4 <- mutate_at(df3, vars("Structures_Destroyed", "Structures_Damaged"), replace_na, 0)

df5 <- mutate(df4, struc_impact = Structures_Damaged + Structures_Destroyed)

#mess with time
library(lubridate)

df6 <- mutate(df5, interv = interval(Start_Date, Controlled_Date), dur = as.duration(interv), days = as.numeric(dur, "days"))

#### Introduction to piping ####

socal.fires <- data %>% 
  select(County_Unit:Controlled_Date, Total_Acres_Burned, Cause:Structures_Damaged) %>% 
  filter(County_Unit %in% c("SANTA BARBARA", "VENTURA", "LOS ANGELES", "SAN DIEGO", "ORANGE", "VENTURA/SANTA BARBARA") & Total_Acres_Burned >= 500) %>%
  arrange(desc(Total_Acres_Burned)) %>% 
  mutate_at(vars("Structures_Destroyed", "Structures_Damaged"), replace_na, 0) %>% 
  mutate(struc_impact = Structures_Damaged + Structures_Destroyed, interv = interval(Start_Date, Controlled_Date), dur = as.duration(interv), days = as.numeric(dur, "days"))

#### first graphs with ggplot ####

# we're going to make a graph of acres burned in the south coast from 2013-18 with the color dependent on which county we're showing
# three things you must tell R to make a grpah in ggplot: 
# 1 - that you're using ggplot 
# 2 - what data you're using (including what should be x and what should be y)
# 3 - what type of graph you want to ceate
# everything after that is extra to make it pretty 

ggplot(socal.fires, aes(x = Start_Date, y = Total_Acres_Burned)) + 
  geom_point(aes(color = County_Unit)) +
  ggtitle("CA South Coast Major Fires \n2014 - 2018") + 
  labs(x = "", y = "Total Acres Burned", color = "County") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_grid(rows = "County_Unit", scales = "free")
  
plot.data <- socal.fires %>% 
  rename(county = County_Unit,
         acres = Total_Acres_Burned,
         start = Start_Date, 
         end = Controlled_Date) %>% 
  mutate(year = year(start),
         county = ifelse(county == "VENTURA/SANTA BARBARA", "VENTURA", county)) 

incidents <- plot.data %>%  
  group_by(county, year) %>% 
  tally()
  ungroup()
  
incidents.plot <- incidents %>%  
  ggplot(aes(x = year, y = n)) +
  geom_point(aes(color = county)) +
  geom_line(aes(color = county)) + 
  geom_point(color = "blue") +
  geom_line(color = "blue") +
  labs(title = "CA South Coast Major Fire Inciddents \n 2014-2018", x = "", y = "Incidents", color = "County") + 
  theme_bw() + 
  facet_grid(rows = "county", scales = "free") + 
  guides(color = F)


all_incidents <- plot.data %>% 
  group_by(year) %>% 
  tally()%>% 
  ungroup()

all_incidents.plot <- all_incidents %>%  
  ggplot(aes(x = year, y = n)) +
  geom_point(color = "blue") +
  geom_line(color = "blue") +
  labs(title = "CA South Coast Major Fire Inciddents \n 2014-2018", x = "", y = "Incidents") + 
  theme_bw()

##### Save Data and Plots ####

saveRDS(socal.fires, file = "Input_Data/week1/output data/socal_fires_data.rds")
write_csv(socal.fires, "Input_Data/week1/output data/socal_fires_data.csv")  

ggsave(filename = "Fire_Incidents", all_incidents.plot, device = "jpeg", "Input_Data/week1/output data/" )


## new change
