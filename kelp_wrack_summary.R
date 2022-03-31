###################################################################################
#                                                                                ##
# kelp_wrack_summary.R                                                           ##
# Data are current as of 2022-03-29                                              ##
# Data source: David Bilderback                                                  ##
# R code prepared by Ross Whippo                                                 ##
# Last updated 2022-03-29                                                        ##
#                                                                                ##
###################################################################################

# SUMMARY:
# Data on kelp wrack morphometrics collected by David Bilderback in Bandon, OR, USA
# from 2018 - present. 


# Required Files (check that script is loading latest version):
# mature_sporophyte_data.csv
# young_sporophyte_data.csv
# noaa_portorford_buoydata_2018.csv

# Associated Scripts:
# none

# TO DO

# 1) For each of the four years, how many days was the beach transect surveyed for 
#    young sporophytes?
#  
# 2) For each of the four years how many days were young sporophytes found on the 
#    transect?
#
# 3) Is there a relationship between stipe length and bulb diameter?
#
# 4) Is there a relationship between bulb diameter and widest blade width?
#
# 5) What is the mean diameter of holdfasts?
#
# 6) For each year, how many sporophytes were singular, and how many were clusters?
#
# 7) What is the count and percentage of each substrate type?
#
# 8) How common were cospecies? Numbers and percentages.
#
# 9) For each year, when does recruitment occur? (as function of % of 0.1-10cm
#    length stipes throughout the year?)
#
# 10) Create separate csv's of groups (0,1/clustered, not clustered)

# TO DO
# pull out northerly winds (they are south of the bed)
# use WSPD instead of GST
# summarise macro-epiphytes
# 180cm stipe length cutoff for epiphyte cover (what is actual cutoff for epibionts?)
# when do epibionts first appear after winter?
# size structure of sporopytes with unique epibiont communities?
# reliable way to age kelp?
# chamical or physical blocks to epibionts?
# temp treatments on kelp growth? SJI?
# don't differentiate epibiont TYPE, just PRESENCE/ABSENCE, related to stipe length
# color change of kelp has FA change?
# add epiphyte status to stipe_length by date figure
# power analysis - how much do you need to collect to get same answers? Stratified?
# subset these data sets and find minimum needed
# how can this be published???? 

###################################################################################
# TABLE OF CONTENTS                                                               #
#                                                                                 #
# RECENT CHANGES TO SCRIPT                                                        #
# LOAD PACKAGES                                                                   #
# READ IN AND PREPARE DATA                                                        #
# MANIPULATE DATA                                                                 #
#                                                                                 #
###################################################################################

###################################################################################
# RECENT CHANGES TO SCRIPT                                                        #
###################################################################################

# 2020-02-13 Script created
# 2022-03-29 Data updated, initial questions identified

###################################################################################
# LOAD PACKAGES                                                                   #
###################################################################################

library(tidyverse)
library(viridis)
library(lme4)
library(lubridate)

# options(max.print = 9999)


###################################################################################
# READ IN AND PREPARE DATA                                                        #
###################################################################################

# read in mature dataset
mature_sporophyte <- read_csv("Data/mature_sporophyte_data.csv", 
                              col_types = cols(date = col_character(),
                                               absd_sori_num = col_number(), 
                                               blade_width = col_number(), 
                                               bulb_diam = col_number(), 
                                               notes = col_character(), 
                                               re_sori_num = col_number(), 
                                               sori_num = col_number(), 
                                               stipe_length = col_number()))

# remove erroneous date (need to fix)
mature_sporophyte <- mature_sporophyte %>%
  filter(date != "1/3/1900") 
mature_sporophyte$date <- mature_sporophyte$date %>%
  recode("7/12/2008" = "7/12/2018")
# turn into date format
mature_sporophyte$date <- as.Date(mature_sporophyte$date, format = "%m/%d/%Y") 

# fix bulb diameter mistake
mature_sporophyte <- mature_sporophyte %>% 
  mutate(bulb_diam = replace(bulb_diam, bulb_diam == 72, 7.2))
                                
# read in young dataset
young_sporophyte <- read_csv("Data/young_sporophyte_data.csv", 
                                  col_types = cols(Date = col_character(), 
                                                   Bulb = col_double(), 
                                                   Blade = col_double(), 
                                                   `Ho Fa` = col_double(), 
                                                   Single = col_character(), 
                                                   Group = col_character()))

# fix erroneous date
young_sporophyte$Date <- young_sporophyte$Date %>%
  recode("11/28/2028" = "11/28/2018",
         "10/9/2028" = "10/9/2018",
         "6/17/2028" = "6/17/2018",
         "56/2/2021" = "5/26/2021")
young_sporophyte$Date <- mdy(young_sporophyte$Date)
# fix capitalization of missing values
young_sporophyte$Stipe <- young_sporophyte$Stipe %>%
  recode("None" = "none")

# import NOAA weather data
noaa_portorford_buoydata_2018 <- read_table("Data/noaa_portorford_buoydata_2018.csv", 
                                             col_types = cols(ATMP = col_number(), 
                                                              GST = col_number(), 
                                                              PRES = col_number(), 
                                                              WSPD = col_number(), 
                                                              WTMP = col_number()))

# remove first metadata row
noaa_portorford_buoydata_2018 <- noaa_portorford_buoydata_2018[-1,]

# create date/time column
noaa_portorford_buoydata_2018 <- noaa_portorford_buoydata_2018 %>%
  unite("date", YY:DD, sep = "-", remove = FALSE)
noaa_portorford_buoydata_2018$date <- as.Date(noaa_portorford_buoydata_2018$date)


###################################################################################
# NOAA WEATHER DATA                                                               #
###################################################################################

# daily means of wind
daily_wind_means <- noaa_portorford_buoydata_2018 %>%
  group_by(date) %>%
  summarise_at(c("WSPD", "GST"), mean, na.rm = TRUE)

# extract mean gusts > 10 m/s
gust_mean_10m <- daily_wind_means %>%
  filter(GST > 10)

# extract mean gusts > 15 m/s
gust_mean_15m <- daily_wind_means %>%
  filter(GST > 15)

ggplot(daily_wind_means, aes(x = date, y = GST)) +
  geom_point()

## test do high wrack days follow windy days?

# tally number of obs per day
wrack_obs <- mature_sporophyte %>%
  group_by(date) %>%
  tally()

# make date a character vector in new winds dataset for joining
daily_wind_means_char <- daily_wind_means

# reduce number of days to match weather data
wrack_obs_reduced <- wrack_obs %>%
  filter(date %in% daily_wind_means_char$date)


# join to wind data
wrack_by_wind <- left_join(wrack_obs, daily_wind_means_char, by = "date")

# quick plot of GST by wrack
ggplot(wrack_by_wind, aes(x = GST, y = n)) +
  geom_point()

# linear model wind gust means 
lm(n ~ GST, data = wrack_by_wind)

# test by offset dataset
wrack_by_wind_offset <- mutate(wrack_by_wind_value, offset = lag(n = 5, n, order_by = date))

# plot of offset
ggplot(wrack_by_wind_offset, aes(x = GST, y = offset)) +
  geom_point()


### Use selection model to determine if gust, wind, have effect on number of wrack observations

wrack_offset <- wrack_by_wind
wrack_offset1 <- mutate(wrack_by_wind, offset1 = lead(n = 1, WSPD, order_by = date))
wrack_offset2 <- mutate(wrack_by_wind, offset2 = lead(n = 2, WSPD, order_by = date))
wrack_offset3 <- mutate(wrack_by_wind, offset3 = lead(n = 3, WSPD, order_by = date))
wrack_offset4 <- mutate(wrack_by_wind, offset4 = lead(n = 4, WSPD, order_by = date))
wrack_offset5 <- mutate(wrack_by_wind, offset5 = lead(n = 5, WSPD, order_by = date))
wrack_offset6 <- mutate(wrack_by_wind, offset6 = lead(n = 6, WSPD, order_by = date))
wrack_offset7 <- mutate(wrack_by_wind, offset7 = lead(n = 7, WSPD, order_by = date))
total_wrack_offset <- wrack_offset
total_wrack_offset$offset1 <- wrack_offset1$offset1
total_wrack_offset$offset2 <- wrack_offset2$offset2
total_wrack_offset$offset3 <- wrack_offset3$offset3
total_wrack_offset$offset4 <- wrack_offset4$offset4
total_wrack_offset$offset5 <- wrack_offset5$offset5
total_wrack_offset$offset6 <- wrack_offset6$offset6
total_wrack_offset$offset7 <- wrack_offset7$offset7

wrack_fit <- lm(n ~ WSPD + offset1 + offset2 + offset3 + offset4 + offset5 + offset6 + offset7, data = total_wrack_offset)
summary(wrack_fit)

###################################################################################
# VISUALIZATIONS                                                                  #
###################################################################################

# mature date by blade width (same as Dave's w/ extras)
ggplot(mature_sporophyte, aes(x = date, y = blade_width, color = sori_num, size = stipe_length)) +
  geom_point() +
  theme_minimal() +
  annotate("rect", xmin = as.Date("2018-01-01"), xmax = as.Date("2018-06-23"), 
           ymin = -Inf, ymax = Inf, 
           alpha = .1) +
  annotate("rect", xmin = as.Date("2018-12-23"), xmax = as.Date("2019-06-23"), 
           ymin = -Inf, ymax = Inf, 
           alpha = .1) +
  annotate("rect", xmin = as.Date("2019-12-23"), xmax = as.Date("2020-04-27"), 
           ymin = -Inf, ymax = Inf, 
           alpha = .1) +
  scale_colour_viridis() +
  geom_vline(data = gust_mean_15m, aes(xintercept = as.numeric(date)),
             color = "red") +
  labs(color = "Number of sori", size = "Stipe length (cm)") + # give the legend the name 
  xlab("Date") + # label for x axis
  ylab("Blade width (cm)") # label for y axis


# young date by stipe length
ggplot(young_sporophyte, aes(x = date, y = stipe, color = group)) +
  geom_point() +
  scale_colour_viridis(option = "D", discrete = TRUE, begin = 0.8, end = 0.2)

# mature stipe by bladeunite(data, col, ..., sep = "_"
ggplot(mature_sporophyte, aes(x = stipe_length, y = blade_width, color = sori_num, size = bulb_diam)) +
  geom_point() +
  scale_colour_viridis()

# mature stipe blade width by sori number
ggplot(mature_sporophyte, aes(x = as.factor(sori_num), y = blade_width)) +
  geom_boxplot() +
  scale_colour_viridis()

blade_sori_fit <- lm(sori_num ~ blade_width, data = mature_sporophyte)
summary(blade_sori_fit)

# mature stipe length by bulb diameter
ggplot(mature_sporophyte, aes(x = stipe_length, y = bulb_diam)) +
  geom_point() +
  geom_smooth(method='gam', formula = y ~ s(log(x)))


  ############### NEW VIZUALIZTIONS 2022-03-29

# 1) For each of the four years, how many days was the beach transect surveyed for 
#    young sporophytes?

young_sporophyte_q1 <- young_sporophyte %>%
  separate(Date, c("year", "month", "day"), sep = "-") %>%
  unite("month-day", month:day, remove = FALSE) %>% 
  select(year, "month-day") %>%
  group_by(year) %>%
  summarise(days = length(unique(`month-day`)))

ggplot(young_sporophyte_q1, aes(x = year, y = days)) +
  geom_col() +
  theme_minimal() +
  labs(y = "Number of Days", x = "Year") +
  geom_text(aes(label = days), vjust = -0.5)

# 2) For each of the four years how many days were young sporophytes found on the 
#    transect?

young_sporophyte_q2 <- young_sporophyte  %>%
  separate(Date, c("year", "month", "day"), sep = "-") %>%
  unite("month-day", month:day, remove = FALSE) %>% 
  filter(Stipe == "none") %>%
  select(year, "month-day") %>%
  group_by(year) %>%
  summarise(days = length(unique(`month-day`))) %>%
  full_join(young_sporophyte_q1, by = "year") %>%
  mutate(sporophytes_found = days.y - days.x)
  


#####
#<<<<<<<<<<<<<<<<<<<<<<<<<<END OF SCRIPT>>>>>>>>>>>>>>>>>>>>>>>>#