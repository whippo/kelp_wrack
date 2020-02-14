###################################################################################
#                                                                                ##
# kelp_wrack_summary.R                                                           ##
# Data are current as of 2020-02-11                                              ##
# Data source: David Bilderback                                                  ##
# R code prepared by Ross Whippo                                                 ##
# Last updated 2020-02-13                                                        ##
#                                                                                ##
###################################################################################

# SUMMARY:
# Data on kelp wrack morphometrics collected by David Bilderback in Bandon, OR, USA
# from 2018 - present. 


# Required Files (check that script is loading latest version):
# mature_sporophyte_data.csv
# young_sporophyte_data.csv

# Associated Scripts:
# none

# TO DO

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

###################################################################################
# LOAD PACKAGES                                                                   #
###################################################################################

library(tidyverse)
library(readr)
library(viridis)


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
mature_sporophyte$date <- as.Date(mature_sporophyte$date, format = "%m/%d/%Y")
                                 
                                
# read in young dataset
young_sporophyte <- read_csv("Data/young_sporophyte_data.csv", 
                                  col_types = cols(date = col_character(),
                                                   stipe = col_number(),
                                                   blade = col_number(), 
                                                   bulb = col_number(), 
                                                   date = col_date(format = "%m/%d/%Y"), 
                                                   group = col_character(), 
                                                   ho_fa = col_number(), 
                                                   single = col_character()))

# remove erroneous date (need to fix)
young_sporophyte <- young_sporophyte %>%
  filter(date != "5/8/2013") %>%
  filter(date != "8/20/2013")
young_sporophyte$date <- young_sporophyte$date %>%
  recode("11/28/2028" = "11/28/2018",
         "10/9/2028" = "10/9/2018",
         "6/17/2028" = "6/17/2018")
young_sporophyte$date <- as.Date(young_sporophyte$date, format = "%m/%d/%Y")

# import NOAA weather data
noaa_portorford_buoydata_2018 <- read_table2("Data/noaa_portorford_buoydata_2018.csv", 
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
# VISUALIZATIONS                                                                  #
###################################################################################

# mature date by blade width (same as Dave's w/ extras)
ggplot(mature_sporophyte, aes(x = date, y = blade_width, color = sori_num, size = stipe_length)) +
  geom_point() +
  scale_colour_viridis() +
  geom_vline(data = gust_mean_15m, aes(xintercept = as.numeric(date)),
             color = "red")


# young date by stipe length
ggplot(young_sporophyte, aes(x = date, y = stipe, color = group)) +
  geom_point() +
  scale_colour_viridis(option = "D", discrete = TRUE, begin = 0.8, end = 0.2)

# mature stipe by bladeunite(data, col, ..., sep = "_"
ggplot(mature_sporophyte, aes(x = stipe_length, y = blade_width, color = sori_num, size = bulb_diam)) +
  geom_point() +
  scale_colour_viridis()

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



############### SUBSECTION HERE

#####
#<<<<<<<<<<<<<<<<<<<<<<<<<<END OF SCRIPT>>>>>>>>>>>>>>>>>>>>>>>>#