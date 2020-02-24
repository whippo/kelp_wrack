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
library(lme4)


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

# fix bulb diameter mistake
mature_sporophyte <- mature_sporophyte %>% 
  mutate(bulb_diam = replace(bulb_diam, bulb_diam == 72, 7.2))
                                
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

# reduce number of days to match weather data
wrack_obs_reduced <- wrack_obs %>%
  filter(date %in% daily_wind_means_char$date)

# make date a character vector in new winds dataset for joining
daily_wind_means_char <- daily_wind_means

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




  ############### SUBSECTION HERE

#####
#<<<<<<<<<<<<<<<<<<<<<<<<<<END OF SCRIPT>>>>>>>>>>>>>>>>>>>>>>>>#