
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                                                             ##
# Kelp Wrack Analyses                                                         ##
# Script created 2020-02-13                                                   ##
# Data source: David Bilderback                                               ##
# R code prepared by Ross Whippo                                              ##
# Last updated 2024-01-31                                                     ##
#                                                                             ##
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# SUMMARY:
# Data on kelp wrack morphometrics collected by David Bilderback in Bandon, OR, USA
# from 2018 - 2022. Script modified from kelp_wrack_summary.R


# Required Files (check that script is loading latest version):
# mature_sporophyte_data.csv
# young_sporophyte_data.csv
# noaa_portorford_buoydata_2018.csv

# Associated Scripts:
# kelp_wrack_summary.R

# TO DO

# PER KATRIN: 
#     - is there a correlation between substrate type and reproductive state?
#     - is there a correlation between wind direction/speed and wrack?
#     - are there unexpected pulses of recruitment?

# PER KAIN 1987: stipe growth in <25cm ~2mm/day, 25-100cm ~5mm/day
# PER DOBKOWSKI: juv <30-40, adult >200
# PER NICHOLSON: stipe groth in juv ~2.5cm/day


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# TABLE OF CONTENTS                                                         ####
#                                                                              +
# RECENT CHANGES TO SCRIPT                                                     +
# LOAD PACKAGES                                                                +
# READ IN AND PREPARE DATA                                                     +
# MANIPULATE DATA                                                              +
#                                                                              +
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# RECENT CHANGES TO SCRIPT                                                  ####
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# 2024-01-31 Script modified from previous wrack analysis

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# LOAD PACKAGES                                                             ####
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(tidyverse)
library(viridis)
library(lme4)
library(lubridate)
library(ggpubr)

options(max.print = 9999)


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# READ IN AND PREPARE DATA                                                  ####
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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

# remove year so months overlap

mature_sporophyte_monthly <- mature_sporophyte %>%
  mutate(Date = as.Date.POSIXct(date, "%Y-%M-%D")) %>%
  mutate(year = as.character(year(Date))) %>%
  mutate(monthday = format(as.Date.POSIXct(Date, "%Y-%M-%D"), "2020-%m-%d")) %>%
  mutate(monthday = as.Date.POSIXct(monthday, "%Y-%M-%D", tz = "PT")) %>%
  arrange(monthday)

mature_sporophyte_cohort <- mature_sporophyte_monthly %>%
  mutate(cohort = case_when(
    stipe_length <= 150 ~ "35-150",
    stipe_length <= 1000 ~ "151-1000",
    stipe_length > 1000 ~ ">1000")
  ) %>%
  filter(!is.na(cohort))

                                
# read in young dataset
young_sporophyte <- read_csv("Data/young_sporophyte_data.csv", 
                                  col_types = cols(Date = col_character(), 
                                                   Stipe = col_number(),
                                                   Bulb = col_number(), 
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
# fix bulb diameter and blade width typos
young_sporophyte['Bulb'][young_sporophyte['Bulb'] == 53.0] <- 5.3
young_sporophyte['Bulb'][young_sporophyte['Bulb'] == 30.0] <- 3.0
young_sporophyte['Blade'][young_sporophyte['Blade'] == 22.0] <- 2.2
# fix substrate typos
young_sporophyte$Subst <- young_sporophyte$Subst %>%
  recode("NL" = "Nl",
         "Ni" = "Nl",
         "MY" = "My",
         "R" = "Ro",
         "B" = "Ba",
         "Bs" = "Ba",
         "Bo" = "Ba",
         "Bl" = "Nl",
         "LS" = "Ls",
         "SS" = "Ss",
         "Hy" = "Hf",
         "Mp" = "Mu",
         "Dl" = "Do")


# fix CoSp1 typos
young_sporophyte$CoSp1 <- young_sporophyte$CoSp1 %>%
  recode("LS" = "Ls",
         "Po" = "Pl",
         "Cor" = "Co",
         "Io" = "Is",
         "Tu" = "none",
         "ls" = "Ls",
         "Ba" = "Bo",
         "By" = "Bo")
         
# remove year so months overlap

young_sporophyte_monthly <- young_sporophyte %>%
  mutate(Date = as.Date.POSIXct(Date, "%Y-%M-%D")) %>%
  mutate(year = as.character(year(Date))) %>%
  mutate(monthday = format(as.Date.POSIXct(Date, "%Y-%M-%D"), 
                           "2020-%m-%d")) %>%
  mutate(monthday = as.Date.POSIXct(monthday, "%Y-%M-%D", 
                                    tz = "PT")) %>%
  arrange(monthday)

young_sporophyte_cohort <- young_sporophyte_monthly %>%
  mutate(cohort = case_when(
    Stipe <= 15 ~ "0.1-15",
    Stipe <= 20 ~ "15.1-20",
    Stipe <= 100 ~ "20.1-100")
  ) %>%
  filter(!is.na(cohort))

# combine mature and young datasets
young_step1 <- young_sporophyte %>%
  select(Date, Stipe) %>%
  filter(Stipe < 160)

mature_step1 <- mature_sporophyte %>%
  select(date, stipe_length) %>%
  filter(stipe_length >= 160) %>%
  rename(Date = date,
         Stipe = stipe_length)

# 160cm chosen as cutoff for datasets, overlapping values
full_sporophyte <- young_step1 %>%
  bind_rows(mature_step1) %>%
  filter(!is.na(Stipe))

all_sporophyte_monthly <- full_sporophyte %>%
  mutate(Date = as.Date.POSIXct(Date, "%Y-%M-%D")) %>%
  mutate(year = as.character(year(Date))) %>%
  mutate(monthday = format(as.Date.POSIXct(Date, "%Y-%M-%D"), 
                           "2020-%m-%d")) %>%
  mutate(monthday = as.Date.POSIXct(monthday, "%Y-%M-%D", 
                                    tz = "PT")) %>%
  arrange(monthday)

all_sporophyte_cohort <- all_sporophyte_monthly %>%
  mutate(cohort = case_when(
    Stipe <= 40 ~ "0.1-40",
    Stipe < 200 ~ "40.1-199.9",
    Stipe >= 200 ~ "200+")
  ) %>%
  filter(!is.na(cohort)) %>%
  mutate(cohort = factor(cohort, levels = c("0.1-40",
                                            "40.1-199.9",
                                            "200+")))

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# RECRUITMENT FIGURES                                                       ####
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

all_sporophyte_cohort %>%
  ggplot() +
  geom_point(aes(x = monthday, y = log10(Stipe+1)/185, col = cohort), alpha = 0.4) +
  scale_y_continuous(
    name = "Log + 1 Stipe Length", 
    lim = c(0,0.018), 
    n.breaks = 6, 
    labels = c("0", "1", "2", "3", "4"), 
    sec.axis = sec_axis(trans = ~.*1, "Stipe Count Density Distribution")) +
  scale_colour_viridis(option = "D", 
                       discrete = TRUE, 
                       begin = 0.9, 
                       end = 0.2) +
  scale_x_date(name = "Month", 
               limits = as.Date(c('2020-01-01', '2020-12-31'), 
                                format="%Y-%M-%D"), 
               date_breaks = "1 month",
               date_labels = ("%b")) +
  geom_density(aes(monthday, group = cohort), 
               inherit.aes = FALSE, size = 1.5, col = "black") +
  geom_density(aes(monthday, color = cohort), 
               inherit.aes = FALSE, size = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = -0.9), 
        axis.title.y.right = element_text(vjust = 2)) +
  annotate(geom = "text", x = as.Date('2020-11-20'), y = 0.018, label = "n = 9743")

# by year
all_sporophyte_cohort %>%
  filter(year %in% c("2018",
                     "2019",
                     "2020",
                     "2021")) %>%
  ggplot() +
  geom_point(aes(x = monthday, y = log10(Stipe+1)/11, col = cohort), alpha = 0.4) +
  scale_y_continuous(
    name = "Log + 1 Stipe Length", 
    #lim = c(0,0.3), 
   # n.breaks = 6, 
  #  labels = c("0", "1", "2", "3", "4"), 
    sec.axis = sec_axis(trans = ~.*1, "Stipe Count Density Distribution")) +
  scale_colour_viridis(option = "D", 
                       discrete = TRUE, 
                       begin = 0.9, 
                       end = 0.2) +
  scale_x_date(name = "Month", 
               limits = as.Date(c('2020-01-01', '2020-12-31'), 
                                format="%Y-%M-%D"), 
               date_breaks = "1 month",
               date_labels = ("%b")) +
  geom_density(aes(monthday, group = cohort), 
               inherit.aes = FALSE, size = 1.5, col = "black") +
  geom_density(aes(monthday, color = cohort), 
               inherit.aes = FALSE, size = 1) +
  theme_bw() +
  theme(axis.title.y.right = element_text(vjust = 2)) +
  facet_wrap(.~year, ncol = 1, scales = "free_y") 


# small only
all_sporophyte_cohort %>%
  filter(cohort == "0.1-40") %>%
  ggplot() +
  geom_point(aes(x = monthday, y = Stipe/2700, col = cohort), alpha = 0.4) +
  scale_y_continuous(
    name = "Stipe Length (cm)", 
   # lim = c(0,0.015), 
   #  n.breaks = 6, 
   #   labels = c("0", "10", "20", "30", "40"), 
    sec.axis = sec_axis(trans = ~.*1, "Stipe Count Density Distribution")) +
  scale_colour_viridis(option = "D", 
                       discrete = TRUE, 
                       begin = 0.9, 
                       end = 0.2) +
  scale_x_date(name = "Month", 
               limits = as.Date(c('2020-01-01', '2020-12-31'), 
                                format="%Y-%M-%D"), 
               date_breaks = "1 month",
               date_labels = ("%b")) +
  geom_density(aes(monthday), 
               inherit.aes = FALSE, size = 1.5, col = "black") +
  geom_density(aes(monthday, col = cohort), 
               inherit.aes = FALSE, size = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = -0.9), 
        axis.title.y.right = element_text(vjust = 2)) +
  annotate(geom = "text", x = as.Date('2020-11-20'), y = 0.015, label = "n = 5211") 

# juvenile only
all_sporophyte_cohort %>%
  filter(cohort == "0.1-40") %>%
  ggplot() +
  scale_colour_viridis(option = "D", 
                       discrete = TRUE, 
                       begin = 0.9, 
                       end = 0.2) +
  scale_x_date(name = "Month", 
               limits = as.Date(c('2020-01-01', '2020-12-31'), 
                                format="%Y-%M-%D"), 
               date_breaks = "1 month",
               date_labels = ("%b")) +
  stat_density(aes(x = monthday), 
              size = 1.5, col = "black") +
  stat_density(aes(x = monthday, col = cohort), 
              size = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = -0.9), 
        axis.title.y.right = element_text(vjust = 2)) +
  facet_wrap(.~year, ncol = 1, scales = "free_y")

# intermediate only
all_sporophyte_cohort %>%
  filter(cohort == "40.1-199.9") %>%
  ggplot() +
  scale_colour_viridis(option = "D", 
                       discrete = TRUE, 
                       begin = 0.55, 
                       end = 0.2) +
  scale_x_date(name = "Month", 
               limits = as.Date(c('2020-01-01', '2020-12-31'), 
                                format="%Y-%M-%D"), 
               date_breaks = "1 month",
               date_labels = ("%b")) +
  stat_density(aes(x = monthday), 
               size = 1.5, col = "black") +
  stat_density(aes(x = monthday, col = cohort), 
               size = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = -0.9), 
        axis.title.y.right = element_text(vjust = 2)) +
  facet_wrap(.~year, ncol = 1, scales = "free_y")

# mature only
all_sporophyte_cohort %>%
  filter(cohort == "200+") %>%
  ggplot() +
  scale_colour_viridis(option = "D", 
                       discrete = TRUE, 
                       begin = 0.25, 
                       end = 0.2) +
  scale_x_date(name = "Month", 
               limits = as.Date(c('2020-01-01', '2020-12-31'), 
                                format="%Y-%M-%D"), 
               date_breaks = "1 month",
               date_labels = ("%b")) +
  stat_density(aes(x = monthday), 
               size = 1.5, col = "black") +
  stat_density(aes(x = monthday, col = cohort), 
               size = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = -0.9), 
        axis.title.y.right = element_text(vjust = 2)) +
  facet_wrap(.~year, ncol = 1, scales = "free_y")


# sori number
sori_presence <- mature_sporophyte_monthly %>%
  filter(!is.na(sori_num)) %>%
  mutate(sori_present = sori_num + re_sori_num + absd_sori_num) %>%
  mutate(sori_absent = case_when(sori_present == 0 ~ 1,
                                 .default = 0)) %>%
  mutate(sori_present = case_when(sori_absent == 0 ~ 1,
                                  .default = 0)) %>%
  pivot_longer(cols = sori_absent:sori_present, names_to = "sori_state",
               values_to = "value") %>%
  mutate(sori_state = factor(sori_state, levels = c("sori_absent",
                                                    "sori_present"))) %>%
  mutate(month = factor(format(Date,"%b"), levels = c("Jan",
                                                      "Feb",
                                                      "Mar",
                                                      "Apr",
                                                      "May",
                                                      "Jun",
                                                      "Jul",
                                                      "Aug",
                                                      "Sep",
                                                      "Oct",
                                                      "Nov",
                                                      "Dec"))) %>%
  group_by(year, month, sori_state) %>%
  summarise(monthly = sum(value))


ggplot() +
  geom_col(data = sori_presence %>%
             group_by(year, month) %>%
             summarise(monthly = sum(monthly)) %>%
             group_by(month) %>%
             summarise(monthly = mean(monthly)),
           aes(x = month, y = monthly), alpha = 0.2) +
  geom_boxplot(data = sori_presence, aes(x = month, y = monthly, fill = sori_state)) +
  scale_fill_viridis(name = "Reproductive State", 
                     labels = c("No Sori", "Sori Present/Abscissed"),
                     option = "A", 
                     discrete = TRUE, 
                     begin = 0.9, 
                     end = 0.4) +
  theme_bw() +
  labs(x = "Month", y = "Count") 


mature_sporophyte_monthly %>%
  filter(!is.na(sori_num)) %>%
  mutate(sori_only = case_when(sori_num > 0 ~ 1,
                               .default = 0)) %>%
  mutate(absd_only = case_when(absd_sori_num > 0 ~ 1,
                               .default = 0)) %>%
  mutate(soriabsd = sori_only + absd_only) %>%
  mutate(soriabsd = soriabsd - 1) %>%
  mutate(soriabsd = case_when(soriabsd == -1 ~ 0,
                              .default = soriabsd)) %>%
  mutate(sori_only = case_when(soriabsd == 1 ~ 0,
                               .default = sori_only)) %>%
  mutate(absd_only = case_when(soriabsd == 1 ~ 0,
                               .default = absd_only)) %>%
  mutate(sum_sori = sori_only + absd_only + soriabsd) %>%
  mutate(no_sori = case_when(sum_sori > 0 ~ 0,
                             .default = 1)) %>%
  select(date, Date, year, monthday, sori_only, absd_only, soriabsd, no_sori) %>%
  pivot_longer(cols = sori_only:no_sori, names_to = "sori_state",
               values_to = "value") %>%
  mutate(sori_state = factor(sori_state, levels = c("no_sori",
                             "sori_only",
                             "absd_only",
                             "soriabsd"))) %>%
  mutate(month = factor(format(Date,"%b"), levels = c("Jan",
                                                         "Feb",
                                                         "Mar",
                                                         "Apr",
                                                         "May",
                                                         "Jun",
                                                         "Jul",
                                                         "Aug",
                                                         "Sep",
                                                         "Oct",
                                                         "Nov",
                                                         "Dec"))) %>%
  group_by(year, month, sori_state) %>%
  summarise(monthly = sum(value)) %>%
  ggplot() +
  geom_boxplot(aes(x = month, y = monthly, fill = sori_state)) +
  scale_fill_viridis(option = "A", 
                       discrete = TRUE, 
                       begin = 0.9, 
                       end = 0.2) +
  theme_bw() 
 

  




## THIS FINALLY WORKED! STIPE LENGTH DOTPLOT AND DENSITY DISTRIBUTION IN ONE FIGURE!

# color coded by length

young_sporophyte_cohort %>%
  filter(is.na(Stipe) | Stipe <= 100) %>%
  ggplot() +
  geom_point(aes(x = monthday, y = Stipe/5600, col = cohort), alpha = 0.4) +
  scale_y_continuous(
    name = "Stipe Length (cm)", lim = c(0,0.018), n.breaks = 6, labels = c("0", "25", "50", "75", "100"), 
    sec.axis = sec_axis(trans = ~.*1, "Stipe Count Density Distribution")) +
  scale_colour_viridis(option = "A", discrete = TRUE, begin = 0.8, end = 0.2, breaks = c("0.1-15", "15.1-20", "20.1-100")) +
  scale_x_date(name = "Month", limits = as.Date(c('2020-05-01', '2020-10-01'), format="%Y-%M-%D"), 
               date_breaks = "1 month",
               date_labels = ("%b")) +
  geom_density(aes(monthday, group = cohort), inherit.aes = FALSE, size = 1.5, col = "black") +
  geom_density(aes(monthday, color = cohort), inherit.aes = FALSE, size = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = -0.9), 
        axis.title.y.right = element_text(vjust = 2)) +
  annotate(geom = "text", x = as.Date('2020-11-15'), y = 0.015, label = "n = 4228")


## THIS FINALLY WORKED! STIPE LENGTH DOTPLOT AND DENSITY DISTRIBUTION IN ONE FIGURE!

# mature sporophyte

mature_sporophyte_cohort %>%
  filter(!is.na(stipe_length) | stipe_length > 35) %>%
  ggplot() +
  geom_point(aes(x = monthday, 
                 y = stipe_length/200000, 
                 col = cohort), alpha = 0.4) +
  scale_y_continuous(
    name = "Stipe Length (cm)", 
    lim = c(0,0.01), n.breaks = 3, 
    labels = c("0", "1000", "2000"), 
    sec.axis = sec_axis(trans = ~.*1, "Stipe Count Density Distribution")) +
  scale_colour_viridis(option = "A", 
                       discrete = TRUE, 
                       begin = 0.8, 
                       end = 0.2, 
                       breaks = c("35-150", "151-1000", ">1000")) +
  scale_x_date(name = "Month", 
               date_breaks = "1 month",
               date_labels = ("%b")) +
  geom_density(aes(monthday, group = cohort), 
               inherit.aes = FALSE, size = 1.5, col = "black") +
  geom_density(aes(monthday, color = cohort), 
               inherit.aes = FALSE, size = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = -0.9), 
        axis.title.y.right = element_text(vjust = 2)) +
  annotate(geom = "text", 
           x = as.Date('2020-11-15'), 
           y = 0.0095, 
           label = "n = 3351")

####
#<<<<<<<<<<<<<<<<<<<<<<<<<<END OF SCRIPT>>>>>>>>>>>>>>>>>>>>>>>>#

# SCRATCH PAD ####
     
# subset and write single and grouped individuals
clustered <- young_sporophyte %>%
  filter(Single == 0) 
# write_csv(clustered, "Data/clustered_young_sporophytes.csv")
singles <- young_sporophyte %>%
  filter(Single == 1)
# write_csv(singles, "Data/single_young_sporophytes.csv")    



# substrate type dataset
substrateData <- read_csv("Data/substrateCodes.csv")

# split columns add none value
substrate <- substrateData %>%
  separate(`Subst=Substrate`, into = c("Subst", "Substrate"), sep = "=") %>%
  bind_rows(c(Subst = "none", Substrate = "none"))



# cospecies dataset
cospeciesData <- read_csv("Data/cospeciesCodes.csv")

# split columns add none and duplicate columns
cospecies <- cospeciesData %>%
  separate(`Cosp = Cospecies`, into = c("Cosp", "Cospecies"), sep =  " = ") %>%
  bind_rows(c(Cosp = "none", Cospecies = "none")) %>%
  mutate(CoSp1 = Cosp,
         CoSp2 = Cosp,
         CoSp3 = Cosp)



# import NOAA weather data
noaa_portorford_buoydata_2018 <- read_tsv("Data/noaa_portorford_buoydata_2018.csv", 
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

SST_2018 <- read_csv("Data/NOAA_SST_2018_PortOrford_PORO3_9431647.csv", 
                     col_types = cols(WSPD = col_number()), 
                     skip = 1)
SST_2018 <- read_csv("Data/NOAA_SST_2019_PortOrford_PORO3_9431647.csv", 
                     col_types = cols(WSPD = col_number()), 
                     skip = 1)
SST_2018 <- read_csv("Data/NOAA_SST_2020_PortOrford_PORO3_9431647.csv", 
                     col_types = cols(WSPD = col_number()), 
                     skip = 1)
SST_2018 <- read_csv("Data/NOAA_SST_2021_PortOrford_PORO3_9431647.csv", 
                     col_types = cols(WSPD = col_number()), 
                     skip = 1)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# NOAA WEATHER DATA                                                         ####
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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





#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# VISUALIZATIONS                                                            ####
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


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
ggplot(young_sporophyte, aes(x = Date, y = Stipe, color = Group)) +
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

  
# 3) Is there a relationship between stipe length and bulb diameter?
#

# NEED TO TEST BETTER MODEL -> LOG? SATURATING?


fit1 <- lm(Bulb ~ Stipe, data = young_sporophyte)
fit2 <- lm(Bulb~poly(Stipe,2,raw=TRUE), data=young_sporophyte)
fit3 <- lm(Bulb~poly(Stipe,3,raw=TRUE), data=young_sporophyte)
fit4 <- lm(Bulb~poly(Stipe,4,raw=TRUE), data=young_sporophyte)
fit5 <- lm(Bulb~poly(Stipe,5,raw=TRUE), data=young_sporophyte)

summary(fit1)$adj.r.squared
summary(fit2)$adj.r.squared
summary(fit3)$adj.r.squared
summary(fit4)$adj.r.squared
summary(fit5)$adj.r.squared

AIC(fit1, fit2, fit3, fit4, fit5)

summary(fit1)

ggplot(young_sporophyte, aes(x = Stipe, y = Bulb, na.rm = TRUE)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal() +
  labs(y = "Bulb diameter (cm)", x = "Stipe length (cm)", title = "Relationship between stipe length and bulb diameter (young sporophytes)") + 
  annotate("text", x = 150, y=9, label = "Adj. R-squared = 0.62, p < 0.0001")


# 4) Is there a relationship between bulb diameter and widest blade width?
#

fit1 <- lm(Bulb ~ Blade, data = young_sporophyte)
fit2 <- lm(Bulb~poly(Blade,2,raw=TRUE), data=young_sporophyte)
fit3 <- lm(Bulb~poly(Blade,3,raw=TRUE), data=young_sporophyte)
fit4 <- lm(Bulb~poly(Blade,4,raw=TRUE), data=young_sporophyte)
fit5 <- lm(Bulb~poly(Blade,5,raw=TRUE), data=young_sporophyte)

summary(fit1)$adj.r.squared
summary(fit2)$adj.r.squared
summary(fit3)$adj.r.squared
summary(fit4)$adj.r.squared
summary(fit5)$adj.r.squared

AIC(fit1, fit2, fit3, fit4, fit5)

summary(fit1)

ggplot(young_sporophyte, aes(x = Blade, y = Bulb, na.rm = TRUE)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal() +
  labs(y = "Bulb diameter (cm)", x = "Blade width (cm)", title = "Relationship between blade width and bulb diameter (young sporophytes)") + 
  annotate("text", x = 6, y=12, label = "Adj. R-squared = 0.16, p < 0.0001")


# 5) What is the mean diameter of holdfasts?
#



mean(young_sporophyte$`Ho Fa`, na.rm = TRUE)
range(young_sporophyte$`Ho Fa`, na.rm = TRUE)
median(young_sporophyte$`Ho Fa`, na.rm = TRUE)


ggplot(young_sporophyte, aes(y = `Ho Fa`)) +
  geom_boxplot() +
  theme_minimal() +
  labs(y = "Holdfast diameter (cm)", title = "Range of holdfast diameters") +
  annotate("text", x = 0.2, y = 11, label = "mean = 3.00; median = 2.60; range = 0.3 - 12.00")

# 6) For each year, how many sporophytes were singular, and how many were clusters?
#

# NEED TO ADD ACTUAL VALUES TO PLOT

young_sporophyte_q3 <- young_sporophyte %>%
  separate(Date, c("year", "month", "day"), sep = "-") %>%
  unite("month-day", month:day, remove = FALSE) %>% 
  filter(Single %in% c("1", "0")) %>%
  mutate(Single = case_when(
      Single == "1" ~ "single",
      Single == "0" ~ "multiple")) %>%
  group_by(year)

young_sporophyte_q3 %>%
  group_by(year, Single) %>%
  summarise(length(year))

ggplot(young_sporophyte_q3, aes(x = year, label = Single)) +
  geom_bar(aes(fill = Single)) +
  theme_minimal() +
  scale_fill_viridis(discrete = TRUE, option = "D", begin = 0.3, end = 0.7) +
  labs(fill = "Sporophyte grouping", x = "Year", y = "Number of individuals")  
  

# 7) What is the count and percentage of each substrate type?
#

# THIS COULD BE A GREAT FIGURE IF ALL OTHER DATA ADDED IN - INFOGRAPHIC STYLE

young_sporophyte_q4 <- young_sporophyte %>%
  mutate(Subst = na_if(Subst, "nd")) %>%
  mutate(Subst = na_if(Subst, "N A")) %>%
  mutate(Subst = na_if(Subst, "Na")) %>%
  filter(!is.na(Subst)) %>%
  left_join(substrate, by = "Subst")

substrate_totals <- young_sporophyte_q4 %>%
  group_by(Substrate) %>%
  summarise('total' = length(Substrate)) %>%
  arrange(desc(total)) %>%
  mutate(percent = (total/sum(total))*100)

ggplot(young_sporophyte_q4 %>%
         group_by(Substrate) %>%
         summarise(count = length(Substrate)) %>%
         filter(count > 20) %>%
         arrange(desc(count)), aes(fct_reorder(Substrate, count), count)) +
  geom_col() +
  theme_minimal() +
  labs(x = "Substrate type", "Number of occurrences") +
  coord_flip()

# 8) How common were cospecies? Numbers and percentages.
#

# primary

CoSp1_data <- cospecies %>%
  select(CoSp1, Cospecies)

young_sporophyte_q8 <- young_sporophyte %>%
  mutate(CoSp1 = na_if(CoSp1, "nd")) %>%
  mutate(CoSp1 = na_if(CoSp1, "N A")) %>%
  mutate(CoSp1 = na_if(CoSp1, "Na")) %>%
  filter(!is.na(CoSp1)) %>%
  left_join(CoSp1_data, by = "CoSp1")

CoSp1_totals <- young_sporophyte_q8 %>%
  group_by(Cospecies) %>%
  summarise('total' = length(Cospecies)) %>%
  arrange(desc(total)) %>%
  mutate(percent = (total/sum(total))*100)

ggplot(young_sporophyte_q8 %>%
         group_by(Cospecies) %>%
         summarise(count = length(Cospecies)) %>%
         filter(count > 100) %>%
         arrange(desc(count)), aes(fct_reorder(Cospecies, count), count)) +
  geom_col() +
  theme_minimal() +
  labs(x = "Primary cospecies", "Number of occurrences") +
  coord_flip()

# percent with multiples

young_sporophyte_q8_1 <- young_sporophyte %>%
  mutate(CoSp1 = na_if(CoSp1, "nd")) %>%
  mutate(CoSp1 = na_if(CoSp1, "N A")) %>%
  mutate(CoSp1 = na_if(CoSp1, "Na")) %>%
  mutate(CoSp2 = na_if(CoSp2, "nd")) %>%
  mutate(CoSp2 = na_if(CoSp2, "N A")) %>%
  mutate(CoSp2 = na_if(CoSp2, "Na")) %>%
  mutate(CoSp3 = na_if(CoSp3, "nd")) %>%
  mutate(CoSp3 = na_if(CoSp3, "N A")) %>%
  mutate(CoSp3 = na_if(CoSp3, "Na")) %>%
  mutate(CoSp1 = ifelse(CoSp1 == "none", "0", "1")) %>%
  mutate(CoSp2 = ifelse(CoSp2 == "none", "0", "1")) %>%  
  mutate(CoSp3 = ifelse(CoSp3 == "none", "0", "1")) %>%
  mutate_at(c("CoSp1", "CoSp2", "CoSp3"), as.numeric) %>%
  mutate(Total_cospecies = CoSp1 + CoSp2 + CoSp3) %>%
  filter(!is.na(Total_cospecies))

cospecies_totals <- young_sporophyte_q8_1 %>%
  group_by(Total_cospecies) %>%
  summarise('total' = length(Total_cospecies)) %>%
  arrange(desc(total)) %>%
  mutate(percent = (total/sum(total))*100)

ggplot(cospecies_totals, aes(x = Total_cospecies, y = total)) +
  geom_col() +
  theme_minimal() +
  labs(x = "Total number of cospecies", y = "Frequency of occurrence")  +
  geom_text(aes(label = total), vjust = -0.5)

# 9) For each year, when does recruitment occur? (as function of % of 0.1-10cm
#    length stipes throughout the year?)

# HOW LONG TILL STIPES REACH 10CM FROM SETTLEMENT?

# THIS IS SUPER SEXY AND POTENTIAL THE CRUX OF THIS DATASET

young_sporophyte_q9 <- young_sporophyte %>%
  complete(Date = seq.Date(as.Date("2018-01-01"), as.Date("2021-12-31"), by="day")) %>%
  separate(Date, c("year", "month", "day"), sep = "-", remove = FALSE) %>%
  unite("month-day", month:day, remove = FALSE) %>% 
  mutate(Stipe = ifelse(is.na(Stipe), 0.5, Stipe)) %>%
  filter(Stipe < 10.1)

ggplot(young_sporophyte_q9, aes(x = Date, y = Stipe)) +
  geom_col() +
  theme_minimal() +
  labs(x = "Date", y = "Count of stipes < 10.1 cm") +
  scale_x_date(breaks = "1 month", labels = date_format("%b")) +
  facet_wrap(~year, nrow = 4, scales = "free_x")


#####
#<<<<<<<<<<<<<<<<<<<<<<<<<<END OF SCRIPT>>>>>>>>>>>>>>>>>>>>>>>>#

##### SCRATCH PAD

# FIGURE for determining seasonality of recruitment events 

# young date by stipe length as timeline
ggplot(young_sporophyte, aes(x = Date, y = Stipe, color = Group)) +
  geom_point() +
  scale_colour_viridis(option = "D", discrete = TRUE, begin = 0.8, end = 0.2)





ggplot(young_sporophyte_monthly) +
  geom_point(aes(x = monthday, y = Stipe, col = year)) +
  geom_point(aes(x = monthday), stat = "count") +
  geom_smooth(aes(x = monthday, group = 1), stat = "count", color = alpha("black", 0.3)) +
  scale_y_continuous(
    name = "Stipe Length (cm)", lim = c(0,350),
    sec.axis = sec_axis(trans = ~.*1, "Stipe Count")) +
  scale_colour_viridis(option = "D", discrete = TRUE, begin = 0.8, end = 0.2) +
  scale_x_date(name = "Month", breaks = "1 month",
               date_labels = ("%b")) +
  theme(axis.text.x = element_text(hjust = -0.2)) +
  facet_wrap(.~year, ncol = 1)



# how many observations on each day? All singletons removed
datecounts <- as_tibble(table(young_sporophyte_monthly['Date'])) %>%
  mutate(monthday = format(as.Date.POSIXct(Date, "%Y-%M-%D"), "2020-%m-%d")) %>%
  mutate(n = n - 1) %>%
  mutate(fakedates = as.Date.POSIXct(monthday, "%Y-%M-%D"))
  

ggplot(young_sporophyte_monthly) +
  geom_point(aes(x = monthday, y = Stipe, col = year)) +
  geom_point(aes(x = monthday), stat = "count") +
  geom_smooth(aes(x = monthday, group = 1), stat = "count", color = alpha("black", 0.3)) +
  scale_y_continuous(
    name = "Stipe Length (cm)", lim = c(0,350),
    sec.axis = sec_axis(trans = ~.*1, "Stipe Count")) +
  scale_colour_viridis(option = "D", discrete = TRUE, begin = 0.8, end = 0.2) +
  scale_x_date(name = "Month", limits = as.Date(c('2020-05-01', '2020-10-01'), format="%Y-%M-%D"), 
               date_breaks = "3 month",
               date_labels = ("%b")) +
  theme(axis.text.x = element_text(hjust = -0.9))



datecounts %>%
  mutate(numberdates = as.numeric(fakedates)) %>%
  summarise(min(numberdates)) # 90 days from jan 1 - mar 31

datecounts %>%
  filter(n != 0) %>%
  mutate(logcount = log10(n)) %>%
  group_by(fakedates) %>%
  summarise(meanlogcount = mean(logcount)) 
datecounts %>%
  filter(n > 1) %>%
ggplot(aes(x = as.numeric(fakedates) - 18352, y = n)) +
  geom_point() +
  stat_smooth(method = "gam", formula = y ~ -1 + I(x^3) + x) 
  #scale_x_date(name = "Month", limits = as.Date(c('2020-05-01', '2020-10-01'), format="%Y-%M-%D"), 
   #            date_breaks = "1 month",
    #           date_labels = ("%b")) +
 # theme(axis.text.x = element_text(hjust = -0.9))

datecounts %>%
  filter(n > 0) %>%
  ggplot(aes(x = as.numeric(fakedates) - 18352, y = n)) +
  geom_point() +
  xlim(0,125) +
  ylim(0,350) +
  stat_smooth(method = "gam", formula = y ~ -1 + I(x^2) + x) 

datecounts %>%
  filter(n > 5) %>%
  ggplot(aes(x = as.numeric(fakedates) - 18352, y = n)) +
  geom_point() +
  xlim(125,260) +
  ylim(0,350) +
  stat_smooth(method = "gam", formula = y ~ -1 + I(x^2) + x) 



ggplot(young_sporophyte_monthly, aes(x = monthday, y = Stipe, color = Group)) +
  geom_count() +
  scale_colour_viridis(option = "D", discrete = TRUE, begin = 0.8, end = 0.2)

ggplot(young_sporophyte_monthly, aes(x = monthday, y = Stipe, color = Group)) +
  geom_line() +
  scale_colour_viridis(option = "D", discrete = TRUE, begin = 0.8, end = 0.2)



# Young sporophyte counts by size class 

ys15 <- young_sporophyte_monthly %>%
  filter(0.1 <= Stipe & Stipe < 15.1)

ys20 <- young_sporophyte_monthly %>%
  filter(0.1 <= Stipe & Stipe < 20.1)

ys25 <- young_sporophyte_monthly %>%
  filter(0.1 <= Stipe & Stipe < 25.1)

ggplot(ys15, aes(x = monthday, y = Stipe)) +
  geom_point() +
  scale_colour_viridis(option = "D", discrete = TRUE, begin = 0.8, end = 0.2) +
  scale_x_date(name = "Month", breaks = "1 month",
               date_labels = ("%b")) +
  theme(axis.text.x = element_text(hjust = -0.2)) +
  facet_wrap(.~year, ncol = 1)

ggplot(ys20, aes(x = monthday, y = Stipe)) +
  geom_point() +
  scale_colour_viridis(option = "D", discrete = TRUE, begin = 0.8, end = 0.2) +
  scale_x_date(name = "Month", breaks = "1 month",
               date_labels = ("%b")) +
  theme(axis.text.x = element_text(hjust = -0.2)) +
  facet_wrap(.~year, ncol = 1)

ggplot(ys25, aes(x = monthday, y = Stipe)) +
  geom_point() +
  scale_colour_viridis(option = "D", discrete = TRUE, begin = 0.8, end = 0.2) +
  scale_x_date(name = "Month", breaks = "1 month",
               date_labels = ("%b")) +
  theme(axis.text.x = element_text(hjust = -0.2)) +
  facet_wrap(.~year, ncol = 1)






# stipe length and density by month
all_sporophyte_cohort %>%
  ggplot() +
  geom_point(aes(x = monthday, y = log(Stipe + 1)/475, col = cohort), alpha = 0.4) +
  scale_y_continuous(
    name = "Log stipe length", lim = c(0,0.017), n.breaks = 10, labels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "10"), 
    sec.axis = sec_axis(trans = ~.*1, "Stipe Count Density Distribution")) +
  scale_colour_viridis(option = "D", discrete = TRUE, begin = 0.9, end = 0.1, breaks = c("0.1-15", "15.1-20", "20.1-25", "25+")) +
  scale_x_date(name = "Month", 
               date_breaks = "1 month",
               date_labels = ("%b")) +
  geom_density(aes(monthday, group = cohort), inherit.aes = FALSE, size = 1.8, col = "white") +
  geom_density(aes(monthday, group = cohort), inherit.aes = FALSE, size = 1.5, col = "black") +
  geom_density(aes(monthday, color = cohort), inherit.aes = FALSE, size = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = -0.9), 
        axis.title.y.right = element_text(vjust = 2)) +
  annotate(geom = "text", x = as.Date('2020-11-15'), y = 0.017, label = "n = 7794")

# color coded by year

young_sporophyte_cohort %>%
  filter(is.na(Stipe) | Stipe <= 25) %>%
  ggplot() +
  geom_point(aes(x = monthday, y = Stipe/2000, col = year), alpha = 0.4) +
  scale_y_continuous(
    name = "Stipe Length (cm)", lim = c(0,0.013), n.breaks = 6, labels = c("0", "5", "10", "15", "20", "25"),
    sec.axis = sec_axis(trans = ~.*1, "Stipe Count Density Distribution")) +
  scale_colour_viridis(option = "D", discrete = TRUE, begin = 0.8, end = 0.2) +
  scale_x_date(name = "Month", limits = as.Date(c('2020-05-01', '2020-10-01'), format="%Y-%M-%D"), 
               date_breaks = "1 month",
               date_labels = ("%b")) +
  geom_density(aes(monthday), inherit.aes = FALSE, size = 1) +
  # scale_y_discrete(lim = c(0,0.026), labels = c("0", "4"))
  theme_bw() +
  theme(axis.text.x = element_text(hjust = -0.9), 
        axis.title.y.right = element_text(vjust = 2)) +
  annotate(geom = "text", x = as.Date('2020-01-20'), y = 0.012, label = "n = 6960")
  
# single stipes color coded by year

young_sporophyte_cohort %>%
  filter(is.na(Stipe) | Stipe <= 25) %>%
  filter(Group == 0) %>%
  ggplot() +
  geom_point(aes(x = monthday, y = Stipe/2000, col = year), alpha = 0.4) +
  scale_y_continuous(
    name = "Stipe Length (cm)", lim = c(0,0.013), n.breaks = 6, labels = c("0", "5", "10", "15", "20", "25"),
    sec.axis = sec_axis(trans = ~.*1, "Stipe Count Density Distribution")) +
  scale_colour_viridis(option = "D", discrete = TRUE, begin = 0.8, end = 0.2) +
  scale_x_date(name = "Month", limits = as.Date(c('2020-05-01', '2020-10-01'), format="%Y-%M-%D"), 
               date_breaks = "1 month",
               date_labels = ("%b")) +
  geom_density(aes(monthday), inherit.aes = FALSE, size = 1) +
  # scale_y_discrete(lim = c(0,0.026), labels = c("0", "4"))
  theme_bw() +
  theme(axis.text.x = element_text(hjust = -0.9), 
        axis.title.y.right = element_text(vjust = 2)) +
  annotate(geom = "text", x = as.Date('2020-04-15'), y = 0.013, label = "n = 1720")

# get number of singles in this dataset:
singles <- young_sporophyte_cohort %>%
  filter(is.na(Stipe) | Stipe <= 25) %>%
  filter(Group == 0)


# separate out 'primary recruits' that have secondary recruits on them, and secondary recruits on holdfasts

young_sporophyte_cohort %>%
  filter(Single == 0 & Subst != "Nl") %>%
  ggplot() +
  geom_point(aes(x = monthday, y = Stipe/10000, col = year), alpha = 0.4) +
  scale_y_continuous(
    name = "Stipe Length (cm)", lim = c(0,0.025), n.breaks = 6, labels = c("0", "50", "100", "150", "200", "250"),
    sec.axis = sec_axis(trans = ~.*1, "Stipe Count Density Distribution")) +
  scale_colour_viridis(option = "D", discrete = TRUE, begin = 0.8, end = 0.2) +
  scale_x_date(name = "Month", limits = as.Date(c('2020-05-01', '2020-10-01'), format="%Y-%M-%D"), 
               date_breaks = "1 month",
               date_labels = ("%b")) +
  geom_density(aes(monthday), inherit.aes = FALSE, size = 1) +
  # scale_y_discrete(lim = c(0,0.026), labels = c("0", "4"))
  theme_bw() +
  theme(axis.text.x = element_text(hjust = -0.9), 
        axis.title.y.right = element_text(vjust = 2)) +
  annotate(geom = "text", x = as.Date('2020-04-15'), y = 0.025, label = "n = 806") 

lengthsubs <- young_sporophyte_cohort %>%
  filter(Single == 0 & Subst != "Nl")

young_sporophyte_cohort %>%
  filter(Single == 0 & Subst == "Nl") %>%
  ggplot() +
  geom_point(aes(x = monthday, y = Stipe/7000, col = year), alpha = 0.4) +
  scale_y_continuous(
    name = "Stipe Length (cm)", lim = c(0,0.03), 
    sec.axis = sec_axis(trans = ~.*1, "Stipe Count Density Distribution")) +
  scale_colour_viridis(option = "D", discrete = TRUE, begin = 0.8, end = 0.2) +
  scale_x_date(name = "Month", limits = as.Date(c('2020-05-01', '2020-10-01'), format="%Y-%M-%D"), 
               date_breaks = "1 month",
               date_labels = ("%b")) +
  geom_density(aes(monthday), inherit.aes = FALSE, size = 1) +
  # scale_y_discrete(lim = c(0,0.026), labels = c("0", "4"))
  theme_bw() +
  theme(axis.text.x = element_text(hjust = -0.9), 
        axis.title.y.right = element_text(vjust = 2)) +
  annotate(geom = "text", x = as.Date('2020-04-15'), y = 0.025, label = "n = 806") 

lengthsubs <- young_sporophyte_cohort %>%
  filter(Single == 0 & Subst != "Nl")

lengthsecs <- young_sporophyte_cohort %>%
  filter(Single == 0 & Subst == "Nl")





ggplot(young_sporophyte, aes(x = Stipe, y = Bulb, na.rm = TRUE)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal() +
  labs(y = "Bulb diameter (cm)", x = "Stipe length (cm)", title = "Relationship between stipe length and bulb diameter (young sporophytes)") + 
  annotate("text", x = 150, y=9, label = "Adj. R-squared = 0.62, p < 0.0001")

ggplot(SBgenerated) +
  geom_point(aes(x = S, y = B))



addvectors <- read.csv("Data/test.csv", header = FALSE)
x <- addvectors[,4]
y <- addvectors[,5]
x + y
sum(x)

addvectors <- read.csv("Data/test2.csv", header = FALSE)
x <- addvectors[,4]
y <- addvectors[,5]
x + y
sum(x)

upwellingMonth %>%
  ggplot(aes(x = as.character(month), y = sqrt(ektrx))) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank()) +
  labs(x = "Month", y = "Upwelling")

SST_temp <- SST_all %>%
  filter(WTMP %notin% c(999, 99)) %>%
  mutate(month = month.abb[MM]) 
SST_temp <- SST_temp %>%
  mutate(month = factor(month, levels = c("Jan",
                                          "Feb",
                                          "Mar",
                                          "Apr",
                                          "May",
                                          "Jun",
                                          "Jul",
                                          "Aug",
                                          "Sep",
                                          "Oct",
                                          "Nov",
                                          "Dec"))) 
SST_temp %>%
  ggplot(aes(x = month, y = WTMP, group = MM)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  theme_bw() +
  labs(x = "Month", y = "Water Temp (C)") 


dataset_young <- young_sporophyte %>%
  mutate(dataset = "young")

dataset_mature <- mature_sporophyte %>%
  mutate(dataset = "mature")

# join both datasets into one
dataset_sporophyte <- dataset_young %>%
  bind_rows(dataset_mature)

# remove year so months overlap
all_dataset_monthly <- dataset_sporophyte %>%
  mutate(Date = as.Date.POSIXct(Date, "%Y-%M-%D")) %>%
  mutate(year = as.character(year(Date))) %>%
  mutate(monthday = format(as.Date.POSIXct(Date, "%Y-%M-%D"), "2020-%m-%d")) %>%
  mutate(monthday = as.Date.POSIXct(monthday, "%Y-%M-%D", tz = "PT")) %>%
  arrange(monthday)


# add cohort values
all_dataset_cohort <- all_dataset_monthly %>%
  mutate(cohort = case_when(
    Stipe <= 15 ~ "0.1-15",
    Stipe <= 20 ~ "15.1-20",
    Stipe <= 25 ~ "20.1-25",
    Stipe > 25 ~ "25+")
  ) %>%
  filter(!is.na(cohort))

# stipe length and density by month
all_dataset_cohort %>%
  ggplot() +
  geom_vline(xintercept = as.Date('2018-12-31')) +
  geom_vline(xintercept = as.Date('2019-12-31')) +
  geom_vline(xintercept = as.Date('2020-12-31')) +  
  geom_point(aes(x = Date, y = log(Stipe + 1)/3000, col = dataset), alpha = 0.4) +
  scale_y_continuous(
    name = "Log stipe length", lim = c(0,0.003), n.breaks = 7, labels = c("0", "", "3", "", "6", "", "8"), 
    sec.axis = sec_axis(trans = ~.*1, "Stipe Count Density Distribution")) +
  scale_colour_viridis(option = "D", discrete = TRUE, begin = 0.1, end = 0.9) +
  scale_x_date(name = "Year", 
               date_breaks = "1 year",
               date_labels = ("%Y")) +
  geom_density(aes(Date, group = dataset), inherit.aes = FALSE, size = 1.8, col = "white") +
  geom_density(aes(Date, group = dataset), inherit.aes = FALSE, size = 1.5, col = "black") +
  geom_density(aes(Date, color = dataset), inherit.aes = FALSE, size = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = -1.8), 
        axis.title.y.right = element_text(vjust = 2)) +
  annotate(geom = "text", x = as.Date('2021-06-15'), y = 0.003, label = "n = 7794")


###### use SST_all from markdown to relate max sustained wind speed to wrack observations

# extract daily mean wind speeds (all directions)
dailyWSP <- SST_all %>%
  unite(Date, `#YY`, MM, DD, sep = "-") %>%
  mutate(Date = ymd(Date), .keep = "all") %>%
  group_by(Date) %>%
  mutate(dailySp = mean(WSPD)) %>%
  select(Date, dailySp) %>%
  distinct(Date, dailySp)

# join daily wind speeds to sporophyte data
all_sporophyte_wind <- all_sporophyte %>%
  left_join(dailyWSP, join_by(Date))

# plot of wind speed and sporophyte observations
all_sporophyte_wind %>%
  filter(!is.na(Stipe)) %>%
  group_by(Date, dailySp) %>% 
  tally() %>%
  ggplot(aes(x = dailySp, y = log(n + 1))) +
  geom_point() +
  geom_smooth(method = "lm")

# extract daily mean wind speeds (onshore only)
dailyWSP_onshore <- SST_all %>%
  unite(Date, `#YY`, MM, DD, sep = "-") %>%
  mutate(Date = ymd(Date), .keep = "all") %>%
  filter(WDIR %in% c(1:179)) %>%
  group_by(Date) %>%
  mutate(dailySp_onshore = mean(WSPD)) %>%
  select(Date, dailySp) %>%
  distinct(Date, dailySp)

# join daily wind speeds to sporophyte data
all_sporophyte_wind <- all_sporophyte %>%
  left_join(dailyWSP_onshore, join_by(Date))

# plot of wind speed and sporophyte observations
all_sporophyte_wind %>%
  filter(!is.na(Stipe)) %>%
  group_by(Date, dailySp) %>% 
  tally() %>%
  ggplot(aes(x = dailySp, y = log(n + 1))) +
  geom_point() +
  geom_smooth(method = "lm")

# plot all by wind direction

dailyDir <- SST_all %>%
  unite(Date, `#YY`, MM, DD, sep = "-") %>%
  mutate(Date = ymd(Date), .keep = "all") %>%
  group_by(Date) %>%
  mutate(dailyDegree = mean(WDIR)) %>%
  select(Date, dailyDegree) %>%
  distinct(Date, dailyDegree)
  
# join daily wind speeds to sporophyte data
all_sporophyte_wind <- all_sporophyte %>%
  left_join(dailyDir, join_by(Date))

# plot of wind speed and sporophyte observations
all_sporophyte_wind %>%
  filter(!is.na(Stipe), dailyDegree < 360) %>%
  group_by(Date, dailyDegree) %>% 
  tally() %>%
  ggplot(aes(x = dailyDegree, y = log(n + 1))) +
  geom_point() +
  geom_smooth(method = "lm")

# add one day offset
all_sporophyte_wind %>%
  mutate(minus1offset = Date %m+% days(1)) %>%
  select(minus1offset, dailyDegree) 
all_sporophyte_wind %>%
  filter(!is.na(Stipe), dailyDegree < 360) %>%
  group_by(Date) %>% 
  tally() %>%
  mutate(minus1offset = Date %m+% days(1)) %>%
  left_join(offset1) %>%
  ggplot(aes(x = dailyDegree, y = log(n + 1))) +
  geom_point() +
  geom_smooth(method = "lm")

# add two day offset
all_sporophyte_wind %>%
  mutate(minus2offset = Date %m-% days(2)) %>%
  filter(!is.na(Stipe), dailyDegree < 359) %>%
  group_by(minus2offset, dailyDegree) %>% 
  tally() %>%
  ggplot(aes(x = dailyDegree, y = log(n + 1))) +
  geom_point() +
  geom_smooth(method = "lm")


###### use upwelling from markdown to relate upwelling index to wrack observations

dailyUpwelling <- upwellingTrans %>%
  mutate(Date = date) %>%
  group_by(Date) %>%
  summarise(meanx = mean(ektrx), meany = mean(ektry))

# join daily upwelling indices to sporophyte observations
all_sporophyte_upwell <- all_sporophyte %>%
  left_join(dailyUpwelling)

# plot of upwelling and sporophyte observations (x componenet)
all_sporophyte_upwell %>%
  filter(!is.na(Stipe)) %>%
  group_by(Date, meanx) %>% 
  tally() %>%
  ggplot(aes(x = meanx, y = log(n + 1))) +
  geom_point() +
  geom_smooth(method = "lm")


# plot of upwelling and sporophyte observations (y component)
all_sporophyte_upwell %>%
  filter(!is.na(Stipe)) %>%
  group_by(Date, meany) %>% 
  tally() %>%
  ggplot(aes(x = meany, y = log(n + 1))) +
  geom_point() +
  geom_smooth(method = "lm")

# create new variable combining x and y
all_sporophyte_upwell %>%
  filter(!is.na(Stipe)) %>%
  mutate(index = meanx + meany) %>%
  group_by(Date, index) %>%
  tally() %>%
  ggplot(aes(x = index, y = log(n + 1))) +
  geom_point() +
  geom_smooth(method = "lm")


################### Quick graph of holdfast snapping


holdfastSnapping <- read_csv("Data/holdfastSnapping.csv", 
                             col_types = cols(WhHf = col_number(), 
                                              ChWhHf = col_number(), BuWhHf = col_number(), 
                                              `Sn Hf` = col_number(), `Sn St` = col_number(), 
                                              Ur = col_number()))

snapLong <- holdfastSnapping %>%
  select(-`WhHf Sub`) %>%
  pivot_longer(WhHf:Ur, names_to = "breakage", values_to = "value")

ggplot(snapLong) +
  geom_col(aes(x = breakage, y = value, fill = breakage)) +
  scale_fill_viridis(discrete = TRUE, option = "D", begin = 0.1, end = 0.8) +
  theme_bw()









