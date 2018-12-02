# extract transition dates from JR
# weekly annotations to be ingested by the
# format_csv() function of phenor
library(zoo)

df <- readRDS("./data/jungle_rhythms_weekly_annotations.rds")

# first, find the date on which an event happens
transition_dates <- df %>%
  filter(phenophase == "flowers") %>%
  mutate(date = as.Date(sprintf("%s-%02d-1",year, week),
                                                "%Y-%W-%u")) %>%
  group_by(id, year, genus, species) %>%
  arrange(date) %>%
  mutate(value = ifelse(value == 0, NA, 1)) %>%
  mutate(value = na.approx(value, maxgap = 2, na.rm = FALSE)) %>%
  mutate(value = ifelse(is.na(value), 0, 1)) %>%
  mutate(diff_value = ifelse(c(value - lag(value)) == 1,
                             1, 0)) %>%
  summarize(date = date[which(diff_value == 1)[1]]) %>%
  na.omit() %>%
  ungroup()

# then reformat those taking the year - date line into
# consideration
transition_dates <- transition_dates %>%
  mutate(doy = as.numeric(format(date, "%j"))) %>%
  mutate(year = ifelse(doy > 180, year + 1, year),
         doy = ifelse(doy > 180, doy - 365, doy),
         site = paste(id,genus,species,sep="_"))

# approximate location of Yangambi
transition_dates$lat <- 0.804593
transition_dates$lon <- 24.452605

test <- transition_dates %>%
  filter(genus == "Anthonotha")

hist(test$doy)
#
# my_request <- list(class="em",
#                 dataset="era20cm",
#                 date="19440101/19440201/19440301/19440401/19440501/19440601/19440701/19440801/19440901/19441001/19441101/19441201",
#                 expver=1,
#                 levtype="sfc",
#                 number=0,
#                 param=167.128,
#                 step=0,
#                 stream="edmm",
#                 time="00:00:00/03:00:00/06:00:00/09:00:00/12:00:00/15:00:00/18:00:00/21:00:00",
#                 type="fc",
#                 format = "netcdf",
#                 target="output.nc")
#
# wf_request(
#   email = "koenhufkens@gmail.com",
#   transfer = TRUE,
#   path = "~",
#   request = my_request)
#
#
