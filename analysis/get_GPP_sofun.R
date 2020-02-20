# GPP Beni
library(raster)
library(tidyverse)

sites <- read.table("data-raw/yangambi_sites.csv",
                    header = TRUE,
                    sep = ",")

# Define lat / lon projection.
lat_lon <- CRS("+init=epsg:4326")

# Read in the coordinates and assign them a projection, # otherwise they remain just 'numbers'
locations <- SpatialPoints(cbind(sites$lon,sites$lat), lat_lon)

# gpp sofun
gpp_sofun <- stack("~/Downloads/s0_fapar3g_v2_global.d.gpp.nc")

# extract modelled gpp values
gpp <- extract(gpp_sofun, locations)

GPP_sofun <- data.frame(value = colMeans(gpp))
GPP_sofun$date <- as.Date(rownames(GPP_sofun), "X%Y.%m.%d")

# create mean values by DOY
GPP_sofun_s <- GPP_sofun %>%
  mutate(doy = as.numeric(format(as.Date(date), "%j"))) %>%
  group_by(doy) %>%
  summarize(GPP = mean(value, na.rm = TRUE))

# plot EVI by site
p4 <- ggplot(GPP_sofun_s, aes(doy, GPP)) +
  geom_point(col = "black") +
  geom_smooth(method = "loess", span = 0.3, col = "black") +
  ylab(bquote("sofun GPP (g C/" ~ m^2~"d)")) +
  xlab("") +
  theme_bw()

print(p4)
