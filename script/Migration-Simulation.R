
library(tidyverse)
library(sf)
library(ebirdst)
library(terra)
library(viridis)

rm(list=ls())

sf_use_s2(FALSE)

setwd("C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Shorebirds/Shorebird-Migration-Sim/")

# -----------------------------------------------------
# Load data
# -----------------------------------------------------

proj <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
lcc <- "+proj=lcc +lon_0=-90 +lat_1=33 +lat_2=45"

if (!file.exists("data/sites.RData")){
  
  load("AdamCSmithCWS-Shorebird_Migration_Trends-cdf9e9a/data/full_observation_dataset.Rdata")
  
  sites <- ssData %>%
    dplyr::select(DecimalLatitude,DecimalLongitude,YearCollected) %>%
    unique() %>%
    group_by(DecimalLatitude,DecimalLongitude) %>%
    rename(lat = DecimalLatitude, lon = DecimalLongitude) %>%
    summarize(nyr = n()) %>%
    subset(nyr >= 10)
  
  sites <- st_as_sf(sites, coords = c("lon", "lat"), crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  
  save(sites, file = "data/sites.RData")
}

load(file = "data/sites.RData")

#load BCR boundaries
study_region <- read_sf("../Arctic-PRISM-Spatial-Model/data/Spatial_Covariates/BCR/BCR_Terrestrial_master.shp") %>%
  subset(WATER == 3) %>%
  st_union() %>%
  st_transform(proj)

sites <- st_transform(sites,proj, remove = FALSE)
sites$x <- as.data.frame(st_coordinates(sites))$X
sites$y <- as.data.frame(st_coordinates(sites))$Y

# -----------------------------------------------------
# Create a large number of potential stopover locations within the migration path
# -----------------------------------------------------

bbox <- st_bbox(study_region)

loc <- data.frame(y = runif(20000,bbox$ymin,bbox$ymax),
                  x = runif(20000,bbox$xmin,bbox$xmax)) %>%
  st_as_sf(coords = c("x", "y"), crs = proj, remove = FALSE) %>%
  bind_rows(sites %>% mutate(PRISM = 1))

# -----------------------------------------------------
# Load ebird raster
# -----------------------------------------------------

species <- subset(ebirdst_runs, common_name == "Red Knot")
ebirdst_download_status(species$species_code,pattern = "_mean_27km_")

ebirdSDM <- load_raster(species$species_code, product = "abundance", 
                        period = "seasonal", metric = "mean", 
                        resolution = "27km")

ebird_mig <- ebirdSDM$prebreeding_migration %>% terra::project(proj) %>% crop(vect(study_region), mask = TRUE)
values(ebird_mig)[values(ebird_mig) == 0] <- NA
ebird_nb <- ebirdSDM$nonbreeding %>% terra::project(proj) %>% crop(vect(study_region), mask = TRUE)
values(ebird_nb)[values(ebird_nb) == 0] <- NA
ebird_br <- ebirdSDM$breeding %>% terra::project(proj) %>% crop(vect(study_region), mask = TRUE)
values(ebird_br)[values(ebird_br) == 0] <- NA

# Only use upper 80% of values for ranges
values(ebird_mig)[values(ebird_mig) < quantile(values(ebird_mig),0.20,na.rm = TRUE)] <- NA
values(ebird_br)[values(ebird_br) < quantile(values(ebird_br),0.50,na.rm = TRUE)] <- NA
values(ebird_nb)[values(ebird_nb) < quantile(values(ebird_nb),0.20,na.rm = TRUE)] <- NA

loc$mig <- extract(ebird_mig,vect(loc %>% st_transform(crs(ebird_mig))))[,2]
loc$nb <- extract(ebird_nb,vect(loc %>% st_transform(crs(ebird_mig))))[,2]
loc$br <- extract(ebird_br,vect(loc %>% st_transform(crs(ebird_mig))))[,2]

loc <- subset(loc, !is.na(mig) | !is.na(nb) | !is.na(br))
loc$id <- 1:nrow(loc)

# Nonbreeding distribution
ggplot(loc)+
  geom_sf(data = study_region, col = "gray85", fill = "gray90")+
  geom_sf(aes(col = nb), size = 5)+
  scale_color_gradientn(colors = viridis(10),trans = "log10", na.value = "transparent")+
  theme_bw()

# Migration routes
ggplot(loc)+
  geom_sf(data = study_region, col = "gray85", fill = "gray90")+
  geom_sf(aes(col = mig), size = 5)+
  scale_color_gradientn(colors = viridis(10),trans = "log10", na.value = "transparent")+
  theme_bw()

# Breeding distribution
ggplot(loc)+
  geom_sf(data = study_region, col = "gray85", fill = "gray90")+
  geom_sf(aes(col = br), size = 5)+
  scale_color_gradientn(colors = viridis(10),trans = "log10", na.value = "transparent")+
  theme_bw()

# Monitoring Sites
ggplot(subset(loc, PRISM == 1))+
  geom_sf(data = study_region, col = "gray85", fill = "gray90")+
  geom_sf()+
  scale_color_gradientn(colors = viridis(10),trans = "log10", na.value = "transparent")+
  theme_bw()

# -----------------------------------------------------
# Construct transition probability matrix for each potential site
# -----------------------------------------------------

# Distances between sites
dists <- st_distance(loc,loc)

# Transition matrix
tmat <- matrix(0,nrow=nrow(loc),ncol = nrow(loc))

for (s in 1:nrow(loc)){
  
  sdat  <- loc[s,]
  sdist <- as.numeric(dists[s,]) # distances to other sites
  sdist[sdist==0] <- NA
  
  # sites to move to
  to <- which(sdist >= 500000 & sdist <= 1000000 & loc$y > sdat$y)
  
  # which sites are breeding locations
  if (length(to) == 0){
    
    # If this is a breeding site, migrants stay here
    if (!is.na(sdat$br)) to <- sdat$id
    
    # If this is a nonbreeding site, move birds to nearest available location
    if (is.na(sdat$br)) to <- which.min(sdist)
    
  }
  
  ggplot()+
    geom_sf(data = study_region, col = "gray80",fill = "gray90")+
    geom_sf(data = loc, col = "gray50")+
    geom_sf(data = st_buffer(st_transform(sdat,lcc),1000000), fill = "transparent")+
    geom_sf(data = sdat, col = "black", size = 3)+
    geom_sf(data = loc[to,], aes(col = mig), size = 2)+
    scale_colour_gradientn(colors = viridis(10), trans = "log10", na.value = "transparent")+
    theme_bw()
  
  # Proportion of birds from the cell arriving at new cells
  tprob <- loc$mig[to] / sum(loc$mig[to],na.rm = TRUE)
  tprob[is.na(tprob)] <- 0
  
  if (sum(tprob) == 0) tprob <- rep(1,length(tprob))/length(tprob)
  
  # Final adjustments
  tprob <- tprob/sum(tprob)
  tmat[to,s] <- tprob 
  
}

# -----------------------------------------------------
# Simulate spatially variable population trends across nonbreeding range
# -----------------------------------------------------

# Trends at each location
loc$trend = -(loc$x - mean(loc$x))/sd(loc$x) * 0.3

# Nonbreeding locations only
nb_loc <- subset(loc, !is.na(nb))

trend_lim <- max(abs(nb_loc$trend))

# Plot nonbreeding trends
ggplot()+
  geom_sf(data = study_region)+
  geom_sf(data = nb_loc, aes(col = trend))+
  scale_color_gradientn(colors = c("red","white","blue"), limits = c(-trend_lim,trend_lim))

# -----------------------------------------------------
# Simulate migration each year
# -----------------------------------------------------

# Initial abundance at each location
prop <- loc$nb
prop[is.na(prop)] <- 0
prop <- prop/sum(prop)
loc$N0 <- 100000*prop

# Simulate migration dynamics and observed data
simdat <- data.frame()

nyr <- 10
for (yr in 1:nyr){
  
  # Abundance at each location
  N <- exp(log(loc$N0) + loc$trend * (yr-1))
  
  simdat <- rbind(simdat, loc %>% mutate(year = yr, week = 0, start = N, end = N))
  for (wk in 1:16){
    
    if (wk == 1){
      wkdat <- loc %>% mutate(year = yr, week = wk, start = N)
      wkdat$start[is.na(wkdat$start)] <- 0
      
    } else{
      
      wkdat$week <- wkdat$week + 1
      wkdat$start <- wkdat$end
    }
    
    wkdat$end <- (tmat %*% wkdat$start)[,1]
    
    simdat <- rbind(simdat, wkdat)
    print("Migration Underway!")
    
  } # weeks
  
}

# -----------------------------------------------------
# Plot observed counts across weeks during a focal year
# -----------------------------------------------------

ggplot()+
  geom_sf(data = study_region, col = "gray90", fill = "gray95")+
  geom_sf(data = subset(subset(simdat, year == 1)), aes(col = end))+
  scale_color_gradientn(colors = viridis(10), limits = c(1,max(simdat$end)),na.value = "transparent", trans = "log10")+
  theme_bw()+
  facet_wrap(week~.)

# -----------------------------------------------------
# Observed annual counts at each location
# -----------------------------------------------------

obsdat <- simdat %>%
  group_by(id,year) %>%
  summarize(count = sum(end),
            nb = mean(nb),
            br = mean(br),
            mig = mean(mig),
            PRISM = mean(PRISM),
            true_trend = mean(trend),
            x = mean(x),
            y = mean(y))

# Remove sites that did not detect any birds
nbirds_detected <- obsdat %>%
  as.data.frame() %>%
  group_by(id,year) %>%
  summarize(mean_count = mean(count))

obsdat <- subset(obsdat, id %in% subset(nbirds_detected, mean_count > 0)$id)

# -----------------------------------------------------
# Observed trends at each location
# -----------------------------------------------------

site_trends <- data.frame()
for (s in unique(obsdat$id)){
  
  sdat <- subset(as.data.frame(obsdat), id == s)
  
  trend <- lm(log(count)~year, data = sdat)
  
  site_trends <- rbind(site_trends, data.frame(id = s,observed_trend = trend$coefficients[2]))
  
}

# Summarize data at each station
obsdat <- full_join(obsdat,site_trends)

# Observed trends at nonbreeding locations
trend_lim <- max(abs(as.data.frame(obsdat)[,c("observed_trend")]))

# nonbreeding sites only
ggplot()+
  geom_sf(data = study_region)+
  geom_sf(data = subset(obsdat, !is.na(nb)), aes(col = observed_trend))+
  scale_color_gradientn(colors = c("red","white","blue"), limits = c(-trend_lim,trend_lim))

# Migration
ggplot()+
  geom_sf(data = study_region)+
  geom_sf(data = subset(obsdat,!is.na(mig) & mig>0), aes(col = observed_trend))+
  scale_color_gradientn(colors = c("red","white","blue"), limits = c(-trend_lim,trend_lim))

# Breeding sites
ggplot()+
  geom_sf(data = study_region)+
  geom_sf(data = subset(obsdat,!is.na(br)), aes(col = observed_trend))+
  scale_color_gradientn(colors = c("red","white","blue"), limits = c(-trend_lim,trend_lim))

# All sites
ggplot()+
  geom_sf(data = study_region)+
  geom_sf(data = obsdat, aes(col = observed_trend))+
  scale_color_gradientn(colors = c("red","white","blue"), limits = c(-trend_lim,trend_lim))


# ggplot()+
#   geom_sf(data = study_region)+
#   geom_sf(data = subset(obsdat, !is.na(nb)), aes(col = true_trend))+
#   scale_color_gradientn(colors = c("red","white","blue"), limits = c(-trend_lim,trend_lim))

# -----------------------------------------------------
# Estimate of population trend based on migration locations,
# compared to truth
# -----------------------------------------------------

observed_total <- simdat %>%
  as.data.frame() %>%
  subset(!is.na(mig) & mig>0) %>%
  group_by(year) %>%
  summarize(count = sum(end))
observed_total$percent_change <- 100*(observed_total$count - observed_total$count[1])/observed_total$count[1]

PRISM_total <- simdat %>%
  as.data.frame() %>%
  subset(PRISM == 1) %>%
  group_by(year) %>%
  summarize(count = sum(end))
PRISM_total$percent_change <- 100*(PRISM_total$count - PRISM_total$count[1])/PRISM_total$count[1]

true_total <- simdat %>%
  as.data.frame() %>%
  subset(!is.na(nb) & week == 1) %>%
  group_by(year) %>%
  summarize(count = sum(end))
true_total$percent_change <- 100*(true_total$count - true_total$count[1])/true_total$count[1]

ggplot()+
  geom_line(data = true_total, aes(x = year, y = percent_change, col = "True"), linewidth = 2)+
  geom_line(data = observed_total, aes(x = year, y = percent_change, col = "Estimated"), linewidth = 2)+
  geom_line(data = PRISM_total, aes(x = year, y = percent_change, col = "PRISM"), linewidth = 2)+
  ylab("Percent Change")+
  scale_color_manual(values = c("dodgerblue","orangered","black"), name = "")+
  theme_bw()
