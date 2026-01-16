# Part 1

# (1) setup ----

library(terra) 
library(tidyterra) 
library(dplyr)
library(ggplot2)

# (2) create AOI ----

## ARP ----
# the Area of Interest (AOI) for this case study is the Arapaho-Roosevelt National Forest (ARP)
### load & process ----
NF_CONUS_vect <- vect("S_USA.FSCommonNames.shp")

plot(NF_CONUS_vect)

# see unique names 
names(NF_CONUS_vect)
unique(NF_CONUS_vect$COMMONNAME)

# select for just ARP 
ARP_vect <- NF_CONUS_vect %>%
  filter(COMMONNAME == "Arapaho and Roosevelt National Forests")

plot(ARP_vect)

# project 
ARP_vect <- project(ARP_vect,"EPSG:5070")

# calc area
expanse(ARP_vect) # 6975245280 m^2
6975245280/4046.86 # 4046.86 m/acre = 1723619 acres

expanse(ARP_vect, unit = "ha")

#### write & read ----
writeVector(ARP_vect, "ARP_vect.shp")
ARP_vect <- vect("ARP_vect.shp")

## SRME ----
# Southern Rocky Mountain Ecoregion
### load & process ----
EPA_ecoregions <- vect("us_eco_l3.shp") # 1250 geoms

# see unique names 
names(EPA_ecoregions)
unique(EPA_ecoregions$US_L3NAME)

# select for just SRME
SRME_s_rockies <- EPA_ecoregions %>% 
  filter(US_L3NAME == "Southern Rockies")
# has 4 separate polygons

# aggregate them together
SRME_aggregated <- terra::aggregate(SRME_s_rockies)
plot(SRME_aggregated)
# has 1 geom & 0 attributes (they are lost after aggregate)

# project 
SRME_vect <- project(SRME_aggregated, "EPSG:5070")

#### write & read ----
writeVector(SRME_vect, "SRME_vect.shp")
SRME_vect <- vect("SRME_vect.shp")

# (3) pre-process data ----

## QMD ----
# this is using quadratic mean diameter (QMD) from TreeMap 2022
QMD_CONUS <- rast("TreeMap2022_CONUS_QMD.tif")
# already in 5070
plot(QMD_CONUS)

### crop and mask ----
SRME_QMD_rast <- crop(QMD_CONUS, SRME_vect, mask=TRUE)
plot(SRME_QMD_rast)
plot(is.na(SRME_QMD_rast))

#### write & read ----
writeRaster(SRME_QMD_rast, "SRME_QMD_rast.tif")
SRME_QMD_rast <- rast("SRME_QMD_rast.tif")

# reclassify with ifel()
SRME_QMD_filt_rast <- ifel(
  SRME_QMD_rast >= 5, 5, NA 
)
# if >= 5 inches, reclassify to 5
# if < 5 inches, reclassify to NA

### stats ----
# entire SRME = X cells 
# all areas with QMD values
global(SRME_QMD_rast, fun = "notNA") # 100858399 cells
(5697616/X)*100 # X % of SRME has QMD values

# areas with QMD > 5 inches
global(SRME_QMD_filt_rast, fun = "notNA") # 74599513
(74599513/X)*100 # X % of SRME has trees > 5 in QMD

### viz ----
# see classified values
plot(SRME_QMD_filt_rast, col = "darkgreen")
polys(SRME_vect, col = "black", alpha=0.01, lwd=1.5)

#### write & read ----
writeRaster(SRME_QMD_filt_rast, "SRME_QMD_filt_rast.tif")
SRME_QMD_filt_rast <- rast("SRME_QMD_filt_rast.tif")


# ** left off here ----
# we want to download the 2024 EVH data first

## EVH ----
# using existing vegetation height (EVH) from LANDFIRE
# these values are not continuous
# also the veg height has an offset added
# e.g. value 103 = tree height of 3 meters

EVH_CONUS <- rast("LC23_EVH_240.tif")
crs(EVH_CONUS) # 5070
res(EVH_CONUS) # 30 30

### crop / mask ----
EVH_SRME <- crop(EVH_CONUS, SRME_vect, mask=TRUE)

### adjust values ----
# define conversion factor
meters_to_feet_factor <- 3.28084

# reclassify with ifel()
SRME_EVH_filt_rast <- ifel(
  # condition 1: it is dominant veg type trees? (values 100-199)
  EVH_SRME >= 100 & EVH_SRME < 200,
  # if TRUE, 
  # condition 2: is it > 10 ft tall? 
  ifel(
    (EVH_SRME - 100) * meters_to_feet_factor > 10, # subtract offset, convert units, filter
    10, # if TRUE, reclassify to 10
    NA # if FALSE, reclassify to NA
  ),
  NA # if not a tree value (condition 1 = FALSE), reclassify to NA
)

### stats ----
# entire ARP = 7776004 cells

# all veg area
global(EVH_SRME >= 100, fun = "sum", na.rm = TRUE) # 7004697 cells
(7004697/7776004)*100 # 90.08093 % of ARP is vegetated 

# all tree area
global(EVH_SRME >= 100 & EVH_SRME < 200, fun = "sum", na.rm = TRUE) # 5324379 cells
(5324379/7776004)*100 # 68.47192 % of ARP has trees 

# trees > 10 ft area
global(SRME_EVH_filt_rast == 100, fun = "sum", na.rm = TRUE) # 5219760 cells
(5219760/7776004)*100 # 67.12651 % of ARP has trees > 10 ft

### viz ----
plot(SRME_EVH_filt_rast, col = "forestgreen")
polys(SRME_vect, col = "black", alpha=0.01, lwd=1.5)

#### write & read ----
writeRaster(SRME_EVH_filt_rast, "SRME_EVH_filt_rast.tif")
SRME_EVH_filt_rast <- rast("SRME_EVH_filt_rast.tif")



## slope ----
# this slope raster is generated using  
# digital elevation model (DEM) tiles, downloaded from The National Map (USGS)
# they are 1 Arc Sec
# these tiles have GEOGCRS NAD83, but are not yet projected

### load & process DEMs ----
#### CO ----
DEM_n41_w106 <- rast("USGS_1_n41w106_20230314.tif")
DEM_n41_w107 <- rast("USGS_1_n41w107_20230314.tif")
DEM_n40_w106 <- rast("USGS_1_n40w106_20230602.tif")
DEM_n40_w107 <- rast("USGS_1_n40w107_20220216.tif")


#### WY ----


#### NM ----



#### UT ----




### combine ----
# mosaic tiles together
SRME_DEM_unprojected <- mosaic(DEM_n41_w106, DEM_n41_w107, DEM_n40_w106, DEM_n40_w107, 
                  fun="first")
# project
SRME_DEM_projected <- project(SRME_DEM, "EPSG:5070")

# crop and mask the DEM to the extent of ARP 
SRME_DEM_rast <- crop(SRME_DEM_projected, SRME_vect, mask=TRUE)
plot(SRME_DEM_rast) # min = 1470.285 , max = 4393.409 (meters)

#### write & read ----
writeRaster(SRME_DEM_rast, "SRME_DEM_rast.tif")
SRME_DEM_rast <- rast("SRME_DEM_rast.tif")


### calc slope ----
SRME_slope_rast = terrain(SRME_DEM_rast, v="slope", unit="degrees")
plot(SRME_slope_rast)

### adjust values ----
minmax(SRME_slope_rast) 
# min = 0, max = 72.59397 
# but the max we want to include is 24 degrees
# and we want 0-24 degree slope to become 0-1 score (normalize)

# make all values > 24 degrees NA, leave other values as-is
SRME_slope_filt_rast <- ifel(SRME_slope_rast > 24, NA, ARP_slope)

### viz ----
plot(SRME_slope_filt_rast)
polys(SRME_vect, col = "black", alpha=0.01, lwd=2)

plot(is.na(SRME_slope_filt_rast))

### stats ----
global(SRME_slope_filt_rast, fun = "notNA") # 7433981 cells 

# entire ARP = 7776004 cells 
# (7433981/7776004)*100 # 95.60156 % remaining after 24* filter 

#### write & read ----
writeRaster(SRME_slope_filt_rast, "SRME_slope_filt_rast.tif")
SRME_slope_filt_rast <- rast("SRME_slope_filt_rast.tif")



## road ----
### inport ----
#### CO----
# import CO roads shapefile
# downloaded from The National Map
roads_CONUS <- vect("Trans_RoadSegment_0.shp")
plot(roads_CONUS)
crs(roads_CONUS) # EPSG 4269

road_df <- as.data.frame(roads_CONUS)
# could filter by road type, we did not

# project, crop & mask 
roads_CONUS = project(roads_CONUS, "EPSG:5070")
crs(roads_CONUS) # EPSG 5070

roads_ARP = crop(roads_CONUS, ARP_vect)
plot(roads_ARP)

#### WY ----


#### NM ----


#### UT ----

### combine ----



### rasterize ----
SRME_road_rast <- rasterize(roads_SRME, ARP_risk_score_rast , touches=TRUE)
plot(ARP_road_rast, col="blue") # all values = 1
plot(is.na(ARP_road_rast)) # values not 1 are NA
# TBH, the raster does not look nearly as contiguous as the road lines from the .shp
# but when I open the .tif in Arc, it looks fine 
# I think it is too much for R studio to render with plot()

#### write & read file ----
writeRaster(ARP_road_rast, "ARP_road_rast.tif")
ARP_road_rast <- rast("ARP_road_rast.tif") 

### distance ----
# we will calculate the distance to nearest road for each raster cell (pixel)
ARP_road_dist_rast <- distance(ARP_road_rast, unit="m", method="haversine") 
plot(ARP_road_dist_rast)
# cell values = distance to nearest road (in meters)

#### write & read file ----
writeRaster(ARP_road_dist_rast, "ARP_road_dist_rast.tif")
ARP_road_dist_rast <- rast("ARP_road_dist_rast.tif")

### adjust values ----
minmax(ARP_road_dist_rast) 
# min = 0, max = 37416.17 
# but the max we want to include is 917.3261 meters (0.57 miles)
# and we want 0-917 m distance to become 0-1 score (normalize)

# make NA all values > 917.3261 meters, leave other values as-is
road_filt_rast <- ifel(ARP_road_dist_rast > 917.3261, NA, ARP_road_dist_rast)
plot(road_filtered)

### crop ----
# need to crop again bc the road distance buffer goes a bit outside of the SRME
SRME_road_filt_rast = crop(road_filt_rast, SRME_vect, mask = TRUE)

### viz ----
plot(ARP_road_filt_rast)
polys(ARP_vect, col = "black", alpha=0.01, lwd=1.5)
plot(is.na(ARP_road_filt_rast))

### stats ----
global(ARP_road_filt_rast, fun = "notNA") # 5213776 cells 

# entire ARP = 7776004 cells 
(5213776/7776004)*100 # 67.04955 % remaining  

#### write & read ----
writeRaster(ARP_road_filt_rast, "ARP_road_filt_rast.tif")
ARP_road_filt_rast <- rast("ARP_road_filt_rast.tif")






