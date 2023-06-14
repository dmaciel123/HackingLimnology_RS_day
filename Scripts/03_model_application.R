# Using RSTAC to get Sentinel-2 image and apply the algorithm developed

require(data.table)
require(dplyr)
require(rstac)
require(terra)
require(mapview)
require(httr)
require(Metrics)
require(geodrawr)
require(svDialogs)
require(rstac)
require(wesanderson)
require(randomForest)
library(rasterVis)
require(RColorBrewer)

# What is STAC? https://stacspec.org/en

# What we need to get the images?

# 1) Collection and STAC provider (from where we gonna get that?)
# 2) Location and/or image name/date (which date and where we want that data?)


# We gonna use the Microsoft Planetary Computer STAC to get Sentinel-2 images

stac_obj <- stac('https://planetarycomputer.microsoft.com/api/stac/v1/')


dates = c("2021-08-10/2021-09-21")
CLOUD_COVER = 20

BBOX = c(-55.87, -2.02, -55.57, -2.01)

it_obj <- stac_obj %>%
  stac_search(collections = "sentinel-2-l2a",
              bbox = BBOX,
              datetime = dates) %>%
  get_request() %>% 
  items_filter(`eo:cloud_cover` < as.numeric(CLOUD_COVER)) %>%
  items_sign(sign_fn = sign_planetary_computer())  


print(it_obj)


crop_pt = ext(600000, 700000,9740000, 9800000)


blue <- paste0("/vsicurl/", it_obj$features[[1]]$assets$B02$href) %>% rast() %>% crop(crop_pt)
green <- paste0("/vsicurl/", it_obj$features[[1]]$assets$B03$href) %>% rast() %>% crop(crop_pt)
red <- paste0("/vsicurl/", it_obj$features[[1]]$assets$B04$href) %>% rast() %>% crop(crop_pt)


# Plot to see if it ok

stack_img = c(blue, green, red)
plotRGB(stack_img, r = 3, g = 2, b =1, stretch = 'lin')


x705 <- paste0("/vsicurl/", it_obj$features[[1]]$assets$B05$href) %>% rast() %>% crop(crop_pt)
x850 <- paste0("/vsicurl/", it_obj$features[[1]]$assets$B8A$href) %>% rast() %>% crop(crop_pt)

# Water Mask 

SCL <- paste0("/vsicurl/", it_obj$features[[1]]$assets$SCL$href) %>% rast() %>% crop(crop_pt)
SCL[SCL[[1]] != 6] = NA

# Reprojecting nir/swir to match the 20m spatial res

blue = project(blue, x705)
green = project(green, x705)
red = project(red, x705)


# full stack

img.full = c(blue, green, red, x705, x850) %>% mask(SCL)

names(img.full) = c('x490', 'x560', 'x660', 'x705', 'x850')

# Calculate index

img.full$NDCI = (img.full$x705-img.full$x660)/(img.full$x705+img.full$x660)+1


## Random forest prediction


# load model

rf.chla = readRDS('Outputs/rf_chla.R')
emp.chla = readRDS('Outputs/emp_chla.R')

rf.chla.pred = predict(img.full, rf.chla)
emp.chla.pred = exp(predict(img.full, emp.chla))


models = c(rf.chla.pred, emp.chla.pred)
names(models) = c("Random Forest", "Empirical NDCI")

colr <- colorRampPalette(rev(brewer.pal(11, 'RdBu')))

levelplot(raster::stack(emp.chla.pred),col.regions = viridis::viridis, maxpixels = 1e6, main = "Chlorophyll-a Concentration (ug/L)",
          colorkey=list(
            space='bottom',                   # plot legend at bottom
            labels=list(at=seq(from = 0, to = 1000, by = 100), font=4)      # legend ticks and labels 
          )) 

par(mfrow=c(1,2))
plot(rf.chla.pred)
plot(emp.chla.pred)





