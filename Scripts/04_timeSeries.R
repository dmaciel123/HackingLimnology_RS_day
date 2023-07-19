require(rstac)
require(gdalcubes)
require(terra)
require(RColorBrewer)
require(rasterVis)
require(animation)

s_obj <- stac("https://planetarycomputer.microsoft.com/api/stac/v1/")

## Select date interval (2021 - 2022)

datas = c(paste('2021-01-01', '2021-12-30', sep = '/'))

CLOUD_COVER = 20

BBOX = c(-46.638311, -23.766166, -46.638311, -23.766166)

it_obj <- s_obj %>%
  stac_search(collections = "sentinel-2-l2a",
              bbox = BBOX, #xmin, ymin, xmax, ymax
              datetime = datas) %>%
  get_request() %>% 
  items_filter(`eo:cloud_cover` < as.numeric(CLOUD_COVER)) %>%
  items_filter(`s2:mgrs_tile` == "23KLP") %>%
  
  items_sign(sign_fn = sign_planetary_computer())


# Select the assets

assets = c("B02","B03", "B04", "B05", 'SCL')

col = stac_image_collection(it_obj$features, asset_names = assets)


v = cube_view(srs = "EPSG:32723", extent = list(t0 = '2021-01-01', t1 = '2021-12-31',
                                                left = 315000, right = 360000, 
                                                bottom = 7350000, top = 7380000),
              dx = 50, dy = 50, dt = "P5D", aggregation = 'min', resampling = 'average')

gdalcubes_options(parallel = 4)



# Colors

colr <- colorRampPalette(rev(brewer.pal(11, 'RdBu')))

## Save TIFS

RES = raster_cube(col, v) %>% #, mask = image_mask("SCL", values = 6, invert = T)) %>% 
  select_bands(assets) %>% 
  apply_pixel("(B05-B04)/(B05+B04)+1", "NDCI", keep_bands = T) %>% 
  write_tif(dir = 'Outputs/TimeSeries')


## Load results

files_saved = list.files(path = 'Outputs/TimeSeries/',pattern = '.tif$',
                         full.names = T)

s2.col = rast(files_saved[10])


## Plot the first one

plotRGB(s2.col, r = 3, g = 2, b= 1, stretch = 'lin')



# Apply the RF algorithm to the time series

# Load model

RF.model = readRDS("Outputs/rf_chla.R")

for(i in 1:length(files_saved)) {
  
  validity = rast(files_saved[i])[[1]] %>% values %>% na.omit()
  
  if(dim(validity)[1] != 0) {
    
    
  stacked = rast(files_saved[i])
    
  # Mask based on SCL (Maybe it is not the best idea... )
  
  stacked[stacked$SCL != 6] = NA
  
  names(stacked) = c('x490', 'x560', 'x660', 'x705',"SCL", 'NDCI')
  
  Chla = predict(stacked, RF.model)
    
  date = strsplit(files_saved[i], split = '-')[[1]][c(2,3)]
  
  save_name = paste('Outputs/TimeSeries/RandomForest/Chla', '2021', date[1], date[2], sep = "_")
  
  writeRaster(Chla, save_name)
    
  print(i)
  }
  

  if(dim(validity)[1] == 0) {
    unlink(files_saved[i])
    
  }
  
}


file_chla_names = list.files("Outputs/TimeSeries/RandomForest/", pattern = '.tif$', full.names = T) 

chla.res = file_chla_names %>% rast()

# Dates: 

dates = list.files("Outputs/TimeSeries/RandomForest/")  %>% 
            gsub(pattern = 'Chla_', replacement = '') %>% 
            gsub(pattern = '.tif', replacement = '')
  
names(chla.res) = dates

colr <- colorRampPalette(rev(brewer.pal(11, 'RdBu')))


levelplot(raster::stack(chla.res),col.regions = colr, maxpixels = 1e6, main = "Chlorophyll-a Concentration (ug/L)",
          colorkey=list(
            space='bottom',                   # plot legend at bottom
            labels=list(at=seq(from = 0, to = 1000, by = 100), font=4)      # legend ticks and labels 
          ), names = dates) 




saveGIF({
  for(i in c(1:dim(chla.res)[3])){
    
    l <- levelplot(raster::stack(chla.res[[i]]),col.regions = colr, maxpixels = 1e6, main = dates[i],
                   colorkey=list(
                     space='bottom',                   # plot legend at bottom
                     labels=list(at=seq(from = 0, to = 1000, by = 100), font=4)      # legend ticks and labels 
                   )) 
    
    plot(l)
  }
}, interval=0.4, movie.name="animation.gif")




