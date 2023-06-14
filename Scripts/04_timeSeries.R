require(rstac)
require(gdalcubes)
require(terra)
require(RColorBrewer)

s_obj <- stac("https://planetarycomputer.microsoft.com/api/stac/v1/")

datas = c(paste('2021-06-01', '2021-06-30', sep = '/'))


###  Série com NDVI 


datas = c(paste('2021-01-01', '2021-12-30', sep = '/'))

CLOUD_COVER = 20

BBOX = c(-46.638311, -23.766166, -46.638311, -23.766166)

it_obj <- s_obj %>%
  stac_search(collections = "sentinel-2-l2a",
              bbox = BBOX, #xmin, ymin, xmax, ymax
              datetime = datas) %>%
  get_request() %>% 
  items_sign(sign_fn = sign_planetary_computer())


# Select the assets

assets = c("B02","B03", "B04", "B05", 'SCL')

col = stac_image_collection(it_obj$features, asset_names = assets)


v = cube_view(srs = "EPSG:32723", extent = list(t0 = '2021-07-01', t1 = '2021-07-31',
                                                left = 315000, right = 360000, 
                                                bottom = 7350000, top = 7380000),
              dx = 50, dy = 50, dt = "P5D", aggregation = 'min', resampling = 'average')

gdalcubes_options(parallel = 4)



# Colors

colr <- colorRampPalette(rev(brewer.pal(11, 'RdBu')))

## animation

RES = raster_cube(col, v) %>% #, mask = image_mask("SCL", values = 6, invert = T)) %>% 
  select_bands(assets) %>% 
  apply_pixel("(B05-B04)/(B05+B04)+1", "NDCI") %>% 
  apply_pixel(names="Chla", 
              FUN=function(x) {
                exp(-1.063 + 2.127*x[1]+1.63*x[1]^2)
              }) %>%
  gdalcubes::animate(zlim = c(0, 500), 
                     key.pos = 1,col=colr,
                     save_as = 'anim2.gif', fps = 1, na.color = '#FFFFFF')




RESrgb = raster_cube(col, v) %>% 
  select_bands(c("B02","B03", "B04")) %>% 
  gdalcubes::animate(rgb=3:1,zlim = c(0,1000),
                     save_as = 'animRGB2.gif', fps = 1, na.color = '#FFFFFF')




raster_cube(col, v, mask = image_mask("SCL", values = 6, invert = T)) %>% 
  select_bands(assets) %>% 
  apply_pixel(names="Chla", 
              FUN=function(x) {
                
                require(randomForest)
                RF_FINAL = readRDS('Outputs/rf_chla.R')
                
                x490 = x['B02']
                x560 = x['B03']
                x660 = x['B04']
                x705 = x['B05']
                x850 = x['B08']
                
                NDCI = (x['B05']-x['B04'])/(x['B05']+x['B04'])

                STACKED = data.frame(x490, x560, x660, x705, x850, NDCI)
                
                print(STACKED)
                names(STACKED) = c('x490', 'x560', 'x660', 'x705', 'x850', 'NDCI')
                
                predict(RF_FINAL, STACKED)
              
            }) %>%
  
  write_tif('teste5.tif')
  
  


