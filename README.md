# Remote Sensing applied to freshwater ecosystems studies 

In this workshop, we gonna learn how to use Remote Sensing data applied to aquatic sciences. We gonna use in situ dataset available in the GLORIA dataset (Lehmann et al. 2023) to generate a Chorophyll-a concentration empirical model based on Normalized Difference Chlorophyll-a Index (NDCI) and Random Forests. Therefore, with a calibrated algorithm, we gonna apply the developed models to Sentinel-2/MSI Surface Reflectance data atmospherically corrected using Sen2Cor and available in the Miscrosoft Planetary Computer STAC platform. 

# Required software 

For running the scripts, we reccomend the attendes to install R and RSTUDIO. 

R could be downloaded here: 

RSTUDIO could be downloaded here:

We also encourage attendes to create a Microsoft Planetary Computer and Google Earth Engine accounts. Althought we're not gonna use these platforms directly, thei could be both used to process and work with satellite big data. Advantages of MPC is that they allow programming in R :)

# Required packages 

During the workshop, we will present the necessary packages. However, we encourage attendes to install it in advance. For that, just run in your R console:

```r

# Required packages

packages = c('data.table','dplyr','terra','mapview','httr','Metrics','geodrawr','
             svDialogs','rstac','wesanderson','PerformanceAnalytics',
             'ggpmisc','gdalcubes','Metrics','randomForest','rasterVis','RColorBrewer')

install.packages(packages)

```

# Other intallations and reccomendations

We don't need to download and/or install other softwares and download any other data. We gonna use R to download, organize and process the in situ and satellite data. 

