## Pre-processing of GLORIA data to predict Chl-a and TSS ####

# loading require packages

require(data.table)
require(dplyr)
require(terra)
require(mapview)
require(httr)
require(Metrics)
require(geodrawr)
require(svDialogs)
require(rstac)
require(wesanderson)
require(PerformanceAnalytics)
library(ggpmisc)
require(gdalcubes)
require(Metrics)
require(randomForest)
library(rasterVis)
require(RColorBrewer)

## Configure the Project Directories 

dir.create("Data")
dir.create("Outputs")
dir.create("Scripts")


###### Download GLORIA data ##########

URL = 'https://download.pangaea.de/dataset/948492/files/GLORIA-2022.zip'

# Before download, let`s set timeout to 200s (sometimes PANGAEA download is slow)

options(timeout=300)

# If the directory with files doesn't exist, let's download GLORIA data.

if(dir.exists('Data/GLORIA_2022/') == FALSE) {
  
  # Download
  download.file(URL, 'Data/GLORIA.zip')
  
  # Extract
  unzip(zipfile = 'Data/GLORIA.zip', exdir = 'Data')
  
}


##### Analyzing GLORIA data #######

meta_and_lab = fread("Data/GLORIA_2022/GLORIA_meta_and_lab.csv")
rrs = fread("Data/GLORIA_2022/GLORIA_Rrs.csv")

head(meta_and_lab)
head(rrs)

##### Plot for different concentrations #######

# High Chl-a

meta_and_lab[meta_and_lab$Chla > 1000, 'GLORIA_ID']

matplot(t(select(rrs, paste("Rrs_", 400:900, sep = ''))[rrs$GLORIA_ID == 'GID_7403',]), ylim = c(0,0.06),
        x= c(400:900), pch = 20, xlab = '', ylab = '', type = 'l')

# High TSS

meta_and_lab[meta_and_lab$TSS > 1000, 'GLORIA_ID']

matplot(t(select(rrs, paste("Rrs_", 400:900, sep = ''))[rrs$GLORIA_ID == 'GID_1805',]), ylim = c(0,0.06),
        x= c(400:900), pch = 20, xlab = '', ylab = '', type = 'l')


# High aCDOM

meta_and_lab[meta_and_lab$aCDOM440 > 15, 'GLORIA_ID']

matplot(t(select(rrs, paste("Rrs_", 400:900, sep = ''))[rrs$GLORIA_ID == 'GID_2468',]), ylim = c(0,0.0005),
        x= c(400:900), pch = 20, xlab = '', ylab = '', type = 'l')


# Band simulation ####

devtools::install_github("dmaciel123/BandSimulation")

require(bandSimulation)

spectra_formated = select(rrs, paste("Rrs_", 400:900, sep = '')) %>% t()

head(spectra_formated[1:10,1:10])

MSI_sim = msi_simulation(spectra = spectra_formated, 
                         point_name = rrs$GLORIA_ID)


#It simulates for Sentinel-2A and Sentinel-2B and gives the results in a list.
# Let's select only Sentinel-2A.

MSI = MSI_sim$s2a[,-1] %>% t() %>% data.frame()

head(MSI[,1:9])


# Add names to a collumn
MSI$GLORIA_ID = row.names(MSI)

head(MSI[,1:9])

# Change band names

names(MSI) = c('x440', "x490", 'x560', 'x660', "x705", 'x740', 'x780', 'x850', 'x865', "GLORIA_ID")


selection = filter(rrs, GLORIA_ID == 'GID_207')
selection.s = filter(MSI, GLORIA_ID == 'GID_207')
meta.s = filter(meta_and_lab, GLORIA_ID == 'GID_207')

##### Plot example #####

matplot(t(select(selection, paste("Rrs_", 400:900, sep = ''))[,]), ylim = c(0,0.05),
        x= c(400:900), pch = 20, xlab = '', ylab = '')

par(new=T)

matplot(t(selection.s[,-10]), x= c(440,490,560,660,705,740,780,842,860), pch = 20,
        ylim = c(0,0.05), xlim = c(400,900), col = 'red', cex = 2, xlab = 'Wavelength (nm)', 
        ylab = 'Rrs (sr-1)')

legend('topleft', legend = c(paste("Chl-a = ", meta.s$Chla),
                             paste("Secchi = ", meta.s$Secchi_depth)))

## Merge dataset and prepare to export ########

## Merge with water quality, lat long (By GLORIA_ID)

merged = merge(select(meta_and_lab, c('GLORIA_ID', 'Chla', "TSS", "Latitude", "Longitude")),
               MSI, by = "GLORIA_ID")

head(merged)


###### Index calculation and NA remove #######

# We want to model TSS and Chla concentration based on RF

merged = merged[is.na(merged$Chla) == FALSE, ]
merged = merged[is.na(merged$TSS) == FALSE, ]


merged = merged[(merged$Chla < 1000 & merged$Chla > 1), ]

merged = merged[(merged$TSS < 1000 & merged$TSS > 1), ]

#Index calculations

merged$NDCI = (merged$x705-merged$x660)/(merged$x705+merged$x660)+1
merged$Blue_red = (merged$x490-merged$x660)/(merged$x490+merged$x660)
merged$green_blue = (merged$x560-merged$x490)/(merged$x560+merged$x490)
merged$nir_red = (merged$x850-merged$x660)/(merged$x850+merged$x660)


# Check dimension

dim(merged)

####### Vizualize #######

summary(merged)

chart.Correlation(log(select(merged, c("Chla", 
                                   'TSS', 
                                   "NDCI", 
                                   'Blue_red',
                                   'green_blue', 
                                   'nir_red'))))


ggplot(merged, aes(x = log(NDCI+1), y = log(Chla))) +
  geom_point(color = 'black') +
  scale_y_continuous(limits = c(-2,10)) +
  scale_x_continuous(limits = c(0,1)) +
  theme_bw() + 
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "R2"))) +
  theme(panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text=element_text(family = "Tahoma"),
        axis.title = element_text(face="bold", size = 15),
        axis.text.x = element_text(colour="black", size = 15),
        axis.text.y = element_text(colour="black", size = 15),
        axis.line = element_line(size=2, colour = "black"),
        strip.text = element_text(size=15)) +
  theme(plot.margin = unit(c(4,4,4,4), "lines"))



vector = vect(merged, 
              geom = c('Longitude', 'Latitude'), 
              "EPSG:4326")


# plot map

mapview(sf::st_as_sf(vector),  zcol = 'Chla')


## Saving results

write.table(merged, file = 'Outputs/sentinel2_simulated_filtered.csv', row.names = F)






