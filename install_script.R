#############################################################
## Modelos de prediccion de habitat en ciencias forestales ##
#############################################################

###########################################################
# Script para instalacion de los paquetes de R necesarios #
###########################################################

# Actualizar a la ultima version de R y R-Studio

install.packages(c(
'rgbif', 'raster', 'rgdal', 'spatstat', 'spatialEco', 'spThin', 'usdm', 'ggplot2', 'corrplot', 'magrittr', 'dismo', 'rgeos', 'biomod2', 'RColorBrewer', 'dplyr', 'sp', 'rasterVis', 'rJava'), dependencies = T)

# Check paquetes instalados
library(rgbif)
library(raster)
library(rgdal)
library(spatstat)
library(spatialEco)
library(spThin)
library(usdm)
library(ggplot2)
library(corrplot)
library(magrittr)
library(dismo)
library(rgeos)
library(biomod2)
library(RColorBrewer)
library(dplyr)
library(sp)
library(rasterVis)
library(rJava)