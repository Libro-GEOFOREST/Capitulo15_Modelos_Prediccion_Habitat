#############################################################
## Modelos de prediccion de habitat en ciencias forestales ##
#############################################################

#####################
# FILTERING RECORDS #
#####################

# Cargar librerias basicas en el espacio de trabajo
###################################################
library(raster)
library(rgdal)
library(spatstat) # Analisis de patrones de puntos espaciales, ajuste de modelos, simulacion, tests, ...
library(spatialEco) # Analisis espacial y utilidades de modelado

# Definir el directorio de trabajo
##################################
setwd("C:/GitHub_repo")

# Rutas para importar/exportar archivos
#######################################
path2sp <- "./DATA/ResponseVariable/Pb.rgbif.csv" # Path para la especie
path2sp.shp <- "./DATA/ResponseVariable/Pb.shp" # Path para la especie (shp)
path2IP.shp <- "./DATA/shp/PI.shp" # Path para el area de estudio (shp)
path2work<-"./DATA/ResponseVariable/filtering_records/" # Path para los resultados

# Cargar los datos
##################

# DF de la especie
sp.df <- read.csv(path2sp, header=T, sep=',', fill=TRUE, check.names=TRUE, stringsAsFactors=FALSE)
head(sp.df, n=5)

# Eliminar campos innecesarios
sp.df[1] <- NULL # Eliminar el campo X
head(sp.df, n=5)

# Cargar el shp de la especie
sp.shp <- readOGR(path2sp.shp)
plot(sp.shp)

# Cargar el shp del area de estudio
area.shp<-readOGR(path2IP.shp)
plot(area.shp)

# Comprobar la clase de los objetos
class(sp.shp) # SpatialPointsDataFrame
class(area.shp) # SpatialPolygonsDataFrame

# Transformar el formato de los objetos
sp.sp<-as(sp.shp,"SpatialPoints") # A SpatialPoints la especie
plot(sp.sp)

area_poly<-as(area.shp, "SpatialPolygons") # A SpatialPolygons el area
plot(area_poly)

# Plot puntos y area de estudio
plot(area_poly, col="gray80")
points(sp.sp,add=TRUE, col = "red")

# Testar agrupamiento de puntos (clusters) en el dataset de las ocurrencias de la especie
###############################################################

# Calculo del Average Nearest Neighbour Index (NNI)
nni(sp.sp, win = "extent")

# Comprobar que el NNI:
# NNI > 1 = ordering (no hay clusters)
# NNI < 1 = clustering

# Tras comprobar que existe agrupamiento entre puntos, se procede a filtrarlos

# Reduccion (thinning) espacial de registros de ocurrencias para el uso en modelos ecológicos
###########################################################

# Instalar y cargar la libreria spThin

# windows systems require Rtool to compile source files: 
# https://cran.r-project.org/bin/windows/Rtools/
if (!require(devtools))
  install.packages('devtools')
devtools:::install_github('mlammens/spThin')

library(spThin)

# Calcular distancia minima (en km) entre dos puntos
# Para ello se selecciona la $expected.mean.distance después de aplicar el NNI

# Obtener distancia en unidades métricas
0.1102866*111 # 111km ~ 1 grado

# Inspeccionar el DF original
head(sp.df, n = 5)

# Run spThin::thin sobre el DF
thin.sp <-
  thin( loc.data = sp.df, 
        lat.col = "lat", long.col = "lon", 
        spec.col = "Pb", 
        thin.par = 10, # Aproximar a la distancia mínima entre puntos calculada en pasos anteriores (~ 10km)
        reps = 10, 
        locs.thinned.list.return = TRUE, 
        write.files = TRUE, 
        max.files = 3, 
        out.dir = paste(path2work, "pb_thinned"),
        write.log.file = TRUE,
        log.file = paste(path2work, "pb_thinned_log_file.txt"))
        
summaryThin(thin.sp, show = TRUE)

# Plot para visualizar las replicas y obtener el numero optimo de registros
plotThin(thin.sp)

# Pulsar <Enter> en la consola para ver los graficos

###############################################################
###############################################################