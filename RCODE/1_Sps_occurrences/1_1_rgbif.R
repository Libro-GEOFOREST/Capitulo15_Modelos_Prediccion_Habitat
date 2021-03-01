#############################################################
## Modelos de prediccion de habitat en ciencias forestales ##
#############################################################

######################################
# Descarga de ocurrencias desde GBIF #
######################################

# Instalar packages (ver script de instalacion)
###############################################

# Cargar librerias basicas en el espacio de trabajo
###################################################
library(rgbif) # Incluye funciones para buscar datos de biodiversidad en la base de datos GBIF
library(raster) # Lectura, escritura, manipulacion, analisis y modelado de datos espaciales en formato 'raster'
library(rgdal) # Manipulacion de datos 'geoespaciales'

# Definir el directorio de trabajo
##################################
setwd("C:/GitHub_repo")

# Rutas para importar/exportar archivos
#######################################
path2PI.shp <- "./DATA/shp/PI.shp" # Para cargar el shapefile de la Peninsula Iberica
path2Pb.shp <- "./DATA/ResponseVariable/Pb.shp" # Para guardar el shapefile de las ocurrencias de la especie modelo
path2Pb.df <- "./DATA/ResponseVariable/Pb.rgbif.csv" # Para guardar la tabla de las ocurrencias de la especie modelo

# Descarga de ocurrencias de 'Pyrus bourgaeana' en la Penisnula Iberica (SP y PT)
###############################################################
pyrus_ES <- occ_search(scientificName = "Pyrus bourgaeana",
                        hasCoordinate = TRUE, country = "ES",
                        return = "data",limit=2000)

pyrus_PT <- occ_search(scientificName = "Pyrus bourgaeana",
                        hasCoordinate = TRUE, country = "PT",
                        return = "data",limit=2000)

# Nombre de los campos (columnas) del dataframe (DF)
colnames(pyrus_ES$data) # Spain
colnames(pyrus_PT$data) # Portugal

# Seleccionar las columnas de interes
pb.ES <- subset(pyrus_ES$data,
                select = c('species',
                           'decimalLatitude',
                           'decimalLongitude'))

pb.PT <- subset(pyrus_PT$data,
                select = c('species',
                           'decimalLatitude',
                           'decimalLongitude'))

# Unir ambas tablas (Peninsula Iberica)
pb.PI <- rbind.data.frame(pb.ES,pb.PT)

# Comprobar el encabezado del DF
head(pb.PI, n = 5)

# Comprobar estructura y clase del DF
str(pb.PI)
class(pb.PI)

# Transformar a DF (tabla)
pb.PI.DF <- as.data.frame(pb.PI)

# Eliminar posibles sesgos en el DF
###################################

# Eliminar NAs
pb.PI.DF <- pb.PI.DF[complete.cases(pb.PI.DF),]
str(pb.PI.DF)

# Identificar registros duplicados y eliminarlos
dups <- duplicated(pb.PI.DF)
dups
pb.PI.DF <- unique(pb.PI.DF)
dim(pb.PI.DF)

# Eliminar coordinadas 0,0
pb.PI.DF<-subset.data.frame(pb.PI.DF, pb.PI.DF$decimalLongitude!=0 | pb.PI.DF$decimalLatitude!=0)

# Añadir información adicional
##############################

# Añadir el campo Pb = 1 a la tabla
Presence<-pb.PI.DF[,c(1,2,3)]
Presence$Pb<-1
summary(Presence)
head(Presence, n = 5)

# Renombrar los campos de las coordenadas como "lon" y "lat"
names(Presence)[names(Presence) == "decimalLongitude"] <- "lon"
names(Presence)[names(Presence) == "decimalLatitude"] <- "lat"
head(Presence, n = 5)

# Convertir el DF en vector de puntos (shapefile)
Pb <- Presence # Asignamos otro nombre al DF
Pb <- SpatialPoints(Pb[,c('lon', 'lat')]) # Asignamos los nombres de las columnas con las coordenadas
Pb # No tiene sistema de referencia espacial - crs

# Asignamos sistema de referencia al DF (e.g. WGS84)
proj4string(Pb)<- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") # Asignamos el sistema de referencia
Pb

# Cargar shapefiles
###################
PI.shp <- readOGR(path2PI.shp) # Colocar el nombre del path al shp

# Obtenemos el plot para la PI
plot(PI.shp)

# Superponemos los puntos de Pb
plot(Pb, add = TRUE)

# Guardar el shapefile de manera externa
raster::shapefile(Pb, path2Pb.shp, overwrite=TRUE) # Colocar el nombre del path al shp

# Guardar DF como tabla csv (separación por comas)
write.csv(Presence, path2Pb.df, sep=",") # Colocamos el nombre del path para el csv

###############################################################
###############################################################