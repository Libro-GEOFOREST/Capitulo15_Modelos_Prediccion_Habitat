#############################################################
## Modelos de prediccion de habitat en ciencias forestales ##
#############################################################

############################################################
# Descarga de variables (bio)climaticas (Worldclim ~10km2) #
############################################################

# 1.- Cargar librerias basicas en el espacio de trabajo
#######################################################
library(raster)

# 2.- Definir el directorio de trabajo
######################################
setwd("C:/GitHub_repo")

# 3.- Rutas para importar/exportar archivos
###########################################
path2bio<-"./DATA/PredictorVariables/1_WorldClim_R" # Path para guardar layers de las variables bioclimaticas
path2alt<-"./DATA/PredictorVariables/1_WorldClim_R/alt.asc" # Path para guardar layer de elevacion

# 4.- Descarga de variables bioclimaticas (10km2)
#################################################???
clim10km<- getData('worldclim', var = "bio", res = 5)
plot(clim10km)

# Explorar la resolucion y extension de las variables
extent(clim10km)
xres(clim10km)
yres(clim10km)

# Descargar variable elevacion (altitud)
alt <- getData('alt', country='ESP', mask=FALSE)
plot(alt)

# Explorar la resolucion y extension de la variable
extent(alt)
xres(alt)
yres(alt)

# 5.- Igualar la resolucion y extension de las variables
########################################################

# Cortar las bioclimaticas con la elevacion
VAR_PI<-crop(clim10km, alt)
plot(VAR_PI)

# Igualar los rasters
altitud<-resample(alt,VAR_PI[[14]]) # En este caso elegimos cualquiera de las bioclimaticas
dim(altitud)<-c(105,166)

# Comprobar extent de los rasters
extent(VAR_PI)
extent(altitud)

# Comprobar tamaño de pixel de los rasters
xres(VAR_PI)
yres(altitud)
xres(VAR_PI)
yres(altitud)

# 6.- Guardar los ficheros (*asc, *tif, ...) de las variables predictoras 
#############################################################

# 19 bioclimaticas
for (k in 1:19){
  name<-paste(path2bio, k, ".asc", sep="")
  writeRaster(VAR_PI[[k]], name, overwrite=TRUE)
  
}

# Elevacion
writeRaster(altitud, path2alt, overwrite=TRUE)

################################################################