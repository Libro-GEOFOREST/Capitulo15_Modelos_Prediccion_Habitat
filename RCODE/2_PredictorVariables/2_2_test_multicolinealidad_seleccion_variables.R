#############################################################
## Modelos de predicción de hábitat en ciencias forestales ##
#############################################################

##########################################################
# Analisis de multicolinealidad y seleccion de variables #
##########################################################

# 1.- Cargar librerias basicas en el espacio de trabajo
#######################################################
library(raster)
library(rgdal)
library(spatstat)
library(spatialEco)
library(usdm) # Analisis de incertidumbre para modelos de distribucion de especies
library(ggplot2) # Creacion de elegantes y complejos plots graficos

# 2.- Definir el directorio de trabajo
######################################
setwd("C:/GitHub_repo")

# 3.- Rutas para importar/exportar archivos
###########################################
path2layers<-"./DATA/PredictorVariables/1_WorldClim_R" # Path para las variables bioclimaticas
path2hfp<- "./DATA/PredictorVariables/2_test_multicolinealidad_seleccion_variables/2_2_Otras_variables/hfp.asc" # Path para el raster de human footprint
path2ndvi<- "./DATA/PredictorVariables/2_test_multicolinealidad_seleccion_variables/2_2_Otras_variables/ndvi.asc" # Path para el raster de NDVI
path2work1<-"./DATA/PredictorVariables/2_test_multicolinealidad_seleccion_variables/2_1_Resultados/" # Path de resultados
path2work2<-"./DATA/PredictorVariables/2_test_multicolinealidad_seleccion_variables/2_1_Resultados/variables_seleccionadas/" # Path de variables seleccionadas

# 4.- Cargar las variables bioclimaticas y elevacion
####################################################
VAR_PI <- stack(list.files(path=path2layers,pattern='asc',full.names=TRUE))
names(VAR_PI)

# Estadistica descriptiva de las variables
summary(VAR_PI)

# Explorar la estructura de las variables
print(VAR_PI)

# Asignar sistema de coordenadas
crs(VAR_PI) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

# Plots de las variables en el espacio geografico
plot(VAR_PI)

# 5.- Eliminar bio_3, bio_14, bio_15 (Varela_et_al_2015_PLoSONE)
###############################################################
names(VAR_PI)
VAR_PI_DEF <- dropLayer(VAR_PI, c(7, 8, 14))
plot(VAR_PI_DEF)

# 6.- Cargar otras variables (e.g. hfp y ndvi)
##############################################
VAR_hfp <- raster(path2hfp)
plot(VAR_hfp)

VAR_ndvi <- raster(path2ndvi)
plot(VAR_ndvi)

# Explorar la estructura de las variables
print(VAR_PI_DEF)
print(VAR_hfp)
print(VAR_ndvi)

# Comprobar extent de los rasters
extent(VAR_PI_DEF)
extent(VAR_hfp)
extent(VAR_ndvi)

# Comprobar tamaño de pixel de los rasters
xres(VAR_PI_DEF)
xres(VAR_hfp)
xres(VAR_ndvi)

# Asignar sistema de coordenadas
crs(VAR_hfp) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
crs(VAR_ndvi) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

# 7.- Igualar la resolucion y extension de las variables
########################################################

# Igualamos los rasters (por ejemplo con la alt)
hfp_new<-resample(VAR_hfp,VAR_PI_DEF[[1]])
dim(hfp_new)<-c(105,166)
plot(hfp_new)

ndvi_new<-resample(VAR_ndvi,VAR_PI_DEF[[1]])
dim(ndvi_new)<-c(105,166)
plot(ndvi_new)

extent(VAR_PI_DEF[[1]])
extent(hfp_new)
extent(ndvi_new)

xres(VAR_PI_DEF[[1]])
xres(hfp_new)
xres(ndvi_new)

# 8.- Calcular el coeficiente de correlacion entre variables
############################################################

# Cargar todas las variables
VAR_PI <- stack(c(VAR_PI_DEF, hfp_new, ndvi_new))
plot(VAR_PI)

# Transformar los rasters en una tabla
var.df<-as.data.frame(VAR_PI)

# Eliminar valores nulos
var.df<-na.omit(var.df)

# Calcular el coeficiente de correlacion de Spearman (o Pearson) entre variables
var.cor<-cor(var.df, method = c("spearman"))

# Explorar los resultados de la correlacion
print(var.cor)

# Plot de correlaciones
library(corrplot) # Visualizacion grafica de una matriz de correlacion, intervalo de confianza

corrplot_number<-corrplot(var.cor, method = "number") # Plot de los indices numericos de correlacion
corrplot_circle<-corrplot(var.cor, type = "upper", method = "circle") # Plot de los indices de correlacion en circulos

# Transformar las correlaciones en un dataframe
cor.df<-as.data.frame(var.cor)

# Explorar la estructura del dataframe
head(cor.df)

# Exportar el dataframe (solo la diagonal superior)
lower<-var.cor
lower[lower.tri(var.cor, diag=TRUE)]<-""
lower.df<-as.data.frame(lower)
lower.df

# Guardar la tabla
write.csv(lower.df, paste(path2work1, "var_correlations.csv"),sep=",")

# Transformar la matriz de correlacion en un matriz de distancias
var.dist <- abs(as.dist(cor.df))

# Calcular el dendrograma basando en distancias (less distance = more correlation)
var.cluster <- hclust(1-var.dist)

# Plot del dendrograma
plot(var.cluster)

# Seleccionar variables con una correlacion < 0.8
abline(h=0.2, lty=2, lwd=2, col="red")

# Explorar el plot y seleccionar las variables no correlacionadas
var.selected <- c("bio5", "bio7", "bio8", "bio11", "bio12", "alt", "hfp", "ndvi")

# Construir una nueva tabla con las variables seleccionadas
var.df2<-var.df[ , var.selected]

# Calcular el coeficiente de correlacion de Spearman (o Pearson) entre variables seleccionadas
var.cor<-cor(var.df2, method = c("spearman"))

# Plot de la tabla de correlacion
corrplot_circle<-corrplot(var.cor, type = "upper", method = "number")

# Podria haber variables que serian combinacion lineal de otras

# 7.- Calcular el Variance Inflation Factor (VIF)
#################################################

# VIF - I. Calcular VIF basandonos en las funciones `vifcor´ y `vifstep´
VIF1 <- vifcor(VAR_PI_DEF, th=0.8) # Encuentra un par de variables que tengan la correlacion lineal maxima (mayor que `th´), y excluye aquella que tenga mayor VIF
VIF1

VIF2 <- vifstep(VAR_PI_DEF, th=3) # Calcula VIF para todas las variables, excluye una con VIF mas alto (mayor que `th´), repite el procedimiento hasta que no queden variables con VIF mayor que `th´
VIF2

VIF3 <- exclude(VAR_PI_DEF,VIF2) # Excluye variables especificadas en un objeto VIF
VIF3

# Plot de las variables con baja colinealidad
plot(VIF3)

# Aplicamos el VIF sobre las variables menos correlacionadas
# VIF - II
result.vif<-vif(var.df2) #  Otra manera de calcular el VIF con la funcion `vif´
result.vif

# Eliminar las de mayor valor de VIF
# bio11
var.df2$bio11<-NULL
result.vif<-vif(var.df2)
result.vif

# bio5
var.df2$bio5<-NULL
result.vif<-vif(var.df2)
result.vif

# Representar un boxplot horizontal basico para el VIF
p <- ggplot(data=result.vif, aes(x=Variables, y=VIF), fill = activity) +
  ggtitle("VIF") +
  geom_bar(stat="identity") +
  coord_flip() +
  theme(
    plot.title = element_text(size=18, face="bold"),
    axis.title=element_text(size=18, face="bold"),
    legend.title=element_text(size=14, face="bold"),
    legend.text=element_text(size=12),
    axis.text=element_text(size=14, colour="black"))

p

# Obtener los nombres de las variables
var.selected<-names(var.df2)
var.selected

# Unir todos los rasters en un unico objeto
var.def<-brick(VAR_PI[[var.selected]])

# Plot variables
plot(var.def)

# Guardar las variables seleccionadas para calibrar los modelos
writeRaster(x=var.def, path2work2, names(var.def), bylayer=TRUE, format="ascii", overwrite=TRUE)

###############################################################